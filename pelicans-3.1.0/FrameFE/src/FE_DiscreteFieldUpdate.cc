/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#include <FE_DiscreteFieldUpdate.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Variable.hh>
#include <doubleVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>

FE_DiscreteFieldUpdate const*
FE_DiscreteFieldUpdate::PROTOTYPE = new FE_DiscreteFieldUpdate() ;

//---------------------------------------------------------------------------
FE_DiscreteFieldUpdate:: FE_DiscreteFieldUpdate( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_DiscreteFieldUpdate" )
   , FIELD( 0 )
   , FIELD_LEVEL( PEL::bad_index() )
   , cFE( 0 )
   , VALUE_TYPE( FE_DiscreteFieldUpdate::from_analytic )
   , VALUE( 0 )
   , COMPO( 0 )
   , COMPO_LEVEL( PEL::bad_index() )
   , PARAMETER( 0 )
{
}

//---------------------------------------------------------------------------
FE_DiscreteFieldUpdate*
FE_DiscreteFieldUpdate:: create_replica(
                             PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_DiscreteFieldUpdate* result =
      new FE_DiscreteFieldUpdate( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_DiscreteFieldUpdate:: FE_DiscreteFieldUpdate(
                             PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , FIELD(
      dom->set_of_discrete_fields()->item( exp->string_data( "field_name" ) ) )
   , FIELD_LEVEL( (size_t) exp->int_data( "field_level" ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , VALUE_TYPE( FE_DiscreteFieldUpdate::from_analytic )
   , VALUE( 0 )
   , COMPO( 0 )
   , COMPO_LEVEL( PEL::bad_index() )
   , PARAMETER( 0 )
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: FE_DiscreteFieldUpdate" ) ;

   check_field_storage_depth( FIELD, FIELD_LEVEL ) ;

   cFE->require_field_calculation( FIELD, PDE_LocalFE::node ) ;

   // Parallel computation : halo color has to be considered :
   if( cFE->is_excluded( GE_Color::halo_color() ) )
   {
      cFE->include_color( GE_Color::halo_color() ) ;
   }

   // Value strategy :
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "DOFs_values" ) ;
   std::string const& type = ee->string_data( "type" ) ;
   if( type == "from_analytic" )
   {
      VALUE_TYPE = FE_DiscreteFieldUpdate::from_analytic ;
      VALUE = ee->abstract_data( this, "value" ) ;
   }
   else if( type == "from_field_composition" )
   {
      VALUE_TYPE = FE_DiscreteFieldUpdate::from_composition ;
      COMPO_LEVEL = (size_t) ee->int_data( "fields_level" ) ;
      COMPO = dom->set_of_field_compositions()->item(
                       ee->string_data( "field_composition_name" ) ) ;
      COMPO->start_variable_iterator() ;
      for( ; COMPO->valid_variable() ; COMPO->go_next_variable() )
      {
         PDE_DiscreteField const* ff = COMPO->variable() ;
         cFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
         if( COMPO_LEVEL >= ff->storage_depth() )
         {
            PEL_Error::object()->raise_plain(
                "*** FE_DiscreteFieldUpdate<"+FIELD->name()+"> error\n"
                   "    bad \"fields_level\" regarding the storage depth\n"
                   "    of the field \""+ff->name()+"\"" ) ;
         }
      }
   }
   else if( type == "from_parameter" )
   {
      VALUE_TYPE = FE_DiscreteFieldUpdate::from_parameter ;
      PARAMETER = prms->item( ee->string_data( "parameter_name" ) ) ;
      PARAMETER->transfer_cell_calculation_requirements(
                                           cFE, FE_Parameter::Val ) ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         ee, "type",
         "   - \"from_analytic\"\n"
         "   - \"from_field_composition\"\n"
         "   - \"from_parameter\"\n" ) ;
   }
   ee->destroy() ; ee = 0 ;
}

//---------------------------------------------------------------------------
FE_DiscreteFieldUpdate:: ~FE_DiscreteFieldUpdate( void )
//---------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//---------------------------------------------------------------------------
void
FE_DiscreteFieldUpdate:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: do_before_time_stepping" ) ;

   FE_OneStepIteration::do_before_time_stepping( t_it ) ;

   // Verify datas :
   {
      if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_analytic )
      {
         cFE->start() ;
         cFE->set_calculation_point( cFE->polyhedron()->center() ) ;
         PEL_Context const* ct = context_at_pt( t_it, cFE ) ;
         if( !VALUE->value_can_be_evaluated( ct ) )
         {
            PEL_Error::object()->raise_plain(
               "*** FE_DiscreteFieldUpdate<"+FIELD->name()+" error :\n"
               "    the expression \"value\" cannot be evaluated\n"
               "    (check the variables)" ) ;
         }
         if( VALUE->data_type() != PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_plain(
               "*** FE_DiscreteFieldUpdate<"+FIELD->name()+" error :\n"
               "    the expression \"value\" has not the good type\n"
               "    (doubleVector is expected)" ) ;
         }
         doubleVector const& val = VALUE->to_double_vector( ct ) ;
         if( val.size() != FIELD->nb_components() )
         {
            PEL_Error::object()->raise_plain(
               "*** FE_DiscreteFieldUpdate<"+FIELD->name()+" error :\n"
               "    the expression \"value\" has not the good dimension" ) ;
         }
      }
   }
   
   if( t_it->time() == t_it->initial_time() ) // Not a restoration
   {
      start_total_timer(
         "FE_DiscreteFieldUpdate<"+FIELD->name()+">:: do_before_time_stepping" ) ;
      set_field_values( 0, t_it ) ;
      for( size_t level = 1 ; level<FIELD->storage_depth() ; ++level )
      {
         FIELD->copy_DOFs_value( 0, level ) ;
      }
      stop_total_timer() ;
   }
}

//---------------------------------------------------------------------------
void
FE_DiscreteFieldUpdate:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
    
   start_total_timer(
      "FE_DiscreteFieldUpdate<"+FIELD->name()+">:: do_one_inner_iteration" ) ;
   set_field_values( FIELD_LEVEL, t_it ) ;
   stop_total_timer() ;
}

//----------------------------------------------------------------------------
void
FE_DiscreteFieldUpdate:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: print" ) ;

   FE_OneStepIteration::print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   os << s << "field_name : \""
      << FIELD->name()  << "\"" << std::endl ;
   os << s << "field_level : "
      << FIELD_LEVEL  << std::endl ;
   os << s << "field_value : " ;
   if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_analytic )
   {
      os << "\"from_analytic\"" << std::endl ;
   }
   else if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_composition )
   {
      os << "\"from_field_composition\"" << std::endl ;
      os << s << "   field_composition_name : \""
         << COMPO->name() << "\"" << std::endl ;
      os << s << "   fields_level : " << COMPO_LEVEL << std::endl ;
   }
   else if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_parameter )
   {
      os << "\"from_parameter\"" << std::endl ;
      os << s << "   parameter_name : \""
         << PARAMETER->name() << "\"" << std::endl ;
   }
   else
   {
      os << "\"undefined\"" << std::endl ;
   }
}

//---------------------------------------------------------------------------
void
FE_DiscreteFieldUpdate:: set_field_values( size_t f_level,
                                           FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: set_field_values" ) ;
   PEL_CHECK( t_it != 0 ) ;
   size_t const nb_comps = FIELD->nb_components() ;
   doubleVector field_value( nb_comps ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      for( size_t i_node=0 ; i_node<cFE->nb_local_nodes( FIELD ) ; ++i_node )
      {
         size_t const n = cFE->global_node( FIELD, i_node ) ;
         cFE->set_calculation_point(
                            cFE->local_node_location( FIELD, i_node ) ) ;
         compute_field_value_at_pt( t_it, cFE, field_value ) ;
         for( size_t ic=0 ; ic<nb_comps ; ++ic )
         {
            FIELD->set_DOF_value( f_level, n, field_value(ic), ic ) ;
         }
      }
   }
}

//---------------------------------------------------------------------------
void
FE_DiscreteFieldUpdate:: compute_field_value_at_pt( FE_TimeIterator const* t_it,
                                                    PDE_LocalFEcell const* fe,
                                                    doubleVector& result )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: compute_field_value_at_pt" ) ;
   PEL_CHECK( t_it != 0 ) ;
   PEL_CHECK( fe != 0 ) ;
   PEL_CHECK( fe->is_valid() ) ;
   PEL_CHECK( fe->calculation_point() != 0 ) ;
   PEL_CHECK( result.size() == FIELD->nb_components() ) ;
   
   if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_analytic )
   {
      PEL_Context const* ct = context_at_pt( t_it, fe ) ;
      result = VALUE->to_double_vector( ct ) ;
   }
   else if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_composition )
   {
      COMPO->start_variable_iterator() ;
      for( ; COMPO->valid_variable() ; COMPO->go_next_variable() )
      {
         PDE_DiscreteField const* f = COMPO->variable() ;
         for( size_t i=0 ; i<f->nb_components() ; ++i )
         {
            double const val = fe->value_at_pt( f, COMPO_LEVEL, i ) ;
            COMPO->set_variable_value( f, i, val ) ;
         }
      }
      COMPO->compute() ;
      for( size_t ic=0 ; ic<result.size() ; ++ic )
      {
         result(ic) = COMPO->value( ic ) ;
      }
   }
   else if( VALUE_TYPE == FE_DiscreteFieldUpdate::from_parameter )
   {
      for( size_t ic=0 ; ic<result.size() ; ++ic )
      {
         result(ic) = PARAMETER->cell_value_at_pt( t_it, fe, ic ) ;
      }
   }
}

//---------------------------------------------------------------------------
PEL_Context const*
FE_DiscreteFieldUpdate:: context_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DiscreteFieldUpdate:: context_at_pt" ) ;
   PEL_CHECK( t_it != 0 ) ;
   PEL_CHECK( fe != 0 ) ;
   PEL_CHECK( fe->is_valid() ) ;
   PEL_CHECK( fe->calculation_point() != 0 ) ;

   static PEL_ContextSimple* result = 0 ;
   static PEL_Double* time = 0 ;
   static PEL_DoubleVector* coords = 0 ;
   if( result == 0 )
   {
      result = PEL_ContextSimple::create( PEL_Root::object() ) ;
      coords = PEL_DoubleVector::create( result, doubleVector( 0 ) ) ;
      result->extend( PEL_Variable::object( "DV_X" ), coords ) ;
      time = PEL_Double::create( result, 0. ) ;
      result->extend( PEL_Variable::object( "DS_T" ), time ) ;
   }
   coords->set( fe->calculation_point()->coordinate_vector() ) ;
   if( t_it->is_started() )
   {
      time->set( t_it->time() ) ;
   }
   else
   {
      time->set( t_it->initial_time() ) ;
   }
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->has_variable( PEL_Variable::object( "DV_X" ) ) ) ;
   PEL_CHECK_POST( result->has_variable( PEL_Variable::object( "DS_T" ) ) ) ;  
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
FE_DiscreteFieldUpdate:: invariant( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( FE_OneStepIteration::invariant() ) ;
   return( true ) ;
}
