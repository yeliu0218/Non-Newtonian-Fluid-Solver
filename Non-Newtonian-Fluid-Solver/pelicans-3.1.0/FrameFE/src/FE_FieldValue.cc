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

#include <FE_FieldValue.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_PointInGridFE.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <FE.hh>
#include <FE_FieldCompositionParameter.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::ios_base ; using std::setprecision ; using std::setw ;
using std::endl ;

FE_FieldValue const* FE_FieldValue::PROTOTYPE = new FE_FieldValue() ;

struct FE_FieldValue_ERROR
{
   static void n0( std::vector<GE_Point const*> const& pts,
                   boolVector const& pt_found ) ;
   static void n1( std::string const& field_name ) ;
   static void n2( PEL_ModuleExplorer const* exp,
                   std::string const& entry_name, size_t nb_dims ) ;
   static void n3( PEL_ModuleExplorer const* exp,
                   std::string const& entry_name, size_t nb_dims ) ;
   static void n4( PDE_DiscreteField const* ff, 
                   std::string const& module_name,
                   size_t comp ) ;
   static void n5( size_t the_size, size_t the_comp ) ;
   static void n6( std::string const& file_name ) ;
   static void n7( PEL_ModuleExplorer const* exp,
                   std::string const& entry_name, int nb_points ) ;
   static void n8( PDE_DiscreteField const* ff ) ;
   static void n9( std::string const& param_name ) ;
   static void n10( FE_Parameter const* p, 
                    std::string const& module_name,
                    size_t comp ) ;
   static void n11( std::string const& compo_name ) ;
   static void n12( FE_Parameter const* p, 
                    std::string const& module_name,
                    size_t comp ) ;
} ;

size_t const CELL_VALUES = 0 ;
size_t const AT_PTS = 1 ;

//----------------------------------------------------------------------
FE_FieldValue:: FE_FieldValue( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_FieldValue" )
   , NB_SP_DIMS( PEL::bad_index() )
   , REMOVE_PTS( false )
   , POINTS( std::vector<GE_Point const*>( 0 ) )
   , ABS_PTS( doubleVector( 0 ) )
   , PIG( 0 )
   , cFE( 0 )
   , PT_CELL( size_t_vector( 0 ) )
   , NB_VALUES( PEL::bad_index() )
   , REDO( false )
   , FIELDS( 0 )
   , FIELD_COMPS( size_t_vector( 0 ) )
   , PARAMS( 0 )
   , PARAM_TYPES( size_t_vector( 0 ) )
   , PARAM_COMPS( size_t_vector( 0 ) )
   , EXPRS( 0 )
   , EXPR_COMPS( size_t_vector( 0 ) )
   , EXPR_DIMS( size_t_vector( 0 ) )
   , COORDS( 0 )
   , TT( 0 )
   , SAVING( one_file )
   , BANNER( false )
   , OFILENAME()
   , RSNAME()
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( doubleVector( 0 ) )
{
}

//----------------------------------------------------------------------
FE_FieldValue*
FE_FieldValue:: create_replica( PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_FieldValue* result = new FE_FieldValue( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_FieldValue:: FE_FieldValue( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , NB_SP_DIMS( dom->nb_space_dimensions() )
   , REMOVE_PTS( false )
   , POINTS( std::vector<GE_Point const*>( 0 ) )
   , ABS_PTS( doubleVector( 0 ) )
   , PIG( PDE_PointInGridFE::create( this, dom ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , PT_CELL( size_t_vector( 0 ) )
   , NB_VALUES( PEL::bad_index() )
   , REDO( dom->adapter_CHARMS() != 0 )
   , FIELDS( 0 )
   , FIELD_COMPS( size_t_vector( 0 ) )
   , PARAMS( 0 )
   , PARAM_TYPES( size_t_vector( 0 ) )
   , PARAM_COMPS( size_t_vector( 0 ) )
   , EXPRS( 0 )
   , EXPR_COMPS( size_t_vector( 0 ) )
   , EXPR_DIMS( size_t_vector( 0 ) )
   , COORDS( 0 )
   , TT( 0 )
   , SAVING( one_file )
   , BANNER( false )
   , OFILENAME()
   , RSNAME()
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( doubleVector( 0 ) )
{
   PEL_LABEL( "FE_FieldValue:: FE_FieldValue" ) ;

   if( cFE->is_excluded(  GE_Color::halo_color() ) )
   {
      cFE->include_color( GE_Color::halo_color() ) ;
   }
   set_fields_and_expressions( dom, prms, exp ) ;
   set_points( dom, exp ) ;
   set_post_processing( dom, exp ) ;
}

//----------------------------------------------------------------------
FE_FieldValue:: ~FE_FieldValue( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_FieldValue:: do_before_time_stepping( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   // Search for NB_VALUES:
   {
      cFE->start() ;
      cFE->set_calculation_point( cFE->polyhedron()->center() ) ;
      NB_VALUES = 0 ;
      for( size_t i=0; i<FIELDS.size(); ++i )
      {
         PDE_DiscreteField const* F =  FIELDS[i] ;
         if( FIELD_COMPS(i) == PEL::bad_index() )
         {
            NB_VALUES += F->nb_components() ;
         }
         else
         {
            NB_VALUES += 1 ;
         }
      }
      for( size_t i=0; i<PARAMS.size(); ++i )
      {
         FE_Parameter const* p = PARAMS[i] ;
         if( PARAM_COMPS(i) == PEL::bad_index() )
         {
            NB_VALUES += p->nb_components() ;
         }
         else
         {
            NB_VALUES += 1 ;
         }
      }
      if( !EXPRS.empty() )
      {
         cFE->start() ;
         set_context( cFE->polyhedron()->center(), t_it->time() ) ;
         for( size_t i=0 ; i<EXPRS.size() ; ++i )
         {
            doubleVector const& value = EXPRS[i]->to_double_vector() ;
            EXPR_DIMS(i) = value.size() ;
            if( EXPR_COMPS(i) == PEL::bad_index() )
            {
               NB_VALUES += value.size() ;
            }
            else
            {
               NB_VALUES += 1 ;
               if( EXPR_COMPS(i) >= value.size() )
               {
                  FE_FieldValue_ERROR::n5( value.size(), EXPR_COMPS(i) ) ;
               }
            }
         }
      }   
   }

   // Initialize post_processing:
   if( SAVING == one_file )
   {
      if( PEL_Exec::communicator()->rank() == 0 )
      {
         initialize_one_file() ;
      }
   }
   if( SAVING_TIMES.size() != 0 )
   {
      if( !t_it->table_of_times_is_valid( SAVING_TIMES ) )
      {
         t_it->raise_invalid_table_of_times( "FE_FieldValue", "saving_times", 
                                             SAVING_TIMES ) ;
      }
      SAVING_NUMBER = 0 ;
      if( SAVING_TIMES(0) == t_it->time() )
      {
         save_field_value( t_it, 0 ) ;
      }
      NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
   }
}

//----------------------------------------------------------------------
void
FE_FieldValue:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------
void
FE_FieldValue:: do_after_time_adaptation( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: do_after_time_adaptation" ) ;
   PEL_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;
      
   if( SAVING_TIMES.size() != 0 )
   {
      if( FE_TimeIterator::greater_or_equal( t_it->time(),
                                             NEXT_SAVING_TIME ) )
      {
         start_total_timer( "FE_FieldValue:: do_after_time_adaptation" ) ;
         
         save_field_value( t_it, 0 ) ;
         NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
         
         stop_total_timer() ;      
      }
   }
}

//----------------------------------------------------------------------
void
FE_FieldValue:: save_other_than_time_and_fields( 
                                            FE_TimeIterator const* t_it,
                                            PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   if( SAVING_TIMES.size() == 0 )
   {
      start_total_timer( "FE_FieldValue:: save_other_than_time_and_fields" ) ;
      
      save_field_value( t_it, rs ) ;
      
      stop_total_timer() ;
   }      
}

//----------------------------------------------------------------------
void
FE_FieldValue:: save_field_value( FE_TimeIterator const* t_it,
                                  PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: save_field_value" ) ;
   PEL_CHECK( IMPLIES( SAVING == result_saver, rs != 0 ) ) ;
   
   doubleArray2D values( POINTS.size(), NB_VALUES ) ;
   
   compute_values( t_it, values ) ;

   if( SAVING == one_file )
   {
      if( PEL_Exec::communicator()->rank() == 0 )
      {
         save_values_in_one_file( t_it->time(), values ) ;
      }
   }
   else if( SAVING == separated_files )
   {
      if( PEL_Exec::communicator()->rank() == 0 )
      {
         size_t const cycle =
            ( rs != 0 ? rs->cycle_number() : ++SAVING_NUMBER ) ;
         save_values_in_separated_files( cycle, t_it->time(), values ) ;
      }
   }
   else if( SAVING == result_saver )
   {
      rs->save_variable( values, RSNAME ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: set_fields_and_expressions( PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: set_fields_and_expressions" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( prms != 0 ) ;
   PEL_CHECK( exp != 0 ) ;

   // Fields to be saved:
   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "fields" ) ;
      sexp->start_module_iterator() ; 
      for( ; sexp->is_valid_module(); sexp->go_next_module() )
      { 
         PEL_ModuleExplorer const* ssexp = sexp->create_subexplorer( 0 ) ;
         std::string field_name = ssexp->string_data( "name" ) ;
         if( !dom->set_of_discrete_fields()->has( field_name ) )
            FE_FieldValue_ERROR::n1( field_name ) ;
         PDE_DiscreteField const* F = 
            dom->set_of_discrete_fields()->item( field_name ) ;
         cFE->require_field_calculation( F, PDE_LocalFE::N ) ;
         FIELDS.push_back( F ) ;
         size_t comp = PEL::bad_index() ;
         if( ssexp->has_entry( "component" ) )
         {
            comp = ssexp->int_data( "component" ) ;
            if( comp >=  F->nb_components() )
            {
               FE_FieldValue_ERROR::n4( F, ssexp->name(), comp ) ;
            }
         }
         FIELD_COMPS.append( comp ) ;
         ssexp->destroy() ; ssexp = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;

      // Check field discretization:
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         for( size_t i=0 ; i<FIELDS.size() ; ++i )
         {
            PDE_DiscreteField const* F = FIELDS[i] ;
            if( cFE->nb_local_nodes( F ) == 0 )
            {
               FE_FieldValue_ERROR::n8( F ) ;
            }
         }
      }
   }

   // Parameters:
   if( exp->has_module( "parameters" ) )
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "parameters" ) ;
      sexp->start_module_iterator() ; 
      for( ; sexp->is_valid_module(); sexp->go_next_module() )
      { 
         PEL_ModuleExplorer const* ssexp = sexp->create_subexplorer( 0 ) ;
         std::string param_name = ssexp->string_data( "name" ) ;
         if( !prms->has( param_name ) )
            FE_FieldValue_ERROR::n9( param_name ) ;
         FE_Parameter* p = prms->item( param_name ) ;
         p->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
         PARAMS.push_back( p ) ;
         size_t comp = PEL::bad_index() ;
         if( ssexp->has_entry( "component" ) )
         {
            comp = ssexp->int_data( "component" ) ;
            if( comp >=  p->nb_components() )
            {
               FE_FieldValue_ERROR::n10( p, ssexp->name(), comp ) ;
            }
         }
         PARAM_COMPS.append( comp ) ;
         std::string const& t = ssexp->string_data( "type" ) ;
         size_t type = AT_PTS ;
         if( t == "cell_values" )
         {
            type = CELL_VALUES ;
         }
         else if( t == "at_points" )
         {
            type = AT_PTS ;
         }
         else
         {
            PEL_Error::object()->raise_bad_data_value(
               sexp, "type",
               "   - \"cell_values\"\n"
               "   - \"at_points\"" ) ;
         }
         PARAM_TYPES.append( type ) ;
         
         ssexp->destroy() ; ssexp = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;      
   }

   // Field compositions:
   if( exp->has_module( "field_compositions" ) )
   {
      PDE_SetOfFieldCompositions const* sfcs =
                                          dom->set_of_field_compositions() ;
      PEL_ModuleExplorer* sexp =
                        exp->create_subexplorer( 0, "field_compositions" ) ;
      sexp->start_module_iterator() ; 
      for( ; sexp->is_valid_module(); sexp->go_next_module() )
      { 
         PEL_ModuleExplorer const* ssexp = sexp->create_subexplorer( 0 ) ;
         std::string param_name = ssexp->string_data( "name" ) ;
         if( !sfcs->has( param_name ) )
            FE_FieldValue_ERROR::n11( param_name ) ;
         size_t const field_levels = 0 ;
         FE_Parameter* p = FE_FieldCompositionParameter::create(
             this, param_name, dom, sfcs->item( param_name ), field_levels ) ;
         p->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
         PARAMS.push_back( p ) ;
         size_t comp = PEL::bad_index() ;
         if( ssexp->has_entry( "component" ) )
         {
            comp = ssexp->int_data( "component" ) ;
            if( comp >=  p->nb_components() )
            {
               FE_FieldValue_ERROR::n12( p, ssexp->name(), comp ) ;
            }
         }
         PARAM_COMPS.append( comp ) ;
         PARAM_TYPES.append( AT_PTS ) ;
         
         ssexp->destroy() ; ssexp = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;      
   }
   
   // And expressions:
   if( exp->has_module( "expressions" ) )
   {
      PEL_ContextSimple* ct = PEL_ContextSimple::create( this ) ;
      COORDS = PEL_DoubleVector::create( ct, doubleVector( 0 ) ) ;
      ct->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;   
      TT = PEL_Double::create( ct, 0.0 ) ;
      ct->extend( PEL_Variable::object( "DS_T" ), TT ) ;   

      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "expressions" ) ;
      sexp->start_module_iterator() ; 
      for( ; sexp->is_valid_module(); sexp->go_next_module() )
      { 
         PEL_ModuleExplorer const* ssexp = sexp->create_subexplorer( 0 ) ;
         PEL_DataWithContext* val = ssexp->abstract_data( this, "value", ct ) ;
         if( !val->value_can_be_evaluated() )
         {
            PEL_Error::object()->raise_not_evaluable(
                                   exp, "value", val->undefined_variables() ) ;
         }
         if( val->data_type() != PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
                        exp, "value", PEL_Data::DoubleVector ) ;
         }
         EXPRS.push_back( val ) ;
         size_t comp = PEL::bad_index() ;
         if( ssexp->has_entry( "component" ) )
         {
            comp = ssexp->int_data( "component" ) ;
         }
         EXPR_COMPS.append( comp ) ;
         EXPR_DIMS.append( PEL::bad_index() ) ;
         ssexp->destroy() ; ssexp = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;
   }
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: set_points( PDE_DomainAndFields const* dom,
                            PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: set_points" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( exp != 0 ) ;
   
   PEL_ModuleExplorer const* sexp =
      exp->create_subexplorer( 0, "points_definition" ) ;
   std::string const& type = sexp->string_data( "type" ) ;
   if( type == "list_of_points" )
   {
      doubleArray2D const& pts = sexp->doubleArray2D_data( "points" ) ;
      if( pts.index_bound( 1 ) != NB_SP_DIMS )
      {
         FE_FieldValue_ERROR::n2( sexp, "points", NB_SP_DIMS ) ;
      }
      doubleVector pt( NB_SP_DIMS ) ;
      for( size_t i=0 ; i<pts.index_bound( 0 ) ; ++i )
      {
         for( size_t ic=0 ; ic<NB_SP_DIMS ; ++ic )
         {
            pt(ic) = pts( i, ic ) ;
         }
         POINTS.push_back( GE_Point::create( this, pt ) ) ;
      }   
   }
   else if( type == "cutline" || type == "regular_cutline" )
   {
      doubleVector const& V0 = sexp->doubleVector_data( "first_endpoint" ) ;
      if( V0.size() != NB_SP_DIMS )
      {
         FE_FieldValue_ERROR::n3( sexp, "first_endpoint", NB_SP_DIMS ) ;
      }
      doubleVector const& V1 = sexp->doubleVector_data( "second_endpoint" ) ;
      if( V1.size() != NB_SP_DIMS )
      {
         FE_FieldValue_ERROR::n3( sexp, "second_endpoint", NB_SP_DIMS ) ;
      }
      if( type == "cutline" )
      {
         ABS_PTS = sexp->doubleVector_data( "curvilinear_abscissae" ) ;
      }
      else
      {
         int const n = sexp->int_data( "number_of_points" ) ;
         if( n<2 )
         {
            FE_FieldValue_ERROR::n7( sexp, "number_of_points", n ) ;
         }
         size_t const nb_pts = (size_t) n ;
         ABS_PTS.re_initialize( nb_pts ) ;
         for( size_t i=0 ; i<nb_pts ; ++i )
         {
            ABS_PTS(i) = double(i)/(double)(nb_pts-1) ;
         }
      }
      doubleVector components( NB_SP_DIMS ) ;
      for( size_t i=0 ; i<ABS_PTS.size() ; ++i )
      {
         double abs = ABS_PTS(i) ;
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
         {
            components( d ) = (1.-abs)*V0( d ) + abs*V1( d ) ;
         }
         GE_Point const* pt = GE_Point::create( this, components ) ;
         POINTS.push_back( pt ) ;
      }
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         sexp, "type",
         " - \"list_of_points\"\n"
         " - \"cutline\"\n"
         " - \"regular_cutline\"" ) ;
   }
   REMOVE_PTS = false ;
   if( sexp->has_entry( "ignore_exterior_points" ) )
   {
      REMOVE_PTS = sexp->bool_data( "ignore_exterior_points" ) ;
      sexp->set_default( "ignore_exterior_points", "false" ) ;
   }
   sexp->destroy() ; sexp = 0 ;
     
   PT_CELL.re_initialize( POINTS.size() ) ;
   PT_CELL.set( PEL::bad_index() ) ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: set_post_processing( PDE_DomainAndFields const* dom,
                                     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: set_post_processing" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( exp != 0 ) ;

   PEL_ModuleExplorer const* sexp =
      exp->create_subexplorer( 0, "post_processing" ) ;
   std::string const& type = sexp->string_data( "type" ) ;
   if( type == "one_file" )
   {
      SAVING = FE_FieldValue::one_file ;
      OFILENAME = sexp->string_data( "file_name" )+".txt" ;
      BANNER = true ;
      if( sexp->has_entry( "banner" ) )
      {
         BANNER = sexp->bool_data( "banner" ) ;
      }
      if( sexp->has_entry( "saving_times" ) )
      {
         SAVING_TIMES = sexp->doubleVector_data( "saving_times" ) ;
      }
   }
   else if( type == "separated_files" )
   {
      SAVING = FE_FieldValue::separated_files ;
      OFILENAME = sexp->string_data( "file_basename" ) ;
      BANNER = true ;
      if( sexp->has_entry( "banner" ) )
      {
         BANNER = sexp->bool_data( "banner" ) ;
      }
      if( sexp->has_entry( "saving_times" ) )
      {
         SAVING_TIMES = sexp->doubleVector_data( "saving_times" ) ;
      }
   }
   else if( type == "result_saver" )
   {
      SAVING = FE_FieldValue::result_saver ;
      RSNAME = sexp->string_data( "variable_name" ) ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         sexp, "type",
         " - \"one_file\"\n"
         " - \"separated_files\"\n"
         " - \"result_saver\"" ) ;
   }
   sexp->destroy() ; sexp = 0 ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: set_context( GE_Point const* pt, double time )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: set_context" ) ;
   PEL_CHECK( pt != 0 ) ;

   TT->set( time ) ;
   COORDS->set( pt->coordinate_vector() ) ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: compute_values( FE_TimeIterator const* t_it,
                                doubleArray2D& values )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: compute_values" ) ;
   PEL_CHECK( t_it != 0 ) ;
   PEL_CHECK( values.index_bound(0) == POINTS.size() ) ;
   PEL_CHECK( values.index_bound(1) == NB_VALUES ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   size_t const rank =  com->rank() ;
   size_t const last =  com->nb_ranks()-1 ;
   
   boolVector pt_found( POINTS.size() ) ;

   // Search for "in grid" values
   {
      if( rank != 0 )
      {
         com->receive( rank-1, pt_found ) ;
         com->receive( rank-1, values ) ;
      }
      else
      {
         pt_found.set( false ) ;
         values.set( PEL::bad_double() ) ;
      }

      for( size_t i=0 ; i<POINTS.size() ; ++i )
      {
         if( POINTS[i] == 0 ) pt_found( i ) = true ;
         if( !pt_found( i ) && PT_CELL(i) != PEL::bad_index() )
         {
            compute_values_at_point( i, t_it, pt_found(i), values ) ;
         }
      }

      if( rank < last )
      {
         com->send( rank+1, pt_found ) ;
         com->send( rank+1, values ) ;
         com->receive( last, pt_found ) ;
         com->receive( last, values ) ;
      }
      else
      {
         for( size_t i=0 ; i<last ; ++i )
         {
            com->send( i, pt_found ) ;
            com->send( i, values ) ;
         }
      }
   }

   // Second search for no found points (for parallel)
   {
      if( rank != 0 )
      {
         com->receive( rank-1, pt_found ) ;
         com->receive( rank-1, values ) ;
      }
      
      for( size_t i=0 ; i<POINTS.size() ; ++i )
      {
         if( !pt_found( i ) )
         {
            compute_values_at_point( i, t_it, pt_found(i), values ) ;
         }
      }

      if( rank < last )
      {
         com->send( rank+1, pt_found ) ;
         com->send( rank+1, values ) ;
         com->receive( last, pt_found ) ;
         com->receive( last, values ) ;
      }
      else
      {
         for( size_t i=0 ; i<last ; ++i )
         {
            com->send( i, pt_found ) ;
            com->send( i, values ) ;
         }
      }
   }

   // Verify that all points were found
   if( rank == 0 )
   {
      bool found = true ;
      for( size_t i=0 ; found && i<pt_found.size() ; ++i )
      {
         found = pt_found(i) ;
      }
      if( !found )
      {
         FE_FieldValue_ERROR::n0( POINTS, pt_found ) ;
         if( !REMOVE_PTS )
         {
            PEL_Error::object()->raise_plain(
               "*** FE_FieldValue:\n"
               "    you may use the optional entry\n"
               "       \"ignore_exterior_points\" = true\n"
               "    so that the above points, that were not found\n"
               "    within the grid, will be ignored" ) ;
         }
         else
         {
            for( size_t i=0 ; i<pt_found.size() ; ++i )
            {
               if( !pt_found(i) ) POINTS[i] = (GE_Point const*) 0 ;
            }
         }
      }
   }
   
   if( REDO ) PT_CELL.set( PEL::bad_index() ) ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: compute_values_at_point( size_t pt_index,
                                         FE_TimeIterator const* t_it,
                                         bool& found, doubleArray2D& values )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: compute_values_at_point" ) ;
   PEL_CHECK( pt_index<POINTS.size() ) ;
   PEL_CHECK( values.index_bound(0) == POINTS.size() ) ;
   PEL_CHECK( values.index_bound(1) == NB_VALUES ) ;
   PEL_CHECK( POINTS[pt_index] != 0 ) ;
   PEL_CHECK( t_it != 0 ) ;

   size_t const f_level = 0 ;
   GE_Point const* pt = POINTS[pt_index] ;

   if( PT_CELL( pt_index ) == PEL::bad_index() )
   {
      if( !cFE->is_valid() ) cFE->start() ;
   }
   else
   {
      cFE->go_i_th( PT_CELL( pt_index ) ) ;
   }
   found = PIG->is_in_grid( pt, cFE ) ;
   if( found )
   {
      PT_CELL( pt_index ) = cFE->mesh_id() ;
      cFE->set_calculation_point( pt ) ;
      size_t ival = 0 ;
      for( size_t i=0 ; i<FIELDS.size() ; ++i )
      {
         PDE_DiscreteField const* F = FIELDS[i] ;
         if( FIELD_COMPS(i) == PEL::bad_index() )
         {
            for( size_t ic=0 ; ic<F->nb_components() ; ++ic )
            {
               values( pt_index, ival++) = cFE->value_at_pt( F, f_level, ic ) ;
            }
         }
         else
         {
            size_t const ic = FIELD_COMPS(i) ;
            values( pt_index, ival++) = cFE->value_at_pt( F, f_level, ic ) ;
         }
      }
      for( size_t i=0 ; i<PARAMS.size() ; ++i )
      {
         FE_Parameter const* p = PARAMS[i] ;
         size_t const type = PARAM_TYPES(i) ;
         if( PARAM_COMPS(i) == PEL::bad_index() )
         {
            for( size_t ic=0 ; ic<p->nb_components() ; ++ic )
            {
               if( type == AT_PTS )
               {
                  values( pt_index, ival++) =
                                        p->cell_value_at_pt( t_it, cFE, ic ) ;
               }
               else if( type == CELL_VALUES )
               {
                  values( pt_index, ival++) = p->cell_value( t_it, cFE, ic ) ;
               }
            }
         }
         else
         {
            size_t const ic = PARAM_COMPS(i) ;
            if( type == AT_PTS )
            {
               values( pt_index, ival++) =
                                        p->cell_value_at_pt( t_it, cFE, ic ) ;
            }
            else if( type == CELL_VALUES )
            {
               values( pt_index, ival++) = p->cell_value( t_it, cFE, ic ) ;
            }
         }
      }
      if( !EXPRS.empty() )
      {
         set_context( pt, t_it->time() ) ;
         for( size_t i=0 ; i<EXPRS.size() ; ++i )
         {
            doubleVector const& value = EXPRS[i]->to_double_vector() ;
            if( EXPR_COMPS(i) == PEL::bad_index() )
            {
               for( size_t ic=0 ; ic<EXPR_DIMS(i) ; ++ic )
               {
                  values( pt_index, ival++) = value(ic) ;
               }
            }
            else
            {
               size_t const ic = EXPR_COMPS(i) ;
               values( pt_index, ival++) = value(ic) ;
            }
         }
      }
   }
   else
   {
      PT_CELL( pt_index ) = PEL::bad_index() ;
   }
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: initialize_one_file( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: initialize_one_file" ) ;
   PEL_CHECK( SAVING == one_file ) ;

   std::ofstream file( OFILENAME.c_str(),
                       std::ios::out | std::ios::trunc ) ;
   if( !file )
   {
      FE_FieldValue_ERROR::n6( OFILENAME ) ;
   }

   if( BANNER )
   {
      file << "#" << std::endl ;
      file << "# FE_FieldValue generated file" << std::endl ;
      file << "#" << std::endl ;;
      file << "#    Field values computed at:" << std::endl ;
      for( size_t i=0 ; i<POINTS.size() ; ++i )
      {
         file << "#        pt"<< i << " = " ;
         POINTS.at(i)->print( file, 0 ) ;
         if( ABS_PTS.size() != 0 )
         {
            file << "   -> absolute coordinate : " << ABS_PTS(i) ;
         }
         file << std::endl ;
      }
      file << "#" << std::endl ;
      file << "#        time #" ;
      for( size_t i=0 ; i<POINTS.size() ; ++i )
      {
         file << "### " ;
         for( size_t j=0 ; j<FIELDS.size() ; ++j )
         {
            size_t const nb_cols =
               ( FIELD_COMPS(j) == PEL::bad_index() ?
                    FIELDS[j]->nb_components() : 1 ) ;
            print_name( file, FIELDS[j]->name(), nb_cols ) ;
         }
         for( size_t j=0 ; j<PARAMS.size() ; ++j )
         {
            size_t const nb_cols =
               ( PARAM_COMPS(j) == PEL::bad_index() ?
                    PARAMS[j]->nb_components() : 1 ) ;
            print_name( file, PARAMS[j]->name(), nb_cols ) ;
         }
         for( size_t j=0 ; j<EXPRS.size() ; ++j )
         {
            size_t const nb_cols =
               ( EXPR_COMPS(j) == PEL::bad_index() ?
                    EXPR_DIMS(j) : 1 ) ;
            print_name( file, "expression", nb_cols ) ;
         }
      }
      file << std::endl ;
   }
   file.close() ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: save_values_in_one_file(
                           double time, doubleArray2D const& values ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: save_values_in_one_file" ) ;
   PEL_CHECK( SAVING == one_file ) ;
   PEL_CHECK( values.index_bound(0) == POINTS.size() ) ;
   PEL_CHECK( values.index_bound(1) == NB_VALUES ) ;
   
   std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
   if( !file )
   {
      FE_FieldValue_ERROR::n6( OFILENAME ) ;
   }
   file << "   "
        << std::setprecision( 4 )
        << std::setiosflags( std::ios::scientific )
        << setw( 5 )
        << time
        << std::resetiosflags( std::ios::scientific ) ;
   for( size_t i=0 ; i<values.index_bound(0) ; ++i )
   {
      if( POINTS[i] != 0 )
      {
         file << "    " ;
         for( size_t j=0 ; j<values.index_bound(1) ; ++j )
         {
            print_value( file, values(i,j) ) ;
         }
      }
   }
   file << std::endl ;
   file.close() ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: save_values_in_separated_files(
           size_t i_cycle, double time, doubleArray2D const& values ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldValue:: save_values_in_separated_files" ) ;
   PEL_CHECK( SAVING == separated_files ) ;
   PEL_CHECK( values.index_bound(0) == POINTS.size() ) ;
   PEL_CHECK( values.index_bound(1) == NB_VALUES ) ;

   PEL_ASSERT( i_cycle<99999 ) ;

   std::ostringstream tmp ;
   tmp << i_cycle ;
   std::string nb_string = tmp.str() ;
   std::string filename = OFILENAME+"_00000" ;
   filename.replace( filename.length()-nb_string.length(), 
                     nb_string.length(), nb_string ) ;
   filename += ".txt";
   
   std::ofstream file( filename.c_str(), std::ios::out | std::ios::trunc ) ;
   if( !file )
   {
      FE_FieldValue_ERROR::n6( filename ) ;
   }

   if( BANNER )
   {
      file << "#" << std::endl ;
      file << "# FE_FieldValue generated file" << std::endl ;
      file << "#" << std::endl ;;
      file << "#    Field values computed at time: " << time << std::endl ;
      file << "#" << std::endl ;

      file << "# " ;
      print_name( file, "coordinates", NB_SP_DIMS ) ;
      if( ABS_PTS.size() != 0 )
      {
         print_name( file, "curve", 1 )  ;
      }
      for( size_t j=0 ; j<FIELDS.size() ; ++j )
      {
         size_t const nb_cols =
            ( FIELD_COMPS(j) == PEL::bad_index() ?
                 FIELDS[j]->nb_components() : 1 ) ;
         print_name( file, FIELDS[j]->name(), nb_cols ) ;
      }
      for( size_t j=0 ; j<PARAMS.size() ; ++j )
      {
         size_t const nb_cols =
            ( PARAM_COMPS(j) == PEL::bad_index() ?
                    PARAMS[j]->nb_components() : 1 ) ;
         print_name( file, PARAMS[j]->name(), nb_cols ) ;
      }
      for( size_t j=0 ; j<EXPRS.size() ; ++j )
      {
         std::string s = "expression" ;
         size_t const nb_cols =
            ( EXPR_COMPS(j) == PEL::bad_index() ?
                 EXPR_DIMS(j) : 1 ) ;
         print_name( file, "expression", nb_cols ) ;
      }
      file << std::endl ;
   }
   
   for( size_t i=0 ; i<values.index_bound(0) ; ++i )
   {
      GE_Point const* pt = POINTS[i] ;
      if( pt != 0 )
      {
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
         {
            print_value( file, pt->coordinate( d ) ) ;
         }
         if( ABS_PTS.size() != 0 )
         {
            print_value( file, ABS_PTS(i) ) ;
         }
         for( size_t j=0 ; j<values.index_bound(1) ; ++j )
         {
            print_value( file, values(i,j) ) ;
         }
         file << std::endl ;
      }
   }
   file << std::endl ;
   file.close() ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: print_name( std::ostream& os, std::string name,
                            size_t nb_col ) const
//-------------------------------------------------------------------------
{
   std::string s = name ;
   if( s.size()>11*nb_col )
   {
      s.erase( 11*nb_col, s.size() ) ;
   }
   os << std::setw( 13*nb_col-2 ) << s << " #" ;
}

//-------------------------------------------------------------------------
void
FE_FieldValue:: print_value( std::ostream& os, double value ) const
//-------------------------------------------------------------------------
{
   std::ios_base::fmtflags original_flags = os.flags() ;
   os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision( 5 ) ;
   os << " " << setw( 12 ) << value ;
   os << std::setprecision(p) ;
   os.flags( original_flags ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n0( std::vector<GE_Point const*> const& pts,
                          boolVector const& pt_found  )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue: the following points were not found" << endl ;
   msg << "                   within the grid" 
       << std::endl ;
   for( size_t i=0 ; i<pts.size() ; ++i )
   {
      if( !pt_found(i) )
      {
         msg << "       - " ;
         pts[i]->print( msg, 0 ) ;
         msg << std::endl ;
      }
   }
   msg << std::endl ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n1( std::string const& field_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the field : \"" <<  field_name << "\" is unknown" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n2( PEL_ModuleExplorer const* exp,
                          std::string const& entry_name, size_t nb_dims )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "    The second dimension should be of size " << nb_dims 
       << "    (i.e. the number of space dimensions)" ;
   PEL_Error::object()->raise_bad_data_value( exp, entry_name, msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n3( PEL_ModuleExplorer const* exp,
                          std::string const& entry_name, size_t nb_dims )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "    The size should be " << nb_dims 
       << "    (i.e. the number of space dimensions)" ;
   PEL_Error::object()->raise_bad_data_value( exp, entry_name, msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n4( PDE_DiscreteField const* ff,
                          std::string const& module_name,
                          size_t comp )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the entry of keyword : \"component\"" << std::endl ;
   msg << "    in MODULE : \"" << module_name << "\"" << std::endl ; 
   msg << "    has the value : " << comp << std::endl ;
   msg << "    whereas the number of components of field \"" << ff->name()
       << "\" is : " << ff->nb_components() ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n5( size_t the_size, size_t the_comp )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    a computation is requested for the component : "
       << the_comp << std::endl ;
   msg << "    of a DoubleVector whose size is only : " << the_size ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n6( std::string const& file_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    Unable to open file \""
       << file_name << "\" for writing" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n7( PEL_ModuleExplorer const* exp,
                          std::string const& entry_name, int nb_points )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "    A value greater than 1 is expected."<< std::endl ;
   if( nb_points == 1 )
   {
      msg << "    Use the \"points_definition\" option for "
          << "single points. " << std::endl ;
   }
   PEL_Error::object()->raise_bad_data_value( exp, entry_name, msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n8( PDE_DiscreteField const* ff )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the field of name \"" << ff->name() << "\"" << std::endl ;
   msg << "    is not defined on all the cells of the grid" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n9( std::string const& param_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the parameter : \"" <<  param_name << "\" is unknown" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n10( FE_Parameter const* p,
                           std::string const& module_name,
                           size_t comp )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the entry of keyword : \"component\"" << std::endl ;
   msg << "    in MODULE : \"" << module_name << "\"" << std::endl ; 
   msg << "    has the value : " << comp << std::endl ;
   msg << "    whereas the number of components of parameter \"" << p->name()
       << "\" is : " << p->nb_components() ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n11( std::string const& compo_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the field composition : \"" <<  compo_name << "\" is unknown" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_FieldValue_ERROR:: n12( FE_Parameter const* p,
                           std::string const& module_name,
                           size_t comp )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_FieldValue error:" << std::endl ;
   msg << "    the entry of keyword : \"component\"" << std::endl ;
   msg << "    in MODULE : \"" << module_name << "\"" << std::endl ; 
   msg << "    has the value : " << comp << std::endl ;
   msg << "    whereas the number of components of field composition \"" << p->name()
       << "\" is : " << p->nb_components() ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
