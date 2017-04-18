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

#include <FE_ParameterSaver.hh>

#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfFieldCompositions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <doubleArray2D.hh>

#include <iostream>

using std::string ;

FE_ParameterSaver const* 
FE_ParameterSaver::PROTOTYPE = new FE_ParameterSaver() ;

size_t const CELL_VALUES = 0 ;
size_t const AT_CELL_CENTERS = 1 ;
size_t const AT_VERTICES = 2 ;

//---------------------------------------------------------------------------
FE_ParameterSaver:: FE_ParameterSaver( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_ParameterSaver")
   , PARAMS( 0 )
   , NAMES( 0 )
   , TYPES( 0 )
   , cFE( 0 )
   , EPS_DBL( PEL::bad_double() )
   , MIN_DBL( PEL::bad_double() )
{
}

//---------------------------------------------------------------------------
FE_OneStepIteration*
FE_ParameterSaver:: create_replica( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ParameterSaver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_ParameterSaver* result =
                       new FE_ParameterSaver( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_ParameterSaver:: FE_ParameterSaver( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , PARAMS( PEL_Vector::create( this, 0 ) )
   , NAMES( 0 )
   , TYPES( 0 )
   , cFE( dom->create_LocalFEcell( this ) )
   , EPS_DBL( exp->has_entry( "dbl_epsilon" ) ?
                 exp->double_data( "dbl_epsilon" ) : 1.E-4 )
   , MIN_DBL( exp->has_entry( "dbl_minimum" ) ?
                 exp->double_data( "dbl_minimum" ) : 1.E-8 )
{
   PEL_LABEL( "FE_ParameterSaver:: FE_ParameterSaver" ) ;

   cFE->include_color( GE_Color::halo_color() ) ;

   PEL_ModuleExplorer* pexp = exp->create_subexplorer( 0, "parameters" ) ;
   pexp->start_module_iterator() ;
   for( ; pexp->is_valid_module() ; pexp->go_next_module() )
   {
      PEL_ModuleExplorer* ee = pexp->create_subexplorer( 0 ) ;
      
      FE_Parameter* param =
                      prms->item( ee->string_data( "parameter_name" ) ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      PARAMS->append( param ) ;

      NAMES.append( ee->string_data( "entry_name" ) ) ;

      size_t type = CELL_VALUES ;
      std::string const& t = ee->string_data( "type" ) ;
      if( t == "cell_values" )
      {
         type = CELL_VALUES ;
      }
      else if( t == "at_cell_centers" )
      {
         type = AT_CELL_CENTERS ;
      }
      else if( t == "at_vertices" )
      {
         type = AT_VERTICES ;
      }
     
      else
      {
         PEL_Error::object()->raise_bad_data_value(
            ee, "type",
            "   - \"cell_values\"\n"
            "   - \"at_cell_centers\"\n"
            "   - \"at_vertices\"" ) ;
      }
      TYPES.append( type ) ;
      
      ee->destroy() ; ee = 0 ;
   }
   pexp->destroy() ; pexp = 0 ;
}

//---------------------------------------------------------------------------
FE_ParameterSaver:: ~FE_ParameterSaver( void )
//---------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//---------------------------------------------------------------------------
void
FE_ParameterSaver:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ParameterSaver:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
FE_ParameterSaver:: save_other_than_time_and_fields(
                           FE_TimeIterator const* t_it, PDE_ResultSaver* rs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ParameterSaver:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   // Field savings should have been desactivated...
   if( rs->grid_is_saved() )
   {
      doubleArray2D values( 0, 0 ) ;
      doubleVector default_values( 0 ) ;
   
      size_t const nb_cells = cFE->nb_meshes() ;
      for( size_t i=0 ; i<PARAMS->index_limit() ; ++i )
      {
         FE_Parameter* param = static_cast<FE_Parameter*>( PARAMS->at(i) ) ;
         size_t nbcs = param->nb_components() ;
         size_t type = TYPES(i) ;
         PDE_ResultSaver::SavingLocation location = 
               ( type == AT_VERTICES ? PDE_ResultSaver::AtVertices :
                 PDE_ResultSaver::AtCellCenters ) ;
         string mesg =
                "*** FE_ParameterSaver : error saving \""+param->name()+"\"" ;
         rs->prepare_for_field_saving( location, NAMES(i), nbcs, 
                                       values, default_values ) ;
         size_t im = 0 ;
         for( cFE->start() ; cFE->is_valid() ; cFE->go_next(), ++im )
         {
            PEL_CHECK( im < nb_cells ) ;
            if( type == CELL_VALUES )
            {
               for( size_t ic=0 ; ic<nbcs ; ++ic )
               {
                  values( ic, im ) = param->cell_value( t_it, cFE, ic ) ;
               }
            }
            else if( type == AT_CELL_CENTERS )
            {
               cFE->set_calculation_point( cFE->polyhedron()->center() ) ;
               for( size_t ic=0 ; ic<nbcs ; ++ic )
               {
                  values( ic, im ) = param->cell_value_at_pt( t_it, cFE, ic ) ;
               }
            }
            else
            {
               GE_Mpolyhedron const* poly = cFE->polyhedron() ;
               for( size_t v=0 ; v<poly->nb_vertices() ; ++v )
               {
                  GE_Point const* pt = poly->vertex( v ) ;
                  size_t const iv = rs->vertex_index( pt ) ;
                  cFE->set_calculation_point( pt ) ;
                  for( size_t ic=0 ; ic<nbcs ; ++ic )
                  {
                     double val = param->cell_value_at_pt( t_it, cFE, ic ) ;
                  
                     // the value at a vertex obtained from different cells
                     // should be the same
                     PDE_ResultSaver::check_value_consistency_at_vertex(
                         mesg, pt, values( ic, iv ), val, EPS_DBL, MIN_DBL ) ;
                     
                     if( ( default_values( ic ) == 
                                   PDE_ResultSaver::undefined_value() ) ||
                         ( val > default_values( ic ) ) )
                     {
                        default_values( ic ) = val ;
                     }
                     
                     values( ic, iv ) = val ;
                  }
               }
            }
         }
         rs->save_field( values, default_values ) ;
      }
   }
}

//---------------------------------------------------------------------------
void
FE_ParameterSaver:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ParameterSaver:::: print" ) ;   

   FE_OneStepIteration::print( os, indent_width ) ;
   std::string const s( indent_width+3, ' ' ) ;
   for( size_t i=0 ; i<PARAMS->index_limit() ; ++i )
   {
      FE_Parameter const* param =
                     static_cast<FE_Parameter const*>( PARAMS->at(i) ) ;
      os << s << "saving \"" << param->name()
         <<"\" (postprocessing entry: \"" << NAMES(i) << "\")" << std::endl ;
   }
}
