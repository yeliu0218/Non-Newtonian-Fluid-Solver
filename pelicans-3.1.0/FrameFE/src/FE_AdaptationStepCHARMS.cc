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

#include <FE_AdaptationStepCHARMS.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SegmentPolyhedron_INT.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE_AdaptationIndicator.hh>
#include <FE_TimeIterator.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

FE_AdaptationStepCHARMS const* 
FE_AdaptationStepCHARMS::PROTOTYPE = new FE_AdaptationStepCHARMS() ;

struct FE_AdaptationStepCHARMS_ERROR
{
   static void n0( std::string const& module_name ) ;
   static void n1( void ) ;
} ;

//----------------------------------------------------------------------
FE_AdaptationStepCHARMS:: FE_AdaptationStepCHARMS( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_AdaptationStepCHARMS" )
{
}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMS*
FE_AdaptationStepCHARMS:: create_replica( PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_AdaptationStepCHARMS* result = 
                      new FE_AdaptationStepCHARMS( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMS:: FE_AdaptationStepCHARMS( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , DOM( dom )
   , DFS( dom->set_of_discrete_fields() )
   , EXP( 0 )
   , DA( 0 )
   , INDICS( PEL_List::create( this ) )
   , INDICS_IT( 0 )
   , KEEP_ADAPTING( false )
   , SAVER( dom->result_saver() )
   , BUILD_INIT( exp->has_entry( "build_initial_refinement" ) ?
                 exp->bool_data( "build_initial_refinement" ) :
                 true ) 
   , NB_ITER_MAX_0( exp->int_data( "nb_iterations_max_before_time_stepping" ) )
   , NB_ITER_MAX_1( exp->int_data( "nb_iterations_max_during_time_stepping" ) )
   , STOP_CV( exp->has_entry( "stop_after_nb_iterations_max" ) ?
              exp->bool_data( "stop_after_nb_iterations_max" ) :
              false )
   , CHECK_FACES( exp->has_entry( "check_faces_consistency" ) ? 
                  exp->bool_data( "check_faces_consistency" ) :
                  false )
   , CHECK_BD_CELLS( exp->has_entry( "check_covering_of_cells_boundary" ) ? 
                     exp->bool_data( "check_covering_of_cells_boundary" ) :
                     false )
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: FE_AdaptationStepCHARMS" ) ;
   
   if( exp->has_module( "list_of_PDE_DiscreteField" ) )
   {
      EXP = exp->create_subexplorer( this, "list_of_PDE_DiscreteField" ) ;
      EXP->start_module_iterator() ;
      for( ; EXP->is_valid_module() ; EXP->go_next_module() )
      {
         PEL_ModuleExplorer* se = EXP->create_subexplorer( 0 ) ;
         PDE_DiscreteField const* sf = DFS->item( 
                                            se->string_data( "current" ) ) ;
         PDE_DiscreteField const* tf = DFS->item( 
                                            se->string_data( "explicit" ) ) ;

         if( ( sf->nb_nodes() != tf->nb_nodes() ) ||
             ( sf->nb_components() != tf->nb_components() ) ||
             ( sf->storage_depth() != tf->storage_depth() ) )
         {
            FE_AdaptationStepCHARMS_ERROR::n0( se->absolute_path_name() ) ;
         }
         se->destroy() ; se = 0 ;
      }
   }
   
   if( NB_ITER_MAX_0 != 0 || NB_ITER_MAX_1 != 0 )
   {
      DA = dom->adapter_CHARMS() ;
      if( exp->has_module( "FE_AdaptationIndicator" ) )
      {
         if( exp->has_module( "list_of_FE_AdaptationIndicator" ) )
                                       FE_AdaptationStepCHARMS_ERROR::n1() ;
         PEL_ModuleExplorer* se = 
                   exp->create_subexplorer( 0, "FE_AdaptationIndicator" ) ;
         size_t verb = verbose_level() ;
         FE_AdaptationIndicator* indic = 
                      FE_AdaptationIndicator::make( DA, DOM, prms, se , verb ) ;
         se->destroy() ; se = 0 ;
         INDICS->append( indic ) ;
         DA->append_indicator( indic ) ;
      }
      else if( exp->has_module( "list_of_FE_AdaptationIndicator" ) )
      {
         if( exp->has_module( "FE_AdaptationIndicator" ) )
                               FE_AdaptationStepCHARMS_ERROR::n1() ;
         PEL_ModuleExplorer* se =
             exp->create_subexplorer( 0, "list_of_FE_AdaptationIndicator" ) ;
         size_t verb = verbose_level() ;
         se->start_module_iterator() ;
         for( ; se->is_valid_module() ; se->go_next_module() )
         {
            PEL_ModuleExplorer* sse = se->create_subexplorer( se ) ;
            FE_AdaptationIndicator* indic =
               FE_AdaptationIndicator::make( DA, DOM, prms, sse, verb ) ;
            INDICS->append( indic ) ;
            DA->append_indicator( indic ) ;
         }
         se->destroy() ; se = 0 ;
      }
   }
   if( INDICS != 0 ) INDICS_IT = INDICS->create_iterator( INDICS ) ;      

}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMS:: ~FE_AdaptationStepCHARMS( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_before_time_stepping( 
                                         FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "FE_AdaptationStepCHARMS:: do_before_time_stepping" ) ;

   // do_saving( t_it ) ; // for an initial save before refining

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   if( DA != 0 )
   {      
      if( BUILD_INIT )
      {
         DA->reset() ;
         
         size_t nb_iter = 0 ;
         do
         {
            nb_iter++ ;
            if( nb_iter > NB_ITER_MAX_0 )
            {
               if( STOP_CV ) 
               {
                  PEL_Error::object()->raise_plain(
                  "*** FE_AdaptationStepCHARMS :  "
                  "    convergence failure before time stepping" ) ;
               }
               else
               {
                  PEL::out() << indent() << endl
                       << "      max nb iterations reached before convergence"
                       << endl << endl ;
                  break ;
               }
            }
            PEL::out() << indent() << "   adaptation... " ;
            give_time_to_indicators( t_it ) ;
            DA->adapt() ;
            KEEP_ADAPTING = com->boolean_or( DA->something_changed() ) ;
            if( !KEEP_ADAPTING )
            {
               PEL::out() << "nothing done" << endl ;
            }
            else
            {
               PEL::out() << endl ;
               DOM->apply_requests_of_DOFs_values_modules( true ) ;
               DA->print_statistics( PEL::out(), indent().size()+6 ) ;

               // do_saving( t_it ) ; // for saving after each refinement
            }

         } while( KEEP_ADAPTING ) ;
      }

      if( EXP != 0 )
      {
         EXP->start_module_iterator() ;
         for( ; EXP->is_valid_module() ; EXP->go_next_module() )
         {
            PEL_ModuleExplorer* se = EXP->create_subexplorer( 0 ) ;
            PDE_DiscreteField* tf = DFS->item( 
                                         se->string_data( "explicit" ) ) ;
         
            DA->add_excluded_field( tf ) ;
         
            se->destroy() ; se = 0 ;
         }
      }
   }

   // DOM->print_grid( cout, 3 ) ;
   
   if( CHECK_FACES    )  check_faces_consistency() ;
   if( CHECK_BD_CELLS )  check_covering_of_cells_boundary() ;
   
   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_before_inner_iterations_stage( 
                                          FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( 
               "FE_AdaptationStepCHARMS:: do_before_inner_iterations_stage" ) ;

   if( DA != 0 ) DA->reset() ;
   KEEP_ADAPTING = false ;

   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_AdaptationStepCHARMS:: do_one_inner_iteration" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   if( DA != 0 )
   {
      size_t nb_iter = 0 ;
      do
      {
         nb_iter++ ;
         if( nb_iter > NB_ITER_MAX_1 )
         {
            if( STOP_CV ) 
            {
               PEL_Error::object()->raise_plain(
                  "*** FE_AdaptationStepCHARMS :  "
                  "    convergence failure at initialization" ) ;
            }
            else
            {
               PEL::out() << indent() 
                    << "   max nb iterations reached before convergence"
                    << endl ;
               break ;
            }
         }
         PEL::out() << indent() << "   adaptation... " ;
         give_time_to_indicators( t_it ) ;
         DA->adapt() ;
         KEEP_ADAPTING = com->boolean_or( DA->something_changed() ) ;
         if( !KEEP_ADAPTING )
         {
            PEL::out() << "nothing done" << endl ;
         }
         else
         {
            PEL::out() << endl ;
            DA->print_statistics( PEL::out(), indent().size()+6 ) ;
                     
         // do_saving( t_it ) ; // for saving after each refinement
         }
      } while( KEEP_ADAPTING ) ;
   }
   else
   {
      KEEP_ADAPTING = false ;
   }

   // DOM->print_grid( cout, 3 ) ;
   
   if( CHECK_FACES    )  check_faces_consistency() ;
   if( CHECK_BD_CELLS )  check_covering_of_cells_boundary() ;

   stop_total_timer() ;
}

//------------------------------------------------------------------------
bool
FE_AdaptationStepCHARMS:: inner_iterations_are_completed( 
                                           FE_TimeIterator const* t_it ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: inner_iterations_are_completed" ) ;
   PEL_CHECK_PRE( inner_iterations_are_completed_PRE( t_it ) ) ;

   start_total_timer( 
               "FE_AdaptationStepCHARMS:: inner_iterations_are_completed" ) ;

   bool result = !KEEP_ADAPTING ;

   stop_total_timer() ;

   return( result ) ;
}

//------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_additional_savings( FE_TimeIterator const* t_it,
                                                 PDE_ResultSaver* rs )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_additional_savings" ) ;
   PEL_CHECK_PRE( do_additional_savings_PRE( t_it, rs ) ) ;

   start_total_timer( "FE_AdaptationStepCHARMS:: do_additional_savings" ) ;

   if( DA != 0 )
   {
      if( t_it->is_started() ) 
      {
         if( verbose_level() >= 2 )
         {
            PEL::out() << indent() << "   save grid..." << endl ;
         }
         rs->save_grid() ;
      }
   }

   stop_total_timer() ;
}

//-------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_after_inner_iterations_stage(
                                         FE_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_after_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( 
               "FE_AdaptationStepCHARMS:: do_after_inner_iterations_stage" ) ;

   if( EXP != 0 )
   {
      EXP->start_module_iterator() ;
      for( ; EXP->is_valid_module() ; EXP->go_next_module() )
      {
         PEL_ModuleExplorer* se = EXP->create_subexplorer( 0 ) ;
         PDE_DiscreteField* sf = DFS->item( se->string_data( "current" ) ) ;
         PDE_DiscreteField* tf = DFS->item( se->string_data( "explicit" ) ) ;
         if( verbose_level() >= 2 )
         {
            PEL::out() << indent() << "   copy " 
                       << "\"" << sf->name() << "\" -> " 
                       << "\"" << tf->name() << "\"..." << endl ;
         }
         tf->set( sf ) ;
         se->destroy() ; se = 0 ;
      }
   }

   if( DA != 0 )
   {
      if( verbose_level() >= 2 )
         PEL::out() << indent() << "   unsplit of the meshes " << endl ;
      DA->unsplit_meshes() ;
      
      if( CHECK_FACES    )  check_faces_consistency() ;
      if( CHECK_BD_CELLS )  check_covering_of_cells_boundary() ;
      
      if( verbose_level() >= 2 )
         DA->print_statistics( PEL::out(), indent().length()+9 ) ;
   }

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: print" ) ;

   FE_OneStepIteration:: print( os, indent_width ) ;

   std::string const s( indent_width+3, ' ' ) ;

   if( DA==0 )
   {
      os << s << "no adaptation" << endl ;
   }

   if( EXP != 0 ) 
   {
      EXP->start_module_iterator() ;
      for( ; EXP->is_valid_module() ; EXP->go_next_module() )
      {
         PEL_ModuleExplorer* se = EXP->create_subexplorer( 0 ) ;     
         PDE_DiscreteField* sf = DFS->item( se->string_data( "explicit" ) ) ;
         PDE_DiscreteField* tf = DFS->item( se->string_data( "current" ) ) ;
         os << s << "copy "
            << "\"" << sf->name() << "\" -> "
            << "\"" << tf->name() << "\"..." << std::endl ;
         se->destroy() ; se = 0 ;
      }
   }
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: give_time_to_indicators( 
                                       FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: give_time_to_indicators" ) ;

   if( INDICS_IT != 0 )
   {
      INDICS_IT->start() ;
      for( ; INDICS_IT->is_valid() ; INDICS_IT->go_next() )
      {
         FE_AdaptationIndicator* indic =
             static_cast< FE_AdaptationIndicator* >( INDICS_IT->item() ) ;
         indic->set_time_iterator( t_it ) ;
      }
   }
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: do_saving( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: do_saving" ) ;

   SAVER->start_cycle() ;
   PEL::out() << endl << "   +++ SAVE FOR POSTPROCESSING" 
              << " *** CYCLE = " << SAVER->cycle_number()
              << " *** TIME = "  << t_it->time() 
              << " ++++++" << endl << endl ;
   SAVER->save_grid() ;
   SAVER->save_fields( 0 ) ;
   SAVER->terminate_cycle() ;
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: check_faces_consistency( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: check_faces_consistency" ) ;
   
   PDE_CursorFEside* s_fe = DOM->create_CursorFEside( 0 ) ;
   
   for( s_fe->start() ; s_fe->is_valid() ; s_fe->go_next() )
   {
      size_t s_level = s_fe->refinement_level() ;
      
      PDE_LocalFEcell const* c_fe_0 = s_fe->adjacent_localFEcell( 0 ) ;
      size_t c0_level = c_fe_0->refinement_level() ;
      
      PDE_LocalFEcell const* c_fe_1 = s_fe->adjacent_localFEcell( 1 ) ;
      size_t c1_level = c_fe_1->refinement_level() ;
      
      bool ok = ( ( (s_level == c0_level) && (s_level >= c1_level) ) ||
                  ( (s_level == c1_level) && (s_level >= c0_level) ) ) ;
      if( !ok )
      {
         PEL::out() << "********" << endl ;
         s_fe->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         c_fe_0->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         c_fe_1->print_current_mesh( PEL::out(), 3 ) ;
         PEL_Error::object()->raise_plain( "error in refinement levels" ) ;
      }
   }
   s_fe->destroy() ; s_fe=0 ;
   
   PDE_LocalFEcell*  c_fe = DOM->create_LocalFEcell( 0 ) ;
   PDE_LocalFEbound* b_fe = DOM->create_LocalFEbound( 0 ) ;
   for( b_fe->start() ; b_fe->is_valid() ; b_fe->go_next() )
   {
      c_fe->go_i_th( b_fe->adjacent_cell_id() ) ;      
      bool ok = ( b_fe->refinement_level() == c_fe->refinement_level() ) ;
      if( !ok )
      {
         PEL::out() << "********" << endl ;
         b_fe->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         c_fe->print_current_mesh( PEL::out(), 3 ) ;
         PEL_Error::object()->raise_plain( "error in refinement levels" ) ;
      }
   }
   c_fe->destroy() ;
   b_fe->destroy() ;
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepCHARMS:: check_covering_of_cells_boundary( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMS:: check_covering_of_cells_boundary" ) ;
   
   if( DOM->nb_space_dimensions() != 2 )
   {
      std::string mesg = "*** FE_AdaptationStepCHARMS: \n" ;
      mesg+= "   internal check of faces only provided in 2D" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
   
   size_t dim = 2 ;
   size_t nb_points = 50 ;
   
   GE_SegmentPolyhedron_INT* sp = 
    GE_SegmentPolyhedron_INT::make( 0, "GE_SegmentPolyhedron1D_INT", dim, 0 ) ;
   
   GE_Point* pt = GE_Point::create( 0, dim ) ;
   PDE_LocalFEcell*  c_fe = DOM->create_LocalFEcell( 0 ) ;
   c_fe->include_color( GE_Color::halo_color() ) ;
   PDE_LocalFEbound* b_fe = DOM->create_LocalFEbound( 0 ) ;
   b_fe->include_color( GE_Color::halo_color() ) ;
   PDE_CursorFEside* s_fe = DOM->create_CursorFEside( 0 ) ;
   s_fe->include_color( GE_Color::halo_color() ) ;
   for( c_fe->start() ; c_fe->is_valid() ; c_fe->go_next() )
   {
      size_t_vector const& sid = c_fe->adjacent_side_ids() ;
      size_t_vector const& bid = c_fe->adjacent_bound_ids() ;

      GE_Mpolyhedron const* poly = c_fe->polyhedron() ;
      GE_Point const* center = poly->center() ;
      double x_c = center->coordinate( 0 ) ;
      double y_c = center->coordinate( 1 ) ;
      double rr = poly->inter_vertices_maximum_distance() ;
      for( size_t i=0 ; i<nb_points ; ++i )
      {
         double theta = 2.0*PEL::pi()*((double)i)/((double)nb_points) ;
         pt->set_coordinate( 0, x_c + rr*PEL::cos( theta ) ) ;
         pt->set_coordinate( 1, y_c + rr*PEL::sin( theta ) ) ;
         size_t nb_inter = 0 ;
         for( size_t j=0 ; j<sid.size() ; ++j )
         {
            s_fe->go_i_th( sid( j ) ) ;
            sp->check_intersection( center, pt, s_fe->polyhedron() ) ;
            if( sp->one_single_intersection() )
            {
               nb_inter++ ;
            }
         }
         for( size_t j=0 ; j<bid.size() ; ++j )
         {
            b_fe->go_i_th( bid( j ) ) ;
            sp->check_intersection( center, pt, b_fe->polyhedron() ) ;
            if( sp->one_single_intersection() )
            {
               nb_inter++ ;
            }
         }
         if( nb_inter!=1 && nb_inter!=2 )
         {
            PEL::out() << "******** nb_inter = " << nb_inter << endl ;
            c_fe->print( PEL::out(), 3 ) ;
            // poly->print( PEL::out(), 3 ) ;
            PEL::out() << "segment : " ;
            center->print( PEL::out(), 0 ) ; PEL::out() << " -> " ;
            pt->print( PEL::out(), 0 ) ; PEL::out() << endl ;
            for( size_t j=0 ; j<sid.size() ; ++j )
            {
               s_fe->go_i_th( sid( j ) ) ;
               s_fe->polyhedron()->print( PEL::out(), 6 ) ;
            }
            for( size_t j=0 ; j<bid.size() ; ++j )
            {
               b_fe->go_i_th( bid( j ) ) ;
               b_fe->polyhedron()->print( PEL::out(), 6 ) ;
            }
            PEL_Error::object()->raise_plain( "error" ) ;
         }
      }
      
   }
   c_fe->destroy() ;
   s_fe->destroy() ;
   b_fe->destroy() ;
   pt->destroy() ;
   sp->destroy() ;
}

//internal--------------------------------------------------------------
void
FE_AdaptationStepCHARMS_ERROR:: n0( std::string const& module_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_AdaptationStepCHARMS error:" << endl ;
   msg << "    within module \"" << module_name << "\"" << endl ;
   msg << "    the entries of keyword \"explicit\" and \"current\"" << endl ;
   msg << "    should represent fields with exactly the same characteristics" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_AdaptationStepCHARMS_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** FE_AdaptationStepCHARMS error:" << endl ;
   msg << "    at most one of the following modules is accepted:" << endl ;
   msg << "       MODULE FE_AdaptationIndicator" << endl ;
   msg << "       MODULE list_of_FE_AdaptationIndicator" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
