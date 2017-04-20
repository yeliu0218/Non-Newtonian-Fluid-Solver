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

#include <FE_RefinementAppli.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>
#include <PEL_Variable.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_SegmentPolyhedron_INT.hh>

#include <PDE.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>

using std::endl ;
using std::cout ;
using std::map ;
using std::set ;
using std::string ;
using std::vector ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE_RefinementAppli const*
FE_RefinementAppli::PROTOTYPE = new FE_RefinementAppli() ;

//----------------------------------------------------------------------
FE_RefinementAppli:: FE_RefinementAppli( void )
//----------------------------------------------------------------------
   : PEL_Application( "FE_RefinementAppli" ) 
{
}

//----------------------------------------------------------------------
FE_RefinementAppli*
FE_RefinementAppli:: create_replica( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   FE_RefinementAppli* result = new FE_RefinementAppli( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
FE_RefinementAppli:: FE_RefinementAppli( PEL_Object* a_owner,
                     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , DOM( 0 )
   , UU( 0 )
   , L_UU( PEL::bad_index() )
   , UU_EXP( 0 )
   , L_UU_EXP( PEL::bad_index() )
   , TIME_IT( 0 )
   , PRMS( 0 )
   , BCs( 0 )
   , SAVER( 0 )
   , NB_ITER_MAX( 100 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( 0 )
   , cFE( 0 )
   , bFE( 0 )
   , ALPHA( 0.0 )
   , LAMBDA( 0.0 )
   , FORCE_TERM( 0 )
   , AA( 0 )
   , PRM_CHECK( 0 )
   , UU_CHECK( 0 )
   , EPS_CHECK( PEL::bad_double() )
   , MIN_CHECK( PEL::bad_double() )
{
   PEL_LABEL( "FE_RefinementAppli:: FE_RefinementAppli" ) ;

   PEL_ModuleExplorer* ee=exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;

   DOM = PDE_DomainAndFields::create( this, ee, PEL_Exec::communicator() ) ;
      
   UU = DOM->set_of_discrete_fields()->item( "uu" ) ;
   UU_EXP = DOM->set_of_discrete_fields()->item( "uu_explicit" ) ;

   BCs = DOM->set_of_boundary_conditions() ;
   SAVER = DOM->result_saver() ;

   cFE = DOM->create_LocalFEcell( this ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( UU_EXP, PDE_LocalFE::N ) ;

   bFE = DOM->create_LocalFEbound( this ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;

   FE::set_geometry( FE::cartesian ) ;

   ee->destroy() ;
   ee = exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   PRMS = FE_SetOfParameters::create( this, DOM, ee ) ;

   ee->destroy() ;
   ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   TIME_IT = FE_TimeIterator::create( this, ee ) ;

   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU, true ) ;
   NMB = PDE_SystemNumbering::create( this, uu_link ) ;
   
   ee->destroy() ;
   ee = exp->create_subexplorer( 0, "FE_RefinementAppli" ) ;
   read_one_step_iteration_data( DOM, PRMS, ee ) ;
   ee->destroy() ;
}

//----------------------------------------------------------------------
FE_RefinementAppli:: ~FE_RefinementAppli( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_RefinementAppli:: read_one_step_iteration_data( PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: read_one_step_iteration_data" ) ;

   L_UU = exp->int_data( "level_of_current" ) ;
   L_UU_EXP = exp->int_data( "level_of_explicit" ) ;

   FORCE_TERM = PRMS->item( exp->string_data( "param_force" ) ) ;
   FORCE_TERM->transfer_cell_calculation_requirements( cFE, 
                                                       FE_Parameter::Val ) ;
   ALPHA = exp->double_data( "alpha" ) ;
   LAMBDA = exp->double_data( "lambda" ) ;

   if( exp->has_entry( "param_advective_field" ) )
   {
      AA = prms->item( exp->string_data( "param_advective_field" ) ) ;
      AA->transfer_cell_calculation_requirements( cFE,
                                                  FE_Parameter::Val ) ;
   }

   if( exp->has_module( "check_prolongation" ) )
   {
      PEL_ModuleExplorer const* ee = 
                         exp->create_subexplorer( 0, "check_prolongation" ) ;
      PRM_CHECK = prms->item( ee->string_data( "param_expected" ) ) ;
      PRM_CHECK->transfer_cell_calculation_requirements( cFE, 
                                                         FE_Parameter::Val ) ;
      UU_CHECK = DOM->set_of_discrete_fields()->item( 
                                          ee->string_data( "field_test" ) ) ;
      cFE->require_field_calculation( UU_CHECK, PDE_LocalFE::N ) ;
      EPS_CHECK = ee->double_data( "dbl_epsilon" ) ;
      MIN_CHECK = ee->double_data( "dbl_minimum" ) ;
      ee->destroy() ;
   }

   if( exp->has_entry( "nb_iterations_max_during_time_stepping" ) )
   {
      NB_ITER_MAX = exp->int_data( "nb_iterations_max_during_time_stepping" ) ;
   }
   
   QRP = GE_QRprovider::object( 
                        exp->string_data( "quadrature_rule_provider" ) ) ;

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;
   
   X_LOC = LA_SeqVector::create( this, 0 ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
void
FE_RefinementAppli:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: run" ) ;
   
   PEL_Communicator const* com = PEL_Exec::communicator() ;

   PEL::out() << endl << "++++++ TIME STEPPING " 
        << " *** INITIAL TIME = " << TIME_IT->initial_time() 
        << " *** FINAL TIME = "   << TIME_IT->final_time() 
        << " ++++++" << endl ;

   do_saving() ;

   PDE_AdapterCHARMS* da = DOM->adapter_CHARMS() ;

   if( da != 0 )
   {
      da->add_excluded_field( UU_EXP ) ;
      if( UU_CHECK != 0 ) da->add_excluded_field( UU_CHECK ) ;
   }

   bool transient = ( TIME_IT->initial_time() + TIME_IT->time_step() < 
                      TIME_IT->final_time() ) ;

   if( transient ) do_before_time_stepping( TIME_IT, da ) ;

   double time = 0.0 ;
   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      PEL::out() << endl << "++++++ ITERATION = "
                         << TIME_IT->iteration_number() 
                 << " *** TIME = "      << TIME_IT->time() 
                 << " *** TIME STEP = " << TIME_IT->time_step() 
                 << " ++++++" << endl << endl ;
      time = TIME_IT->time() ;

      size_t nb_iter = 0 ;
      bool keep_adapting = false ;

      if( da != 0 ) da->reset() ;

      do
      {
         ++nb_iter ;

         do_one_inner_iteration( TIME_IT ) ;

         do_saving() ;

         if( da != 0 )
         {
            PEL::out() << "      adaptation... " ;
            da->adapt() ;
            keep_adapting = com->boolean_or( da->something_changed() ) ;
            if( keep_adapting && PRM_CHECK!=0 )
            {
               bool ok = check_consistency( TIME_IT )  ;
               if( !ok ) do_saving() ;
               PEL_ASSERT( ok ) ;
            }
            if( !keep_adapting )
            {
               PEL::out() << "nothing done" << endl ;
            }
            else
            {
               check_faces_consistency() ;
               check_covering_of_cells_boundary() ;
               
               PEL::out() << endl ;
               da->print_statistics( PEL::out(), 9 ) ;
            }
         }
      }
      while( keep_adapting && nb_iter < NB_ITER_MAX ) ;

      PEL::out() << "      copy \"" << UU->name()
                 << "\" -> \"" << UU_EXP->name() << "\"..." << endl ;
      UU_EXP->set( UU ) ;
      if( da != 0 )
      {
         PEL::out() << "      unsplit of the meshes " << endl ;
         da->unsplit_meshes() ;
         
         check_faces_consistency() ;
         check_covering_of_cells_boundary() ;
         
         da->print_statistics( PEL::out(), 9 ) ;
         if( da->something_changed() )
         {
            do_saving() ;
         }
      }
      
      // DOM->print_grid( cout, 0 ) ;
   }
}

//---------------------------------------------------------------------------
void
FE_RefinementAppli:: do_saving( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: do_saving" ) ;

   SAVER->start_cycle() ;
   PEL::out() << endl << "   +++ SAVE FOR POSTPROCESSING" 
              << " *** CYCLE = " << SAVER->cycle_number()
              << " *** TIME = "  << TIME_IT->time() 
              << " ++++++" << endl << endl ;
   SAVER->save_grid() ;
   SAVER->save_fields( L_UU ) ;
   SAVER->terminate_cycle() ;
}

//---------------------------------------------------------------------------
void
FE_RefinementAppli:: iterate( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: iterate" ) ;

   cout << endl ;
   cout << "====================================" << endl ;
   cout << " ADAPTATION RESULTS " << endl ;
   cout << "====================================" << endl << endl ;

   cout << "********* CELLS ********************" << endl ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cout << "------------------------------" << endl ;
      cFE->print( cout, 0 ) ;
   }

   cout << "********* BOUNDS ********************" << endl ;
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      bFE->set_row_and_col_fields( UU, UU ) ;
      cout << "------------------------------" << endl ;
      cout << "current mesh" << endl ;
      bFE->print_current_mesh( cout, 3 ) ;
      cout << "Local discretization of fields" << endl ;
      bFE->print_local_discretization_of_fields( cout, 3 ) ;
      cout << "Integration points" << endl ;
      bFE->start_IP_iterator( QRP ) ;
      for( ; bFE->valid_IP() ; bFE->go_next_IP() )
      {
         bFE->print_current_IP( cout, 3 ) ;
         bFE->print_values_at_current_IP( cout, 3 ) ;
      }
   }

   cout << "**********************************" << endl ;
   cout << "NMB->nb_global_unknowns()=" << NMB->nb_global_unknowns() << endl ;
   cout << "**********************************" << endl ;
   UU->print( cout, 0 ) ;
   cout << "**********************************" << endl ;
}

//---------------------------------------------------------------------------
bool
FE_RefinementAppli:: check_consistency( FE_TimeIterator const* t_it ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: check_consistency" ) ;

   bool result = true ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double theo = PRM_CHECK->cell_value_at_IP( t_it, cFE ) ;
         double xx = cFE->value_at_IP( UU_CHECK, 0 ) ;
         bool eq = PEL::double_equality( theo, xx, EPS_CHECK, MIN_CHECK ) ;
         if( !eq ) display_error( cFE, theo, xx ) ;
         result = result && eq ;
      }
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE_RefinementAppli:: do_before_time_stepping( FE_TimeIterator const* t_it,
                                              PDE_AdapterCHARMS* da )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: do_before_time_stepping" ) ;

   if( da == 0 ) return ;
   // --------------------
   
   da->reset() ;
   
   bool keep_adapting = false ;
   do
   {
      PEL::out() << "   adaptation... " ;
      da->adapt() ;
      keep_adapting = da->something_changed() ;
      if( !keep_adapting )
      {
         PEL::out() << "nothing done" << endl ;
      }
      else
      {
         PEL::out() << endl ;
         da->print_statistics( PEL::out(), 6 ) ;
         DOM->apply_requests_of_DOFs_values_modules( true ) ;
         do_saving() ;
      }
   } while( keep_adapting ) ;

   UU_EXP->set( UU ) ;
   da->unsplit_meshes() ;
   if( da->something_changed() )
   {
      do_saving() ;
   }
}

//---------------------------------------------------------------------------
void
FE_RefinementAppli:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: do_one_inner_iteration" ) ;
   
   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector aa( nb_dims ) ;
   double dt = t_it->time_step() ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         FE::add_row_col_S( ELEMENT_EQ, cFE, ALPHA/dt ) ;

         FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, LAMBDA ) ;

         double source = ALPHA * cFE->value_at_IP( UU_EXP, L_UU_EXP )/dt ;
         if( FORCE_TERM != 0 )
         {  
            source += FORCE_TERM->cell_value_at_IP( t_it, cFE ) ;
         }
         FE::add_row( ELEMENT_EQ, cFE, source ) ;

         if( AA != 0 )
         {
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               aa( d ) = AA->cell_value_at_IP( t_it, cFE, d ) ;
            }
            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;
         }
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      bFE->set_row_and_col_fields( UU, UU ) ;

      ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 1,
                              bFE->col_field_node_connectivity(), 1 ) ;

      GE_Color const* color = bFE->color() ;
      if( BCs->has_BC( color, UU ) )
      {
	      PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;
         string const& type = ee->string_data( "type" ) ;
         if( type == "robin" )
         {
            double hh = ee->double_data( "prefactor" ) ;
            double uu_out = ee->double_data( "field_out_value" ) ;
            bFE->start_IP_iterator( QRP ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {
               FE::add_row_col_S( ELEMENT_EQ, bFE, hh ) ;
               FE::add_row( ELEMENT_EQ, bFE, hh*uu_out ) ;
            }
         }
         else if( type == "Neumann" )
         {
            double duu = ee->double_data( "derivative_value" ) ;
            bFE->start_IP_iterator( QRP ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {
               FE::add_row( ELEMENT_EQ, bFE, duu ) ;
            }
         }
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }

   PEL::out() << "      resolution... " ;
   A->synchronize() ;
   F->synchronize() ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;

   PEL::out() << " update of \"" << UU->name() << "(" << L_UU << ")\"..." 
              << endl ;
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, X_LOC ) ;
   UU->update_free_DOFs_value( L_UU, X_LOC, NMB->link() ) ;
}

//---------------------------------------------------------------------------
void
FE_RefinementAppli:: check_faces_consistency( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: check_faces_consistency" ) ;
   
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
FE_RefinementAppli:: check_covering_of_cells_boundary( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementAppli:: check_covering_of_cells_boundary" ) ;
   PEL_ASSERT( DOM->nb_space_dimensions() == 2 ) ;
   
   size_t dim = 2 ;
   size_t nb_points = 30 ;
   
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
            DOM->print_grid( cout, 0 ) ;
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

//-----------------------------------------------------------------------
void FE_RefinementAppli:: display_error( PDE_LocalFEcell const* fe, 
                                         double theo, double xx ) const
//-----------------------------------------------------------------------
{
   std::cout << "cell : " << fe->mesh_id() 
             << " theo " << theo << " xx " << xx << std::endl ;
}



