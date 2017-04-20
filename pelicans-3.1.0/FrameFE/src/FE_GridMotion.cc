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

#include <FE_GridMotion.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridMover.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_LocalBCsBuilder.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;
using std::ostringstream ;

FE_GridMotion const* FE_GridMotion::PROTOTYPE = new FE_GridMotion() ;

//---------------------------------------------------------------------------
FE_GridMotion:: FE_GridMotion( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_GridMotion" )
{
}

//---------------------------------------------------------------------------
FE_GridMotion*
FE_GridMotion:: create_replica( PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridMotion:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_GridMotion* result = new FE_GridMotion( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_GridMotion:: FE_GridMotion( PEL_Object* a_owner,
			       PDE_DomainAndFields const* dom,
			       FE_SetOfParameters const* prms,
			       PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , CC( dom->set_of_discrete_fields()->item( 
                                      exp->string_data( "grid_velocity" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , ALPHA( 0 )
   , GRID_MOVER( PDE_GridMover::create( this, dom ) )
   , BCs( dom->set_of_boundary_conditions() )
   , LOCAL_BC( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "FE_GridMotion:: FE_GridMotion" ) ;

   check_field_storage_depth( CC, L_UPDATE ) ;

   if( exp->has_entry( "param_alpha" ) )
   {
      ALPHA = prms->item( exp->string_data( "param_alpha" ) ) ;
      check_param_nb_components( ALPHA, "param_alpha", 1 ) ;
      ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   }

   if( ALPHA != 0 || FE::geometry() == FE::axisymmetrical )
      cFE->require_field_calculation( CC, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( CC, PDE_LocalFE::dN ) ;

   PDE_LinkDOF2Unknown* cc_link = PDE_LinkDOF2Unknown::create( 0, CC,
                                          "sequence_of_the_components",
                                          true ) ;
   NMB = PDE_SystemNumbering::create( this, cc_link ) ;  

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   U = A->create_vector( this ) ;

   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   if( exp->has_entry( "boundary_conditions_types" ) )
   {
      stringVector const& bt = 
                     exp->stringVector_data( "boundary_conditions_types" ) ;
      string const& fname = exp->string_data( "grid_velocity" ) ;
      LOCAL_BC = FE_LocalBCsBuilder::create( this, dom, fname, bt, prms ) ;
      LOCAL_BC->transfer_calculation_requirements( bFE ) ;
   }
}

//---------------------------------------------------------------------------
FE_GridMotion:: ~FE_GridMotion( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_GridMotion:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridMotion:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_GridMotion:: do_one_inner_iteration" ) ;

   NMB->reset() ;

   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   U->re_initialize( n_glob, n_loc ) ;
   
   U_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( U ) ;
   
   start_assembling_timer() ;
   // ---------------------
   
   size_t nb_dims = CC->nb_components() ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( CC, CC ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_dims,
                              cFE->col_field_node_connectivity(), nb_dims ) ;

      for( cFE->start_IP_iterator( QRP ) ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         if( ALPHA != 0 )
         {
            double xx = ALPHA->cell_value_at_IP( t_it, cFE, 0 ) ;
            FE::add_row_col_S( ELEMENT_EQ, cFE, xx ) ;
         }
         FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }

   if( LOCAL_BC != 0 )
   {
      for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
      {
         GE_Color const* color = bFE->color() ;
         if( BCs->has_BC( color, CC ) )
         {
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, CC ) ;
            LOCAL_BC->set_current_BC_type( ee->string_data( "type" ) ) ;
            if( LOCAL_BC->current_BC_type_is_ok() )
            {
               bFE->set_row_and_col_fields( CC, CC ) ;
               ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 
                                       nb_dims,
                                       bFE->col_field_node_connectivity(), 
                                       nb_dims ) ;

               LOCAL_BC->build_current_BC( ELEMENT_EQ, bFE, t_it, QRP ) ;

               PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
            }
         }
      }
   }
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------
   
   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_initial_guess_nonzero( false ) ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, U ) ;
   SOLVER->unset_matrix() ;
   
   stop_solving_timer() ;
   
   if( ! SOLVER->solution_is_achieved() )
   {
      PEL_Error::object()->display_info(
         "*** FE_GridMotion error\n"
         "    No convergence of the linear solver" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      LA_Scatter const* sca = NMB->scatter() ;
      sca->get( U, U_LOC ) ;
      CC->update_free_DOFs_value( L_UPDATE, U_LOC, NMB->link() ) ;
   }

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void FE_GridMotion:: do_after_inner_iterations_stage( 
                                               FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridMotion:: do_after_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "FE_GridMotion:: do_after_inner_iterations_stage" ) ;

   if( verbose_level() >= 2 ) PEL::out() << indent() << "   move grid..." << endl ;
   GRID_MOVER->move_grid( CC, L_UPDATE, t_it->time_step() ) ;

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
FE_GridMotion:: save_other_than_time_and_fields( 
						    FE_TimeIterator const* it,
						    PDE_ResultSaver* rs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridMotion:: save_other_than_time_and_fields" ) ;

   start_total_timer( "FE_GridMotion::save_other_than_time_and_fields" ) ;

   if( verbose_level() >= 2 ) PEL::out() << indent() << "   save grid..." << endl ;
   if( it->is_started() ) rs->save_grid() ;

   stop_total_timer() ;
}
