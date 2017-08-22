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

#include <LA_PenaltyProjection.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Timer.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

LA_PenaltyProjection const* 
LA_PenaltyProjection:: PROTOTYPE = new LA_PenaltyProjection() ;

//----------------------------------------------------------------------
LA_PenaltyProjection:: LA_PenaltyProjection( void )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( "LA_PenaltyProjection" )
{
}

//----------------------------------------------------------------------
LA_PenaltyProjection*
LA_PenaltyProjection:: create_replica( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   LA_PenaltyProjection* result = new LA_PenaltyProjection( a_owner, exp ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PenaltyProjection:: LA_PenaltyProjection( PEL_Object* a_owner,
                                             PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( a_owner, exp )
   , A( 0 )
   , B( 0 )
   , L( 0 )
   , M( 0 )
   , F( 0 )
   , G( 0 )
   , S( 0 )
   , K( 0 )
   , BtS( 0 )
   , P0( 0 )
   , SOLVER_A( 0 )
   , SOLVER_L( 0 )
   , SOLVER_M( 0 )
   , RR( exp->double_data( "augmentation_parameter" ) )
   , PROJP( false )
   , TIMER_RR( 0 )
{
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "solver_A" ) ;
   SOLVER_A = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = exp->create_subexplorer( 0, "solver_L" ) ;
   SOLVER_L = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = exp->create_subexplorer( 0, "solver_Mv" ) ;
   SOLVER_M = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   if( exp->has_entry( "H1_projection_of_explicit_pressure" ) )
   {
      PROJP = exp->bool_data( "H1_projection_of_explicit_pressure" ) ;
   }
   
   if( RR != 0 )
   {
      TIMER_RR = PEL_Timer::create( this ) ;
   }
}

//----------------------------------------------------------------------
LA_PenaltyProjection:: ~LA_PenaltyProjection( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_matrix_prototype_sub( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_matrix_prototype_sub" ) ;
   
   PEL_ASSERT( BtS == 0 ) ;
   BtS = mat->create_matrix( this ) ;
   
   PEL_ASSERT( P0 == 0 ) ;
   P0 = mat->create_vector( this ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: re_initialize_internals_sub( size_t nv_glob, 
                                                    size_t np_glob,
                                                    size_t nv_loc, 
                                                    size_t np_loc,
                                                    size_t& nv_loc_final,
                                                    size_t& np_loc_final )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: re_initialize_internals_sub" ) ;
   
   BtS->re_initialize( nv_glob, np_glob, nv_loc, np_loc ) ;
   P0->re_initialize( np_glob, np_loc ) ;
   
   np_loc_final = P0->nb_local_rows() ;
   nv_loc_final = BtS->nb_local_rows() ;
   
   PEL_ASSERT( BtS->nb_local_cols()  == np_loc_final ) ;
   PEL_ASSERT( P0->nb_local_rows()   == np_loc_final ) ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: dtinv_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: S_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_S( LA_Vector* a_S )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_S" ) ;
   PEL_CHECK_PRE( set_S_PRE( a_S ) ) ;
   
   S = a_S ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: L_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_L( LA_Matrix* a_L )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_L" ) ;
   PEL_CHECK_PRE( set_L_PRE( a_L ) ) ;
   
   L = a_L ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: K_is_required( void ) const
//----------------------------------------------------------------------
{
   return( PROJP ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_K( LA_Vector* a_K )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_K" ) ;
   PEL_CHECK_PRE( set_K_PRE( a_K ) ) ;
   
   K = a_K ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: MV_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_MV( LA_Matrix* a_MV )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_MV" ) ;
   PEL_CHECK_PRE( set_MV_PRE( a_MV ) ) ;
   
   M = a_MV ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                                       LA_Vector* a_F, LA_Vector* a_G,
                                       LA_Matrix* a_C )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: set_system_sub" ) ;
   
   check_zero_C( a_C ) ;
   
   if( dtinv()==PEL::bad_double() || M==0 || L==0 || S==0 || ( K==0 && PROJP ) ) 
      raise_invalid_usage() ;
   
   A = a_A ;
   B = a_B ;
   F = a_F ;
   G = a_G ;
   
   augment_system( RR ) ;
   
   SOLVER_A->set_matrix( A ) ;
   SOLVER_L->set_matrix( L ) ;
   SOLVER_M->set_matrix( M ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: unset_system_sub( void )
//----------------------------------------------------------------------
{
   SOLVER_A->unset_matrix() ;
   SOLVER_L->unset_matrix() ;
   SOLVER_M->unset_matrix() ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                              bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: estimate_unknowns_sub" ) ;
   
   if( has_init_P )
   {
      P0->set( P ) ;
   }
   else
   {
      P0->nullify() ;
   }
   
   if( PROJP )
   {
      start_substep( "   0. projection of explicit pressure" ) ;
      //******************************************************
      // Solve L.P0 = K
      
      if( !inner_solve( SOLVER_L, K, P0, false ) ) return ; //??? false
   }
   
   start_substep( "   1. velocity prediction" ) ;
   //******************************************
   // F -= BT.P0
   // Solve A.U = F
   
   B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
   if( !inner_solve( SOLVER_A, F, U, has_init_U ) ) return ;
   
   start_substep( "   2. pressure increment" ) ;
   //*****************************************
   // G = B.U - G
   // Solve L.P = G
   
   B->multiply_vec_then_add( U, G, 1.0, -1.0 ) ;   
   if( !inner_solve( SOLVER_L, G, P, false ) ) return ;
   
   start_substep( "   3. velocity correction" ) ;
   //******************************************
   // F = M.U
   // F -= BT.PHI
   // Solve M.U = F
   
   M->multiply_vec_then_add( U, F, 1.0, 0.0 ) ;
   B->tr_multiply_vec_then_add( P, F, -1.0, 1.0 ) ;
   if( !inner_solve( SOLVER_M, F, U, true ) ) return ;
      
   start_substep( "   4. pressure update" ) ;
   //**************************************
   // P = P0 + beta*PHI/DT
   
   P->scale( dtinv() ) ;
   P->sum( P0 ) ;
         
   if( RR != 0.0 )
   {
      // G = S*G
      // P = P + RR*G
      
      G->set_as_v_product( G, S ) ;            
      P->sum( G, RR ) ;
   }

   if( verbose_level() > 0 ) PEL::out() << std::endl ;   
   notify_success( true ) ;
}

//----------------------------------------------------------------------
bool
LA_PenaltyProjection:: inner_solve( LA_Solver* solver, 
                                    LA_Vector* rhs, 
                                    LA_Vector* sol,
                                    bool sol_is_init )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: inner_solve" ) ;
   
   solver->set_initial_guess_nonzero( sol_is_init ) ;
   solver->solve( rhs, sol ) ;

   bool result = true ;
   if( !solver->solution_is_achieved() )
   {
      result = false ;
   }
   else
   {
      if( verbose_level() > 0 && solver->is_iterative() )
         PEL::out() << ": " << solver->nb_iterations_achieved()
                    << " iterations" ;
      if( verbose_level() > 0 ) PEL::out() << std::endl ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: print_times( size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: print_times" ) ;

   std::string space( indent_width, ' ' ) ;
   if( TIMER_RR != 0 )
   {
      PEL::out() << space << "LA_PenaltyProjection" << endl ;
      PEL::out() << space << "   augmentation: " ;
      TIMER_RR->print( PEL::out(), 0 ) ;
      PEL::out() << endl ;
   }
}

//----------------------------------------------------------------------
void
LA_PenaltyProjection:: augment_system( double augmentation_parameter )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: augment_system" ) ;
   
   if( augmentation_parameter > 0.0 )
   {
      // Bt = transpose( B )
      BtS->nullify() ;
      BtS->add_tMat( B ) ;
      BtS->synchronize() ;
      
      // BtS = Bt.S
      BtS->scale_as_mat_diag_mat( S ) ;
      
      // BtS = r.Bt.S
      BtS->scale( augmentation_parameter ) ;
   
      // A += BBT.B
      A->add_Mat_Mat( BtS, B, 1.0 ) ;
      A->synchronize() ;
   
      // F += BtS.G
      BtS->multiply_vec_then_add( G, F, 1.0, 1.0 ) ;
   }
}

//-----------------------------------------------------------------------
void
LA_PenaltyProjection:: start_substep( std::string const& title ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PenaltyProjection:: start_substep" ) ;

   if( verbose_level() > 0 ) PEL::out() << indent() << "   " << title ;
}

