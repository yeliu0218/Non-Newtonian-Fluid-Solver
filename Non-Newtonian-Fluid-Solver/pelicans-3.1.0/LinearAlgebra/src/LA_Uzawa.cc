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

#include <LA_Uzawa.hh>

#include <PEL.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_MemoryTracer.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

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

LA_Uzawa const* LA_Uzawa:: PROTOTYPE = new LA_Uzawa() ;

//----------------------------------------------------------------------
LA_Uzawa:: LA_Uzawa( void )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( "LA_Uzawa" )
{
}

//----------------------------------------------------------------------
LA_Uzawa*
LA_Uzawa:: create_replica( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   LA_Uzawa* result = new LA_Uzawa( a_owner, exp ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Uzawa:: LA_Uzawa( PEL_Object* a_owner,
                     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( a_owner, exp )
   , A( 0 )
   , B( 0 )
   , C( 0 )
   , F( 0 )
   , G( 0 )
   , S( 0 )
   , BtS( 0 )
   , ImrC( 0 )
   , F0( 0 )
   , PRES( 0 )
   , P0( 0 )
   , DU( 0 )
   , H( 0 )
   , SOLVER_U( 0 )
   , SOLVER_P( 0 )
   , RR( exp->double_data( "augmentation_parameter" ) )
   , RHO( PEL::bad_double() )
   , TOL_VELO( exp->double_data( "tolerance_on_velocity_increment" ) )
   , TOL_DIV( exp->double_data( "tolerance_on_divergence" ) )
   , MAXITS( 10000 )
{
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "solver_A" ) ;
   SOLVER_U = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   if( exp->has_module( "solver_C" ) )
   {
      ee = exp->create_subexplorer( 0, "solver_C" ) ;
      SOLVER_P = LA_Solver::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   
   exp->test_data( "augmentation_parameter", "augmentation_parameter>0.0" ) ;
   
   if( exp->has_entry( "descent_parameter" ) )
   {
      RHO = exp->double_data( "descent_parameter" ) ;
      exp->test_data( "descent_parameter", "descent_parameter>0.0" ) ;
   }
   else
   {
      RHO = RR ;
   }
   
   if( exp->has_entry( "nb_iterations_max" ) )
   {
      MAXITS = exp->int_data( "nb_iterations_max" ) ;
      exp->test_data( "nb_iterations_max", "nb_iterations_max>0" ) ;
   }
}

//----------------------------------------------------------------------
LA_Uzawa:: ~LA_Uzawa( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_Uzawa:: set_matrix_prototype_sub( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: set_matrix_prototype_sub" ) ;
   
   //??? this matrix my be useless
   PEL_ASSERT( BtS == 0 ) ;
   BtS = mat->create_matrix( this ) ;
   
   PEL_ASSERT( F0 == 0 ) ;
   F0 = mat->create_vector( this ) ;
   
   PEL_ASSERT( P0 == 0 ) ;
   P0 = mat->create_vector( this ) ;
   
   PEL_ASSERT( DU == 0 ) ;
   DU = mat->create_vector( this ) ;
   
   PEL_ASSERT( PRES == 0 ) ;
   PRES = mat->create_vector( this ) ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: re_initialize_internals_sub( size_t nv_glob, 
                                        size_t np_glob,
                                        size_t nv_loc, 
                                        size_t np_loc,
                                        size_t& nv_loc_final,
                                        size_t& np_loc_final )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: re_initialize_internals_sub" ) ;
   
   BtS->re_initialize( nv_glob, np_glob, nv_loc, np_loc ) ;
   F0->re_initialize( nv_glob, nv_loc ) ;
   DU->re_initialize( nv_glob, nv_loc ) ;
   P0->re_initialize( np_glob, np_loc ) ;
   PRES->re_initialize( np_glob, np_loc ) ;
   
   np_loc_final = P0->nb_local_rows() ;
   nv_loc_final = DU->nb_local_rows() ;
   
   PEL_ASSERT( BtS->nb_local_rows()  == nv_loc_final ) ;
   PEL_ASSERT( BtS->nb_local_cols()  == np_loc_final ) ;
   PEL_ASSERT( F0->nb_local_rows()   == nv_loc_final ) ;
   PEL_ASSERT( PRES->nb_local_rows() == np_loc_final ) ;
}

//----------------------------------------------------------------------
bool
LA_Uzawa:: S_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: set_S( LA_Vector* a_S )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: set_S" ) ;
   PEL_CHECK_PRE( set_S_PRE( a_S ) ) ;
   
   S = a_S ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                           LA_Vector* a_F, LA_Vector* a_G,
                           LA_Matrix* a_C )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: set_system_sub" ) ;
   
   if( S == 0 ) raise_invalid_usage() ;
   
   A = a_A ;
   B = a_B ;
   C = a_C ;
   F = a_F ;
   G = a_G ;
   
   if( verbose_level() >= 1 )
   {
      size_t nbF = F->nb_rows() ;
      size_t nbG = G->nb_rows() ;
      size_t nbA = A->nb_stored_items() ;
      size_t nbB = B->nb_stored_items() ;
      PEL::out() << indent() << "   nb velocity unknowns: " << nbF << endl ;
      PEL::out() << indent() << "   nb pressure unknowns: " << nbG << endl ;
      PEL::out() << indent() << "      A : " << nbA << " stored items"
                           << " (" << ((double)nbA)/nbF << " per line)"<< endl ;
      PEL::out() << indent() << "      B : " << nbB << " stored items"
                           << " (" << ((double)nbB)/nbG << " per line)"<< endl ;
   }

   if( C != 0 )
   {
      if( ImrC == 0 ) ImrC = C->create_matrix( this ) ;
      ImrC->re_initialize( C->nb_rows(), C->nb_cols(), 
                           C->nb_local_rows(), C->nb_local_cols() ) ;
      if( H == 0 ) H = C->create_vector( this ) ;
      H->re_initialize( C->nb_rows(), C->nb_local_rows() ) ;
   }

   augment_system() ;
   
   if( verbose_level() >= 1 )
   {
      size_t nbF = F->nb_rows() ;
      size_t nbA = A->nb_stored_items() ;
      PEL::out() << indent() << "      A augmented : " << nbA << " stored items"
                           << " (" << ((double)nbA)/nbF << " per line)"
                           << endl ;
   }
   
   if( verbose_level() >= 1 )
   {
      PEL::out() << indent() << "   set_matrix: " ;
      PEL_MemoryTracer::display_memory( PEL::out(),
                                        PEL_MemoryTracer::used_memory() ) ;
      PEL::out() << " -> " ;
   }

   SOLVER_U->set_matrix( A ) ;
   
   if( verbose_level() >= 1 )
   {
      PEL_MemoryTracer::display_memory( PEL::out(),
                                        PEL_MemoryTracer::used_memory() ) ;
      PEL::out() << endl ;
   }
   
   if( C != 0 )
   {
      PEL_ASSERT( SOLVER_P != 0 ) ;
      SOLVER_P->set_matrix( ImrC ) ;
   }
}

//----------------------------------------------------------------------
void
LA_Uzawa:: unset_system_sub( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: unset_system_sub" ) ;
   
   SOLVER_U->unset_matrix() ;
   if( SOLVER_P != 0 ) SOLVER_P->unset_matrix() ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                  bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: estimate_unknowns_sub" ) ;
   
   //??? these two functions may be merged
   if( C == 0 ) 
      estimate_unknowns_sub_1( has_init_U, U, has_init_P, P ) ;
   else 
      estimate_unknowns_sub_2( has_init_U, U, has_init_P, P ) ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: estimate_unknowns_sub_1( bool has_init_U, LA_Vector* U, 
                                    bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: estimate_unknowns_sub_1" ) ;
   
   if( has_init_P )
   {
      P0->set( P ) ;
   }
   else
   {
      P0->nullify() ;
   }
   
   if( !has_init_U ) U->nullify() ;
   
   double uvar = PEL::max_double() ;
   double delta_P = PEL::max_double() ;

   // F0 = F
   F0->set( F ) ;
   
   size_t nb_it = PEL::bad_index() ;
   
   size_t it = 0 ;   
   do
   {
      ++it ;
      
      // rhsA = F
      F->set( F0 ) ;
      
      // rhsA -= BT.P
      B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      
      double u_old_norm = U->two_norm() ;

      // DU = U
      DU->set( U ) ;

      // U = A-1.F
      if( it == 1 )
         SOLVER_U->set_initial_guess_nonzero( has_init_U ) ;
      else
         SOLVER_U->set_initial_guess_nonzero( true ) ;
      SOLVER_U->solve( F, U ) ;
      
      if( ! SOLVER_U->solution_is_achieved() ) return ; // <---

      if( SOLVER_U->is_iterative() ) nb_it = SOLVER_U->nb_iterations_achieved() ;

      // DU = U - DU
      DU->scale( -1.0 ) ;
      DU->sum( U ) ;
     
      uvar = DU->two_norm() ;
      
      if( u_old_norm > 1.e-12 ) uvar /= u_old_norm ;
      bool cv1 = ( uvar < TOL_VELO ) ;

      // RES = -G
      PRES->set( G ) ;

      // RES = B.U - G
      B->multiply_vec_then_add( U, PRES, 1.0, -1.0 ) ;
                                     
      // RES = INVP*RES
      PRES->set_as_v_product( PRES, S ) ;
      
      // P0 += rho.RES
      P0->sum( PRES, RHO ) ;
         
      delta_P = PRES->two_norm() ;
      
      bool cv2 = ( delta_P < TOL_DIV ) ;
      
      if( cv1 && cv2 )
      {
         // P = P0
         P->set( P0 ) ;
         
         if( verbose_level()==1 ) print_end_errors( it, delta_P, uvar ) ;
         
         notify_success( true ) ;
         return ; // <---
      }
      if( verbose_level() > 1 ) print_errors( it, delta_P, uvar, nb_it ) ;

   } while( it < MAXITS ) ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: estimate_unknowns_sub_2( bool has_init_U, LA_Vector* U, 
                                    bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: estimate_unknowns_sub_2" ) ;
   
   PEL_ASSERT( SOLVER_P != 0 ) ;
      
   if( has_init_P )
   {
      P0->set( P ) ;
   }
   else
   {
      P0->nullify() ;
   }
   
   if( !has_init_U ) U->nullify() ;
   
   double uvar = PEL::max_double() ;
   double delta_P = PEL::max_double() ;
   
   F0->set( F ) ; // F0 = F
   
   size_t nb_it_1 = PEL::bad_index() ;
   size_t nb_it_2 = PEL::bad_index() ;
   
   size_t it = 0 ;
   do
   {
      ++it ;

      //  F = F0
      //  F = F - (BT).P0
      //  PRES = C.P0
      //  F = F - (BtS).PRES
      F->set( F0 ) ;
      B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      C->multiply_vec_then_add( P0, PRES, 1.0, 0.0 ) ;
      BtS->multiply_vec_then_add( PRES, F, -1.0, 1.0 ) ;

      //  DU = U
      double u_old_norm = U->two_norm() ;
      DU->set( U ) ;
      
      //  computation of U solution of (A).U = F
      if( it == 1 )
         SOLVER_U->set_initial_guess_nonzero( has_init_U ) ;
      else
         SOLVER_U->set_initial_guess_nonzero( true ) ;
      SOLVER_U->solve( F, U ) ;

      if( ! SOLVER_U->solution_is_achieved() ) return ; // <---

      if( SOLVER_U->is_iterative() ) nb_it_1 = SOLVER_U->nb_iterations_achieved() ;
      
      //  DU = U - DU
      DU->scale( -1.0 ) ;
      DU->sum( U ) ;
      
      uvar = DU->two_norm() ;
      
      if( u_old_norm > 1.e-12 ) uvar /= u_old_norm ;
      bool cv1 = ( uvar < TOL_VELO ) ;

      //  PRES = - G
      //  PRES = PRES + (B).U
      PRES->set( G ) ;
      
      B->multiply_vec_then_add( U, PRES, 1.0, -1.0 ) ;

      //  H = PRES
      //  H = S.H
      //  H = P0 + d.H
      H->set_as_v_product( PRES, S ) ;
      
      H->scale( RHO ) ;
      H->sum( P0, 1.0 ) ;

      //  computation of P0 solution of (ImrC).P0 = H     
      if( it == 1 )
         SOLVER_P->set_initial_guess_nonzero( has_init_U ) ;
      else
         SOLVER_P->set_initial_guess_nonzero( true ) ;
      SOLVER_P->solve( H, P0 ) ;

      if( ! SOLVER_P->solution_is_achieved() ) return ; // <---
      
      if( SOLVER_P->is_iterative() ) nb_it_2 = SOLVER_P->nb_iterations_achieved() ;

      //  PRES = PRES + (C).P0
      //  PRES = (MpInv).PRES
      C->multiply_vec_then_add( P0, PRES, 1.0, 1.0 ) ;
      
      PRES->set_as_v_product( PRES, S ) ;
       
      delta_P = PRES->two_norm() ;
      
      bool cv2 = ( delta_P < TOL_DIV ) ;
      
      if( cv1 && cv2 )
      {
         //  P = P0
         P->set( P0 ) ;
         
         if( verbose_level()==1 ) print_end_errors_2( it, delta_P, uvar ) ;
         
         notify_success( true ) ;
         return ; // <---
         
      }
      if( verbose_level() > 1 )
         print_errors_2( it, delta_P, uvar, nb_it_1, nb_it_2 ) ;

   } while( it < MAXITS ) ;
}

//----------------------------------------------------------------------
void
LA_Uzawa:: augment_system( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Uzawa:: augment_system" ) ;
      
   // BtS = transpose( B )
   BtS->nullify() ; //??? my be useless
   BtS->add_tMat( B ) ;
   BtS->synchronize() ;
      
   // BtS = Bt.S
   BtS->scale_as_mat_diag_mat( S ) ;
      
   // BtS = r.Bt.S
   BtS->scale( RR ) ;
   
   // A += r.Bt.S.B
   A->add_Mat_Mat( BtS, B, 1.0 ) ;
   
   // F += r.Bt.S.G
   BtS->multiply_vec_then_add( G, F, 1.0, 1.0 ) ;
   
   if( C != 0 )
   {
      //  ImrC = -rho.C
      ImrC->add_Mat( C, -RHO ) ;
      ImrC->synchronize() ;
         
      //  ImrC = -rho.S.C
      ImrC->scale_as_diag_mat_mat( S ) ;
         
      //  ImrC = Id - rho.S.C
      ImrC->start_local_modifs() ;
      size_t i = ImrC->row_distribution()->first_local_index() ;
      for( ; i<ImrC->row_distribution()->local_index_limit() ; ++i )
      {
         ImrC->add_to_item( i, i, 1.0 ) ;
      }
      ImrC->stop_local_modifs() ;
   }
}

//-----------------------------------------------------------------------
void 
LA_Uzawa:: print_errors( size_t n, double err1, double err2, size_t nb_it )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( n==1 )
   {
      PEL::out() << indent()
                 << setw( 3+15 ) << "S[BU-G]" 
                 << setw( 15 )   << "DU/U" ;
      if( nb_it != PEL::bad_index() ) PEL::out() << setw( 9 ) << "AM_its" ;
      PEL::out() << endl ;
   }
   PEL::out() << indent()
              << setw( 3 ) << n 
              << setprecision( 6 ) << setw( 15 ) << err1
              << setprecision( 6 ) << setw( 15 ) << err2 ;
   if( nb_it != PEL::bad_index() ) PEL::out() << setw( 9 ) << nb_it ;
   PEL::out() << endl ;

   PEL::out().flags( original_flags ) ;
}

//-----------------------------------------------------------------------
void 
LA_Uzawa:: print_end_errors( size_t n, double err1, double err2 )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   PEL::out() << indent()
              << n << " iterations" 
              << "   S[BU-G] = " << setprecision( 6 ) << err1
              << "   DU/U = "    << setprecision( 6 ) << err2 
              << endl ;

   PEL::out().flags( original_flags ) ;
}

//??? print_errors and print_errors_2 should be merged
//-----------------------------------------------------------------------
void 
LA_Uzawa:: print_errors_2( size_t n,
                           double err1, double err2,
                           size_t nb_it_1, size_t nb_it_2 )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( n==1 )
   {
      PEL::out() << indent()
           << setw( 3+15 ) << "[BU+CP-G]/Mp" 
           << setw( 15 )   << "DU/U" ;
      if( nb_it_1 != PEL::bad_index() ) PEL::out() << setw( 9 ) << "AM_its" ;
      if( nb_it_2 != PEL::bad_index() ) PEL::out() << setw( 9 ) << "IC_its" ;
      PEL::out() << endl ;
   }
   PEL::out() << indent() 
        << setw( 3 ) << n 
        << setprecision( 6 ) << setw( 15 ) << err1
        << setprecision( 6 ) << setw( 15 ) << err2 ;
   if( nb_it_1 != PEL::bad_index() ) PEL::out() << setw( 9 ) << nb_it_1 ;
   if( nb_it_2 != PEL::bad_index() ) PEL::out() << setw( 9 ) << nb_it_2 ;
   PEL::out() << endl ;

   PEL::out().flags( original_flags ) ;
}

//??? print_end_errors and print_end_errors_2 should be merged
//-----------------------------------------------------------------------
void 
LA_Uzawa:: print_end_errors_2( size_t n, double err1, double err2 )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   PEL::out() << indent() 
        << n << " iterations" 
        << "   [BU-G]/Mp = " << setprecision( 6 ) << err1 //????? C ????
        << "   DU/U = "      << setprecision( 6 ) << err2 
        << endl ;

   PEL::out().flags( original_flags ) ;
}


