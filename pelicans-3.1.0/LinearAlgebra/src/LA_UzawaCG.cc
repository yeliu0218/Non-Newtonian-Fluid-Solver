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

#include <LA_UzawaCG.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Solver.hh>
#include <LA_UzawaPreconditioner.hh>
#include <LA_Vector.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

LA_UzawaCG const* LA_UzawaCG:: PROTOTYPE = new LA_UzawaCG() ;

struct LA_UzawaCG_ERROR
{
   static void n0( void ) ;
   static void n1( size_t i_iter ) ;
   static void n3( void ) ;
} ;

//----------------------------------------------------------------------
LA_UzawaCG:: LA_UzawaCG( void )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( "LA_UzawaCG" )
{
}

//----------------------------------------------------------------------
LA_UzawaCG*
LA_UzawaCG:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   LA_UzawaCG* result = new LA_UzawaCG( a_owner, exp ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_UzawaCG:: LA_UzawaCG( PEL_Object* a_owner,
                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( a_owner, exp )
   , A( 0 )
   , B( 0 )
   , C( 0 )
   , F( 0 )
   , G( 0 )
   , S( 0 )
   , PREC( 0 )
   , SOLVER_A( 0 )
   , TOL( exp->double_data( "tolerance" ) )
   , MAXITS( 10000 )
{
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "solver_A" ) ;
   SOLVER_A = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = exp->create_subexplorer( 0, "LA_UzawaPreconditioner" ) ;
   PREC = LA_UzawaPreconditioner::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   if( exp->has_entry( "nb_iterations_max" ) )
   {
      MAXITS = exp->int_data( "nb_iterations_max" ) ;
      exp->test_data( "nb_iterations_max", "nb_iterations_max>0" ) ;
   }
}

//----------------------------------------------------------------------
LA_UzawaCG:: ~LA_UzawaCG( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: set_matrix_prototype_sub( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: re_initialize_internals_sub( size_t nv_glob, 
                                          size_t np_glob,
                                          size_t nv_loc, 
                                          size_t np_loc,
                                          size_t& nv_loc_final, 
                                          size_t& np_loc_final  )
//----------------------------------------------------------------------
{
   //??? completely false
   nv_loc_final = nv_loc ;
   np_loc_final = np_loc ;
}

//----------------------------------------------------------------------
bool
LA_UzawaCG:: S_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: set_S( LA_Vector* a_S )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: set_S" ) ;
   PEL_CHECK_PRE( set_S_PRE( a_S ) ) ;
   
   S = a_S ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                             LA_Vector* a_F, LA_Vector* a_G,
                             LA_Matrix* a_C )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: set_system_sub" ) ;
   
   if( S == 0 ) raise_invalid_usage() ;
   
   A = a_A ;
   B = a_B ;
   C = a_C ;
   F = a_F ;
   G = a_G ;
   
   SOLVER_A->set_matrix( A ) ;
   
   PREC->build( A, B, C ) ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: unset_system_sub( void )
//----------------------------------------------------------------------
{
   SOLVER_A->unset_matrix() ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                    bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: estimate_unknowns_sub" ) ;
   
   PEL_ASSERT( C != 0 ) ; //??? this case should be taken into account

   size_t iter = 0 ;
      
   //??? reduce the number of local vectors
   //??? make them attributes
   LA_Vector* r   = G->create_vector( 0 ) ;
   LA_Vector* r0  = G->create_vector( 0 ) ;
   LA_Vector* gg  = G->create_vector( 0 ) ;
   LA_Vector* g0  = G->create_vector( 0 ) ;
   LA_Vector* d   = G->create_vector( 0 ) ;
   LA_Vector* zt  = F->create_vector( 0 ) ;
   LA_Vector* z   = G->create_vector( 0 ) ;
   LA_Vector* residCV = G->create_vector( 0 ) ;

   double rho, rNorm, rg, r0g0 ;

   bool no_error = true ;
   
   //--- INITIALIZATION OF y
   
   if( !has_init_P )
   {
      LA_Vector* vecAux = F->create_vector( 0 ) ;
      LA_Vector* rhsB = G->create_vector( 0 ) ;
      
      SOLVER_A->set_initial_guess_nonzero( true ) ;
      vecAux->set( F ) ;
      SOLVER_A->solve( F, vecAux ) ;
      no_error = SOLVER_A->solution_is_achieved() ;
      if( !no_error )
      {
         P->nullify() ;
         LA_UzawaCG_ERROR::n3() ;
      }
      else
      {
         rhsB->set( G ) ;
         B->multiply_vec_then_add( vecAux, rhsB, 1.0, -1.0 ) ;
         PREC->solve( rhsB, P ) ;
         no_error = PREC->successful_solve() ;
         if( !no_error )
         {
            P->nullify() ;
            LA_UzawaCG_ERROR::n3() ;
         }
      }

      vecAux->destroy() ; vecAux = 0 ;
      rhsB->destroy() ; rhsB = 0 ;
   }

   //--- INITIALIZATION OF x
   LA_Vector* rhsA = F->create_vector( 0 ) ;
   rhsA->set( F ) ;
   B->tr_multiply_vec_then_add( P, rhsA, -1.0, 1.0 ) ;

   SOLVER_A->set_initial_guess_nonzero( false ) ;
   SOLVER_A->solve( rhsA, U ) ;
   if( !SOLVER_A->solution_is_achieved() )
   {
      U->nullify() ;
      no_error = false ;
      LA_UzawaCG_ERROR::n0() ;
   }
      
   //--- COMPUTATION OF THE INITIAL RESIDUAL
   r->set( G ) ;
   B->multiply_vec_then_add( U, r, -1.0, 1.0 ) ;
   C->multiply_vec_then_add( P, r, -1.0, 1.0 ) ;
   
   iter = 0 ;

   residCV->set_as_v_product( r, S ) ;
   rNorm = residCV->max_norm() ;
   bool converged = ( rNorm < TOL ) ;

   if( no_error && verbose_level() > 1 )
   {
      print_iteration( iter, rNorm, SOLVER_A ) ;
   }
   converged = ( rNorm < TOL ) ;

   while( no_error && (!converged) && (iter < MAXITS) )
   {
      iter++ ;
   
      //--- COMPUTATION OF THE DESCENT DIRECTION FOR Y
      PREC->solve( r, gg ) ;
      no_error = PREC->successful_solve() ;
      if( !no_error )
      {
         LA_UzawaCG_ERROR::n3() ;
         break ;
      }
      if( iter==1 )
      {
         d->set( gg ) ;
      }
      else
      {
         r0g0 = r0->dot( g0 ) ;
         rg = r->dot( gg ) ;
         d->scale( rg/r0g0 ) ;
         d->sum( gg ) ;
       }

      //--- COMPUTATION OF THE DESCENT DIRECTION FOR X
      B->tr_multiply_vec_then_add( d, rhsA ) ;
      SOLVER_A->solve( rhsA, zt ) ;
      if( !SOLVER_A->solution_is_achieved() )
      {
         no_error = false ;
         LA_UzawaCG_ERROR::n1( iter ) ;
         break ;
      }

      //--- COMPUTATION OF THE DESCENT INCREMENT
      B->multiply_vec_then_add( zt, z ) ;
      C->multiply_vec_then_add( d, z, -1.0, 1.0 ) ;
      rho = d->dot( r ) / d->dot( z ) ;

      //--- UPDATE OF X AND Y
      U->sum( zt, rho ) ;
      P->sum( d, -rho ) ;

      //--- UPDATE OF THE RESIDUAL
      r0->set( r ) ;
      g0->set( gg ) ;
      r->sum( z, -rho ) ;

      //--- Convergence test :      
      residCV->set_as_v_product( r, S ) ;
      rNorm = residCV->max_norm() ;
      if( verbose_level() > 1 ) 
      {
         print_iteration( iter, rNorm, SOLVER_A ) ;
      }
      converged = ( rNorm < TOL ) ;
   }

   if( no_error && verbose_level() > 1 ) print_residuals( U, P ) ;
   
   if( no_error && converged )
   {
      notify_success( true ) ;
   }
   
   rhsA->destroy() ; rhsA = 0 ;
   r->destroy()    ; r  = 0   ;
   r0->destroy()   ; r0 = 0   ;
   gg->destroy()   ; gg = 0   ;
   g0->destroy()   ; g0 = 0   ;
   d->destroy()    ; d  = 0   ;
   zt->destroy()   ; zt = 0   ;
   z->destroy()    ; z  = 0   ;
   residCV->destroy() ; residCV = 0 ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: print_iteration( size_t iter, double residual,
                              LA_Solver const* Asolver ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: print_iteration" ) ;
   
   std::ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   std::streamsize p = PEL::out().precision() ;
   PEL::out() << std::setprecision( 7 ) ;
   PEL::out() << "   "
              << "iteration " << std::setw( 3 ) << iter
              << " : residual = " << residual ;
   if( Asolver->is_iterative() )
   {
      PEL::out() << " ( " << std::setw( 4 )
                 << Asolver-> nb_iterations_achieved()
                 << " internal iterations)" ;
   }
   PEL::out() << "\n" ;
   PEL::out() << std::setprecision(p) ;
   PEL::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
LA_UzawaCG:: print_residuals( LA_Vector const* U,
                              LA_Vector const* P ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaCG:: print_residuals" ) ;
   
   LA_Vector* err1 = F->create_vector( 0 ) ;
   err1->set( F ) ;
   A->multiply_vec_then_add( U, err1, -1., 1. ) ;
   B->tr_multiply_vec_then_add( P, err1, -1., 1. ) ;

   LA_Vector* err2=  G->create_vector( 0 ) ;
   err2->set( G ) ;
   B->multiply_vec_then_add( U, err2, -1., 1. ) ;
   C->multiply_vec_then_add( P, err2, -1., 1. ) ;

   std::ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   PEL::out() << std::setprecision( 7 ) ;
   PEL::out() << "            ||A*x+tB*y-F||=" << err1->two_norm() << "\n"
              << "            ||B*x+C*y-g||=" << err2->two_norm() << "\n" ;
   PEL::out().flags( original_flags ) ;
   err1->destroy() ; err1 = 0 ;
   err2->destroy() ; err2 = 0 ;
}

//internal--------------------------------------------------------------
void 
LA_UzawaCG_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_UzawaCG:" << std::endl ;
   mesg << "***   internal solve failure during the" << endl ;
   mesg << "***   computation of the initial residual" << endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
LA_UzawaCG_ERROR:: n1( size_t i_iter )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_UzawaCG:" << std::endl ;
   mesg << "***   internal solve failure" << endl ;
   mesg << "***   iteration number : " << i_iter << endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
LA_UzawaCG_ERROR:: n3( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_UzawaCG:" << std::endl ;
   mesg << "***   preconditioner inversion failure" << endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

