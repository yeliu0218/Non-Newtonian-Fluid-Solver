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

#include <LA_BiCGSTAB_IS.hh>

#include <LA_ConvergenceTest.hh>
#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

LA_BiCGSTAB_IS const* LA_BiCGSTAB_IS::PROTOTYPE = new LA_BiCGSTAB_IS() ;

//----------------------------------------------------------------------
LA_BiCGSTAB_IS:: LA_BiCGSTAB_IS( void )
//----------------------------------------------------------------------
   : LA_IterativeSolver( "LA_BiCGSTAB_IS" )
   , SIZE( 0 )
   , R( 0 )
   , RP( 0 )
   , V( 0 )
   , T( 0 )
   , S( 0 )
   , P( 0 )
{
}

//----------------------------------------------------------------------
LA_BiCGSTAB_IS*
LA_BiCGSTAB_IS:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_BiCGSTAB_IS* result = new LA_BiCGSTAB_IS( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BiCGSTAB_IS:: LA_BiCGSTAB_IS( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, exp )
   , SIZE( 0 )
   , R( 0 )
   , RP( 0 )
   , V( 0 )
   , T( 0 )
   , S( 0 )
   , P( 0 )
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: LA_BiCGSTAB_IS" ) ;
}

//----------------------------------------------------------------------
LA_BiCGSTAB_IS:: LA_BiCGSTAB_IS( PEL_Object* a_owner,
                                 LA_BiCGSTAB_IS const* other )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, other )
   , SIZE( 0 )
   , R( 0 )
   , RP( 0 )
   , V( 0 )
   , T( 0 )
   , S( 0 )
   , P( 0 )
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: LA_BiCGSTAB_IS" ) ;
}

//----------------------------------------------------------------------
LA_BiCGSTAB_IS:: ~LA_BiCGSTAB_IS( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: ~LA_BiCGSTAB_IS" ) ;
}

//----------------------------------------------------------------------
LA_BiCGSTAB_IS*
LA_BiCGSTAB_IS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: create_clone" ) ;
   
   LA_BiCGSTAB_IS* result = new LA_BiCGSTAB_IS( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;   
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_BiCGSTAB_IS:: reset_internals( LA_Vector const* prototype )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: reset_internals" ) ;
   PEL_CHECK( invariant() ) ;
   
   if( ( R == 0 ) ||
       ( R->implementation() != prototype->implementation() ) ||
       ( SIZE != prototype->nb_rows() ) )
   {
      if( R!=0 )
      {
         destroy_possession( R )  ; R  = 0 ;
         destroy_possession( RP ) ; RP = 0 ;
         destroy_possession( V )  ; V  = 0 ;
         destroy_possession( T )  ; T  = 0 ;
         destroy_possession( S )  ; S  = 0 ;
         destroy_possession( P )  ; P  = 0 ;
      }
      R  = prototype->create_vector( this ) ;
      RP = prototype->create_vector( this ) ;
      V  = prototype->create_vector( this ) ;
      T  = prototype->create_vector( this ) ;
      S  = prototype->create_vector( this ) ;
      P  = prototype->create_vector( this ) ;
   }

   SIZE = prototype->nb_rows() ;
}

//----------------------------------------------------------------------
void
LA_BiCGSTAB_IS:: do_solve( LA_Matrix const* A,
                           LA_Vector const* b,
                           LA_Preconditioner* prec,
                           bool zero_initial_guess,
                           LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BiCGSTAB_IS:: do_solve" ) ;
   PEL_CHECK_PRE( do_solve_PRE( A, b, prec, zero_initial_guess, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   reset_internals( b ) ;

   double rho   = PEL::bad_double() ;
   double beta  = PEL::bad_double() ;
   double omega = PEL::bad_double() ;
   double d1 = PEL::bad_double(), d2 = PEL::bad_double() ;

   double rhoold   = 1.0 ;
   double alpha    = 1.0 ;
   double omegaold = 1.0 ;

   bool no_error = true ;

   //--- INITIALIZATION :
   size_t it = 0 ;
   T->set( b ) ;
   if( zero_initial_guess )
   {
      x->set( 0.0 ) ;
   }
   else
   {
      A->multiply_vec_then_add( x , T, -1.0, 1.0 ) ; // t =  b - Ax
   }
   
   apply_prec( it, prec, T, R, no_error ) ;
   if( !no_error )
   {
      set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
      return ; // <---
   }
   
   test_convergence( it, R->two_norm(), A, b, prec, zero_initial_guess ) ;
   if( convergence_achieved() ) return ; // <---

   RP->set( R ) ;
   P->nullify() ;
   V->nullify() ;

   do
   {
      ++it ;

      rho = R->dot( RP ) ;
      if( rho == 0.0 )
      {
         set_converged_reason( it, LA_ConvergenceTest::DivergedBreakdown ) ;
         break ; // <---
      }

      beta = (rho/rhoold) * (alpha/omegaold) ;
      P->sum( V, -omegaold ) ;

      //???? doit etre un VecAYPY
      P->scale( beta ) ;
      P->sum( R, 1.0 ) ;

      A->multiply_vec_then_add( P, T ) ;
      apply_prec( it, prec, T, V, no_error ) ;
      if( !no_error ) break ; // <---

      d1 = V->dot( RP ) ;

      alpha = rho/d1 ;

      //??? c'est un VecWAXPY
      S->set( R ) ; // s <- r
      S->sum( V, -alpha ) ; // s <- s - alpha*v

      A->multiply_vec_then_add( S, R ) ; // r <- A s
      
      apply_prec( it, prec, R, T, no_error ) ; // t <- C^{-1} s
      if( !no_error )
      {
         set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
         break ; // <---
      }

      d1 = S->dot( T ) ;
      d2 = T->dot( T ) ;
      if( d2 == 0.0 )
      {
         // t is 0.  
         // if s is 0, then alpha v == r, and hence alpha p may be 
         //            our solution.  Give it a try? 
         d1 = S->dot( S ) ;
         if( d1 == 0.0 )
         {
            set_converged_reason( it, LA_ConvergenceTest::DivergedBreakdown ) ;
            break ; // <---
         }
         x->sum( P, alpha ) ; // x <- x + alpha*p
         it++ ;
         set_converged_reason( it, LA_ConvergenceTest::ConvergedRtol ) ;
         break ; // <---
      }
      PEL_ASSERT( d2 != 0 ) ;

      omega = d1 / d2 ;
      x->sum( P, alpha ) ;
      x->sum( S, omega ) ;

      //????? doit etre un VecWAXPY
      R->set( S ) ;
      R->sum( T, -omega ) ;

      rhoold   = rho;
      omegaold = omega;

      test_convergence( it, R->two_norm(), 
                        A, b, prec, zero_initial_guess ) ;
      if( convergence_achieved() ) break ; // <---
      
   }
   while( it < max_nb_iterations() ) ;
   
}

//----------------------------------------------------------------------
bool
LA_BiCGSTAB_IS:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( R==0 || 
     ( R->nb_rows() == SIZE && RP!=0 && V!=0 && T!=0 && S!=0 && P!=0 ) ) ;
   PEL_ASSERT( R==0 || 
     ( RP->nb_rows() == SIZE && RP->implementation()==R->implementation() ) ) ;
   PEL_ASSERT( R==0 || 
     ( V->nb_rows() == SIZE && V->implementation()==R->implementation() ) ) ;
   PEL_ASSERT( R==0 || 
     ( T->nb_rows() == SIZE && T->implementation()==R->implementation() ) ) ;
   PEL_ASSERT( R==0 || 
     ( S->nb_rows() == SIZE && S->implementation()==R->implementation() ) ) ;
   PEL_ASSERT( R==0 || 
     ( P->nb_rows() == SIZE && P->implementation()==R->implementation() ) ) ;
   return( true ) ;
}
