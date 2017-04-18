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

#include <LA_CG_IS.hh>

#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

LA_CG_IS const* LA_CG_IS::PROTOTYPE = new LA_CG_IS() ;

//----------------------------------------------------------------------
LA_CG_IS:: LA_CG_IS( void )
//----------------------------------------------------------------------
   : LA_IterativeSolver( "LA_CG_IS" )
   , SIZE( 0 )
   , R( 0 )
   , Z( 0 )
   , P( 0 )
{
}

//----------------------------------------------------------------------
LA_CG_IS*
LA_CG_IS:: create_replica( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CG_IS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_CG_IS* result = new LA_CG_IS( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CG_IS:: LA_CG_IS( PEL_Object* a_owner, PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, exp )
   , SIZE( 0 )
   , R( 0 )
   , Z( 0 )
   , P( 0 )
{
   PEL_LABEL( "LA_CG_IS:: LA_CG_IS" ) ;
   
   set_norm_type( exp ) ;
}

//----------------------------------------------------------------------
LA_CG_IS*
LA_CG_IS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CG_IS:: create_clone" ) ;
   
   LA_CG_IS* result = new LA_CG_IS( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CG_IS:: LA_CG_IS( PEL_Object* a_owner,
                     LA_CG_IS const* other )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, other )
   , SIZE( 0 )
   , R( 0 )
   , Z( 0 )
   , P( 0 )
{
   PEL_LABEL( "LA_CG_IS:: LA_CG_IS" ) ;
}

//----------------------------------------------------------------------
LA_CG_IS:: ~LA_CG_IS( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CG_IS:: ~LA_CG_IS" ) ;
}

//----------------------------------------------------------------------
void
LA_CG_IS:: reset_internals( LA_Vector const* prototype ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CG_IS:: reset_internals" ) ;
   PEL_CHECK( invariant() ) ;
   
   if( ( R == 0 ) ||
       ( R->implementation() != prototype->implementation() ) ||
       ( SIZE != prototype->nb_rows() ) )
   {
      if( R!=0 )
      {
         destroy_possession( R ) ; R = 0 ;
         destroy_possession( Z ) ; Z = 0 ;
         destroy_possession( P ) ; P = 0 ;
      }
      R = prototype->create_vector( this ) ;
      Z = prototype->create_vector( this ) ;
      P = prototype->create_vector( this ) ;
   }
   
   SIZE = prototype->nb_rows() ;
}

//----------------------------------------------------------------------
void
LA_CG_IS:: do_solve( LA_Matrix const* A,
                     LA_Vector const* b,
                     LA_Preconditioner* prec,
                     bool zero_initial_guess,
                     LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CG_IS:: do_solve" ) ;
   PEL_CHECK_PRE( do_solve_PRE( A, b, prec, zero_initial_guess, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   reset_internals( b ) ;
   
   double r_norm  = PEL::bad_double() ;
   double r_z_new = PEL::bad_double() ;
   double r_z_old = 1.0 ;
   double alpha   = PEL::bad_double() ;
   double beta    = PEL::bad_double() ;
   double dpi     = PEL::bad_double() ;
   
   bool no_error = true ;

   //--- INITIALIZATION   
   size_t it = 0 ;
   R->set( b ) ;
   if( zero_initial_guess )
   {
      x->set( 0.0 ) ;
   }
   else
   {
      A->multiply_vec_then_add( x , R, -1.0, 1.0 ) ; // r =  b - Ax      
   }

   apply_prec( it, prec, R, Z, no_error ) ;
   if( !no_error )
   {
      set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
      return ; // <---
   }
   
   r_z_new = Z->dot( R ) ;
   if( norm_type() == LA_ConvergenceTest::Preconditioned )
   {
      r_norm = Z->two_norm() ;
   }
   else if( norm_type() == LA_ConvergenceTest::Unpreconditioned )
   {
      r_norm = R->two_norm() ;
   }
   test_convergence( it, r_norm, A, b, prec, zero_initial_guess ) ;
   if( convergence_achieved() ) return ; // <---

   do
   {
      it++ ;
      if( r_z_new == 0 )
      {
         set_converged_reason( it, LA_ConvergenceTest::ConvergedAtol ) ;
         break ; // <---
      }
      else if( r_z_new < 0.0 )
      {
         set_converged_reason( it, LA_ConvergenceTest::DivergedIndefinitePC ) ;
         break ; // <---
      }
      if( it == 1 )
      {
         P->set( Z ) ;
         beta = 0.0 ;
      }
      else
      {
         beta = r_z_new / r_z_old ;
         P->scale( beta ) ;
         P->sum( Z ) ;
      }
      r_z_old = r_z_new ;
      A->multiply_vec_then_add( P, Z ) ;
      dpi = P->dot( Z ) ;
      if( dpi <= 0.0 ) 
      {
         set_converged_reason( it, LA_ConvergenceTest::DivergedIndefiniteMat ) ;
         break ; // <---
      }
      
      alpha = r_z_new / dpi ;
      
      x->sum( P,  alpha ) ;  // x <- x + alpha p 
      R->sum( Z, -alpha ) ;  // r <- r - alpha z
      
      apply_prec( it, prec, R, Z, no_error ) ;
      if( !no_error )
      {
         set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
         break ; // <---
      }
      r_z_new = R->dot( Z ) ;

      if( norm_type() == LA_ConvergenceTest::Preconditioned )
      {
         r_norm = Z->two_norm() ;
      }
      else if( norm_type() == LA_ConvergenceTest::Unpreconditioned )
      {
         r_norm = R->two_norm() ;
      }
      
      test_convergence( it, r_norm,
                        A, b, prec, zero_initial_guess ) ;
      if( convergence_achieved() ) break ; // <---
   }
   while( it < max_nb_iterations() ) ;
   
}

//----------------------------------------------------------------------
bool
LA_CG_IS:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( R==0 || 
       ( R->nb_rows() == SIZE && Z!=0 && P!=0 ) ) ;
   PEL_ASSERT( R==0 || 
       ( Z->nb_rows() == SIZE && Z->implementation()==R->implementation() ) ) ;
   PEL_ASSERT( R==0 || 
       ( P->nb_rows() == SIZE && P->implementation()==R->implementation() ) ) ;
   return( true ) ;
}
