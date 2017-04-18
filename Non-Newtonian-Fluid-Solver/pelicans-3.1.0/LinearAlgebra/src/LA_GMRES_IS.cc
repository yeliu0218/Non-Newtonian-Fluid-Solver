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
#include <LA_GMRES_IS.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_SeqVector.hh>
#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_DenseMatrix.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>


#include <iostream>
#include <sstream>
using std::endl ;

// happy breakdown tolerance
//    it may not exactly correspond to that of PETSc, but it should...
const double LA_GMRES_IS:: HAPTOL = 1.e-30 ;

LA_GMRES_IS const* LA_GMRES_IS::PROTOTYPE = new LA_GMRES_IS() ;

//----------------------------------------------------------------------
LA_GMRES_IS:: LA_GMRES_IS( void )
//----------------------------------------------------------------------
   : LA_IterativeSolver( "LA_GMRES_IS" )
   , H( 0, 0 )
   , S( 0 )
   , cRot( 0 )
   , sRot( 0 )
   , UPY( 0 )
{
}

//----------------------------------------------------------------------
LA_GMRES_IS*
LA_GMRES_IS:: create_replica( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_GMRES_IS* result = new LA_GMRES_IS( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_GMRES_IS:: LA_GMRES_IS( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
  : LA_IterativeSolver( a_owner, exp )
  , H( 0, 0 )
  , S( 0 )
  , cRot( 0 )
  , sRot( 0 )
  , UPY( 0 )
  , vTemp( 0 )
{
   int restart = exp->int_data( "restart" ) ;
   exp->test_data( "restart", "restart>0" ) ;

   init( restart ) ;
}

//----------------------------------------------------------------------
LA_GMRES_IS:: LA_GMRES_IS( PEL_Object* a_owner,
                             LA_GMRES_IS const* other )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, other )
   , H( 0, 0 )
   , S( 0 )
   , cRot( 0 )
   , sRot( 0 )
   , UPY( 0 )
   , vTemp( 0 )
{
   init( other->RESTART ) ;
}

//----------------------------------------------------------------------
void
LA_GMRES_IS:: init( size_t restart )
//----------------------------------------------------------------------
{
   SIZE = 0 ;
   RESTART = restart ;
   vBasis = 0 ;
   H.re_initialize( restart+1, restart ) ;
   S.re_initialize( restart+1 ) ;
   cRot.re_initialize( restart ) ;
   sRot.re_initialize( restart ) ;
   UPY.re_initialize( restart+1 ) ;
   vBasis = PEL_Vector::create( this, restart+1  ) ;
}

//----------------------------------------------------------------------
LA_GMRES_IS:: ~LA_GMRES_IS( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_GMRES_IS*
LA_GMRES_IS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: create_clone" ) ;
   return( new LA_GMRES_IS( a_owner, this ) ) ;
}

//----------------------------------------------------------------------
void
LA_GMRES_IS:: do_solve( LA_Matrix const* A,
                        LA_Vector const* b,
                        LA_Preconditioner* prec,
                        bool zero_initial_guess,
                        LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: do_solve" ) ;
   PEL_CHECK_PRE( do_solve_PRE( A, b, prec, zero_initial_guess, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   reset_internals( b ) ;

   double xx, rNorm ;
   bool no_error = true ;

   bool zero_init = zero_initial_guess ;
   size_t it = 0 ;
   do
   {
      r->set( b ) ;
      if( zero_init )
      {
         x->set( 0.0 ) ;
      }
      else
      {
         A->multiply_vec_then_add( x, r, -1.0, 1.0 ) ;
      }
      apply_prec( it, prec, r, vTemp, no_error ) ;
      if( !no_error )
      {
         set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
         break ; // <---
      }

      rNorm = vTemp->two_norm() ;
      test_convergence( it, rNorm, A, b, prec, zero_init ) ;
      if( convergence_achieved() ) break ; // <---

      //--- CONSTRUCT THE FIRST VECTOR OF THE OTHONORMAL BASIS OF THE
      //--- KRYLOV SUBSPACE, AND THE CURRENT VALUE OF s
      LA_Vector* vb0 = static_cast<LA_Vector*>( vBasis->at( 0 ) ) ;
      vb0->set( 0.0 ) ;
      vb0->sum( vTemp, 1.0/rNorm ) ;
      S.set( 0.0 ) ;
      S( 0 ) = rNorm ;

      int k=0 ;
      int i ;
      do
      {
         ++it ;

         //--- CONSTRUCT THE ( k+1 )-TH VECTOR OF THE ORTHONORMAL BASIS
         //--- OF THE KRYLOV SUBSPACE, I.E. CONSTRUCT THE k-TH COLUMN
         //--- OF H ORTHONORMAL TO THE PREVIOUS ( K-1 )-TH COLUMNS
         LA_Vector* vbk   = static_cast<LA_Vector*>( vBasis->at( k   ) ) ;
         LA_Vector* vbkp1 = static_cast<LA_Vector*>( vBasis->at( k+1 ) ) ;
         A->multiply_vec_then_add( vbk, Avk ) ;
         apply_prec( it, prec, Avk, vbkp1, no_error ) ;
         if( !no_error )
         {
            set_converged_reason( it, LA_ConvergenceTest::PreconditionerFailure ) ;
            break ; // <---
         }

         // Use modified Gram-Schmidt orthogonalization method
         for( i=0 ; i<=k ; i++ )
         {
            LA_Vector* vbi = static_cast<LA_Vector*>( vBasis->at( i ) ) ;
            H( i, k ) = vbi->dot( vbkp1 ) ;
            vbkp1->sum( vbi, - H( i, k ) ) ;
         }

         double xNorm = vbkp1->two_norm() ;

         //--- APPLY THE PREVIOUS GIVENS ROTATIONS TO THE K-TH COLUMN
         //--- OF H
         for( i=0 ; i<=( k-1 ) ; i++ )
         {
            xx = H( i, k ) ;
            H( i,   k ) = cRot( i )*xx - sRot( i )*H( i+1, k ) ;
            H( i+1, k ) = sRot( i )*xx + cRot( i )*H( i+1, k ) ;
         }

         //--- TEST FOR A POSSIBLE BREAKDOWN
         if( xNorm<HAPTOL )
         {
            update_solution( k, H, S, UPY, vBasis, x ) ;
            r->set( b ) ;
            A->multiply_vec_then_add( x, r, -1.0, 1.0 ) ;
            rNorm = r->two_norm() ;

            test_convergence( it, rNorm, A, b, prec, zero_init ) ;
            no_error = convergence_achieved() ;
            if( no_error )
            {
               set_converged_reason( it,
                             LA_ConvergenceTest::ConvergedHappyBreakdown ) ;
            }
            else
            {
               set_converged_reason( it,
                             LA_ConvergenceTest::DivergedBreakdown ) ;
            }
            break ; // <---
         }

         H( k+1, k ) = xNorm ;
         vbkp1->scale( 1.0/H( k+1, k ) ) ;

         //--- CONSTRUCT THE K-TH GIVENS ROTATION MATRIX
         xx = PEL::sqrt( H( k, k )*H( k, k ) + H( k+1, k )*H( k+1, k ) ) ;
         cRot( k ) = H( k, k )/xx  ;
         sRot( k ) =  -1.0*H( k+1, k )/xx ;

         //--- APPLY THE K-TH GIVENS ROTATION MATRIX TO H SO THAT
         //--- H( k+1, k )=0
         H( k,   k ) = cRot( k )*H( k, k ) - sRot( k )*H( k+1, k )  ;
         H( k+1, k ) = 0.0 ;

         //--- APPLY THE K-TH GIVENS ROTATION MATRIX TO s
         xx = S( k ) ;
         S( k )   = cRot( k )*xx - sRot( k )*S( k+1 ) ;
         S( k+1 ) = sRot( k )*xx + cRot( k )*S( k+1 ) ;

         //--- SET THE NORM OF THE RESIDUAL AT STEP K
         rNorm = PEL::abs( S( k+1 ) ) ;

         //--- TEST ON THE NORM OF THE RESIDUAL
         test_convergence( it, rNorm, A, b, prec, zero_init ) ;
         if( convergence_achieved() )
         {
            update_solution( k, H, S, UPY, vBasis, x ) ;
            break ;
         }
         else
         {
            k++ ;
         }
      }
      while( ( k < (int)RESTART ) && ( it < max_nb_iterations() ) ) ;

      if( !no_error || convergence_achieved() ) break ; // <---

      //--- COMPUTATION OF THE NEW INITIAL VECTOR FOR THE RESTART
      k = RESTART-1 ;
      update_solution( k, H, S, UPY, vBasis, x ) ;

      zero_init = false ;
   }
   while( it < max_nb_iterations() ) ;
}

//----------------------------------------------------------------------
void
LA_GMRES_IS:: reset_internals( LA_Vector const* prototype )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: reset_internals" ) ;
   PEL_CHECK_PRE( prototype != 0 ) ;
   
   size_t n = prototype->nb_rows() ;
   bool vTemp_is_compatible = true ;

   if(  ( vTemp == 0 ) ||
       !( vTemp->row_distribution()->is_compatible( 
                                    prototype->row_distribution() ) ) || 
        ( vTemp->implementation() != prototype->implementation() ) )
   {
      vTemp_is_compatible = false ;

      if( vTemp!=0 )
      {
         destroy_possession( vTemp ) ; vTemp = 0 ;
         destroy_possession( r ) ; r = 0 ;
         destroy_possession( Avk ) ; Avk = 0 ;
      }
      vTemp = prototype->create_vector( this ) ;
      r = prototype->create_vector( this ) ;
      Avk = prototype->create_vector( this ) ;
   }

   if( ( vBasis->at( 0 ) != 0 ) &&
       ( SIZE == n ) &&
       vTemp_is_compatible )
   {
      // No initialization needed
   }
   else
   {
      size_t bsize = RESTART+1 ;
      for( size_t i=0 ; i<bsize ; i++ )
      {
         PEL_Object* old = vBasis->at( i ) ;
         if( old != 0 ) vBasis->destroy_possession( old ) ;

         LA_Vector* v = prototype->create_vector( vBasis ) ;
         v->set( 0.0 ) ;
         vBasis->set_at( i, v ) ;
      }
   }
   SIZE = n ;
   
   PEL_CHECK_POST( vTemp->implementation() == prototype->implementation() ) ;
   PEL_CHECK_POST( vTemp->is_synchronized() ) ;
   PEL_CHECK_POST( vTemp->row_distribution()->is_compatible( 
                                            prototype->row_distribution() ) ) ;
   PEL_CHECK_POST( vTemp->nb_rows() == prototype->nb_rows() ) ;

   PEL_CHECK_POST( Avk->implementation() == prototype->implementation() ) ;
   PEL_CHECK_POST( Avk->is_synchronized() ) ;
   PEL_CHECK_POST( Avk->row_distribution()->is_compatible( 
                                            prototype->row_distribution() ) ) ;
   PEL_CHECK_POST( Avk->nb_rows() == prototype->nb_rows() ) ;

   PEL_CHECK_POST( r->implementation() == prototype->implementation() ) ;
   PEL_CHECK_POST( r->is_synchronized() ) ;
   PEL_CHECK_POST( r->row_distribution()->is_compatible( 
                                            prototype->row_distribution() ) ) ;
   PEL_CHECK_POST( r->nb_rows() == prototype->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_GMRES_IS:: update_solution( int k,
                               doubleArray2D const& Hmat,
                               doubleVector const& s,
                               doubleVector& work,
                               PEL_Vector* vec_basis,
                               LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: update_solution" ) ;

   int i, j ;

   //--- COMPUTATION OF THE COMPONENTS OF THE CURRENT SOLUTION VECTOR
   //--- WITH RESPECT TO THE BASIS ( V0, ...Vk )  OF THE k-TH KRYLOV
   //--- SPACE : RESOLUTION OF A LINEAR SYSTEM WITH UPPER
   //--- TRIANGULAR MATRIX
   for( i=k ; i>=0 ; i-- )
   {
      work( i ) = s( i ) ;
      for( j=i+1 ; j<=k ; j++ )
      {
         work( i ) -= Hmat( i, j ) * work( j ) ;
      }
      work( i ) /= Hmat( i, i ) ;
   }

   //--- COMPUTATION OF THE CURRENT SOLUTION VECTOR
   for( i=0 ; i<=k ; i++ )
   {
      LA_Vector* vbi = static_cast<LA_Vector*>( vec_basis->at( i ) ) ;

      x->sum( vbi, work( i ) ) ;
   }
}

//----------------------------------------------------------------------
void
LA_GMRES_IS:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GMRES_IS:: print_more" ) ;
   std::string s( indent_width, ' ') ;
   os << s << "restart: " << RESTART << std::endl ;
}
