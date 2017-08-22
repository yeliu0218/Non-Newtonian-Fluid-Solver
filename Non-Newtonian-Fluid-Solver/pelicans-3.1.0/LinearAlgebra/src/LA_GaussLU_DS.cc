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

#include <LA_GaussLU_DS.hh>

#include <LA_DenseMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

LA_GaussLU_DS const* LA_GaussLU_DS::PROTOTYPE = new LA_GaussLU_DS() ;

struct LA_GaussLU_DS_ERROR
{
   static void n0( LA_Matrix const* A ) ;
   static void n1( void ) ;
} ;

//----------------------------------------------------------------------
LA_GaussLU_DS:: LA_GaussLU_DS( void )
//----------------------------------------------------------------------
   : LA_Solver( "LA_GaussLU_DS" )
   , PIV_MIN_VAL( PEL::bad_double() )
   , MAT( 0 )
   , PIV( size_t_vector(0) )
   , DET( PEL::bad_double() )
{
}

//----------------------------------------------------------------------
LA_GaussLU_DS*
LA_GaussLU_DS:: create( PEL_Object* a_owner,
                        double const a_pivot_minimal_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: create" ) ;
   PEL_CHECK_PRE( a_pivot_minimal_value>0. ) ;
   
   LA_GaussLU_DS* result = new LA_GaussLU_DS( a_owner,
                                              a_pivot_minimal_value ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_iterative() ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_GaussLU_DS*
LA_GaussLU_DS:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double const a_piv_min_val = exp->double_data( "pivot_minimal_value" ) ;
   if( a_piv_min_val<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "pivot_minimal_value", "a value greater than 0 is expected" ) ;
   }
   
   LA_GaussLU_DS* result = new LA_GaussLU_DS( a_owner, a_piv_min_val ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_GaussLU_DS:: LA_GaussLU_DS( PEL_Object* a_owner,
                               double const a_pivot_minimal_value )
//----------------------------------------------------------------------
   : LA_Solver( a_owner )
   , PIV_MIN_VAL( a_pivot_minimal_value )
   , MAT( LA_DenseMatrix::create( this, 0, 0 ) )
   , PIV( size_t_vector(0) )
   , DET( PEL::bad_double() )
{
   set_iterative( false ) ;
}

//----------------------------------------------------------------------
LA_GaussLU_DS*
LA_GaussLU_DS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: create_clone" ) ;
   
   LA_GaussLU_DS* result = new LA_GaussLU_DS( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
LA_GaussLU_DS:: LA_GaussLU_DS( PEL_Object* a_owner,
                               LA_GaussLU_DS const* other )
//----------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , PIV_MIN_VAL( other->PIV_MIN_VAL )
   , MAT( LA_DenseMatrix::create( this, 0, 0 ) )
   , PIV( size_t_vector(0) )
   , DET( PEL::bad_double() )
{
}

//----------------------------------------------------------------------
LA_GaussLU_DS:: ~LA_GaussLU_DS( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: set_matrix_self( LA_Matrix const* mat,
                                 bool &ok, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   if( dynamic_cast<LA_SeqMatrix const* >( mat ) == 0 )
      LA_GaussLU_DS_ERROR:: n1() ;
   
   size_t const n = mat->nb_rows() ;
   
   MAT->re_initialize( n, n ) ;
   MAT->set( mat ) ;
   DET = PEL::bad_double() ;
   ok = true ;

   if( n>3 )
   {
      PIV.re_initialize( n ) ;
      factorize_LU( MAT, PIV, ok ) ;
   }
   else
   {
      if( n==1 )
      {
         DET = MAT->item( 0, 0 ) ;
      }
      else if( n==2 )
      {
         DET = MAT->item( 0, 0 )*MAT->item( 1, 1 )
                                 -MAT->item( 1, 0 )*MAT->item( 0, 1 ) ;
      }
      else if( n==3 )
      {
         DET = MAT->item( 0, 0 )*( MAT->item( 1, 1 )*MAT->item( 2, 2 )
                                  -MAT->item( 1, 2 )*MAT->item( 2, 1 ) )
             - MAT->item( 0, 1 )*( MAT->item( 1, 0 )*MAT->item( 2, 2 )
                                  -MAT->item( 1, 2 )*MAT->item( 2, 0 ) )
             + MAT->item( 0, 2 )*( MAT->item( 1, 0 )*MAT->item( 2, 1 )
                                  -MAT->item( 1, 1 )*MAT->item( 2, 0 ) ) ;
      }
      if( PEL::abs( DET )<PIV_MIN_VAL )
      {
         if( stop_on_error() ) LA_GaussLU_DS_ERROR:: n0( MAT ) ;
         ok = false ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: unset_matrix_self( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: unset_matrix_self" ) ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: solve_self(  LA_Vector const* b, 
                             LA_Vector* x,
                             size_t &nb_iter,
                             bool &ok )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( b ) != 0 ) ;
   LA_SeqVector const* bb = static_cast<LA_SeqVector const*>( b ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( x ) != 0 ) ;
   LA_SeqVector* bx = static_cast<LA_SeqVector*>( x ) ;
   
   const size_t n = b->nb_rows() ;

   switch( n )
   {
      case 1 :
         gauss1x1( bb, bx ) ;
         break ;
      case 2 :
         gauss2x2( bb, bx ) ;
         break ;
      case 3 :
         gauss3x3( bb, bx ) ;
         break ;
      default :
         solve_LU( MAT, PIV, bb, bx ) ;
         break ;
   }
   nb_iter = 1 ;
   ok = true ;

   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "Direct solver: \"LA_GaussLU_DS\"" << std::endl ;
   os << s << "   pivot_minimal_value = " << PIV_MIN_VAL << std::endl ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: factorize_LU( LA_DenseMatrix* A,
                              size_t_vector& piv, bool &ok ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: factorize_LU" ) ;
   PEL_CHECK( A != 0 ) ;
   PEL_CHECK( A->nb_rows() == A->nb_cols() ) ;
   PEL_CHECK( piv.size() == A->nb_rows() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   const size_t n = A->nb_rows() ;
  
   size_t i,j,k ;
   
   double aMax ;
   double pivInv, c ;
   size_t jMax ;

   // List of the permutations of the lines :
   for( k=0 ; k<n ; k++ )
   {
      piv( k ) =  k ;
   }
   double val ;
   for( k=0 ; k<n-1 ; k++ )
   {
      // Searching for the biggest element on row k
      aMax = PEL::abs( A->item( piv(k), k ) ) ;
      jMax = k ;
      for( i=k+1 ; i<n ; i++ )
      {
         val = PEL::abs( A->item( piv(i), k ) ) ;
         if( val > aMax )
	 {
            aMax = val ;
            jMax = i ;
	 }
      }
      if( jMax!= k )
      {
         size_t temp = piv(k) ;
         piv( k ) =  piv(jMax) ;
         piv( jMax ) = temp ;
      }
      size_t kRow = piv(k) ;

      if( aMax<PIV_MIN_VAL )
      {
         if( stop_on_error() )
         {
            LA_GaussLU_DS_ERROR:: n0( A) ;
         }
         ok = false ;
         return ;
      }
      
      pivInv = 1.0/A->item( kRow, k ) ;
      for( i=k+1 ; i<n ; i++ )
      {
         c = A->item( piv(i), k )*pivInv ;
         A->set_item( piv(i), k, c ) ;
         c = -c ;
         for( j=k+1 ; j<n ; j++ )
	 {
            A->add_to_item( piv(i), j, c*A->item( kRow, j ) ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: solve_LU( LA_SeqMatrix const* LU,
                          size_t_vector const& piv,
                          LA_SeqVector const* b,
                          LA_SeqVector* x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: solve_LU" ) ;
   PEL_CHECK( LU != 0 ) ;
   PEL_CHECK( b != 0 ) ;
   PEL_CHECK( x != 0 ) ;
   PEL_CHECK( LU->nb_rows() == LU->nb_cols() ) ;
   PEL_CHECK( piv.size() == LU->nb_cols() ) ;
   PEL_CHECK( b->nb_rows() == LU->nb_cols() ) ;
   PEL_CHECK( x->nb_rows() == LU->nb_cols() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t n = LU->nb_rows() ;
   
   for( size_t i=0 ; i<n ; i++ )
   {
      double c = 0.0 ;
      for( size_t j=0 ; j<i ; j++ )
      {
         c += LU->item( piv(i), j )*x->item(j) ;
      }
      x->set_item( i, b->item( piv(i) ) - c ) ;
   }
   for( int i=n-1 ; i>=0 ; i-- )
   {
      double c = 0.0 ;
      for( size_t j=i+1 ; j<n ; j++ )
      {
         c += LU->item( piv(i), j )*x->item(j) ;
      }
      x->set_item( i, ( x->item(i) - c )/LU->item( piv(i), i ) ) ;
   }   
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: gauss1x1( LA_SeqVector const* b,
                          LA_SeqVector* x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: gauss1x1" ) ;
   PEL_CHECK( size() == 1 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   x->set_item( 0, b->item( 0 )/DET ) ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: gauss2x2( LA_SeqVector const* b,
                          LA_SeqVector* x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: gauss2x2" ) ;
   PEL_CHECK( size() == 2 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double const x0 =
          b->item( 0 )*MAT->item( 1, 1 )-b->item( 1 )*MAT->item( 0, 1 ) ;
   double const x1 =
          MAT->item( 0, 0 )*b->item( 1 )-MAT->item( 1, 0 )*b->item( 0 ) ;
   x->set_item( 0, x0/DET ) ;
   x->set_item( 1, x1/DET ) ;
}

//----------------------------------------------------------------------
void
LA_GaussLU_DS:: gauss3x3( LA_SeqVector const* b,
                          LA_SeqVector* x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_GaussLU_DS:: gauss3x3" ) ;
   PEL_CHECK( size() == 3 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double const x0 =
      b->item( 0 )*( MAT->item( 1, 1 )*MAT->item( 2, 2 )  
                            -MAT->item( 1, 2 )*MAT->item( 2, 1 ) )
    - MAT->item( 0, 1 )*( b->item( 1 )*MAT->item( 2, 2 )
                            -MAT->item( 1, 2 )*b->item( 2 ) )
    + MAT->item( 0, 2 )*( b->item( 1 )*MAT->item( 2, 1 )
                            -MAT->item( 1, 1 )*b->item( 2 ) ) ;
   double const x1 =
      MAT->item( 0, 0 )*( b->item( 1 )*MAT->item( 2, 2 )
                            -MAT->item( 1, 2 )*b->item( 2 ) )
    - b->item( 0 )*( MAT->item( 1, 0 )*MAT->item( 2, 2 )
                            -MAT->item( 1, 2 )*MAT->item( 2, 0 ) )
    + MAT->item( 0, 2 )*( MAT->item( 1, 0 )*b->item( 2 )
                            -b->item( 1 )*MAT->item( 2, 0 ) ) ;
   double const x2 =
      MAT->item( 0, 0 )*( MAT->item( 1, 1 )*b->item( 2 )
                            -b->item( 1 )*MAT->item( 2, 1 ) )
    - MAT->item( 0, 1 )*( MAT->item( 1, 0 )*b->item( 2 )
                            -b->item( 1 )*MAT->item( 2, 0 ) )
    + b->item( 0 )*( MAT->item( 1, 0 )*MAT->item( 2, 1 )
                            -MAT->item( 1, 1 )*MAT->item( 2, 0 ) ) ;
   x->set_item( 0, x0/DET ) ;
   x->set_item( 1, x1/DET ) ;
   x->set_item( 2, x2/DET ) ;
}

//----------------------------------------------------------------------
bool
LA_GaussLU_DS:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( EQUIVALENT( is_a_prototype(), MAT==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_a_prototype(), !matrix_is_set() ) ) ;
   PEL_ASSERT( IMPLIES( MAT!=0, MAT->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( MAT!=0, MAT->nb_rows() == MAT->nb_cols() ) ) ;
   PEL_ASSERT( IMPLIES( matrix_is_set(), size() == MAT->nb_rows() ) ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
LA_GaussLU_DS_ERROR:: n0( LA_Matrix const* A )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_GaussLU_DS: " << std::endl ;
   mesg << "    null pivot found" ;
   if( A->nb_rows() < 20 )
   {
      mesg << " for matrix : " << std::endl ;
      A->print( mesg, 5 ) ;
   }
   mesg << std::endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
LA_GaussLU_DS_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_GaussLU_DS: " << std::endl ;
   mesg << "    a matrix of type \"LA_SeqMatrix\" is expected" << std::endl ;
   PEL_Error::object()->raise_internal( mesg.str() ) ;
}
