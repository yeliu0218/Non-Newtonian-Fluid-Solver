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

#include <LA_SymmetricMatrix.hh>

#include <fstream>

#include <LA_SeqVector.hh>
#include <LA_DenseMatrix.hh>
#include <LA_SymmetricMatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

LA_SymmetricMatrix const*
LA_SymmetricMatrix:: PROTOTYPE = new LA_SymmetricMatrix() ;

//----------------------------------------------------------------------
LA_SymmetricMatrix*
LA_SymmetricMatrix:: create( PEL_Object* a_owner, size_t dim )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: create" ) ;

   LA_SymmetricMatrix* result = new LA_SymmetricMatrix( a_owner, dim ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == "LA_SymmetricMatrix" ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->nb_rows() == dim ) ;
   PEL_CHECK_POST( result->nb_cols() == dim ) ;
   PEL_CHECK_POST( result->nb_stored_items() == dim*dim ) ;
   PEL_CHECK_POST( result->is_symmetric() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                      FORALL( ( size_t j=0 ; j<result->nb_cols() ; ++j ),
                         result->item(i,j) == 0.0 ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix*
LA_SymmetricMatrix:: create_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: create_matrix" ) ;
   PEL_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;

   LA_SymmetricMatrix* result = new LA_SymmetricMatrix( a_owner, SIZE  ) ;

   PEL_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix*
LA_SymmetricMatrix:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   int size = 0 ;
   if( exp->has_entry( "size" ) )
   {
      size = exp->int_data( "size" ) ;
      exp->test_data( "size", "size > 0" ) ;
   }

   LA_SymmetricMatrix* result = LA_SymmetricMatrix::create( a_owner, size ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix:: LA_SymmetricMatrix( void )
//----------------------------------------------------------------------
   : LA_SeqMatrix( "LA_SymmetricMatrix" )
   , SIZE( 0 )
   , DIAG( 0 )
{
   PEL_LABEL( "LA_SymmetricMatrix:: LA_SymmetricMatrix( prototype )" ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == "LA_SymmetricMatrix" ) ;
   PEL_CHECK_POST( is_symmetric() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( nb_rows() == 0 ) ;
   PEL_CHECK_POST( nb_cols() == 0 ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix:: LA_SymmetricMatrix( PEL_Object* a_owner, size_t dim )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_SymmetricMatrix" )
   , SIZE( dim )
   , DIAG( dim*(dim+1)/2 )
{
   PEL_LABEL( "LA_SymmetricMatrix:: LA_SymmetricMatrix" ) ;

   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( name() == "LA_SymmetricMatrix" ) ;
   PEL_CHECK_POST( is_symmetric() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( nb_rows() == dim ) ;
   PEL_CHECK_POST( nb_cols() == dim ) ;
   PEL_CHECK_POST( nb_stored_items() == dim*dim ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix:: ~LA_SymmetricMatrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: re_initialize_with_global_sizes( size_t a_nb_rows,
                                                      size_t a_nb_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: re_initialize_with_global_sizes" ) ;
   PEL_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows,
                                                       a_nb_cols ) ) ;

   if( a_nb_rows != a_nb_cols )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_SymmetricMatrix error:\n"
         "    re_initialize expects equal number"
         " of rows and columns" ) ;
   }

   SIZE = a_nb_rows ;
   DIAG.re_initialize( SIZE*(SIZE+1)/2 ) ;

   PEL_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows,
                                                         a_nb_cols ) ) ;
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: is_desynchronizable" ) ;

   bool result = false ;

   PEL_CHECK_POST( result == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( SIZE ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   return( SIZE ) ;
}

//----------------------------------------------------------------------
double
LA_SymmetricMatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: item" ) ;
   PEL_CHECK_PRE( item_PRE( i, j ) ) ;

   double result = DIAG( index( i, j ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   return( SIZE*SIZE ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   return( DIAG.size()*sizeof(double) ) ;
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrix:: is_symmetric( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: is_symmetric" ) ;

   bool result = true ;

   PEL_CHECK_POST( is_symmetric_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_SymmetricMatrix:: create_stored_item_iterator(
                                             PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: create_stored_item_iterator" ) ;
   PEL_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;

   LA_MatrixIterator* result =
                   LA_SymmetricMatrixIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( set_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   DIAG( index( i, j ) ) = x ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   DIAG( index( i, j ) ) += x ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: set_stored_items( double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: set_stored_items" ) ;

   DIAG.set( val ) ;
   synchronize() ;

   PEL_CHECK_POST( set_stored_items_POST( val ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;

   LA_SymmetricMatrix const* sA =
                         dynamic_cast<LA_SymmetricMatrix const* >( A ) ;
   if( sA==0 )
   {
      LA_SeqMatrix::set( A, same_pattern ) ;
   }
   else
   {
      DIAG = sA->DIAG ;
   }
   synchronize() ;

   PEL_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;

   for( size_t i=0 ; i<DIAG.size() ; ++i )
   {
      DIAG(i) *= alpha ;
   }

   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: add_Mat(
                   LA_Matrix const* A, double alpha, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: add_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;

   LA_SymmetricMatrix const* SA =
                           dynamic_cast<LA_SymmetricMatrix const*>( A ) ;
   if( SA != 0 )
   {
      for( size_t i=0 ; i<DIAG.size() ; ++i )
      {
         DIAG(i) += SA->DIAG(i)*alpha ;
      }
   }
   else
   {
      PEL_Error::object()->raise_bad_types( this, "add_Mat", A ) ;
   }

   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: add_Mat_tMat( LA_SeqMatrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: add_Mat_tMat( same matrix )" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( A->nb_rows()==nb_rows() ) ;
   PEL_CHECK_PRE( A->is_synchronized() ) ;

   size_t n = A->nb_rows() ;
   size_t p = A->nb_cols() ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<=i ; j++ )
      {
         double x = item(i,j) ;
         for( size_t k=0 ; k<p ; k++ )
         {
            x += alpha*A->item( i, k )*A->item( j, k ) ;
         }
         set_item( i, j, x) ;
      }
   }
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: invert( double& det )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: invert" ) ;
   PEL_CHECK_PRE( nb_rows()>0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t n = nb_rows() ;
   LA_DenseMatrix* tempMat = LA_DenseMatrix::create( 0, n , n ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         tempMat->set_item( i, j, item( i, j ) ) ;
      }
   }
   tempMat->synchronize() ;
   tempMat->invert( det ) ;
   if( det != 0.0 )
   {
      for( size_t i=0 ; i<n ; i++ )
      {
         for( size_t j=i ; j<n ; j++ )
         {
            set_item( i, j, tempMat->item( i, j ) ) ;
         }
      }
   }
   else
   {

      re_initialize_with_global_sizes( 0, 0 ) ;
   }
   synchronize() ;
   tempMat->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
   PEL_CHECK_POST( IMPLIES( det==0, nb_rows()==0 && nb_cols()==0 ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: eigen_reduce( size_t matord,
                                   LA_SeqVector* egval,
                                   LA_DenseMatrix* egvecs,
                                   bool& success ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: eigen_reduce" ) ;
   PEL_CHECK_PRE( matord <= nb_rows() ) ;
   PEL_CHECK_PRE( egval != 0 ) ;
   PEL_CHECK_PRE( egval->nb_rows() == matord ) ;
   PEL_CHECK_PRE( egvecs != 0 ) ;
   PEL_CHECK_PRE( egvecs->nb_rows() == matord ) ;
   PEL_CHECK_PRE( egvecs->nb_cols() == matord ) ;

   LA_SeqVector* w2 = LA_SeqVector::create( 0, matord ) ;
   reduce( matord, egval, w2, egvecs ) ;
   success = eigenV( matord, egval, w2, egvecs ) ;
   w2->destroy() ; w2 = 0 ;
   if( !success )
   {
      egval->set( PEL::bad_double() ) ;
      egvecs->set_stored_items( PEL::bad_double() ) ;
   }
   egval->synchronize() ;
   egvecs->synchronize() ;

   PEL_CHECK_POST( egval->is_synchronized() ) ;
   PEL_CHECK_POST( egvecs->is_synchronized() ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: reduce( size_t matord,
                             LA_SeqVector* diagTDMA,
                             LA_SeqVector* subdTDMA,
                             LA_DenseMatrix* orth ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: reduce" ) ;
   PEL_CHECK_PRE( matord <= nb_rows() ) ;
   PEL_CHECK_PRE( diagTDMA != 0 ) ;
   PEL_CHECK_PRE( diagTDMA->nb_rows() == matord ) ;
   PEL_CHECK_PRE( subdTDMA != 0 ) ;
   PEL_CHECK_PRE( subdTDMA->nb_rows() == matord ) ;
   PEL_CHECK_PRE( orth != 0 ) ;
   PEL_CHECK_PRE( orth->nb_rows() == matord ) ;
   PEL_CHECK_PRE( orth->nb_cols() == matord ) ;

   size_t i ;
   for( i=0 ; i < matord ; i++ )
   {
      for( size_t j=0 ; j < matord ; j++ )
      {
         orth->set_item( i, j, item( i, j ) ) ;
      }
   }
   for( i=matord-1 ; i >= 1 ; i-- )
   {
      size_t l = i-1 ;
      double h = 0. ;
      double scal = 0. ;
      if( l > 0 )
      {
         size_t k ;
         for(  k = 0 ; k <= l ; k++ )
         {
            scal += PEL::abs( orth->item( i, k ) ) ;
         }
         if( scal == 0.0 )
         // skip transformation
         {
            subdTDMA->set_item( i, orth->item( i, l ) ) ;
         }
         else
         {
            for( k= 0 ; k <= l ; k++ )
            {
               orth->set_item( i, k, orth->item( i, k ) / scal ) ;
               h += orth->item( i, k )*orth->item( i, k ) ;
            }
            double f = orth->item( i, l ) ;
            double g = ( f >= 0 ? -PEL::sqrt( h ) : PEL::sqrt( h ) ) ;
            subdTDMA->set_item( i, scal*g ) ;
            h -= f*g ;
            orth->set_item( i, l, f-g ) ;
            f = 0. ;
            size_t j ;
            for(  j = 0 ; j <= l ; j++ )
            {
               orth->set_item( j, i, orth->item( i, j )/h ) ;
               g = 0. ;
               for( k=0 ; k <= j ; k++ )
               {
                  g += orth->item( j, k )*orth->item( i, k ) ;
               }
               for( k=j+1 ; k <= l ; k++ )
               {
                  g += orth->item( k, j )*orth->item( i, k ) ;
               }
               subdTDMA->set_item( j, g/h ) ;
               f += subdTDMA->item( j )*orth->item( i, j ) ;
            }
            double hh = f/( h+h ) ;
            for( j=0 ; j <= l ; j++ )
            {
               f = orth->item( i, j ) ;
               g = subdTDMA->item( j ) - hh*f ;
               subdTDMA->set_item( j, g ) ;
               for( k = 0 ; k <= j ; k++ )
               {
                  orth->set_item( j, k, orth->item( j, k )-
                                     ( f*subdTDMA->item( k ) +
                                       g*orth->item( i, k ) ) ) ;
               }
            }
         }
      }
      else
      {
         subdTDMA->set_item( i, orth->item( i, l ) ) ;
      }
      diagTDMA->set_item( i, h ) ;
   }


   //  Accumulation of transformation matrices
   subdTDMA->set_item( 0, 0. ) ;
   diagTDMA->set_item( 0, 0. ) ;
   for( i=0 ; i < matord ; i++ )
   {
      int l = i-1 ;
      if( diagTDMA->item( i ) != 0. )
      {
         for( int j=0 ; j <= l ; j++ )
         {
            double g = 0. ;
            int k ;
            for( k=0 ; k<= l ; k++ )
               g += orth->item( i, k )*orth->item( k, j ) ;
            for( k = 0 ; k <= l ; k++ )
               orth->set_item( k, j, orth->item( k, j ) -
                                     g*orth->item( k, i ) ) ;
         }
      }
      diagTDMA->set_item( i, orth->item( i, i ) ) ;
      orth->set_item( i, i, 1.0 ) ;
      for( int j=0 ; j <= l ; j++ )
      {
         orth->set_item( i, j, 0. ) ;
         orth->set_item( j, i, 0. ) ;
      }
   }
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrix:: eigenV( size_t matord,
                             LA_SeqVector* diagTDMA,
                             LA_SeqVector* subdTDMA,
                             LA_DenseMatrix* eigv ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: eigenV" ) ;
   PEL_CHECK_PRE( matord <= nb_rows() ) ;
   PEL_CHECK_PRE( diagTDMA != 0 ) ;
   PEL_CHECK_PRE( diagTDMA->nb_rows() == matord ) ;
   PEL_CHECK_PRE( subdTDMA != 0 ) ;
   PEL_CHECK_PRE( subdTDMA->nb_rows() == matord ) ;
   PEL_CHECK_PRE( eigv != 0 ) ;
   PEL_CHECK_PRE( eigv->nb_rows() == matord ) ;
   PEL_CHECK_PRE( eigv->nb_cols() == matord ) ;

   if( matord == 1 )
   {
      return( true ) ;
   }

   size_t i ;
   for( i=1 ; i < matord ; i++ )
   {
      subdTDMA->set_item( i-1, subdTDMA->item( i ) ) ;
   }

   subdTDMA->set_item( matord-1, 0.0 ) ;
   for( size_t l=0 ; l < matord ; l++ )
   {
      size_t iter = 0 ;
      size_t m ;

      do
      {
         for( m=l ; m < matord-1 ; m++ )
         {
            double dd = PEL::abs( diagTDMA->item( m ) )
                       +PEL::abs( diagTDMA->item( m+1 ) ) ;
            if( dd!=0.0 &&
                PEL::abs( subdTDMA->item( m )/dd )< 1.E-14 )
            {
               break ;
            }
         }
         if( m != l )
         {
            if( iter++ == 30 )
            {
               return( false ) ;
            }
            // form shift
            double g= ( diagTDMA->item( l+1 )-diagTDMA->item( l ) )
                     /( 2.0*subdTDMA->item( l ) ) ;
            double r = PEL::sqrt( g*g+1.0 ) ;
            double signr = ( g >= 0.0 ? r : -r ) ;

            g = diagTDMA->item( m ) - diagTDMA->item( l )
               +subdTDMA->item( l )/( g+signr ) ;

            double s = 1.0 ;
            double c = 1.0 ;
            double p = 0. ;

            for( i = m-1 ; i >= l && i<=m-1 ; i-- )
            {
               double f = s*subdTDMA->item( i ) ;
               double b = c*subdTDMA->item( i ) ;
               if( PEL::abs( f ) >= PEL::abs( g ) )
               {
                  c = g/f ;
                  r = PEL::sqrt( c*c+1.0 ) ;
                  subdTDMA->set_item( i+1, f*r ) ;
                  s = 1.0/r ;
                  c *= s ;
               }
               else
               {
                  s = f/g ;
                  r = PEL::sqrt( ( s*s ) + 1.0 ) ;
                  subdTDMA->set_item( i+1, g*r ) ;
                  c = 1.0/r ;
                  s *= c ;
               }
               g = diagTDMA->item( i+1 ) - p ;
               r = ( diagTDMA->item( i )-g )*s + 2.0*c*b ;
               p = s*r ;
               diagTDMA->set_item( i+1, g + p ) ;
               g = c*r - b ;

               // Form eigenvectors
               for( size_t k=0 ; k < matord ; k++ )
               {
                  f = eigv->item( k, i+1 ) ;
                  eigv->set_item( k, i+1,
                                      s*eigv->item( k, i ) + c*f ) ;
                  eigv->set_item( k, i,
                                      c*eigv->item( k, i ) - s*f ) ;
               }
            }

            diagTDMA->set_item( l, diagTDMA->item( l ) - p ) ;
            subdTDMA->set_item( l, g ) ;
            subdTDMA->set_item( m, 0 ) ;
         }

      } while( m != l ) ;
   }

   //  Order eigenvectors and eigenvalues
   for( size_t ii=1 ; ii < matord ; ii++ )
   {
      i = ii-1 ;
      size_t k = i ;
      double p = diagTDMA->item( i ) ;
      size_t j ;
      for( j = ii ; j < matord ; j++ )
      {
         if( PEL::abs( diagTDMA->item( j ) ) < p )
         {
            k = j ;
            p = diagTDMA->item( j ) ;
         }
      }
      if( k > i )
      {
         diagTDMA->set_item( k, diagTDMA->item( i ) ) ;
         diagTDMA->set_item( i, p ) ;
         for( j = 0 ; j < matord ; j++ )
         {
            p = eigv->item( j, i ) ;
            eigv->set_item( j, i, eigv->item( j, k ) ) ;
            eigv->set_item( j, k, p ) ;
         }
      }
   }
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrix:: writeMM( std::string const& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: writeMM" ) ;
   std::ofstream out( file.c_str() ) ;
   if( !out )
   {
      PEL_Error::object()->raise_file_handling( file, "open" ) ;
   }
   out << "%%MatrixMarket matrix array real symmetric" << std::endl ;
   out << nb_rows() << " " << nb_cols() << std::endl ;
   out.precision( 16 ) ;
   for( size_t j=0; j<nb_cols() ; j++ )
   {
      for( size_t i=j ; i<nb_rows() ; i++ )
      {
         out << item( i, j ) << " " ;
      }
      out << std::endl ;
   }
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrix:: index( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix:: index" ) ;
   PEL_CHECK_PRE( i<nb_rows() ) ;
   PEL_CHECK_PRE( j<nb_rows() ) ;

   size_t ii ;
   size_t kk ;
   if( i>j )
   {
      kk = i-j ;
      ii = j ;
   }
   else
   {
      kk = j-i ;
      ii = i ;
   }

   size_t result = ii+SIZE*kk-kk*(kk-1)/2 ;

   PEL_CHECK_POST( result<DIAG.size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_SeqMatrix::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrix:: is_symmetric_POST( bool result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::is_symmetric_POST( result ) ) ;
   PEL_ASSERT( result ) ;
   return( true ) ;
}




