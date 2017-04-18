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

#include <LA_DenseMatrix.hh>

#include <LA_DenseMatrixIterator.hh>
#include <LA_SymmetricMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>

LA_DenseMatrix* LA_DenseMatrix:: PROTOTYPE = new LA_DenseMatrix() ;

//----------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix:: create( PEL_Object* a_owner,
                         size_t a_nb_rows, size_t a_nb_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: create" ) ;
   PEL_CHECK_PRE( EQUIVALENT( a_nb_rows == 0, a_nb_cols == 0 ) ) ;

   LA_DenseMatrix* result =
                  new LA_DenseMatrix( a_owner, a_nb_rows, a_nb_cols ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == "LA_DenseMatrix" ) ;
   PEL_CHECK_POST( !result->is_symmetric() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   PEL_CHECK_POST( !result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->nb_rows() == a_nb_rows ) ;
   PEL_CHECK_POST( result->nb_cols() == a_nb_cols ) ;
   PEL_CHECK_POST( result->nb_stored_items() == a_nb_rows*a_nb_cols ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                      FORALL( ( size_t j=0 ; j<result->nb_cols() ; ++j ),
                         result->item(i,j) == 0.0 ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix:: create( PEL_Object* a_owner, doubleArray2D const& amat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: create" ) ;
   PEL_CHECK_PRE( amat.index_bound(0)!=0 && amat.index_bound(1)!=0 ) ;

   LA_DenseMatrix* result =
      new LA_DenseMatrix( a_owner, amat.index_bound(0), amat.index_bound(1) ) ;

   for( size_t i=0 ; i<result->nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<result->nb_cols() ; ++j )
      {
         result->set_item( i, j, amat( i, j ) ) ;
      }
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == "LA_DenseMatrix" ) ;
   PEL_CHECK_POST( !result->is_symmetric() ) ;
   PEL_CHECK_POST( !result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->nb_rows() == amat.index_bound(0) ) ;
   PEL_CHECK_POST( result->nb_cols() == amat.index_bound(1) ) ;
   PEL_CHECK_POST( result->nb_stored_items() ==
                             amat.index_bound(0)*amat.index_bound(1) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                      FORALL( ( size_t j=0 ; j<result->nb_cols() ; ++j ),
                         result->item(i,j) == amat(i,j) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix:: create_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: create_matrix" ) ;
   PEL_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;

   LA_DenseMatrix* result = new LA_DenseMatrix( a_owner, NB_ROWS, NB_COLS ) ;

   PEL_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   int a_nb_rows = 0 ;
   int a_nb_cols = 0 ;
   if( exp->has_entry( "nb_rows" ) )
   {
      a_nb_rows = exp->int_data( "nb_rows" ) ;
      a_nb_cols = exp->int_data( "nb_cols" ) ;
      if( a_nb_rows <= 0 )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "nb_rows", "A strickly positive value is expected" ) ;
      }
      a_nb_cols = exp->int_data( "nb_cols" ) ;
      if( a_nb_cols <= 0 )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "nb_cols", "A strickly positive value is expected" ) ;
      }
   }

   LA_DenseMatrix* result =
      LA_DenseMatrix::create( a_owner, a_nb_rows, a_nb_cols ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}
//----------------------------------------------------------------------
LA_DenseMatrix:: LA_DenseMatrix( void )
//----------------------------------------------------------------------
   : LA_SeqMatrix( "LA_DenseMatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , SIZE( 0 )
   , MAT( 0 )
{
   PEL_LABEL( "LA_DenseMatrix:: LA_DenseMatrix( prototype )" ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == "LA_DenseMatrix" ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( nb_rows() == 0 ) ;
   PEL_CHECK_POST( nb_cols() == 0 ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrix:: LA_DenseMatrix( PEL_Object* a_owner,
                                 size_t a_nb_rows,
                                 size_t a_nb_cols )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_DenseMatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , SIZE( 0 )
   , MAT( 0 )
{
   PEL_LABEL( "LA_DenseMatrix:: LA_DenseMatrix" ) ;

   re_initialize_with_global_sizes( a_nb_rows, a_nb_cols ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( name() == "LA_DenseMatrix" ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( nb_rows() == a_nb_rows ) ;
   PEL_CHECK_POST( nb_cols() == a_nb_cols ) ;
   PEL_CHECK_POST( nb_stored_items() == a_nb_rows*a_nb_cols ) ;
   PEL_CHECK_POST( !is_desynchronizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::NoDistribution ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrix:: ~LA_DenseMatrix( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: ~LA_DenseMatrix" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( MAT!=0 )
   {
      delete [] MAT[0] ;
      MAT[0] = 0 ;
      delete [] MAT ;
      MAT = 0 ;
   }
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: re_initialize_with_global_sizes( size_t a_nb_rows,
                                                  size_t a_nb_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: re_initialize_with_global_sizes" ) ;
   PEL_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows,
                                                       a_nb_cols ) ) ;

   if( NB_ROWS!=a_nb_rows || NB_COLS!=a_nb_cols )
   {
      if( MAT!=0 )
      {
         delete [] MAT[0] ;
         MAT[0] = 0 ;
         delete [] MAT ;
         MAT = 0 ;
      }
      NB_ROWS = a_nb_rows ;
      NB_COLS = a_nb_cols ;
      SIZE = a_nb_rows*a_nb_cols ;

      if( SIZE > 0 )
      {
         MAT = new double* [NB_ROWS] ;
         PEL_CHECK( MAT != 0 ) ;
         MAT[0] = new double [SIZE] ;
         PEL_CHECK( MAT[0] != 0 ) ;
         for( size_t i=1 ; i<NB_ROWS ; i++ )
         {
            MAT[i] = MAT[i-1]+NB_COLS ;
         }
      }
   }
   if( a_nb_rows>0 )
   {
      double* ptr = MAT[0] ;
      for( size_t i=0 ; i<SIZE ; i++ ) ptr[i] = 0.0 ;
   }

   PEL_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows,
                                                         a_nb_cols ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           FORALL( ( size_t j=0 ; j<nb_cols() ; ++j ),
                                   item(i,j) == 0.0 ) ) ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: allocated_memory" ) ;
   return( nb_stored_items()*sizeof(double) ) ;
}

//----------------------------------------------------------------------
bool
LA_DenseMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: is_desynchronizable" ) ;

   bool result = false ;

   PEL_CHECK_POST( result == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: nb_rows" ) ;
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: nb_cols" ) ;
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
double
LA_DenseMatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: item" ) ;
   PEL_CHECK_PRE( item_PRE( i, j ) ) ;
   return( MAT[i][j] ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: nb_stored_items" ) ;
   return( SIZE ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_DenseMatrix:: create_stored_item_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: create_stored_item_iterator" ) ;
   PEL_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;

   LA_DenseMatrixIterator* result =
                        LA_DenseMatrixIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: set_stored_items( double val  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: set_stored_items" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( SIZE != 0 )
   {
      double* deb = MAT[0] ;
      for( size_t i=0 ; i<SIZE ; ++i ) deb[i] = val ;
   }
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_stored_items_POST(val) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( set_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   MAT[i][j] = x  ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   MAT[i][j] += x  ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_column( size_t j,
                             LA_SeqVector const* column,
                             double alpha,
                             double beta )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_column" ) ;
   PEL_CHECK_PRE( j<nb_cols() ) ;
   PEL_CHECK_PRE( column->nb_rows()==nb_rows() ) ;
   for( size_t i=0 ; i<nb_rows() ; i++ )
   {
      MAT[i][j] = alpha*column->item(i) + beta*MAT[i][j]  ;
   }
   PEL_CHECK_POST( FORMAL(
       FORALL( (size_t i=0;i<nb_rows() ;i++),
               item(i,j)==alpha*column->item(i)+beta*OLD(item(i,j))))) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: set" ) ;
   PEL_CHECK_PRE(set_PRE( A, false )  ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* dA = dynamic_cast<LA_DenseMatrix const* >( A ) ;
   if( dA==0 )
   {
      LA_SeqMatrix::set( A, same_pattern ) ;
   }
   else
   {
      if( SIZE!= 0 && (A!=this) )
      {
         double* dest = MAT[0] ;
         double* src = dA->MAT[0] ;

         for( size_t i=0 ; i<SIZE ; ++i )
            dest[i] = src[i] ;
      }
      synchronize() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST(A) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat( LA_Matrix const* A, double alpha,
                          bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   if( AF != 0 )
   {
      add_Mat_IMP( AF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_Mat(A,alpha,same_pattern) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat( LA_Matrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   if( AF != 0 )
   {
      add_tMat_IMP( AF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_tMat( A, alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                              double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   LA_DenseMatrix const* BF = dynamic_cast<LA_DenseMatrix const*>(B) ;

   if( AF!=0 && BF!=0 )
   {
      add_Mat_Mat_IMP( AF, BF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_Mat_Mat( A, B, alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                               double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat_Mat" ) ;
   PEL_CHECK_PRE( add_tMat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   LA_DenseMatrix const* BF = dynamic_cast<LA_DenseMatrix const*>(B) ;
   if( AF!=0 && BF!=0 )
   {
      add_tMat_Mat_IMP( AF, BF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_tMat_Mat( A, B, alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_tMat_Mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                               double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat_tMat" ) ;
   PEL_CHECK_PRE( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   LA_DenseMatrix const* BF = dynamic_cast<LA_DenseMatrix const*>(B) ;

   if( AF!=0 && BF!=0 )
   {
      add_Mat_tMat_IMP( AF, BF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_Mat_tMat( A, B, alpha ) ;
   }


   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                                double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat_tMat" ) ;
   PEL_CHECK_PRE( add_tMat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DenseMatrix const* AF = dynamic_cast<LA_DenseMatrix const*>(A) ;
   LA_DenseMatrix const* BF = dynamic_cast<LA_DenseMatrix const*>(B) ;
   if( AF!=0 && BF!=0 )
   {
      add_tMat_tMat_IMP( AF, BF, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_tMat_tMat( A, B, alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_tMat_tMat_POST() ) ;
}
//----------------------------------------------------------------------
double
LA_DenseMatrix:: determinant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: determinant" ) ;
   PEL_CHECK_PRE( nb_rows() == nb_cols() ) ;
   PEL_CHECK_PRE( nb_rows() > 0 ) ;
   double result ;

   if( NB_ROWS == 1 )
   {
      result = MAT[0][0] ;
   }
   else if( NB_ROWS == 2 )
   {
      result = MAT[0][0]*MAT[1][1]
          - MAT[0][1]*MAT[1][0] ;
   }
   else if( NB_ROWS == 3 )
   {
      result = MAT[0][0]*(MAT[1][1]*MAT[2][2]-MAT[1][2]*MAT[2][1])
              -MAT[0][1]*(MAT[1][0]*MAT[2][2]-MAT[1][2]*MAT[2][0])
              +MAT[0][2]*(MAT[1][0]*MAT[2][1]-MAT[1][1]*MAT[2][0]) ;
   }
   else
   {
      LA_DenseMatrix* clone = create_matrix( 0 ) ;
      clone->set( this ) ;
      clone->invert( result ) ;
      clone->destroy() ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DenseMatrix:: col_mean( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: col_mean" ) ;
   PEL_CHECK_PRE( nb_rows() >0 ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_cols() ) ;

   LA_SeqVector* tempVect = LA_SeqVector::create( 0, nb_rows() ) ;
   for( size_t i=0 ; i<nb_cols() ; i++ )
   {
      extract_col( i, tempVect ) ;
      result->set_item( i, tempVect->mean() ) ;
   }

   tempVect->destroy() ;
   PEL_CHECK_POST( a_owner == result->owner() ) ;
   PEL_CHECK_POST( nb_cols() == result->nb_rows() ) ;
   return result ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DenseMatrix:: col_standard_deviation( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: col_standard_deviation" ) ;
   PEL_CHECK_PRE( nb_rows()>0 ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_cols() ) ;
   LA_SeqVector* tempVect = LA_SeqVector::create( 0, nb_rows() ) ;
   for( size_t i=0 ; i<nb_cols() ; i++ )
   {
      extract_col( i, tempVect ) ;
      result->set_item( i, tempVect->standard_deviation() ) ;
   }

   tempVect->destroy() ;
   PEL_CHECK_POST( a_owner == result->owner() ) ;
   PEL_CHECK_POST( nb_cols() == result->nb_rows() ) ;
   return result ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix*
LA_DenseMatrix:: col_variance( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: col_variance" ) ;
   PEL_CHECK_PRE( nb_rows()>1 ) ;
   LA_SymmetricMatrix* result = LA_SymmetricMatrix::create( a_owner, nb_cols() ) ;
   const size_t n = nb_rows() ;
   const size_t p = nb_cols() ;
   LA_SeqVector* means = col_mean( 0 ) ;
   double rn = double( n ) ;
   for( size_t j=0 ; j<p ; j++ )
   {
      for( size_t i=j ; i<p ; i++ )
      {
         double x = 0.0 ;
         for( size_t is=0 ; is < n ; is++ )
         {
            x +=  ( MAT[is][i]-means->item( i ) )
                * ( MAT[is][j]-means->item( j ) ) ;
         }
         x /= rn ;
         result->set_item( j, i, x ) ;
      }
   }
   means->destroy() ;
   PEL_CHECK_POST( result->nb_rows() == nb_cols() ) ;
   PEL_CHECK_POST( result->nb_cols() == nb_cols() ) ;
   return result ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix*
LA_DenseMatrix:: col_correlation( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: col_correlation" ) ;
   PEL_CHECK_PRE( nb_rows() > 1 ) ;
   const size_t p = nb_cols() ;

   LA_SeqVector* stdev = col_standard_deviation( 0 ) ;
   LA_SymmetricMatrix* result = col_variance( a_owner ) ;
   for( size_t j=0 ; j<p ; j++ )
   {
      for( size_t i=j ; i < p ; i++ )
      {
         double x = result->item( j, i )/stdev->item( i )
                                         /stdev->item( j ) ;
         result->set_item( j, i, x ) ;
      }
   }
   stdev->destroy() ;
   PEL_CHECK_POST( result->nb_rows() == nb_cols() ) ;
   PEL_CHECK_POST( result->nb_cols() == nb_cols() ) ;
   return result ;
}

//----------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix:: col_reduced( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: col_reduced" ) ;
   PEL_CHECK_PRE( nb_rows()>0 ) ;
   const size_t p = nb_cols() ;
   const size_t n = nb_rows() ;
   LA_DenseMatrix* result = LA_DenseMatrix::create( a_owner, n, p ) ;

   LA_SeqVector* means = col_mean( 0 ) ;
   LA_SeqVector* stdev = col_standard_deviation( 0 ) ;
   for( size_t i=0 ; i < n ; i++ )
   {
      for( size_t j=0 ; j < p ; j++ )
      {
         double x = ( MAT[i][j]-means->item( j ) )
                                           / stdev->item( j ) ;
         result->set_item( i, j, x ) ;
      }
   }
   means->destroy() ;
   stdev->destroy() ;
   PEL_CHECK_POST( result->nb_rows() == nb_rows() ) ;
   PEL_CHECK_POST( result->nb_cols() == nb_cols() ) ;
   return result ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: invert( double& det )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: invert" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( nb_rows()==nb_cols() ) ;
   PEL_CHECK_PRE( nb_rows()>0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   det = 1.0 ;
   size_t_vector ik( nb_rows() ) ;
   size_t_vector jk( nb_rows() ) ;

   size_t k ;
   for( k=0 ; k < nb_rows() ; k++ )
   {
      //   find largest element in rest of matrix
      double max = 0.0 ;
      if( !largelem( max, k, ik, jk ) )
      {
         det *= max ;
      }
      else
      {
         det = 0.0 ;
         break ;
      }
   }
   if( det != 0.0 )
   {
      // restore ordering of matrix
      for( size_t l = 0 ; l < nb_rows() ; l++ )
      {
         k = nb_rows() - l -1 ;
         size_t j = ik( k ) ;
         if( j > k )
         {
            for( size_t i = 0 ; i < nb_rows() ; i++ )
            {
               double sav = item( i, k ) ;
               set_item( i, k, -item( i, j ) ) ;
               set_item( i, j, sav ) ;
            }
         }
         size_t i = jk( k ) ;
         if( i > k )
         {
            for( j = 0 ; j < nb_rows() ; j++ )
            {
               double sav = item( k, j ) ;
               set_item( k, j, -item( i, j ) ) ;
               set_item( i, j, sav ) ;
            }
         }
      }
   }
   else
   {
      re_initialize_with_global_sizes( 0, 0 ) ;
   }
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
   PEL_CHECK_POST( IMPLIES( det==0, nb_rows()==0 && nb_cols()==0 ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat_IMP( LA_DenseMatrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat_IMP" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         MAT[i][j] += alpha * A->MAT[i][j] ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat_IMP( LA_DenseMatrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK( add_tMat_PRE( A, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         MAT[i][j] += alpha * A->MAT[j][i] ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat_Mat_IMP( LA_DenseMatrix const* A,
                                  LA_DenseMatrix const* B,
                                  double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat_Mat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const kMax = A->nb_cols() ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         for( size_t k=0 ; k<kMax ; ++k )
         {
            MAT[i][j] += alpha * A->MAT[i][k] * B->MAT[k][j]  ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat_Mat_IMP( LA_DenseMatrix const* A,
                                   LA_DenseMatrix const* B,
                                   double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat_Mat_IMP" ) ;
   PEL_CHECK( add_tMat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const kMax = A->nb_rows() ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         for( size_t k=0 ; k<kMax ; ++k )
         {
            MAT[i][j] += alpha * A->MAT[k][i] * B->MAT[k][j]  ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( add_tMat_Mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_Mat_tMat_IMP( LA_DenseMatrix const* A,
                                   LA_DenseMatrix const* B,
                                   double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_Mat_tMat_IMP" ) ;
   PEL_CHECK( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const kMax = A->nb_cols() ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         for( size_t k=0 ; k<kMax ; ++k )
         {
            MAT[i][j] += alpha * A->MAT[i][k] * B->MAT[j][k]  ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( add_Mat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrix:: add_tMat_tMat_IMP( LA_DenseMatrix const* A,
                                    LA_DenseMatrix const* B,
                                    double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: add_tMat_tMat_IMP" ) ;
   PEL_CHECK( add_tMat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const kMax = A->nb_rows() ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      for( size_t j=0 ; j<nb_cols() ; ++j )
      {
         for( size_t k=0 ; k<kMax ; ++k )
         {
            MAT[i][j] += alpha * A->MAT[k][i] * B->MAT[j][k] ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( add_tMat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
bool
LA_DenseMatrix:: largelem( double &mxel,
                           size_t k,
                           size_t_vector& ik,
                           size_t_vector& jk )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix:: largelem" ) ;
   PEL_CHECK( k<nb_rows() ) ;
   PEL_CHECK( nb_cols() == nb_rows() ) ;
   PEL_CHECK( invariant() ) ;

   for( size_t i = k ; i < nb_rows() ; i++ )
   {
      for( size_t j = k ; j < nb_rows() ; j++ )
      {
         double xij = item( i, j ) ;
         if( PEL::abs( mxel ) <= PEL::abs( xij ) )
         {
            mxel = xij ;
            ik( k ) = i ;
            jk( k ) = j ;
         }
      }
   }

   //  interchange rows and columns to put max in position ( k, k )
   if( PEL::abs( mxel ) > 0.0 )
   {
      size_t i = ik(k) ;
      if( i < k )
      {
         return( largelem( mxel, k, ik, jk ) ) ;
      }
      if( i > k )
      {
         for( size_t j=0 ; j < nb_rows() ; j++ )
         {
            double sav = item( k, j ) ;
            set_item( k, j, item( i, j ) ) ;
            set_item( i, j, -sav ) ;
         }
      }
      size_t j = jk( k ) ;
      if( j < k )
      {
         return( largelem( mxel, k, ik, jk ) ) ;
      }
      if( j > k )
      {
         for( i = 0 ; i < nb_rows() ; i++ )
         {
            double sav = item( i, k ) ;
            set_item( i, k, item( i, j ) ) ;
            set_item( i, j, -sav ) ;
         }
      }
      // accumulate elements of inverse matrix
      for( i=0 ; i < nb_rows() ; i++ )
      {
         if( i != k )
         {
            set_item( i, k,
                        item( i, k ) / ( -mxel ) ) ;
         }
      }
      for( i = 0 ; i < nb_rows() ; i++ )
      {
         for( j=0 ; j < nb_rows() ; j++ )
         {
            if( i != k && j != k )
            {
               add_to_item( i, j,
                           item( i, k ) * item( k, j ) ) ;
            }
         }
      }
      for( j = 0 ; j < nb_rows() ; j++ )
      {
         if( j != k )
         {
            set_item( k, j, item( k, j ) / mxel ) ;
         }
      }
      set_item( k, k, 1.0/mxel ) ;


      return( false ) ;
   }
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_DenseMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( ( MAT==0 && NB_ROWS==0 && NB_COLS==0 ) ||
               ( MAT!=0 && MAT[0]!=0 && NB_ROWS>0 && NB_COLS>0 ) ) ;
   PEL_ASSERT( SIZE == NB_ROWS*NB_COLS ) ;
   return true ;
}
