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

#include <LA_CRSmatrix.hh>

#include <fstream>
#include <memory.h>

#include <LA_CRSmatrixIterator.hh>
#include <LA_PelMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <sstream>

struct LA_CRSmatrix_ERROR
{
   static void n0( void ) ;
} ;

LA_CRSmatrix const* LA_CRSmatrix::PROTOTYPE = new LA_CRSmatrix() ;

//----------------------------------------------------------------------
LA_CRSmatrix*
LA_CRSmatrix:: create( PEL_Object* a_owner, LA_SeqMatrix const* other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: create" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   PEL_CHECK_PRE( other->is_synchronized() ) ;

   LA_CRSmatrix* result =  new LA_CRSmatrix( a_owner, other ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( ! result->is_symmetric() ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( ! result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;
   PEL_CHECK_POST( result->nb_rows() == other->nb_rows() ) ;
   PEL_CHECK_POST( result->nb_cols() == other->nb_cols() ) ;
   PEL_CHECK_POST( result->nb_stored_items() == other->nb_stored_items() ) ;
   PEL_CHECK_POST( ! result->insertion_mode() ) ;
   PEL_CHECK_POST( ! result->is_sending_in_place() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix*
LA_CRSmatrix:: create_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: create_matrix" ) ;
   PEL_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_CRSmatrix* result = new LA_CRSmatrix( a_owner, NB_ROWS, NB_COLS ) ;
   result->set_insertion_mode( INSERT ) ;

   PEL_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( !result->is_sending_in_place() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix*
LA_CRSmatrix:: create_copy( PEL_Object* a_owner,
                            LA_SeqMatrix const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: create_copy" ) ;
   PEL_CHECK_PRE( create_copy_PRE( a_owner, other ) ) ;

   LA_CRSmatrix* result = new LA_CRSmatrix( a_owner, other ) ;

   PEL_CHECK_POST( create_copy_POST( result, a_owner, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix*
LA_CRSmatrix:: create_replica( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   size_t a_nb_rows = 0 ;
   size_t a_nb_cols = 0 ;
   if( exp->has_entry( "nb_rows" ) ||
       exp->has_entry( "nb_cols" ) )
   {
      a_nb_rows = (size_t) exp->int_data( "nb_rows" ) ;
      a_nb_cols = (size_t) exp->int_data( "nb_cols" ) ;
      exp->test_data( "nb_rows", "nb_rows>=0" ) ;
      exp->test_data( "nb_cols", "nb_cols>=0" ) ;
   }

   LA_CRSmatrix* result = new LA_CRSmatrix( a_owner, a_nb_rows, a_nb_cols ) ;

   if( exp->has_entry( "insertion_mode" ) )
   {
      result->set_insertion_mode( exp->bool_data( "insertion_mode" ) ) ;
      exp->set_default( "insertion_mode", "false" ) ;
   }

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   PEL_CHECK( result->distribution_strategy() == LA::NoDistribution ) ;
   PEL_CHECK( !result->is_sending_in_place() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix:: LA_CRSmatrix( PEL_Object* a_owner,
                             size_t a_nb_rows,
                             size_t a_nb_cols )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_CRSmatrix" )
   , NB_ROWS( a_nb_rows )
   , NB_COLS( a_nb_cols )
   , VALUES( 0 )
   , DIAG( a_nb_rows, PEL::bad_index() )
   , COL( 0 )
   , START( a_nb_rows+1 )
   , NB_ELEMS(0)
   , INSERT( false )
{
   PEL_LABEL( "LA_CRSmatrix:: LA_CRSmatrix( a_nb_rows, a_nb_cols )" ) ;

   for( size_t i=0 ; i<3 ; ++i ) DIMS[i] = 0 ;
   for( size_t i=0 ; i<4 ; ++i ) REQUEST[i] = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == "LA_CRSmatrix" ) ;
   PEL_CHECK_POST( nb_rows() == a_nb_rows ) ;
   PEL_CHECK_POST( nb_cols() == a_nb_cols ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix:: LA_CRSmatrix( PEL_Object* a_owner,
                             LA_SeqMatrix const* other )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_CRSmatrix" )
   , NB_ROWS( other->nb_rows() )
   , NB_COLS( other->nb_cols() )
   , VALUES( 0 )
   , DIAG( other->nb_rows(), PEL::bad_index()  )
   , COL( 0 )
   , START( other->nb_rows() +1 )
   , NB_ELEMS( 0 )
   , INSERT( false )
{
   PEL_LABEL( "LA_CRSmatrix:: LA_CRSmatrix( LA_SeqMatrix )" ) ;

   set( other, false ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == "LA_CRSmatrix" ) ;
   PEL_CHECK_POST( nb_rows() == other->nb_rows() ) ;
   PEL_CHECK_POST( nb_cols() == other->nb_cols() ) ;
   PEL_CHECK_POST( nb_stored_items() == other->nb_stored_items() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
   PEL_CHECK_POST( !insertion_mode() ) ;
   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix:: LA_CRSmatrix( PEL_Object* a_owner,
                             PEL_Communicator const* com,
                             size_t src )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_CRSmatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , VALUES( 0 )
   , DIAG( 0 )
   , COL( 0 )
   , START( 0 )
   , NB_ELEMS(0)
   , INSERT( false )
{
   PEL_LABEL( "LA_CRSmatrix::LA_CRSmatrix( PEL_Communicator ) " ) ;
   PEL_CHECK_PRE( com!=0 ) ;
   PEL_CHECK_PRE( src < com->nb_ranks() && src != com->rank() ) ;

   int dims[3] ;
   com->receive( src, dims, 3 ) ;
   NB_ROWS = dims[0] ;
   NB_COLS = dims[1] ;
   NB_ELEMS = dims[2] ;
   START.re_initialize( NB_ROWS+1 ) ;
   COL.re_initialize( NB_ELEMS ) ;
   VALUES.re_initialize( NB_ELEMS ) ;
   DIAG.re_initialize( NB_ROWS, PEL::bad_index() ) ;

   com->receive( src, const_cast<int*>( START.data() ), NB_ROWS+1 ) ;
   if( NB_ELEMS>0 )
   {
      com->receive( src, const_cast<int*>( COL.data() ), NB_ELEMS ) ;
      com->receive( src, const_cast<double*>( VALUES.data() ), NB_ELEMS ) ;
   }
   size_t i1 = START(0) ;
   for( size_t i=0 ; i<NB_ROWS ; i++ )
   {
      size_t i2 = START(i+1) ;
      for( size_t j=i1 ; j<i2 ; j++ )
      {
         if( COL(j)==(int)i )
         {
            DIAG(i)=j ;
            break ;
         }
      }
      i1 = i2 ;
   }
   synchronize() ;

   for( size_t i=0 ; i<3 ; ++i ) DIMS[i] = 0 ;
   for( size_t i=0 ; i<4 ; ++i ) REQUEST[i] = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == "LA_CRSmatrix" ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix:: LA_CRSmatrix( void )
//----------------------------------------------------------------------
   : LA_SeqMatrix( "LA_CRSmatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , VALUES( 0 )
   , DIAG( 0 )
   , COL( 0 )
   , START( 0 )
   , INSERT( false )
{
   PEL_LABEL( "LA_CRSmatrix:: LA_CRSmatrix( prototype )" ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == "LA_CRSmatrix" ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix:: ~LA_CRSmatrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: re_initialize_with_global_sizes( size_t a_nb_rows,
                                                size_t a_nb_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: re_initialize_with_global_sizes" ) ;
   PEL_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows,
                                                       a_nb_cols ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   NB_ROWS = a_nb_rows ;
   NB_COLS = a_nb_cols ;
   NB_ELEMS = 0 ;
   VALUES.re_initialize( 0 ) ;
   DIAG.re_initialize( NB_ROWS, PEL::bad_index() ) ;
   COL.re_initialize( 0 ) ;
   START.re_initialize( NB_ROWS+1 ) ;
   START.set( 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows,
                                                         a_nb_cols ) ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   return( NB_ELEMS ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   size_t result =
      NB_ELEMS*( sizeof(double) + sizeof(int) ) + (NB_ROWS+1)*sizeof(int) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_CRSmatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: is_desynchronizable" ) ;
   
   bool result = false ;
   
   PEL_CHECK_POST( result == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
double
LA_CRSmatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: item" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( item_PRE( i, j ) ) ;

   double result = 0.0 ;

   size_t const i1 = START(i) ;
   size_t const i2 = START(i+1) ;
   for( size_t ii=i1 ; ii<i2 ; ++ii )
   {
      if( COL(ii) >= (int)j )
      {
         if( COL(ii) == (int)j ) result = VALUES(ii) ;
         break ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: extract_diag( LA_Vector* diag ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: extract_diag" ) ;
   PEL_CHECK_PRE( extract_diag_PRE( diag ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector* >( diag ) != 0 ) ;
   LA_SeqVector* d = static_cast<LA_SeqVector* >( diag ) ;

   double* ptr_d = d->data() ;
   for( size_t i=0 ; i<NB_ROWS ; ++i )
   {
      if( DIAG(i)==PEL::bad_index() )
      {
         ptr_d[i] = 0. ;
      }
      else
      {
         ptr_d[i] = VALUES( DIAG(i) ) ;
      }
   }
   d->synchronize() ;

   PEL_CHECK_POST( extract_diag_POST( diag ) ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_CRSmatrix:: create_stored_item_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: create_stored_item_iterator" ) ;
   PEL_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_CRSmatrixIterator* result =
      LA_CRSmatrixIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_CRSmatrix:: insertion_mode( void ) const
//----------------------------------------------------------------------
{
   return( INSERT ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: set_insertion_mode( bool allowed )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: set_insertion_mode" ) ;
   INSERT = allowed ;
   PEL_CHECK_POST( insertion_mode()==allowed ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: set_stored_items( double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: set_stored_items" ) ;
   PEL_CHECK_INV( invariant() ) ;

   VALUES.set( val ) ;
   if( is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_undef ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_stored_items_POST( val ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: nullify_row( size_t i )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: nullify_row" ) ;
   PEL_CHECK_PRE( nullify_row_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const i1 = START(i) ;
   size_t const i2 = START(i+1) ;
   for( size_t ii=i1 ; ii<i2 ; ++ii )
   {
      VALUES(ii) = 0.0 ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nullify_row_POST( i ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( set_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool done = false ;
   size_t const i1 = START(i) ;
   size_t const i2 = START(i+1) ;
   for( size_t ii=i1 ; ii<i2 ; ++ii )
   {
      if( COL(ii) >= (int)j )
      {
         if( COL(ii) == (int)j )
         {
            VALUES(ii)=x ;
         }
         else
         {
            insert( i, j, x, ii ) ;
         }
         done = true ;
         break ;
      }
   }
   if( !done ) insert( i, j, x, i2 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: insert( size_t i, size_t j, double x, size_t pos )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: insert" ) ;
   PEL_CHECK_PRE( i < nb_rows() ) ;
   PEL_CHECK_PRE( j < nb_cols() ) ;
   PEL_CHECK_PRE( pos <= NB_ELEMS ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( !INSERT )
   {
      LA_CRSmatrix_ERROR::n0() ;
   }

   VALUES.resize(NB_ELEMS+1) ;
   COL.resize(NB_ELEMS+1) ;
   for( size_t k=i+1 ; k<NB_ROWS+1 ; ++k )
   {
      START(k)++ ;
   }
   for( size_t k=i+1 ; k<NB_ROWS ; ++k )
   {
      if( DIAG(k) != PEL::bad_index() ) DIAG(k)++ ;
   }
   for( int jj=(int)NB_ELEMS ; jj>(int)pos ; jj-- )
   {
      VALUES(jj) = VALUES(jj-1) ;
      COL(jj) = COL(jj-1) ;
   }
   VALUES(pos) = x ;
   COL(pos) = j ;

   if( j<i )
   {
      if(  DIAG(i) != PEL::bad_index() ) DIAG(i)++ ;
   }
   else if( i==j )
   {
      DIAG(i) = pos ;
   }

   NB_ELEMS++ ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const i1 = START(i) ;
   size_t const i2 = START(i+1) ;
   bool done = false ;
   for( size_t ii=i1 ; ii<i2 ; ++ii )
   {
      if( COL(ii)>= (int)j )
      {
         if( COL(ii)==(int)j )
         {
            VALUES(ii) += x ;
         }
         else
         {
            insert( i, j, x, ii ) ;
         }
         done = true ;
         break ;
      }
   }
   if( !done ) insert( i, j, x, i2 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_CRSmatrix const* AS = dynamic_cast<LA_CRSmatrix const*>(A) ;
   if( AS==0 )
   {
      size_t N = A->nb_stored_items() ;
      VALUES.resize( N ) ;
      COL.resize( N ) ;
      DIAG.set( PEL::bad_index() ) ;

      LA_MatrixIterator* it = A->create_stored_item_iterator( 0 ) ;
      NB_ELEMS = 0 ;
      for( size_t i=0 ; i<NB_ROWS ; i++ )
      {
         START( i ) = NB_ELEMS ;
         for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
         {
            VALUES( NB_ELEMS ) = it->item() ;
            COL( NB_ELEMS ) = it->col() ;
            if( it->col() == i ) DIAG(i) = NB_ELEMS ;
            NB_ELEMS++ ;
         }
      }
      START( NB_ROWS ) = NB_ELEMS ;
      it->destroy() ; it = 0 ;
      for( size_t i=0 ; i<NB_ROWS ; i++ )
      {
         size_t const i1 = START(i) ;
         size_t const i2 = START(i+1) ;
         bool sorted = true ;
         for( size_t ii=i1+1 ; ii<i2 ; ++ii )
         {
            sorted &= ( COL(ii)>COL(ii-1) ) ;
         }
         if( !sorted )
         {
            for( size_t ii=i1 ; ii<i2 ; ++ii )
            {
               size_t min_idx = ii ;
               for( size_t jj=ii+1 ; jj<i2 ; ++jj )
               {
                  if( COL(jj)<COL(min_idx) ) min_idx = jj ;
               }
               if( min_idx != ii )
               {
                  size_t col = COL(ii) ;
                  double x = VALUES(ii) ;
                  COL(ii) = COL(min_idx) ;
                  VALUES(ii) = VALUES(min_idx) ;
                  COL(min_idx) = col ;
                  VALUES(min_idx) = x ;
                  if( COL(ii)==(int)i ) DIAG(i) = ii ;
                  if( COL(min_idx)==(int)i ) DIAG(i) = min_idx ;
               }
            }
         }
      }
   }
   else
   {
      NB_ROWS = AS->nb_rows() ;
      NB_COLS = AS->nb_cols() ;
      VALUES = AS->VALUES ;
      COL = AS->COL ;
      START = AS->START ;
      NB_ELEMS = AS->NB_ELEMS ;
      DIAG = AS->DIAG ;
      INSERT = AS->INSERT ;
   }
   for( size_t i=0 ; i<3 ; ++i ) DIMS[i] = 0 ;
   for( size_t i=0 ; i<4 ; ++i ) REQUEST[i] = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( A ) ) ;
   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
				      double alpha, double beta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector* >( y ) != 0 ) ;
   LA_SeqVector* by = static_cast<LA_SeqVector* >( y ) ;

   double const* in = bx->data() ;
   double* resu = by->data() ;

   int const * col = COL.data() ;
   double const * v = VALUES.data() ;

   int i2 = START(0) ;
   if( beta==1.0 )
   {
      for( size_t i=0 ; i<NB_ROWS ; ++i )
      {
         int i1 = i2 ;
         i2 = START(i+1) ;
         if( i2!=i1 )
         {
            double sum = 0.0 ;
            for( int j=i1 ; j<i2 ; j++ )
               sum += v[j] * in[ col[j] ] ;

            resu[i] += alpha * sum ;
         }
      }
   }
   else if( beta==0.0 )
   {
      for( size_t i=0 ; i<NB_ROWS ; ++i )
      {
         int i1 = i2 ;
         i2 = START(i+1) ;
         double sum = 0.0 ;
         for( int j=i1 ; j<i2 ; j++ )
            sum += v[j] * in[ col[j] ] ;

         resu[i] = alpha * sum ;
      }
   }
   else
   {
      for( size_t i=0 ; i<NB_ROWS ; ++i )
      {
         int i1 = i2 ;
         i2 = START(i+1) ;
         double sum = 0.0 ;
         for( int j=i1 ; j<i2 ; j++ )
            sum += v[j] * in[ col[j] ] ;

         resu[i] = beta*resu[i] + alpha * sum ;

      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( alpha != 1.0 )
   {
      for( size_t i=0 ; i<NB_ELEMS ; ++i )
      {
         VALUES(i) *= alpha ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: send( PEL_Communicator const* com, size_t dest,
                     bool in_place ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: send" ) ;
   PEL_CHECK_PRE( com != 0 ) ;
   PEL_CHECK_PRE( dest < com->nb_ranks() && dest!=com->rank() ) ;
   PEL_CHECK_INV( invariant() ) ;

   DIMS[0] = NB_ROWS ;
   DIMS[1] = NB_COLS ;
   DIMS[2] = NB_ELEMS ;
   REQUEST[0] = 0 ;

   if( in_place )
   {
      REQUEST[0] = com->Isend( dest, &(DIMS[0]), 3 ) ;
      REQUEST[1] = com->Isend( dest, START.data(), NB_ROWS+1 ) ;
      if( NB_ELEMS>0 )
      {
         REQUEST[2] = com->Isend( dest, COL.data(), NB_ELEMS ) ;
         REQUEST[3] = com->Isend( dest, VALUES.data(), NB_ELEMS ) ;
      }
   }
   else
   {
      com->send( dest, &(DIMS[0]), 3 ) ;
      com->send( dest, START.data(), NB_ROWS+1 ) ;
      if( NB_ELEMS>0 )
      {
         com->send( dest, COL.data(), NB_ELEMS ) ;
         com->send( dest, VALUES.data(), NB_ELEMS ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( EQUIVALENT( in_place, is_sending_in_place() ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: wait_send( PEL_Communicator const* com ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: wait_send" ) ;
   PEL_CHECK_PRE( is_sending_in_place() ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<4 ; i++ )
   {
      com->wait( REQUEST[i] ) ;
      REQUEST[i] = 0 ;
   }

   PEL_CHECK_POST( !is_sending_in_place() ) ;
}

//----------------------------------------------------------------------
bool
LA_CRSmatrix:: is_sending_in_place( void ) const
//----------------------------------------------------------------------
{
   return( REQUEST[0]!=0 ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix*
LA_CRSmatrix:: receive( PEL_Object* a_owner,
                        PEL_Communicator const* com,
                        size_t src )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: receive" ) ;
   PEL_CHECK_PRE( com != 0 ) ;
   PEL_CHECK_PRE( src < com->nb_ranks() && src!=com->rank() ) ;

   LA_CRSmatrix* result =  new LA_CRSmatrix( a_owner, com, src ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_sending_in_place() ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
void
LA_CRSmatrix:: factorize_MILU0( bool modified, double piv_min )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: factorize_MILU0" ) ;
   PEL_CHECK_PRE( factorize_MILU0_PRE( modified, piv_min ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   int const n = (int) NB_ROWS ;

   for( int i=0 ; i<n ; ++i )
   {
      if( DIAG(i)==PEL::bad_index() )
      {
         raise_MILU0_zero_pivot( i ) ;
         if( START(i) < START(i+1) )
         {
            int const ii = START(i+1) - 1 ;
            VALUES(ii) = 1.0 ;
            COL(ii) = i ;
            DIAG(i) = ii ;
         }
         else
         {
            PEL_Error::object()->raise_plain( "Unable to decompose matrix" ) ;
         }

      }
      else
      {
         PEL_ASSERT( COL( DIAG(i) ) == i ) ;
      }
   }

   if( PEL::abs( VALUES(DIAG(0)) )<piv_min )
   {
      nullify_row( 0 ) ;
      VALUES( DIAG( 0 ) ) = 1. ;
      raise_MILU0_zero_pivot( 0 ) ;
   }

   for( int i=1 ; i<n ; ++i )
   {
      double drop = 0.0 ;
      for( int itI=START(i) ; itI<START(i+1) ; ++itI )
      {
         int k = COL(itI) ;
         if( k<i )
         {
            double akk = VALUES( DIAG( k ) ) ;

            double aik = VALUES(itI)/akk ;
            VALUES(itI) = aik ;
            for( int itK=START(k) ; itK<START(k+1) ; ++itK )
            {
               int j = COL(itK) ;
               if( j>k )
               {
                  double akj = VALUES(itK) ;
                  for( int itJ=itI+1 ; itJ<START(i+1) ; ++itJ )
                  {
                     if( COL(itJ) >= j )
                     {
                        if( COL(itJ)==j )
                        {
                           VALUES(itJ) -= aik*akj ;
                        }
                        else
                        {
                           drop -= aik*akj ;
                        }
                        itJ = START(i+1) ;
                     }
                  }
               }
            }
         }
      }
      double aii = VALUES( DIAG( i ) ) ;
      if( modified )
      {
         aii -= drop ;
         VALUES( DIAG( i ) ) = aii ;
      }
      if( PEL::abs( aii )<piv_min )
      {
         nullify_row( i ) ;
         VALUES( DIAG( i ) ) = 1. ;
         raise_MILU0_zero_pivot( i ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( factorize_MILU0_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: solve_LU( LA_SeqVector const* rhs, LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: solve_LU" ) ;
   PEL_CHECK_PRE( solve_LU_PRE( rhs, sol ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const* ptr_rhs = rhs->data() ;
   double* ptr_sol = sol->data() ;

   int const n = (int) NB_ROWS ;
   int const* col = COL.data() ;
   double const* v = VALUES.data() ;

   for( int i=0 ; i<n ; ++i )
   {
      int i1 = START(i) ;
      int i2 = (int) DIAG(i) ;
      double r = ptr_rhs[ i ] ;
      for( int j=i1 ; j<i2 ; ++j )
         r -= v[j] * ptr_sol[ col[j] ] ;
      ptr_sol[ i ] = r ;
   }
   for( int i=n-1 ; i>=0 ; --i )
   {
      int i1 = (int) DIAG(i) ;
      int i2 = START(i+1) ;
      double r = ptr_sol[ i ] ;
      for( int j=i1+1 ; j<i2 ; ++j )
         r -= v[j] * ptr_sol[ col[j] ] ;
      double d = v[i1] ;
      ptr_sol[ i ] = r/d ;
   }
   sol->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( solve_LU_POST( rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: relax( double omega,
                      LA_SeqMatrix::relaxation_mode mode,
                      LA_SeqVector const* omega_inv_diag,
                      LA_SeqVector const* rhs,
                      LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: relax" ) ;
   PEL_CHECK_PRE( relax_PRE( omega, mode, omega_inv_diag, rhs, sol ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const* ptr_rhs = rhs->data() ;
   double const* ptr_omega_inv_diag = omega_inv_diag->data() ;
   double* ptr_sol = sol->data() ;

   int const n = (int) NB_ROWS ;
   int const* start = START.data() ;
   int const* col = COL.data() ;
   double const* v = VALUES.data() ;

   if( mode == forward || mode == symmetric )
   {
      for( int i=0 ; i<n ; ++i )
      {
         double xt = ptr_rhs[i] ;
         int i1 = start[i] ;
         int i2 = start[i+1] ;
         for( int j=i1 ; j<i2 ; ++j )
            xt -= v[j] * ptr_sol[ col[j] ] ;
         ptr_sol[ i ] += xt * ptr_omega_inv_diag[i] ;
      }
   }

   if( mode == backward || mode == symmetric )
   {
      for( int i=n-1 ; i>=0 ; i-- )
      {
         double xt = ptr_rhs[i] ;
         int i1 = start[i] ;
         int i2 = start[i+1] ;
         for( int j=i1 ; j<i2 ; ++j )
            xt -= v[j] * ptr_sol[ col[j] ] ;
         ptr_sol[ i ] += xt * ptr_omega_inv_diag[i] ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( relax_POST( omega, mode, omega_inv_diag, rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: readMM( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: readMM" ) ;
   PEL_CHECK_PRE( readMM_PRE( file ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   // Avoid unefficient (or forbidden) call to insert function:
   LA_PelMatrix* m = LA_PelMatrix::create( 0, 0, 0 ) ;
   m->readMM( file ) ;
   re_initialize_with_global_sizes( m->nb_rows(), m->nb_cols() ) ;
   set( m ) ;
   m->destroy() ; m = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( readMM_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: restore( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: restore" ) ;
   PEL_CHECK_PRE( restore_PRE( exp ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   // Avoid unefficient (or forbidden) call to insert function:
   LA_PelMatrix* m = LA_PelMatrix::create( 0, 0, 0 ) ;
   m->restore( exp ) ;
   re_initialize_with_global_sizes( m->nb_rows(), m->nb_cols() ) ;
   set( m ) ;
   m->destroy() ; m = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( restore_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrix:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix:: print" ) ;

   LA_Matrix::print( os, indent_width ) ;
   std::string s( indent_width, ' ' ) ;
   os << s << "insertion_mode: " << ( INSERT ? "\"true\"" : "\"false\"" )
      << std::endl ;
}

//----------------------------------------------------------------------
bool
LA_CRSmatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( ( NB_ROWS==0 && NB_COLS==0 ) ||
               ( NB_ROWS * NB_COLS !=0 ) ) ;
   PEL_ASSERT( IMPLIES( ! is_a_prototype(),
                        START.size()==NB_ROWS+1 ) ) ;
   PEL_ASSERT( IMPLIES( ! is_a_prototype(),
                        DIAG.size()==NB_ROWS ) ) ;
   PEL_ASSERT( IMPLIES( ! is_a_prototype(),
                        COL.size()==NB_ELEMS ) ) ;
   PEL_ASSERT( IMPLIES( ! is_a_prototype(),
                        VALUES.size()==NB_ELEMS ) ) ;
//    PEL_ASSERT( IMPLIES( ! is_a_prototype(),
//                         FORALL( ( size_t i=0 ; i<NB_ROWS ; ++i ),
//                                 DIAG(i)==PEL::bad_index() ||
//                                              COL(DIAG(i))==(int)i ) ) ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void
LA_CRSmatrix_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_internal(
      "*** LA_CRSmatrix error:\n"
      "    insertion of new items for CRS matrix is not efficient\n"
      "    and is not allowed\n"
      "    (use LA_CRSmatrix::set_insertion_mode(true) to allow\n"
      "     anyway insertion of new items)" ) ;
}
