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

#include <LA_BlockSeqMatrix.hh>

#include <LA_MatrixIterator.hh>
#include <LA_BlockSeqMatrixIterator.hh>
#include <LA_PelMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_Iterator.hh>
#include <intVector.hh>

#include <iostream>

LA_BlockSeqMatrix const* 
LA_BlockSeqMatrix:: PROTOTYPE = new LA_BlockSeqMatrix() ;

//----------------------------------------------------------------------
LA_BlockSeqMatrix*
LA_BlockSeqMatrix:: create( PEL_Object* a_owner,
                            size_t_vector const& a_row_partitioning,
                            size_t_vector const& a_col_partitioning,
                            LA_SeqMatrix const* a_block_proto )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: create" ) ;
   PEL_CHECK_PRE( a_row_partitioning.size()>0 ) ;
   PEL_CHECK_PRE( a_col_partitioning.size()>0 ) ;
   PEL_CHECK_PRE( a_block_proto != 0 ) ;
   PEL_CHECK_PRE( !a_block_proto->is_desynchronizable() ) ;
   PEL_CHECK_PRE( !a_block_proto->is_symmetric() ) ;
   PEL_CHECK_PRE( a_block_proto->is_resizable() ) ;

   LA_BlockSeqMatrix* result =  new LA_BlockSeqMatrix( a_owner,
                                                       a_row_partitioning,
                                                       a_col_partitioning,
                                                       a_block_proto ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == "LA_BlockSeqMatrix" ) ;
   PEL_CHECK_POST( ! result->is_symmetric() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   PEL_CHECK_POST( ! result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->row_partition() == a_row_partitioning ) ;
   PEL_CHECK_POST( result->col_partition() == a_col_partitioning ) ;
   PEL_CHECK_POST( result->nb_stored_items() == 0 ) ;
   PEL_CHECK_POST( result->nb_rows() == a_row_partitioning.sum() ) ;
   PEL_CHECK_POST( result->nb_cols() == a_col_partitioning.sum() ) ;
   PEL_CHECK_POST( result->block_prototype()->name() == a_block_proto->name() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrix*
LA_BlockSeqMatrix:: create_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: create_matrix" ) ;
   PEL_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_BlockSeqMatrix* result = new LA_BlockSeqMatrix(
                  a_owner, ROW_PARTITION, COL_PARTITION, BLOCK_PROTO ) ;

   PEL_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrix*
LA_BlockSeqMatrix:: create_replica( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   intVector const& rp = exp->intVector_data( "row_partitioning" ) ;
   size_t const nr = rp.size() ;
   if( nr == 0 )
   {
      PEL_Error::object()->raise_data_error(
         exp, "row_partitioning",
         "invalid null row-partitioning" ) ;
   }
   size_t_vector row_partitioning( nr ) ;
   for( size_t i=0 ; i<nr ; ++i )
   {
      if( rp( i ) <= 0 )
      {
         PEL_Error::object()->raise_data_error(
            exp, "row_partitioning",
            "invalid null or negative row dimension" ) ;
      }
      row_partitioning( i ) = (size_t) rp( i ) ;
   }

   intVector const& cp = exp->intVector_data( "col_partitioning" ) ;
   size_t const nc = cp.size() ;
   if( nc == 0 )
   {
      PEL_Error::object()->raise_data_error(
         exp, "col_partitioning",
         "invalid null col-partitioning" ) ;
   }
   size_t_vector col_partitioning( nc ) ;
   for( size_t i=0 ; i<nc ; ++i )
   {
      if( cp( i ) <= 0 )
      {
         PEL_Error::object()->raise_data_error(
            exp, "col_partitioning",
            "invalid null or negative column dimension" ) ;
      }
      col_partitioning( i ) = (size_t) cp( i ) ;
   }

   PEL_ModuleExplorer const* sexp =
                         exp->create_subexplorer( 0, "block_prototype" ) ;
   LA_SeqMatrix* a_block_proto = LA_SeqMatrix::make( 0, sexp ) ;
   if( !a_block_proto->is_resizable() ||
       a_block_proto->is_desynchronizable() ||
       a_block_proto->is_symmetric() )
   {
      PEL_Error::object()->raise_module_error(
             sexp, "invalid block prototype matrix" ) ;
   }
   sexp->destroy() ; sexp = 0 ;

   LA_BlockSeqMatrix* result = new LA_BlockSeqMatrix( a_owner,
                                                      row_partitioning,
                                                      col_partitioning,
                                                      a_block_proto ) ;
   a_block_proto->destroy() ; a_block_proto = 0 ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrix:: LA_BlockSeqMatrix( void )
//----------------------------------------------------------------------
   : LA_SeqMatrix( "LA_BlockSeqMatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_PARTITION( 0 )
   , COL_PARTITION( 0 )
   , ROW_BLOCK_2_FIRST_ELEM( 0 )
   , COL_BLOCK_2_FIRST_ELEM( 0 )
   , ROW_ELEM_2_BLOCK( 0 )
   , COL_ELEM_2_BLOCK( 0 )
   , VEC_X( 0 )
   , VEC_Y( 0 )
   , BLOCK_PROTO( 0 )
   , BLOCKS( 0 )
{
   PEL_LABEL( "LA_BlockSeqMatrix:: LA_BlockSeqMatrix( prototype )" ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == "LA_BlockSeqMatrix" ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( !is_desynchronizable() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrix:: LA_BlockSeqMatrix(
                                  PEL_Object* a_owner,
                                  size_t_vector const& a_row_partitioning,
                                  size_t_vector const& a_col_partitioning,
                                  LA_SeqMatrix const* a_block_proto )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_BlockSeqMatrix" )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_PARTITION( a_row_partitioning )
   , COL_PARTITION( a_col_partitioning )
   , ROW_BLOCK_2_FIRST_ELEM( 0 )
   , COL_BLOCK_2_FIRST_ELEM( 0 )
   , ROW_ELEM_2_BLOCK( 0 )
   , COL_ELEM_2_BLOCK( 0 )
   , VEC_X( LA_SeqVector::create( this, 0 ) )
   , VEC_Y( LA_SeqVector::create( this, 0 ) )
   , BLOCK_PROTO( a_block_proto->create_matrix( this ) )
   , BLOCKS( 0 )
{
   PEL_LABEL( "LA_BlockSeqMatrix:: LA_BlockSeqMatrix" ) ;
   PEL_CHECK_PRE( a_row_partitioning.size()>0 ) ;
   PEL_CHECK_PRE( a_col_partitioning.size()>0 ) ;
   PEL_CHECK_PRE( a_block_proto != 0 ) ;
   PEL_CHECK_PRE( !a_block_proto->is_symmetric() ) ;
   PEL_CHECK_PRE( !a_block_proto->is_desynchronizable() ) ;
   PEL_CHECK_PRE( a_block_proto->is_resizable() ) ;

   init() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == "LA_BlockSeqMatrix" ) ;
   PEL_CHECK_POST( !is_symmetric() ) ;
   PEL_CHECK_POST( !is_desynchronizable() ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( row_partition() == a_row_partitioning ) ;
   PEL_CHECK_POST( col_partition() == a_col_partitioning ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
   PEL_CHECK_POST( nb_rows() == a_row_partitioning.sum() ) ;
   PEL_CHECK_POST( nb_cols() == a_col_partitioning.sum() ) ;
   PEL_CHECK_POST( block_prototype()->name() == a_block_proto->name() ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrix:: ~LA_BlockSeqMatrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_block_prototype( LA_SeqMatrix const* a_proto )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_block_prototype" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( a_proto != 0 ) ;
   PEL_CHECK_PRE( a_proto->owner() == this ) ;
   PEL_CHECK_PRE( a_proto->nb_rows() == 0 ) ;
   PEL_CHECK_PRE( a_proto->nb_cols() == 0 ) ;
   PEL_CHECK_PRE( !a_proto->is_desynchronizable() ) ;
   PEL_CHECK_PRE( !a_proto->is_symmetric() ) ;
   PEL_CHECK_PRE( a_proto->is_resizable() ) ;
   PEL_CHECK_INV( invariant() ) ;

   BLOCK_PROTO = a_proto ;

   for( size_t i=0 ; i<BLOCKS->index_limit() ; i++ )
   {
      LA_SeqMatrix const* old_mat =
               static_cast<LA_SeqMatrix const*>( BLOCKS->at( i ) ) ;
      if( old_mat != 0 )
      {
         LA_SeqMatrix* new_mat =
                          BLOCK_PROTO->create_copy( BLOCKS, old_mat ) ;
         new_mat->make_non_resizable() ;
         BLOCKS->set_at( i, new_mat ) ;
         PEL_ASSERT( new_mat->state() == old_mat->state() ) ;
         BLOCKS->destroy_possession( old_mat ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( block_prototype() == a_proto ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t ib=0 ; ib<row_partition().size() ; ++ib ),
      FORALL( ( size_t jb=0 ; jb<col_partition().size() ; ++jb ),
              IMPLIES( has_submatrix( ib, jb ),
                       submatrix( ib, jb )->name() == a_proto->name() ) ) ) ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix const*
LA_BlockSeqMatrix:: block_prototype( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: block_prototype" ) ;

   LA_SeqMatrix const* result = BLOCK_PROTO ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_rows() == 0 ) ;
   PEL_CHECK_POST( result->nb_cols() == 0 ) ;
   PEL_CHECK_POST( !result->is_desynchronizable() ) ;
   PEL_CHECK_POST( !result->is_symmetric() ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: allocated_memory" ) ;

   size_t result = 0 ;

   for( size_t i=0 ; i<BLOCKS->index_limit() ; ++i )
   {
      LA_SeqMatrix const* mat =
                  static_cast<LA_SeqMatrix const* >( BLOCKS->at( i ) ) ;
      if( mat!=0 )
      {
         result += mat->allocated_memory() ;
      }
   }

   size_t const nb_size_t =
             ROW_PARTITION.size()
           + COL_PARTITION.size()
           + ROW_BLOCK_2_FIRST_ELEM.size()
           + COL_BLOCK_2_FIRST_ELEM.size()
           + ROW_ELEM_2_BLOCK.size()
           + ROW_ELEM_2_BLOCK.size() ;
   result += nb_size_t*sizeof(size_t) ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: is_desynchronizable" ) ;

   bool result = false ;

   PEL_CHECK_POST( result == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: nb_stored_items" ) ;
   PEL_CHECK_PRE( nb_stored_items_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = 0 ;
   for( size_t i=0 ; i<BLOCKS->index_limit() ; ++i )
   {
      LA_SeqMatrix const* mat =
           dynamic_cast<LA_SeqMatrix const* >( BLOCKS->at( i ) ) ;
      if( mat!=0 )
      {
         result += mat->nb_stored_items() ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
double
LA_BlockSeqMatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: item" ) ;
   PEL_CHECK_PRE( item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const ib = ROW_ELEM_2_BLOCK( i ) ;
   size_t const jb = COL_ELEM_2_BLOCK( j ) ;

   double result = 0.0 ;
   if( has_submatrix( ib, jb ) )
   {
      LA_SeqMatrix const* mat = submatrix( ib, jb ) ;
      result = mat->item( i-ROW_BLOCK_2_FIRST_ELEM( ib ),
                          j-COL_BLOCK_2_FIRST_ELEM( jb ) ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_BlockSeqMatrix:: row_partition( void ) const
//----------------------------------------------------------------------
{
   return( ROW_PARTITION ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_BlockSeqMatrix:: col_partition( void ) const
//----------------------------------------------------------------------
{
   return( COL_PARTITION ) ;
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrix:: has_submatrix( size_t ib, size_t jb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix::has_submatrix" ) ;
   PEL_CHECK( ib<ROW_PARTITION.size() ) ;
   PEL_CHECK( jb<COL_PARTITION.size() ) ;

   size_t const idx = block_idx( ib, jb ) ;
   return( BLOCKS->at( idx ) != 0 ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix const*
LA_BlockSeqMatrix:: submatrix( size_t ib, size_t jb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: submatrix" ) ;
   PEL_CHECK_PRE( ib < row_partition().size() ) ;
   PEL_CHECK_PRE( jb < col_partition().size() ) ;
   PEL_CHECK_PRE( has_submatrix( ib, jb ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const idx = block_idx( ib, jb ) ;
   LA_SeqMatrix const* result =
              static_cast<LA_SeqMatrix *>( BLOCKS->at( idx ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->nb_rows() == row_partition()(ib) ) ;
   PEL_CHECK_POST( result->nb_cols() == col_partition()(jb) ) ;
   PEL_CHECK_POST( result->name() == block_prototype()->name() ) ;
   PEL_CHECK_POST( !result->is_desynchronizable() ) ;
   PEL_CHECK_POST( !result->is_symmetric() ) ;
   PEL_CHECK_POST( !result->is_resizable() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_BlockSeqMatrix:: create_stored_item_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: create_stored_item_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_BlockSeqMatrixIterator* result =
                       LA_BlockSeqMatrixIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: re_initialize_with_global_sizes( size_t a_nb_rows,
                                                     size_t a_nb_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: re_initialize_with_global_sizes" ) ;
   PEL_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows,
                                                       a_nb_cols ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t_vector rows( 1 ) ;
   size_t_vector cols( 1 ) ;
   rows( 0 ) = a_nb_rows ;
   cols( 0 ) = a_nb_cols ;
   set_partitioning_then_nullify( rows, cols ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows,
                                                         a_nb_cols ) ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_partitioning_then_nullify(
                                     size_t_vector const& row_partitioning,
                                     size_t_vector const& col_partitioning )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_partitioning_then_nullify" ) ;
   PEL_CHECK_PRE( is_resizable() ) ;
   PEL_CHECK_PRE( row_partitioning.size() != 0 ) ;
   PEL_CHECK_PRE( col_partitioning.size() != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   ROW_PARTITION = row_partitioning ;
   COL_PARTITION = col_partitioning ;
   init() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( row_partition() == row_partitioning ) ;
   PEL_CHECK_POST( col_partition() == col_partitioning ) ;
   PEL_CHECK_POST( nb_rows() == row_partitioning.sum() ) ;
   PEL_CHECK_POST( nb_cols() == col_partitioning.sum() ) ;
   PEL_CHECK_POST( nb_stored_items() == 0 ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_stored_items( double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_stored_items" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<BLOCKS->index_limit() ; i++ )
   {
      PEL_Object* obj = BLOCKS->at( i ) ;
      LA_SeqMatrix* mat = static_cast<LA_SeqMatrix* >( obj ) ;
      if( mat!=0 )
      {
         mat->set_stored_items(val) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_stored_items_POST( val ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: nullify_row( size_t i )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: nullify_row" ) ;
   PEL_CHECK_PRE( nullify_row_PRE( i ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "nullify_row" ) ;

   PEL_CHECK_POST( nullify_row_POST( i ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: nullify_item( size_t i, size_t j )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: nullify_item" ) ;
   PEL_CHECK_PRE( i < nb_rows() ) ;
   PEL_CHECK_PRE( j < nb_cols() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const ib = ROW_ELEM_2_BLOCK( i ) ;
   size_t const jb = COL_ELEM_2_BLOCK( j ) ;

   LA_SeqMatrix* m = submat( ib, jb ) ;
   m->set_item( i-ROW_BLOCK_2_FIRST_ELEM( ib ),
                j-COL_BLOCK_2_FIRST_ELEM( jb ), 0. ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( set_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const ib = ROW_ELEM_2_BLOCK( i ) ;
   size_t const jb = COL_ELEM_2_BLOCK( j ) ;
   set_submatrix_item( ib, jb,
                       i-ROW_BLOCK_2_FIRST_ELEM( ib ),
                       j-COL_BLOCK_2_FIRST_ELEM( jb ), x ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const ib = ROW_ELEM_2_BLOCK( i ) ;
   size_t const jb = COL_ELEM_2_BLOCK( j ) ;
   add_to_submatrix_item( ib, jb,
                          i-ROW_BLOCK_2_FIRST_ELEM( ib ),
                          j-COL_BLOCK_2_FIRST_ELEM( jb ), x ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_submatrix_item( size_t ib, size_t jb,
                                        size_t i, size_t j,
                                        double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_submatrix_item" ) ;
   PEL_CHECK_PRE( ib < row_partition().size() ) ;
   PEL_CHECK_PRE( jb < col_partition().size() ) ;
   PEL_CHECK_PRE( i < row_partition()(ib) ) ;
   PEL_CHECK_PRE( j < col_partition()(jb) ) ;
   PEL_CHECK_INV( invariant() ) ;

   submat( ib, jb )->set_item( i, j, x ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( has_submatrix( ib, jb ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: add_to_submatrix_item( size_t ib, size_t jb,
                                           size_t i, size_t j,
                                           double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_to_submatrix_item" ) ;
   PEL_CHECK_PRE( ib < row_partition().size() ) ;
   PEL_CHECK_PRE( jb < col_partition().size() ) ;
   PEL_CHECK_PRE( i < row_partition()(ib) ) ;
   PEL_CHECK_PRE( j < col_partition()(jb) ) ;
   PEL_CHECK_INV( invariant() ) ;

   submat( ib, jb )->add_to_item( i, j, x ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( has_submatrix( ib, jb ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha, double beta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   y->scale( beta ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( x ) != 0 ) ;
   LA_SeqVector const* seq_x = static_cast<LA_SeqVector const*>( x ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( y ) != 0 ) ;
   LA_SeqVector* seq_y = static_cast<LA_SeqVector* >( y ) ;

   double* x_ptr = const_cast<double*>( seq_x->data() ) ;
   double* y_ptr = seq_y->data() ;

   for( size_t jb=0 ; jb<COL_PARTITION.size() ; jb++ )
   {
      VEC_X->set( COL_PARTITION(jb),
                  x_ptr+COL_BLOCK_2_FIRST_ELEM(jb) ) ;
      for( size_t ib=0 ; ib<ROW_PARTITION.size() ; ib++ )
      {
         if( has_submatrix( ib, jb ) )
         {
            VEC_Y->set( ROW_PARTITION(ib),
                        y_ptr+ROW_BLOCK_2_FIRST_ELEM(ib) ) ;
            LA_SeqMatrix const* mat = submatrix( ib, jb ) ;
            mat->multiply_vec_then_add( VEC_X, VEC_Y,
                                        alpha, 1.0 ) ;
         }
      }
   }

   VEC_X->set( 0, (double*) 0 ) ;
   VEC_Y->set( 0, (double*) 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: tr_multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha, double beta )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: tr_multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   y->scale( beta ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( x ) != 0 ) ;
   LA_SeqVector const* seq_x = static_cast<LA_SeqVector const* >( x ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( y ) != 0 ) ;
   LA_SeqVector* seq_y = static_cast<LA_SeqVector* >( y ) ;


   double* x_ptr = const_cast<double*>( seq_x->data() ) ;
   double* y_ptr = seq_y->data() ;

   for( size_t jb=0 ; jb<COL_PARTITION.size() ; jb++ )
   {
      VEC_Y->set( COL_PARTITION(jb),
                  y_ptr+COL_BLOCK_2_FIRST_ELEM(jb) ) ;
      for( size_t ib=0 ; ib<ROW_PARTITION.size() ; ib++ )
      {
         if( has_submatrix( ib, jb ) )
         {
            VEC_X->set( ROW_PARTITION(ib),
                        x_ptr+ROW_BLOCK_2_FIRST_ELEM(ib) ) ;
            LA_SeqMatrix const* mat = submatrix( ib, jb ) ;
            mat->tr_multiply_vec_then_add( VEC_X, VEC_Y,
                                           alpha, 1.0 ) ;
         }
      }
   }

   VEC_X->set( 0, (double*) 0 ) ;
   VEC_Y->set( 0, (double*) 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( alpha == 0.0 )
   {
      nullify() ;
   }
   else if( alpha != 1.0 )
   {
      for( size_t i=0 ; i<BLOCKS->index_limit() ; ++i )

      {
         PEL_Object* obj = BLOCKS->at(i) ;
         if( obj != 0 )
         {
            static_cast<LA_SeqMatrix*>( obj )->scale( alpha ) ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_BlockSeqMatrix const* AS = dynamic_cast<LA_BlockSeqMatrix const*>( A ) ;
   if( AS==0 || ( ROW_PARTITION != AS->ROW_PARTITION )
             || ( COL_PARTITION != AS->COL_PARTITION ) )
   {
      LA_SeqMatrix::set( A, same_pattern ) ;
   }
   else
   {
      set_IMP( AS, same_pattern ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: add_Mat( LA_Matrix const* A, double alpha,
                             bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_BlockSeqMatrix const* AS =
      dynamic_cast<LA_BlockSeqMatrix const*>(A) ;
   if( AS==0 ||
       ( ROW_PARTITION != AS->ROW_PARTITION ) ||
       ( COL_PARTITION != AS->COL_PARTITION ) )
   {
      LA_SeqMatrix::add_Mat( A, alpha, same_pattern ) ;
   }
   else
   {
      add_Mat_IMP( AS, alpha, same_pattern ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_Mat_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_BlockSeqMatrix const* BA =
      dynamic_cast<LA_BlockSeqMatrix const*>( A ) ;
   LA_BlockSeqMatrix const* BB =
      dynamic_cast<LA_BlockSeqMatrix const*>( B ) ;
   if( ( BA!=0 && BB!=0 ) &&
       ( BA->ROW_PARTITION == ROW_PARTITION &&
         BB->COL_PARTITION == COL_PARTITION &&
         BA->COL_PARTITION == BB->ROW_PARTITION ) )
   {
      add_Mat_Mat_IMP( BA, BB, alpha ) ;
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
LA_BlockSeqMatrix:: readMM( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: readMM" ) ;
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
LA_BlockSeqMatrix:: restore( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: restore" ) ;
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
LA_BlockSeqMatrix:: add_Mat_Mat_IMP( LA_BlockSeqMatrix const* A,
                                     LA_BlockSeqMatrix const* B,
                                     double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_Mat_Mat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_PRE( A->row_partition() == row_partition() &&
                  B->col_partition() == col_partition() ) ;
   PEL_CHECK_PRE( A->col_partition() == B->row_partition() )  ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t ib = 0 ; ib<A->ROW_PARTITION.size() ; ++ib )
   {
      for( size_t kb = 0 ; kb<A->COL_PARTITION.size() ; ++kb )
      {
         if( A->has_submatrix( ib, kb ) )
         {
            LA_SeqMatrix const* matA = A->submatrix( ib, kb ) ;
            for( size_t jb = 0 ; jb<COL_PARTITION.size() ; ++jb )
            {
               if( B->has_submatrix( kb, jb ) )
               {
                  LA_SeqMatrix const* matB = B->submatrix( kb, jb ) ;
                  submat( ib, jb )->add_Mat_Mat( matA, matB, alpha ) ;
               }
            }
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: set_IMP( LA_BlockSeqMatrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: set_IMP" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_PRE( A->row_partition() == row_partition() &&
                  A->col_partition() == col_partition() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t nbi = ROW_PARTITION.size() ;
   for( size_t ib = 0 ; ib<nbi ; ++ib )
   {
      PEL_CHECK( ROW_PARTITION( ib ) == A->ROW_PARTITION( ib ) ) ;
      size_t nbj = COL_PARTITION.size() ;
      for( size_t jb = 0 ; jb<nbj ; ++jb )
      {
         PEL_CHECK( COL_PARTITION( jb ) == A->COL_PARTITION( jb ) ) ;
         if( A->has_submatrix( ib, jb ) )
         {
            LA_SeqMatrix const* mat = A->submatrix( ib, jb ) ;
            submat( ib, jb )->set( mat, same_pattern ) ;
         }
         else
         {
            PEL_ASSERT( !same_pattern ) ;
            if( has_submatrix( ib, jb ) )
            {
               size_t const idx = block_idx( ib, jb ) ;
               LA_SeqMatrix* old_mat = submat( ib, jb ) ;
               BLOCKS->destroy_possession( old_mat ) ;
               BLOCKS->set_at( idx, 0 ) ;
            }
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: add_Mat_IMP( LA_BlockSeqMatrix const* A,
                                 double alpha, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: add_Mat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_PRE( A->row_partition() == row_partition() &&
                  A->col_partition() == col_partition() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t nbi = ROW_PARTITION.size() ;
   for( size_t ib = 0 ; ib<nbi ; ++ib )
   {
      PEL_CHECK( ROW_PARTITION( ib ) == A->ROW_PARTITION( ib ) ) ;
      size_t nbj = COL_PARTITION.size() ;
      for( size_t jb = 0 ; jb<nbj ; ++jb )
      {
         PEL_CHECK( COL_PARTITION( jb ) == A->COL_PARTITION( jb ) ) ;
         if( A->has_submatrix( ib, jb ) )

         {
            LA_SeqMatrix const* mat = A->submatrix( ib, jb ) ;
            submat( ib, jb )->add_Mat( mat, alpha, same_pattern ) ;
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_BlockSeqMatrix:: submat( size_t ib, size_t jb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: submat" ) ;
   PEL_CHECK_PRE( ib < row_partition().size() ) ;
   PEL_CHECK_PRE( jb < col_partition().size() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const idx = block_idx( ib, jb ) ;
   LA_SeqMatrix* result =
              static_cast<LA_SeqMatrix *>( BLOCKS->at( idx ) ) ;
   if( result == 0 )
   {
      result = BLOCK_PROTO->create_matrix( BLOCKS ) ;
      result->re_initialize( ROW_PARTITION( ib ), COL_PARTITION( jb ) ) ;
      result->make_non_resizable() ;
      BLOCKS->set_at( idx, result ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->name() == block_prototype()->name() ) ;
   PEL_CHECK_POST( !result->is_symmetric() ) ;
   PEL_CHECK_POST( !result->is_desynchronizable() ) ;
   PEL_CHECK_POST( !result->is_resizable() ) ;
   PEL_CHECK_POST( result->nb_rows() == row_partition()(ib) ) ;
   PEL_CHECK_POST( result->nb_cols() == col_partition()(jb) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrix:: init( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix:: init" ) ;
   PEL_CHECK( ROW_PARTITION.size() > 0 ) ;
   PEL_CHECK( COL_PARTITION.size() > 0 ) ;

   ROW_BLOCK_2_FIRST_ELEM.re_initialize( ROW_PARTITION.size()+1 ) ;
   COL_BLOCK_2_FIRST_ELEM.re_initialize( COL_PARTITION.size()+1 ) ;

   NB_ROWS = 0 ;
   for( size_t i = 0 ; i<ROW_PARTITION.size() ; i++ )
   {
      ROW_BLOCK_2_FIRST_ELEM( i ) = NB_ROWS  ;
      NB_ROWS += ROW_PARTITION(i) ;
   }
   ROW_BLOCK_2_FIRST_ELEM( ROW_PARTITION.size() ) = NB_ROWS ;

   ROW_ELEM_2_BLOCK.re_initialize( NB_ROWS ) ;
   for( size_t i = 0 ; i<ROW_PARTITION.size() ; i++ )
   {
      for( size_t j = ROW_BLOCK_2_FIRST_ELEM( i ) ;
           j<ROW_BLOCK_2_FIRST_ELEM( i+1 ) ; j++ )
      {
         ROW_ELEM_2_BLOCK( j ) = i ;
      }
   }

   NB_COLS = 0 ;
   for( size_t i = 0 ; i<COL_PARTITION.size() ; i++ )
   {
      COL_BLOCK_2_FIRST_ELEM( i ) = NB_COLS ;
      // PEL_ASSERT( COL_PARTITION(i)>0 ) ;
      NB_COLS += COL_PARTITION(i) ;
   }
   COL_BLOCK_2_FIRST_ELEM( COL_PARTITION.size() ) = NB_COLS ;

   COL_ELEM_2_BLOCK.re_initialize( NB_COLS ) ;
   for( size_t i = 0 ; i<COL_PARTITION.size() ; i++ )
   {
      for( size_t j = COL_BLOCK_2_FIRST_ELEM( i ) ;
           j<COL_BLOCK_2_FIRST_ELEM( i+1 ) ; j++ )
      {
         COL_ELEM_2_BLOCK( j ) = i ;
      }
   }

   if( BLOCKS != 0 )
   {
      destroy_possession( BLOCKS ) ;
      BLOCKS = 0 ;
   }
   BLOCKS = PEL_Vector::create(
                      this, ROW_PARTITION.size()*COL_PARTITION.size() ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrix:: block_idx( size_t ib, size_t jb ) const
//----------------------------------------------------------------------
{
   return( ib+jb*ROW_PARTITION.size() ) ;
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ||
               ( ROW_BLOCK_2_FIRST_ELEM.size() == ROW_PARTITION.size()+1 &&
                 COL_BLOCK_2_FIRST_ELEM.size() == COL_PARTITION.size()+1 ) ) ;
   PEL_ASSERT( ROW_ELEM_2_BLOCK.size() == NB_ROWS ) ;
   PEL_ASSERT( COL_ELEM_2_BLOCK.size() == NB_COLS ) ;
   PEL_ASSERT( is_a_prototype() || (BLOCK_PROTO != 0) ) ;
   PEL_ASSERT( is_a_prototype() ||
               ( BLOCKS != 0 &&
                 BLOCKS->index_limit() ==
                 ROW_PARTITION.size()*COL_PARTITION.size() ) ) ;

   return( true ) ;
}

