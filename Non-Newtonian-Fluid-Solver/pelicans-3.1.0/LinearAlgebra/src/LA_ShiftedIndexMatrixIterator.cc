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

#include <LA_ShiftedIndexMatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator*
LA_ShiftedIndexMatrixIterator:: create( PEL_Object* a_owner,
                                        size_t row_shift,
                                        size_t col_shift,
                                        LA_Matrix const* A,
                                        LA_MatrixIterator* internal )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( internal != 0 ) ;
   PEL_CHECK_PRE( internal->owner()==0 ) ;
   PEL_CHECK_PRE( internal->matrix()->nb_rows()+row_shift <= A->nb_rows() ) ;
   PEL_CHECK_PRE( internal->matrix()->nb_cols()+col_shift <= A->nb_cols() ) ;

   LA_ShiftedIndexMatrixIterator* result =
      new LA_ShiftedIndexMatrixIterator( a_owner,
                                         row_shift, col_shift, A,
                                         internal  ) ;
   internal->set_owner( result ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   PEL_CHECK_PRE( internal->owner()==result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator:: LA_ShiftedIndexMatrixIterator(
                                            PEL_Object* a_owner,
                                            size_t row_shift,
                                            size_t col_shift,
                                            LA_Matrix const* A,
                                            LA_MatrixIterator* internal )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , ROW_SHIFT( row_shift )
   , COL_SHIFT( col_shift )
   , INTERNAL( internal )
   , NULLROW(false)
{
}

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator:: ~LA_ShiftedIndexMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_ShiftedIndexMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( !NULLROW && INTERNAL->is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;

   NULLROW = false ;
   INTERNAL->start_all_items() ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   NULLROW = i_row<ROW_SHIFT || i_row>=INTERNAL->nb_rows()+ROW_SHIFT ;
   if( !NULLROW )
   {
      INTERNAL->start_row_items( (size_t)( (int)i_row-(int)ROW_SHIFT ) );
   }
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   INTERNAL->go_next() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_ShiftedIndexMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = INTERNAL->row()+ROW_SHIFT ;

   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result = INTERNAL->col()+COL_SHIFT ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_ShiftedIndexMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   return( INTERNAL->item() ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   INTERNAL->set_item( x ) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }


   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ShiftedIndexMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   INTERNAL->add_to_item( x ) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
