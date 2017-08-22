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

#include <LA_PairOfMatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_PairOfMatrixIterator*
LA_PairOfMatrixIterator:: create( PEL_Object* a_owner,
                                  LA_Matrix const* A,
                                  LA_MatrixIterator* first,
                                  LA_MatrixIterator* second )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( first->owner() == 0 ) ;
   PEL_CHECK_PRE( first->nb_rows() == A->nb_rows() ) ;
   PEL_CHECK_PRE( second->owner()==0 ) ;
   PEL_CHECK_PRE( second->nb_rows() == A->nb_rows() ) ;


   LA_PairOfMatrixIterator* result =
         new LA_PairOfMatrixIterator( a_owner, A, first, second  ) ;
   first->set_owner( result ) ;
   second->set_owner( result ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   PEL_CHECK_PRE( first->owner()==result ) ;
   PEL_CHECK_PRE( second->owner()==result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PairOfMatrixIterator:: LA_PairOfMatrixIterator( PEL_Object* a_owner,
                                                   LA_Matrix const* A,
                                                   LA_MatrixIterator* first,
                                                   LA_MatrixIterator* second )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , NB_ROWS( first->nb_rows() )
   , FIRST( first )
   , SECOND( second )
   , CURR( 0 )
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: LA_PairOfMatrixIterator" ) ;
}

//----------------------------------------------------------------------
LA_PairOfMatrixIterator:: ~LA_PairOfMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_PairOfMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( CURR!=0 && CURR->is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;

   FIRST->start_all_items() ;
   SECOND->start_all_items() ;
   CURR=( FIRST->is_valid() ? FIRST : SECOND ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   FIRST->start_row_items(i_row) ;
   SECOND->start_row_items(i_row) ;
   CURR=( FIRST->is_valid() ? FIRST : SECOND ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   CURR->go_next() ;
   if( !CURR->is_valid() && CURR!=SECOND )
   {
      CURR = SECOND ;
   }
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_PairOfMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = CURR->row() ;

   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result = CURR->col() ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_PairOfMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   return( CURR->item() ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CURR->set_item(x) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PairOfMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CURR->add_to_item(x) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
