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

#include <LA_PelMatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <LA_PelMatrix.hh>

//----------------------------------------------------------------------
LA_PelMatrixIterator*
LA_PelMatrixIterator:: create( PEL_Object* a_owner,
                               LA_PelMatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;

   LA_PelMatrixIterator* result =  new LA_PelMatrixIterator( a_owner, A ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
LA_PelMatrixIterator:: LA_PelMatrixIterator( PEL_Object* a_owner,
                                             LA_PelMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , CUR_ELEM( 0 )
   , I_ROW( 0 )
   , STAY_IN_ROW( false )
{
}

//----------------------------------------------------------------------
LA_PelMatrixIterator:: ~LA_PelMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_PelMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: is_valid" ) ;

   bool result = false ;
   if( MAT->ROW_TABLE != 0 )
   {
      result = ( CUR_ELEM != 0 ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;

   I_ROW = 0 ;
   STAY_IN_ROW = false ;
   CUR_ELEM = 0 ;
   if( MAT->ROW_TABLE != 0 )
   {
      CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
      if( CUR_ELEM == 0 )
      {
         go_next_element() ;
      }
   }
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   I_ROW = i_row ;
   STAY_IN_ROW = true ;
   CUR_ELEM = 0 ;
   if( MAT->ROW_TABLE != 0 )
   {
      CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
   }
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   go_next_element() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_PelMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW ;

   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result =  CUR_ELEM->iCol ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_PelMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   return( CUR_ELEM->xVal ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CUR_ELEM->xVal = x ;
   if( MAT->UNSYNCHRO && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CUR_ELEM->xVal += x ;
   if( MAT->UNSYNCHRO && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: go_next_element( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PelMatrixIterator:: go_next_element" ) ;

   if( CUR_ELEM!=0 )
   {
      CUR_ELEM = CUR_ELEM->next ;
   }
   if( !STAY_IN_ROW )
   {
      while( I_ROW<MAT->nb_rows()-1 && CUR_ELEM==0 )
      {
         I_ROW++ ;
         CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
      }
   }
}
