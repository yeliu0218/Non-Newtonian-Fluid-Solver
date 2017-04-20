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

#include <LA_SymmetricMatrixIterator.hh>

#include <LA_SymmetricMatrix.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

//----------------------------------------------------------------------
LA_SymmetricMatrixIterator*
LA_SymmetricMatrixIterator:: create( PEL_Object* a_owner,
                                     LA_SymmetricMatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;

   LA_SymmetricMatrixIterator* result =
                          new LA_SymmetricMatrixIterator( a_owner, A ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrixIterator:: LA_SymmetricMatrixIterator(
                      PEL_Object* a_owner, LA_SymmetricMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , I_ROW( PEL::bad_index() )
   , I_COL( PEL::bad_index() )
   , STAY_IN_ROW( false )
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: LA_SymmetricMatrixIterator" ) ;
   DIM = MAT->nb_rows() ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrixIterator:: ~LA_SymmetricMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_SymmetricMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( I_ROW<DIM && I_COL<DIM ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;
   
   STAY_IN_ROW = false ;
   I_COL = 0 ;
   I_ROW = 0 ;
   
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   DIM = MAT->nb_rows() ;
   
   STAY_IN_ROW = true ;
   I_COL = 0 ;
   I_ROW = i_row ; 
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   I_COL++ ;
   
   if( !STAY_IN_ROW && ( I_COL==DIM ) )
   {
      I_COL = 0 ;
      I_ROW++ ;
   }
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_SymmetricMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW ;
   
   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SymmetricMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result = I_COL ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_SymmetricMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;
   
   return( MAT->item( I_ROW, I_COL ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SymmetricMatrix* dummy = const_cast<LA_SymmetricMatrix*>( MAT ) ;
   dummy->set_item( I_ROW, I_COL, x ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SymmetricMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SymmetricMatrix* dummy = const_cast<LA_SymmetricMatrix*>( MAT ) ;
   dummy->add_to_item( I_ROW, I_COL, x ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
