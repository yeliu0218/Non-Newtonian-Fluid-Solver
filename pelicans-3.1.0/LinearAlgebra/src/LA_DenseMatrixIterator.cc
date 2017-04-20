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

#include <LA_DenseMatrixIterator.hh>

#include <LA_DenseMatrix.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

//----------------------------------------------------------------------
LA_DenseMatrixIterator*
LA_DenseMatrixIterator:: create( PEL_Object* a_owner,
                               LA_DenseMatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;

   LA_DenseMatrixIterator* result =  new LA_DenseMatrixIterator( a_owner, A ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DenseMatrixIterator:: LA_DenseMatrixIterator( PEL_Object* a_owner,
                                                 LA_DenseMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , I_ROW( 0 )
   , I_COL( 0 )
   , STAY_IN_ROW( false )
   , ROW_DIM( 0 )
   , COL_DIM( 0 )
{
}

//----------------------------------------------------------------------
LA_DenseMatrixIterator:: ~LA_DenseMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_DenseMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( I_ROW<ROW_DIM && I_COL<COL_DIM ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;
   
   ROW_DIM = MAT->nb_rows() ;
   COL_DIM = MAT->nb_cols() ;
   
   STAY_IN_ROW = false ;
   I_COL = 0 ;
   I_ROW = 0 ;
   
}

//----------------------------------------------------------------------
void
LA_DenseMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   ROW_DIM = MAT->nb_rows() ;
   COL_DIM = MAT->nb_cols() ;
   
   STAY_IN_ROW = true ;
   I_COL = 0 ;
   I_ROW = i_row ; 
}

//----------------------------------------------------------------------
void
LA_DenseMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   I_COL++ ;
   
   if( !STAY_IN_ROW && ( I_COL==COL_DIM ) )
   {
      I_COL = 0 ;
      I_ROW++ ;
   }
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_DenseMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW  ;
   
   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DenseMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result = I_COL ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_DenseMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;
   
   return( MAT->MAT[I_ROW][I_COL] ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   MAT->MAT[I_ROW][I_COL] = x ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DenseMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   MAT->MAT[I_ROW][I_COL] += x ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}

