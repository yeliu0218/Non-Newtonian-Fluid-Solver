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

#include <LA_CRSmatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <LA_CRSmatrix.hh>

//----------------------------------------------------------------------
LA_CRSmatrixIterator*
LA_CRSmatrixIterator:: create( PEL_Object* a_owner,
                               LA_CRSmatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;

   LA_CRSmatrixIterator* result =  new LA_CRSmatrixIterator( a_owner, A ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrixIterator:: LA_CRSmatrixIterator( PEL_Object* a_owner,
                                             LA_CRSmatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , VALID( false )
   , I_ROW( 0 )
   , STAY_IN_ROW( false )
   , I_CUR( PEL::bad_index() )
{
}

//----------------------------------------------------------------------
LA_CRSmatrixIterator:: ~LA_CRSmatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_CRSmatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( VALID ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;
   
   STAY_IN_ROW = false ;
   I_CUR = 0 ;
   I_ROW = 0 ;
   VALID = ( I_CUR< (int)MAT->NB_ELEMS ) ;
   if( VALID )
   {
      while( MAT->START(I_ROW+1)<=I_CUR ) I_ROW++ ;
      PEL_CHECK( I_CUR >=  MAT->START(I_ROW) ) ;
      PEL_CHECK( I_CUR <  MAT->START(I_ROW+1) ) ;
   }
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   I_ROW = i_row ;
   STAY_IN_ROW = true ;
   I_CUR = MAT->START(i_row) ;
   VALID = I_CUR < MAT->START(i_row+1) ;        
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   I_CUR++ ;
   if( STAY_IN_ROW )
   {
      VALID = ( I_CUR<MAT->START(I_ROW+1) ) ;
   }
   else
   {
      VALID = ( I_CUR<(int)MAT->NB_ELEMS ) ;
      if( VALID )
      {
         while( MAT->START(I_ROW+1)<=I_CUR ) I_ROW++ ;
         PEL_CHECK( I_CUR >=  MAT->START(I_ROW) ) ;
         PEL_CHECK( I_CUR <  MAT->START(I_ROW+1) ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_CRSmatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW  ;
   
   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;

   size_t result = MAT->COL(I_CUR) ;

   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_CRSmatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   return( MAT->VALUES(I_CUR) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_CRSmatrix* dummy = const_cast<LA_CRSmatrix*>( MAT ) ;
   dummy->VALUES(I_CUR) = x ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_CRSmatrix* dummy = const_cast<LA_CRSmatrix*>( MAT ) ;
   dummy->VALUES(I_CUR) += x ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
