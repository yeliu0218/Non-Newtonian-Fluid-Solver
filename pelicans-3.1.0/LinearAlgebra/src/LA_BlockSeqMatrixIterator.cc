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

#include <LA_BlockSeqMatrixIterator.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <LA.hh>
#include <LA_BlockSeqMatrix.hh>

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator*
LA_BlockSeqMatrixIterator:: create(
                       PEL_Object* a_owner, LA_BlockSeqMatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: create" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   
   LA_BlockSeqMatrixIterator* result =
                   new LA_BlockSeqMatrixIterator( a_owner, A ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->matrix() == A ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator:: LA_BlockSeqMatrixIterator(
                       PEL_Object* a_owner, LA_BlockSeqMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , MAT_IT( 0 )
   , IS_VALID( false )
   , IB( PEL::bad_index() )
   , JB( PEL::bad_index() )
   , STAY_IN_ROW( false )
   , IROW( PEL::bad_index() )
   , ROW_SHIFT( PEL::bad_index() )
   , COL_SHIFT( PEL::bad_index() )
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: LA_BlockSeqMatrixIterator" ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator:: ~LA_BlockSeqMatrixIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: ~LA_BlockSeqMatrixIterator" ) ;
   
   if( MAT_IT!=0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }   
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( IS_VALID ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: start_all_items" ) ;
   PEL_CHECK_PRE( start_all_items_PRE() ) ;
   
   IB = 0 ;
   JB = 0 ;
   if( MAT_IT != 0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }
   STAY_IN_ROW = false ;
   IROW = PEL::bad_index() ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: start_row_items" ) ;
   PEL_CHECK_PRE( start_row_items_PRE( i_row ) ) ;
   
   IB = MAT->ROW_ELEM_2_BLOCK( i_row ) ;
   JB = 0 ;
   if( MAT_IT != 0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }
   STAY_IN_ROW = true ;
   IROW = i_row-MAT->ROW_BLOCK_2_FIRST_ELEM( IB ) ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_BlockSeqMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   PEL_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: nb_rows" ) ;
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: row" ) ;
   PEL_CHECK_PRE( row_PRE() ) ;

   PEL_CHECK( MAT_IT!=0 ) ;
   size_t result = MAT_IT->row()+ROW_SHIFT ;

   PEL_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: col" ) ;
   PEL_CHECK_PRE( col_PRE() ) ;
   
   PEL_CHECK( MAT_IT!=0 ) ;
   size_t result = MAT_IT->col() + COL_SHIFT ;
   
   PEL_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_BlockSeqMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   PEL_CHECK( MAT_IT!=0 ) ;   
   PEL_CHECK( MAT_IT->item()==MAT->item( row(), col() ) ) ;
   
   return( MAT_IT->item() ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: set_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( set_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( MAT_IT!=0 ) ;
   MAT_IT->set_item( x ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: add_to_item" ) ;
   PEL_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   PEL_CHECK_PRE( add_to_item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( MAT_IT!=0 ) ;
   MAT_IT->add_to_item( x ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrixIterator:: next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrixIterator:: next" ) ;
   
   bool result = false ;
   if( MAT_IT!=0 )
   {
      MAT_IT->go_next() ;
      result = MAT_IT->is_valid() ;
      if( !result )
      {
         MAT_IT->destroy() ;
         MAT_IT = 0 ;
         JB++ ;
      }
   }
   for( ; !result && IB<MAT->ROW_PARTITION.size() ; ++IB )
   {
      for( ; !result && JB<MAT->COL_PARTITION.size() ; ++JB )
      {
         if( MAT->has_submatrix( IB, JB ) )
         {
            LA_SeqMatrix* matS = MAT->submat( IB, JB ) ;
            if( !matS->is_synchronized() ) matS->synchronize() ;
            MAT_IT = matS->create_stored_item_iterator( 0 ) ;
            if( STAY_IN_ROW )
            {
               MAT_IT->start_row_items( IROW ) ;
            }
            else
            {
               MAT_IT->start_all_items() ;
            }
            result = MAT_IT->is_valid() ;
            if( !result )
            {
               MAT_IT->destroy() ;
               MAT_IT = 0 ;
            }
            else
            {
               ROW_SHIFT = MAT->ROW_BLOCK_2_FIRST_ELEM( IB ) ;
               COL_SHIFT = MAT->COL_BLOCK_2_FIRST_ELEM( JB ) ;
            }
         }
         if( result  )
         {
            break ;
         }
      }
      if( result || STAY_IN_ROW )
      {
         break ;
      }
      JB = 0 ;
   }
   PEL_CHECK( !result || MAT_IT!=0 ) ;
   return( result ) ;
}
   




