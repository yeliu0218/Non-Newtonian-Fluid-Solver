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

#include <PEL_Vector.hh>

#include <iostream>

#include <PEL_System.hh>
#include <PEL_VectorIterator.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_Vector*
PEL_Vector:: create( PEL_Object* a_owner, size_t size )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: create" ) ;
   PEL_Vector* result =  new PEL_Vector( a_owner, size ) ;

   PEL_CHECK_POST( result != 0  ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->index_limit() == size ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; ++i ),
                            result->at(i) == 0 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Vector:: PEL_Vector( PEL_Object* a_owner, size_t size )
//----------------------------------------------------------------------
   : PEL_Sequence( a_owner )
   , VECTOR( 0 )
   , LENGTH( 0 )
   , NB_ENTRIES( 0 )
   , CAPACITY( 0 )
{
   PEL_LABEL( "PEL_Vector:: PEL_Vector" ) ;
   re_initialize( size ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Vector:: ~PEL_Vector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: ~PEL_Vector" ) ;
   PEL_CHECK_INV( invariant() ) ;

   clear() ;
}

//----------------------------------------------------------------------
PEL_Vector* 
PEL_Vector:: create_clone( PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
//????????? on devrait faire appel au copy constructeur ?????????????
//????????? trop trop lent
{
   PEL_LABEL( "PEL_Vector:: create_clone" ) ;
   PEL_Vector* result = new PEL_Vector( a_owner, 0 ) ;
   result->copy( this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: re_initialize( size_t size ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: re_initialize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( LENGTH!=size )
   {
      clear() ;
      LENGTH = size ;
      CAPACITY = size ;
      if( LENGTH>0 )
      {
         VECTOR = new PEL_Object*[ LENGTH ] ;
      }
      else
      {
         VECTOR=0 ;
      }
   }
   nullify() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( state_id() != OLD( state_id) ) ;
   PEL_CHECK_POST( index_limit() == size ) ;
   PEL_CHECK_POST( count()==0 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                            at(i) == 0 ) ) ;
   
}

//----------------------------------------------------------------------
void 
PEL_Vector:: copy( PEL_Vector const* other ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: copy" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   re_initialize( other->LENGTH ) ;
   NB_ENTRIES = other->count() ;
   for( size_t i=0 ; i<LENGTH ; i++ )
      VECTOR[i] = other->VECTOR[i] ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( count() == other->count() ) ;
   PEL_CHECK_POST( index_limit() == other->index_limit() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == other->at(i) ) ) ;
   PEL_CHECK_POST( state_id() != OLD( state_id ) ) ;
}

//----------------------------------------------------------------------
void 
PEL_Vector:: resize( size_t size )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: resize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   if( size > CAPACITY ) 
   {
      CAPACITY = PEL_System::new_block_size( CAPACITY, size ) ;
      PEL_Object** oldVector = VECTOR ;
      VECTOR = new PEL_Object*[ CAPACITY ] ;
      for( size_t i=LENGTH ; i<CAPACITY ; i++ )
         VECTOR[i] = 0 ;
      if( oldVector != 0 )
      {
         for( size_t i=0 ; i<LENGTH ; i++ )
            VECTOR[i] = oldVector[i] ;
        delete[] oldVector ;
      }
   }
   else if( size < LENGTH )
   {
      for( size_t i=size ; i<LENGTH  ; ++i )
      {
         if( VECTOR[ i ]!=0 )
	 {
	    NB_ENTRIES-- ;
	 }
         VECTOR[ i ] = 0 ;
      }
   }
   LENGTH = size ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( index_limit() == size ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t ii=OLD(index_limit) ; ii<index_limit() ; ++ii ),
              at(ii) == 0 ) ) ;
   PEL_CHECK_POST( state_id() != OLD( state_id ) ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Vector:: index_limit( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: index_limit" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( LENGTH ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: nullify( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: nullify" ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; i++ )
         VECTOR[i] = 0 ;
   }
   NB_ENTRIES = 0 ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( state_id() != OLD( state_id ) ) ;
   PEL_CHECK_POST( index_limit() == OLD( index_limit ) ) ;
   PEL_CHECK_POST( count()==0 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == 0 ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: append( PEL_Object* object )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: append" ) ;
   PEL_CHECK_PRE( append_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t oldSize = LENGTH ;
   resize( oldSize+1 ) ;
   set_at( oldSize, object ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( append_POST( OLD( index_limit ),
                                OLD( count ),
                                object,
                                OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: prepend( PEL_Object* object )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: prepend" ) ;
   PEL_CHECK_PRE( prepend_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( LENGTH == 0 )
   {
      append( object ) ;
   }
   else
   {
      insert_at( 0, object ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( prepend_POST( OLD( index_limit ),
                                 OLD( count ),
                                 object,
                                 OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: set_at( size_t i, PEL_Object* object )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: set_at" ) ;
   PEL_CHECK_PRE( ( object==0 && i < index_limit() ) || set_at_PRE( i, object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   PEL_Object* ith = VECTOR[ i ] ;
   VECTOR[ i ] = object ;
   if( ith==0 )
   {
      NB_ENTRIES++ ;
   }
   if( object==0 )
   {
      NB_ENTRIES-- ;
   }
   report_state_change() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_at_POST( OLD( index_limit ),
                                OLD( count ),
                                i,
                                object,
                                OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: insert_at( size_t i, PEL_Object* object )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: insert_at" ) ;
   PEL_CHECK_PRE( insert_at_PRE( i, object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t size = LENGTH+1 ;
   if( size > CAPACITY )
   {
      CAPACITY = PEL_System::new_block_size( CAPACITY, size ) ;
      PEL_Object** oldVector = VECTOR ;
      VECTOR = new PEL_Object*[ CAPACITY ] ;
      for( size_t j = 0 ; j < i ; ++j )
      {
         VECTOR[ j ] = oldVector[ j ] ;
      }
      for( size_t j = i ; j < LENGTH ; ++j )
      {
         VECTOR[ j+1 ] = oldVector[ j ] ;
      }
      delete [] oldVector ;
   }
   else
   {
      for( size_t j = LENGTH ; j > i ; --j )
      {
         VECTOR[ j ] = VECTOR[ j-1 ] ;
      }
   }
   VECTOR[ i ] = object ;
   NB_ENTRIES++ ;   
   LENGTH++ ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( insert_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   i,
                                   object,
                                   OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Vector:: count( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: count" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( count_POST( NB_ENTRIES ) ) ;
   return( NB_ENTRIES ) ;
}

//----------------------------------------------------------------------
PEL_Object*
PEL_Vector:: item( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: item" ) ;
   PEL_CHECK_PRE( item_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Object* result = 0 ;
   size_t i = 0 ;

   for( i = 0 ; i < LENGTH ; i++ )
   {
      if(  VECTOR[ i ]!=0 && VECTOR[ i ]->is_equal( object ) )
      {
         result = VECTOR[ i ] ;
         break ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( result, object ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Object*
PEL_Vector:: at( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: at" ) ;
   PEL_CHECK_PRE( at_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Object* resu = VECTOR[ i ] ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( at_POST( resu, i ) ) ;

   return( resu ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Vector:: index_of( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: index_of" ) ;
   PEL_CHECK_PRE( index_of_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t ret = badIndex ;
   
   for( size_t i = 0 ; i < LENGTH ; i++ )
   {
      if( VECTOR[ i ]!=0 && matching_items( VECTOR[ i ], object ) )
      {
         ret = i ;
         break ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( index_of_POST( ret, object ) ) ;

   return( ret ) ;
}

//----------------------------------------------------------------------
PEL_VectorIterator*
PEL_Vector:: create_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: create_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_VectorIterator* result = PEL_VectorIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
PEL_Vector:: remove_at( size_t i )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: remove_at" ) ;
   PEL_CHECK_PRE( remove_at_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   if( VECTOR[ i ] != 0 )
   {
      NB_ENTRIES-- ;
   }
   for( size_t j = i ; j < LENGTH - 1 ; ++j )
   {
      VECTOR[ j ] = VECTOR[ j + 1 ] ;
   }
   VECTOR[ LENGTH - 1 ] = 0 ;
  
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   OLD( state_id ) ) ) ;
}

//---------------------------------------------------------------------
void
PEL_Vector:: remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: remove_section" ) ;
   PEL_CHECK_PRE( remove_section_PRE( iFirst, length ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   for( size_t j = iFirst ; j < iFirst+length ; j++ )
   {
      if( VECTOR[ j ]!=0 )
      {
	 NB_ENTRIES-- ;
         VECTOR[ j ] = 0 ;
      }
   }
   for( size_t j = iFirst ; j < LENGTH - length ; j++ )
   {
      VECTOR[ j ] = VECTOR[ j + length ] ;
      VECTOR[ j + length ] = 0 ;
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}

//---------------------------------------------------------------------
void
PEL_Vector:: destroy_items_and_remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: destroy_items_and_remove_section" ) ;
   PEL_CHECK_PRE( destroy_items_and_remove_section_PRE( iFirst, length ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   for( size_t j = iFirst ; j < iFirst+length ; j++ )
   {
      if( VECTOR[ j ]!=0 )
      {
	 NB_ENTRIES-- ;
         VECTOR[ j ]->destroy() ;
         VECTOR[ j ] = 0 ;
      }
   }
   for( size_t j = iFirst ; j < LENGTH - length ; j++ )
   {
      VECTOR[ j ] = VECTOR[ j + length ] ;
      VECTOR[ j + length ] = 0 ;
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Vector:: clear( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Vector:: clear" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( VECTOR!=0 )
   {
      delete [] VECTOR ;
      VECTOR = 0 ;
      LENGTH = 0 ;
      CAPACITY = 0 ;
      NB_ENTRIES = 0 ;
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_Vector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::invariant() ) ;

   PEL_ASSERT( count() <= index_limit() ) ;
   PEL_ASSERT( LENGTH==0 || VECTOR!=0 ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Vector:: remove_at_POST( size_t old_index_limit,
                             size_t old_count,
                             size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::remove_at_POST( old_index_limit,
                                             old_count,
                                             old_state_id  ) ) ;

   PEL_ASSERT( index_limit() == old_index_limit ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Vector:: remove_section_POST( size_t old_index_limit,
                                  size_t old_count,
                                  size_t length,
                                  size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::remove_section_POST( old_index_limit, 
                                                  old_count, 
                                                  length,
                                                  old_state_id ) ) ;

   PEL_ASSERT( index_limit() == old_index_limit ) ;

   return( true ) ;
}

