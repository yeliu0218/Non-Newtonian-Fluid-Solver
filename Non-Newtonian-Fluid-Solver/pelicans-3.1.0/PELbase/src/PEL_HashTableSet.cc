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

#include <PEL_HashTableSet.hh>

#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>
#include <PEL_HashTableSetIterator.hh>
#include <PEL_assertions.hh>

//-------------------------------------------------------------------------
PEL_HashTableSet*
PEL_HashTableSet:: create( PEL_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: create" ) ;
   PEL_HashTableSet* result =  new PEL_HashTableSet( a_owner, size ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_buckets() == size ) ;
   PEL_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_HashTableSet:: PEL_HashTableSet( PEL_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
   : PEL_Set( a_owner ), 
     the_size( size )
{
   PEL_LABEL( "PEL_HashTableSet:: PEL_HashTableSet" ) ;
   hTable = PEL_Vector::create( this, the_size ) ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      hTable->set_at( i, PEL_List::create( this ) ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_HashTableSet:: ~PEL_HashTableSet( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: ~PEL_HashTableSet" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_HashTableSet*
PEL_HashTableSet:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_HashTableSet* result = PEL_HashTableSet::create( a_owner, 1 ) ;
   result->hTable->re_initialize( the_size ) ;
   result->the_size = the_size ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      result->hTable->set_at( i, getList(i)->create_clone( result ) ) ;     
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;

   return result ; 
}



//-------------------------------------------------------------------------
size_t
PEL_HashTableSet:: nb_buckets( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: nb_buckets( void )" ) ;
   return( the_size ) ;
}



//-------------------------------------------------------------------------
void
PEL_HashTableSet:: extend( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: extend" ) ;
   PEL_CHECK_PRE( extend_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( bool, has, has(object) ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   getList(entry)->extend( object ) ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_HashTableSet:: count( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: count" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t ret = 0 ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      ret += getList(i)->count() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( count_POST( ret ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_Object*
PEL_HashTableSet:: item( PEL_Object const* object ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: item" ) ;
   PEL_CHECK_PRE( item_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   PEL_Object* ret = getList(entry)->item( object ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( ret, object ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_HashTableSetIterator*
PEL_HashTableSet:: create_iterator( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: create_iterator" ) ;
   PEL_HashTableSetIterator* result = 
                          PEL_HashTableSetIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
PEL_HashTableSet:: remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: remove" ) ;
   PEL_CHECK_PRE( remove_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   getList(entry)->remove( object ) ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_POST( OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_HashTableSet:: destroy_item_and_remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: destroy_item_and_remove" ) ;
   PEL_Error::object()->raise_not_implemented( this, 
                                               "destroy_item_and_remove" ) ;
}



//-------------------------------------------------------------------------
void
PEL_HashTableSet:: clear( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: clear" ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   for( size_t i = 0 ; i<the_size ; i++ )
   {
       getList(i)->clear() ;
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_HashTableSet:: destroy_items_and_clear( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSet:: destroy_items_and_clear" ) ;
   PEL_Error::object()->raise_not_implemented( this, 
                                               "destroy_items_and_clear" ) ;
}



//-------------------------------------------------------------------------
PEL_List *
PEL_HashTableSet:: getList( size_t entry )
//-------------------------------------------------------------------------
{
   PEL_List * ret = dynamic_cast<PEL_List*>( hTable->at( entry ) ) ;
   PEL_CHECK( ret!=0 ) ;
   return ret ;
}



//-------------------------------------------------------------------------
PEL_List const*
PEL_HashTableSet:: getList( size_t entry ) const
//-------------------------------------------------------------------------
{
   PEL_List const* ret = dynamic_cast<PEL_List const*>( hTable->at( entry ) ) ;
   PEL_CHECK( ret!=0 ) ;
   return ret ;
}



//-------------------------------------------------------------------------
bool
PEL_HashTableSet:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( the_size>0 ) ;
   size_t cpt = 0 ;
   
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      PEL_List const* pList = getList( i ) ;
      PEL_ASSERT( pList!=0 ) ;
      cpt += pList->count() ;
   }
   PEL_ASSERT( cpt==count() ) ;
   return( true ) ;
} 
