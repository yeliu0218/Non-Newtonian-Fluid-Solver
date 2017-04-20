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

#include <PEL_Map.hh>

#include <PEL_Error.hh>
#include <PEL_KeyItemPair.hh>
#include <PEL_MapIterator.hh>
#include <PEL_HashTableSet.hh>
#include <PEL_assertions.hh>

//-------------------------------------------------------------------------
PEL_Map*
PEL_Map:: create( PEL_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: create" ) ;
   PEL_Map* result = new PEL_Map( a_owner, size ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_buckets() == size ) ;
   PEL_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_Map:: PEL_Map( PEL_Object* a_owner, size_t aSize ) 
//-------------------------------------------------------------------------
   : PEL_Container( a_owner ), 
     hList( 0 )
{
   PEL_LABEL( "PEL_Map:: PEL_Map" ) ;
   hList = PEL_HashTableSet::create( this, aSize ) ;

   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_Map:: ~PEL_Map( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: ~PEL_Map" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_Map*
PEL_Map:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: create_clone" ) ;
   return( new PEL_Map( a_owner, this ) ) ;
}



//-------------------------------------------------------------------------
PEL_Map:: PEL_Map( PEL_Object* a_owner, PEL_Map const* other )
//-------------------------------------------------------------------------
   : PEL_Container( a_owner ),
     hList( 0 )
{
   PEL_LABEL( "PEL_Map:: PEL_Map" ) ;
   hList = other->hList->create_clone( this ) ;

   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_Map:: nb_buckets( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: nb_buckets" ) ;
   return( hList->nb_buckets() ) ;
}



//-------------------------------------------------------------------------
void
PEL_Map:: set_item_at( PEL_Object* key, PEL_Object* a_item  )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: set_item_at" ) ;
   PEL_CHECK_PRE( key != 0 ) ;
   PEL_CHECK_PRE( a_item != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Object* object = hList->item( key ) ;
   PEL_KeyItemPair* cont = 0 ;
   PEL_Object* old = 0 ;
   
   if( object != 0 )
   {
      cont = static_cast<PEL_KeyItemPair*>( object ) ;
      old = cont->item() ;
      cont->set_item( a_item ) ;
   }
   else
   {
      cont = PEL_KeyItemPair::create( hList, key, a_item ) ;
      hList->extend( cont ) ;
      old = a_item ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_at( key ) == a_item ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_Map:: count( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: count" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t resu = hList->count() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( count_POST( resu ) ) ;

   return( resu ) ;
}



//-------------------------------------------------------------------------
bool
PEL_Map:: has_key( PEL_Object const* key ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: has_key" ) ;
   PEL_CHECK_PRE( key != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( item_at( key ) != 0 ) ;
}



//-------------------------------------------------------------------------
PEL_Object*
PEL_Map:: item( PEL_Object const* object ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: item" ) ;
   PEL_Error::object()->raise_not_implemented( this, "item" ) ;

   return( 0 ) ;
}



//-------------------------------------------------------------------------
PEL_Object* 
PEL_Map:: item_at( PEL_Object const* key ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: item_at" ) ;
   PEL_CHECK_PRE( key!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Object* result = 0 ;

   PEL_Object* pair = hList->item( key ) ;
   if( pair != 0 )
   {
      result = static_cast<PEL_KeyItemPair*>( pair )->item() ;
   }
   PEL_CHECK( pair==0 || result!=0 ) ;

   PEL_CHECK_POST( result==0 || has_key( key ) ) ; //?????????
   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_MapIterator*
PEL_Map:: create_iterator( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: create_iterator" ) ;
   return( PEL_MapIterator::create( a_owner, this ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_Map:: remove_at( PEL_Object const* key )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: remove_at" ) ;
   PEL_CHECK_PRE( key!=0 ) ;
   PEL_CHECK_PRE( has_key( key ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   hList->remove( key ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !has_key( key ) ) ; 
}



//-------------------------------------------------------------------------
void
PEL_Map:: clear( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Map:: clear" ) ;
   hList->clear() ;

   PEL_CHECK_POST( count() == 0 ) ;
}



//-------------------------------------------------------------------------
bool
PEL_Map:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( hList!=0 ) ;

   return( true ) ;
} 
