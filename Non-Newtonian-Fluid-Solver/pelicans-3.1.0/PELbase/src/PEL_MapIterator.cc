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

#include <PEL_MapIterator.hh>

#include <PEL_Map.hh>
#include <PEL_KeyItemPair.hh>
#include <PEL_HashTableSet.hh>
#include <PEL_ListIterator.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_MapIterator*
PEL_MapIterator:: create( PEL_Object* a_owner,
                          PEL_Map const* map  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: create" ) ;

   PEL_MapIterator* result =  new PEL_MapIterator( a_owner, map ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( EQUIVALENT( map->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_MapIterator:: PEL_MapIterator( PEL_Object* a_owner,
                                   PEL_Map const* hTable )
//----------------------------------------------------------------------
   : PEL_Iterator( a_owner, hTable->hList ), 
     table( hTable ), 
     where( 0 )
{
   PEL_LABEL( "PEL_MapIterator:: PEL_MapIterator" ) ;
   where = hTable->hList->create_iterator( this ) ;
}



//----------------------------------------------------------------------
PEL_MapIterator:: ~PEL_MapIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: ~PEL_MapIterator" ) ;
}



//----------------------------------------------------------------------
bool
PEL_MapIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: is_valid" ) ;
   PEL_CHECK_PRE( is_valid_PRE() ) ;
   
   return( where->is_valid() ) ;
}



//----------------------------------------------------------------------
void
PEL_MapIterator:: start( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: start" ) ;
   
   where->start() ;
   record_container_state_id() ;
   
   PEL_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
PEL_MapIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   
   where->go_next() ;
}



//----------------------------------------------------------------------
PEL_Object*
PEL_MapIterator:: key( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: key" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   PEL_KeyItemPair* cont = static_cast<PEL_KeyItemPair *>( where->item() ) ;
   PEL_Object* result = cont->key() ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_Object*
PEL_MapIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MapIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   PEL_KeyItemPair* cont = static_cast< PEL_KeyItemPair* >( where->item() ) ;
   PEL_Object* result = cont->item() ;

   PEL_CHECK_POST( item_POST( result ) ) ;   
   return( result ) ;
}



