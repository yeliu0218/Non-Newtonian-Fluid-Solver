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

#include <PEL_HashTableSetIterator.hh>

#include <PEL_HashTableSet.hh>
#include <PEL_ListIterator.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_HashTableSetIterator*
PEL_HashTableSetIterator:: create( PEL_Object* a_owner,
                                   PEL_HashTableSet const* a_table  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: create" ) ;

   PEL_HashTableSetIterator* result = 
                               new PEL_HashTableSetIterator( a_owner, a_table ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( EQUIVALENT( a_table->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_HashTableSetIterator:: PEL_HashTableSetIterator( PEL_Object* a_owner,
                                             PEL_HashTableSet const* hTable )
//----------------------------------------------------------------------
   : PEL_Iterator( a_owner, hTable ), 
     table( hTable ), 
     where( 0 ),
     init( true ),
     bucket( 0 )
{
   PEL_LABEL( "PEL_HashTableSetIterator:: PEL_HashTableSetIterator" ) ;
   start() ;
}



//----------------------------------------------------------------------
PEL_HashTableSetIterator:: ~PEL_HashTableSetIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: ~PEL_HashTableSetIterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( where != 0 )
   {
      where->destroy() ;
      where = 0 ;
   }
}



//----------------------------------------------------------------------
bool
PEL_HashTableSetIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: is_valid" ) ;
   PEL_CHECK_PRE( is_valid_PRE() ) ;
   return( where!=0 ) ;
}



//----------------------------------------------------------------------
void
PEL_HashTableSetIterator:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: start" ) ;
   if( where!=0 )
   {
      where->destroy() ;
   }
   where = 0 ;
   
   for( bucket = 0 ; bucket<table->the_size ; bucket++ )
   {
      PEL_List const* list = table->getList( bucket ) ;
      if( list->count() > 0 )
      {
         where = list->create_iterator( 0 ) ;
         PEL_CHECK( where->is_valid() ) ;
         break ;
      }
   }
   record_container_state_id() ;
   PEL_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
PEL_HashTableSetIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   
   if( where!=0 )
   {
      where->go_next() ;
      if( !where->is_valid() )
      {
         where->destroy() ;
         where = 0 ;
         for( bucket++ ; bucket<table->the_size ; bucket++ )
         {
            PEL_List const* list = table->getList( bucket ) ;
            if( list->count() > 0 )
            {
               where = list->create_iterator( 0 ) ;
               PEL_CHECK( where->is_valid() ) ;
               break ;
            }
         }
      }
   }
}



//----------------------------------------------------------------------
PEL_Object*
PEL_HashTableSetIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_HashTableSetIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Object* resu = where->item() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( resu ) ) ;   

   return( resu ) ;
}



