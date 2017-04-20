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

#include <PEL_List.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL_ListItem.hh>
#include <PEL_ListIterator.hh>

using std::ostream ;
using std::endl ;


//----------------------------------------------------------------------
PEL_List*
PEL_List:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: create" ) ;
   PEL_List* result = new PEL_List( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->index_limit() == 0 ) ;

   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_List:: PEL_List( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_Sequence( a_owner ),
     theList(0), 
     theLast(0), 
     nbElem( 0 )
{
   PEL_LABEL( "PEL_List:: PEL_List" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_List:: ~PEL_List( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: ~PEL_List" ) ;
   PEL_CHECK_INV( invariant() ) ;

   clear() ;
}



//-------------------------------------------------------------------------
PEL_List*
PEL_List:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_List* result = create( a_owner ) ;
   PEL_ListItem* ret = theList ;         
   while( ret!=0 )
   {
      result->append( ret->val() ) ;
      ret = ret->next() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;  

   return( result ) ;
}



//-------------------------------------------------------------------------
void      
PEL_List:: copy( PEL_List const* other )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: copy" ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   PEL_CHECK_PRE( other !=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( this != other )
   {
      clear() ;
      PEL_ListItem* ret = other->theList ;
         
      while( ret!=0 )
      {
         append( ret->val() ) ;
         ret = ret->next() ;
      }
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( state_id() != OLD(state_id) ) ;
   PEL_CHECK_POST( index_limit() == other->index_limit() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == other->at(i) ) ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_List:: index_limit( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: index_limit" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( nbElem ) ;
}



//-------------------------------------------------------------------------
void
PEL_List:: append( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: append" ) ;
   PEL_CHECK_PRE( append_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* theNew = new PEL_ListItem( object ) ;
   
   if( nbElem==0 )
   {
      theList = theNew ;
   }
   else
   {
      theLast->doLink( theNew ) ;
   }
   theLast = theNew ;
   nbElem++ ;
   report_state_change() ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( append_POST( OLD( index_limit ),
                                OLD( count ),
                                object,
                                OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_List:: prepend( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: prepend" ) ;
   PEL_CHECK_PRE( prepend_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* theNew = new PEL_ListItem( object ) ;
   PEL_ListItem* old = theList ;

   theList = theNew ;

   if( old==0 )
   {
      theLast = theNew ;
   }
   else
   {
      theNew->doLink( old ) ;
   }
   
   nbElem++ ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( prepend_POST( OLD( index_limit ),
                                 OLD( count ),
                                 object,
                                 OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_List:: set_at( size_t i, PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: set_at" ) ;
   PEL_CHECK_PRE( set_at_PRE( i, object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* it = theItem( i ) ;
   it->replaceValue(object) ;
   report_state_change() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_at_POST( OLD( index_limit ),
                                OLD( count ),
                                i,
                                object,
                                OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_List:: insert_at( size_t i, PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: insert_at" ) ;
   PEL_CHECK_PRE( insert_at_PRE( i, object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* theNew = new PEL_ListItem( object ) ;
   PEL_ListItem* old ;
   
   if( i==0 )
   {
      old = theList ;
      theList = theNew ;
   }
   else
   {
      PEL_ListItem* ptr = theItem( i-1 ) ;
      old = ptr->next() ;
      ptr->doLink( theNew ) ;
   }
   theNew->doLink( old ) ;
   if( old==0 )
   {
      theLast=theNew ;
   }
   nbElem++ ;
   report_state_change() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( insert_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   i,
                                   object,
                                   OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_List:: count( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: count" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( count_POST( nbElem ) ) ;
   return( nbElem ) ;
}



//-------------------------------------------------------------------------
PEL_Object*
PEL_List:: item( PEL_Object const* object ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: item" ) ;
   PEL_CHECK_PRE( item_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_ListItem* ptr = theList ;
   PEL_Object* result = 0 ;
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         result = ptr->val() ;
         break ;
      }
      ptr = ptr->next() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( result, object ) ) ;

   return( result ) ; 
}



//-------------------------------------------------------------------------
PEL_Object*
PEL_List:: at( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: at" ) ;
   PEL_CHECK_PRE( at_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Object* result = theItem(i)->val() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( at_POST( result , i ) ) ;

   return( result ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_List:: index_of( PEL_Object const* object ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: index_of" ) ;
   PEL_CHECK_PRE( index_of_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_ListItem* ptr = theList ;
   size_t ret = badIndex ;
   size_t cpt = 0 ;
   
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         ret = cpt ;
         break ;
      }
      ptr = ptr->next() ;
      cpt++ ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( index_of_POST( ret, object ) ) ;

   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_ListIterator*
PEL_List:: create_iterator( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: create_iterator" ) ;
   PEL_ListIterator* result = PEL_ListIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
PEL_List:: remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: remove" ) ;
   PEL_CHECK_PRE( remove_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* ptr = theList ;
   PEL_ListItem* old = 0 ;
   PEL_Object * ret = 0 ;
   
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         PEL_ListItem* nextItem = ptr->next() ;
         if( old==0 )
         {
            theList = nextItem ;
            if( theLast==ptr )
            {
               theLast=0 ;
            }
         }
         else
         {
            old->doLink( nextItem ) ;
            if( theLast==ptr )
            {
               theLast=old ;
            }
         }
         ret = ptr->val() ;
         delete ptr ;
         break ;
      }
      old = ptr ;
      ptr = ptr->next() ;
   }
   nbElem-- ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_POST( OLD( count ), object, OLD( state_id ) ) ) ;

//   return( ret ) ;
}



//-------------------------------------------------------------------------
void 
PEL_List:: remove_at( size_t i )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: remove_at" ) ;
   PEL_CHECK_PRE( remove_at_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   PEL_ListItem* old ;
   
   if( i==0 )
   {
      old=theList ;
      theList = old->next() ;
      if( old == theLast )
      {
         theLast = 0 ;
      }
   }
   else
   {
      PEL_ListItem* ptr = theItem( i-1 ) ;
      old = ptr->next() ;
      ptr->doLink( old->next() ) ;
      if( old == theLast )
      {
         theLast = ptr ;
      }
   }
   
   delete old ;
   nbElem-- ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   OLD( state_id ) ) ) ;
}



//---------------------------------------------------------------------
void
PEL_List:: remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: remove_section" ) ;
   PEL_CHECK_PRE( remove_section_PRE( iFirst, length ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
      
   PEL_ListItem* old = theItem( iFirst+length-1 ) ;
   PEL_ListItem* end = old->next() ;
   PEL_ListItem* firstToDelete ;
   
   if( iFirst==0 )
   {
      firstToDelete = theList ;
      theList = end ;
      if( old == theLast )
      {
         theLast = 0 ;
      }
   }
   else
   {
      PEL_ListItem* prev = theItem( iFirst-1 ) ;
      firstToDelete = prev->next() ;
      prev->doLink( end ) ;
      if( old == theLast )
      {
         theLast = prev ;
      }

   }
   for( size_t cpt = 0 ; cpt<length ; cpt++ )
   {
      PEL_ListItem* ptr = firstToDelete->next() ;
      delete firstToDelete ;
      firstToDelete = ptr ;
   }
   
   PEL_CHECK( end == firstToDelete ) ;
   
   nbElem -= length ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}



//---------------------------------------------------------------------
void
PEL_List:: destroy_items_and_remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: destroy_items_and_remove_section" ) ;
   PEL_CHECK_PRE( destroy_items_and_remove_section_PRE( iFirst, length ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, index_limit, index_limit() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* deb = theItem( iFirst ) ;
   for( size_t cpt = 0 ; cpt<length ; cpt++ )
   {
      PEL_ListItem* nextItem = deb->next() ;
      deb->val()->destroy() ;
      deb=nextItem ;
   }
   remove_section( iFirst, length ) ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
   PEL_CHECK_POST( index_limit() == OLD( index_limit ) - length ) ;
   PEL_CHECK_POST( count() == OLD( count ) - length ) ;
}



//----------------------------------------------------------------------
void
PEL_List:: destroy_items_and_clear( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: destroy_items_and_clear" ) ;
   PEL_CHECK_PRE( destroy_items_and_remove_section_PRE( 0, index_limit() ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* ptr = theList ;
   while( ptr!=0 )
   {
      PEL_ListItem* nextItem = ptr->next() ;
      ptr->val()->destroy() ;
      ptr = nextItem ;
   }
   clear() ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;

}



//-------------------------------------------------------------------------
void
PEL_List:: clear( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_List:: clear" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   PEL_ListItem* ptr = theList ;
   PEL_ListItem* old ;
   
   while( ptr!=0 )
   {
      old = ptr ;
      ptr = ptr->next() ;
      delete old ;
   }
   theList = 0 ;
   theLast = 0 ;
   nbElem = 0 ;     
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//----------------------------------------------------------------------
// ostream&
// operator<<( ostream& out, PEL_List const& l )
//----------------------------------------------------------------------
// {
//    out << "List : " << endl ;
//    out << "  Nb elms = " << l.si() << endl ;
//    return( out ) ;
// }



//-------------------------------------------------------------------------
PEL_ListItem*
PEL_List:: theItem( size_t n ) const
//-------------------------------------------------------------------------
{
   size_t nb = n ;
   PEL_CHECK( n<nbElem ) ; //??????????????????????????
   PEL_ListItem* ret = theList ;
         
   while( nb>0 )
   {
      nb-- ;
      ret=ret->next() ;
   }
   return( ret ) ;
}



//-------------------------------------------------------------------------
bool
PEL_List:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::invariant() ) ;
   PEL_ASSERT( count() == index_limit() ) ;

   return( true ) ;
}



//------------------------------------------------------------------------
bool
PEL_List:: set_at_POST( size_t old_index_limit,
                        size_t old_count,
                        size_t i,
                        PEL_Object const* object,
                        size_t old_state_id ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::set_at_POST( old_index_limit, 
                                          old_count, 
                                          i, 
                                          object,
                                          old_state_id )  ) ;

   PEL_ASSERT( old_count == count() ) ;
   PEL_ASSERT( state_id() != old_state_id ) ;
   
   return( true ) ;
}



//-------------------------------------------------------------------------
bool
PEL_List:: at_POST( PEL_Object const* result, size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::at_POST( result, i ) ) ;
   PEL_ASSERT( result != 0 ) ;

   return( true ) ;
}


//----------------------------------------------------------------------
bool
PEL_List:: remove_at_POST( size_t old_index_limit,
                           size_t old_count,
                           size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::remove_at_POST( old_index_limit,
                                             old_count,
                                             old_state_id  ) ) ;

   PEL_ASSERT( index_limit() == old_index_limit - 1 ) ;
   PEL_ASSERT( count() == old_count - 1 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_List:: remove_section_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t length,
                                size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Sequence::remove_section_POST( old_index_limit, 
                                                  old_count, 
                                                  length,
                                                  old_state_id ) ) ;

   PEL_ASSERT( index_limit() == old_index_limit - length ) ;
   PEL_ASSERT( count() == old_count - length ) ;

   return( true ) ;
}


