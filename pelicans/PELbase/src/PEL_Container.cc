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

#include <PEL_Container.hh>

#include <iostream>

#include <PEL_Iterator.hh>
#include <PEL_assertions.hh>

using std::endl ;
using std::ostream ;

//-----------------------------------------------------------------------------
PEL_Container:: PEL_Container( PEL_Object* a_owner )
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , MY_STATE( 0 )
{
   PEL_LABEL( "PEL_Container:: PEL_Container" ) ;
}

//----------------------------------------------------------------------
PEL_Container:: ~PEL_Container( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Container:: ~PEL_Container" ) ;
   report_state_change() ;
}

//-------------------------------------------------------------------------
bool
PEL_Container:: matching_items( PEL_Object const* object1,
                                PEL_Object const* object2 ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Container:: matching_items" ) ;
   PEL_CHECK_PRE( matching_items_PRE( object1, object2 ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( object1->is_equal( object2 ) ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------
size_t
PEL_Container:: state_id( void ) const
//-----------------------------------------------------------------------
{
   return MY_STATE ;
}

//-------------------------------------------------------------------------
bool
PEL_Container:: has( PEL_Object const* object ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Container:: has" ) ;
   PEL_CHECK_PRE( object != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   bool result = ( item( object ) != 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0 , count()!=0 ) ) ;
   PEL_CHECK_POST( OLD( state_id )==state_id() ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PEL_Container:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Container:: print" ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   std::string const s( indent_width, ' ' ) ;
   os << s << "Collection of type " << type_name() ;
   if( count()==1 )
   {
      os << " contains " << count() << " item :" << endl ;
      PEL_Iterator* it = create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         PEL_Object const* obj = it->item() ;
         obj->display_info( os, indent_width+3 ) ;
         obj->print( os, indent_width+3 ) ;
         os << endl ;
      }
      it->destroy() ;
   }
   else if( count()!=0 )
   {
      os << " contains " << count() << " items :" << endl ;
      PEL_Iterator* it = create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         PEL_Object const* obj = it->item() ;
         os << s << " - " << endl ;
         obj->display_info( os, indent_width+3 ) ;
         obj->print( os, indent_width+3 ) ;
         os << endl ;
         os << s << "    -----------------------------" << endl ;
      }
      it->destroy() ;
   }
   else
   {
      os << " contains " << count() << " item" << endl ;
   }
   
   PEL_CHECK_POST( OLD( state_id )==state_id() ) ;

}

//----------------------------------------------------------------------
bool
PEL_Container:: create_clone_POST( PEL_Container* result,
                                   PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{ 
   PEL_ASSERT( PEL_Object::create_clone_POST( result, a_owner ) ) ;
   PEL_ASSERT( result->count() == count() ) ;
   PEL_Iterator* it1 = 0 ;
   PEL_Iterator* it2 = 0 ;
   PEL_ASSERT( 
      FORALL( ( it1 = create_iterator( 0 ),
                it2 = result->create_iterator( 0 ) ;
                it1->is_valid() && it2->is_valid() ;
                it1->go_next(), it2->go_next() ),
                it1->item()->has_same_address( it2->item() ) ) ) ;
   
   it1->destroy() ;
   it2->destroy() ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Container:: count_POST( size_t result ) const
//----------------------------------------------------------------------
{
   PEL_Iterator* it = create_iterator(0) ;
   size_t val = 0 ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      if( it->item()!=0 )
      {
         val++ ;
      }
   }
   it->destroy() ;
   PEL_ASSERT( result==val ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Container:: matching_items_PRE( PEL_Object const* object1,
                                    PEL_Object const* object2 ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( (object1 != 0) && (object2 != 0) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Container:: item_PRE( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( object != 0 ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Container:: item_POST( PEL_Object const* result,
                           PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result==0 || 
               ( count() != 0 &&  matching_items( result, object ) ) ) ;

  
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
PEL_Container:: create_iterator_POST( PEL_Iterator* result,
                                      PEL_Object* a_owner ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( EQUIVALENT( count() != 0 , result->is_valid() ) ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------
void
PEL_Container:: report_state_change( void ) 
//-----------------------------------------------------------------------
{
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   MY_STATE ++ ;
   PEL_CHECK_POST( OLD(state_id)!=state_id() ) ;
}

