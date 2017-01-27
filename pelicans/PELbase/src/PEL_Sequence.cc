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

#include <PEL_Sequence.hh>

#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_assertions.hh>

using std::string ;

//----------------------------------------------------------------------
size_t const PEL_Sequence:: badIndex = static_cast<size_t>(~0) ;
//----------------------------------------------------------------------



//-----------------------------------------------------------------------------
PEL_Sequence:: PEL_Sequence( PEL_Object* a_owner )
//-----------------------------------------------------------------------------
   : PEL_Collection( a_owner )
{
   PEL_LABEL( "PEL_Sequence:: PEL_Sequence" ) ;
}



//-----------------------------------------------------------------------------
PEL_Sequence:: ~PEL_Sequence( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: ~PEL_Sequence" ) ;
}



//-------------------------------------------------------------------------
bool
PEL_Sequence:: comparable( PEL_Object const* other ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: comparable" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   // these preconditions are less restrictive than PEL_Object::
   return PEL_Object::comparable( other ) ||
      dynamic_cast<PEL_Sequence const*>( other ) != 0  ;
}


//-------------------------------------------------------------------------
bool
PEL_Sequence:: is_equal( PEL_Object const* other ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: is_equal" ) ;
   
   // these preconditions are less restrictive than PEL_Object::is_equal_PRE()
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = ( three_way_comparison( other ) == 0 ) ;

   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
int
PEL_Sequence:: three_way_comparison( PEL_Object const* other ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: three_way_comparison" ) ;
   
   PEL_CHECK_PRE( three_way_comparison_PRE( other )) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Sequence const* otherList = static_cast<PEL_Sequence const*>( other ) ;

   PEL_Iterator* itSelf = create_iterator( 0 ) ;
   PEL_Iterator* itOther = otherList->create_iterator( 0 ) ;
   int result = 0 ;
   
   while( itSelf->is_valid() && itOther->is_valid() )
   {
      int pRet = itSelf->item()->three_way_comparison( itOther->item() ) ;
      if( pRet!=0 )
      {
         result = pRet ;
         break ;
      }
      itSelf->go_next() ;
      itOther->go_next() ;
   }
   
   if( result==0 )
   {
      if( itSelf->is_valid() && !itOther->is_valid() )
      {
         result = 1 ;
      }
      else if( !itSelf->is_valid() && itOther->is_valid() )
      {
         result = -1 ;
      }
   }

   itSelf->destroy() ;
   itOther->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_Sequence:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: hash_code" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t result = 0 ;
   PEL_Iterator* itSelf = create_iterator( 0 ) ;
   for( ; itSelf->is_valid() ; itSelf->go_next() )
   {
      result += itSelf->item()->hash_code() ;
   }
   itSelf->destroy() ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
PEL_Sequence:: remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: remove" ) ;
   
   PEL_CHECK_PRE( remove_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t id = index_of( object ) ;
   PEL_ASSERT( id != badIndex ) ;
   remove_at( id ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_POST( OLD( count ), object, OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_Sequence:: sort( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: sort" ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( count()>0 )
   {
      PEL_BalancedBinaryTree* bin=PEL_BalancedBinaryTree::create(0) ;
      PEL_Iterator* it=create_iterator(0) ;
      for(it->start();it->is_valid();it->go_next())
      {
         bin->extend( it->item() ) ;
      }
      clear() ;
      it->destroy() ;
      it=bin->create_iterator(0) ;
      for(it->start();it->is_valid();it->go_next())
      {
         append( it->item() ) ;
      }
      bin->destroy() ;
      it->destroy() ;
   }
   
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( count()==OLD(count) ) ;
   PEL_CHECK_POST( state_id()!=OLD(state_id) ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=1 ; i<count() ; ++i ),
              at(i-1)->three_way_comparison(at(i))<=0 ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_Sequence:: extend( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: extend" ) ;
   
   PEL_CHECK_PRE( extend_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( bool, has, has( object ) ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( !has( object ) )
   {
      append( object ) ;
   }
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_Sequence:: destroy_item_and_remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: destroy_item_and_remove" ) ;
   
   PEL_CHECK_PRE( object!=0 ) ;
   PEL_CHECK_PRE( has( object ) ) ;
   PEL_CHECK_PRE( object->owner()==0 ) ;
   
   //???????????????????? temps calcul
   PEL_Object* obj = item( object ) ;
   remove( object ) ;
   obj->destroy() ;
}



//---------------------------------------------------------------------
void
PEL_Sequence:: destroy_item_and_remove_at( size_t i )
//---------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: destroy_item_and_remove_at" ) ;
   
   PEL_CHECK_PRE( i<index_limit() ) ;
   PEL_CHECK_PRE( at( i )!=0 ) ;
   PEL_CHECK_PRE( at( i )->owner()==0 ) ;
   
   //?????????????????? temps calcul
   PEL_Object* obj = at( i ) ;
   PEL_CHECK( obj!=0 ) ;
   obj->destroy() ;
   remove_at( i ) ;

   PEL_CHECK_POST( at(i)==0 ) ;
}



//----------------------------------------------------------------------
void
PEL_Sequence:: destroy_items_and_clear( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: destroy_items_and_clear" ) ;
   
   PEL_CHECK_PRE( destroy_items_and_remove_section_PRE( 0, index_limit() ) ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   destroy_items_and_remove_section( 0, index_limit() ) ;
   clear() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: invariant" ) ;
   
   PEL_ASSERT( PEL_Collection::invariant() ) ;
   PEL_ASSERT( count() <= index_limit() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: create_clone_POST( PEL_Sequence* result,
                                  PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: create_clone_POST" ) ;
   
   PEL_ASSERT( PEL_Container::create_clone_POST( result, a_owner ) ) ;
   PEL_ASSERT( result->index_limit() == index_limit() ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                       result->at(i) == at(i) ) ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: append_PRE( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: append_PRE" ) ;
   
   PEL_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: append_POST( size_t old_index_limit,
                            size_t old_count,
                            PEL_Object const* object,
                            size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: append_POST" ) ;

//   PEL_ASSERT( has( object ) ) ;
   PEL_ASSERT( count() == old_count+1 ) ;
   PEL_ASSERT( index_limit() == old_index_limit+1 ) ;
   PEL_ASSERT( at(index_limit()-1) == object ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;
   
   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: prepend_PRE( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: prepend_PRE" ) ;
   
   PEL_ASSERT( object!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: prepend_POST( size_t old_index_limit,
                             size_t old_count,
                             PEL_Object const* object,
                             size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: prepend_POST" ) ;

//   PEL_ASSERT( has( object ) ) ;
   PEL_ASSERT( count() == old_count+1 ) ;
   PEL_ASSERT( index_limit() == old_index_limit+1 ) ;
   PEL_ASSERT( at(0) == object ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: set_at_PRE( size_t i,
                           PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: set_at_PRE" ) ;
   
   PEL_ASSERT( i < index_limit() ) ;
   PEL_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: set_at_POST( size_t old_index_limit,
                            size_t old_count,
                            size_t i,
                            PEL_Object const* object,
                            size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: set_at_POST" ) ;

//   PEL_ASSERT( has( object ) ) ;
   PEL_ASSERT( old_index_limit == index_limit() ) ;
   PEL_ASSERT( old_count <= count() ) ;
   PEL_ASSERT( at(i) == object ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: insert_at_PRE( size_t i,
                              PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: insert_at_PRE" ) ;
   
   PEL_ASSERT( i<index_limit() ) ;
   PEL_ASSERT( object!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: insert_at_POST( size_t old_index_limit,
                               size_t old_count,
                               size_t i,
                               PEL_Object const* object,
                               size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: insert_at_POST" ) ;

//   PEL_ASSERT( has( object ) ) ;
   PEL_ASSERT( index_limit() == old_index_limit+1 ) ;
   PEL_ASSERT( count() == old_count+1 ) ;
   PEL_ASSERT( at(i) == object ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: remove_at_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: remove_at_PRE" ) ;
   
   PEL_ASSERT( i < index_limit() ) ;
   PEL_ASSERT( at(i)!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: remove_at_POST( size_t old_index_limit,
                               size_t old_count,
                               size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: remove_at_POST" ) ;
   
   PEL_ASSERT( count() <= old_count-1 ) ;
   PEL_ASSERT( index_limit() <= old_index_limit ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: remove_section_PRE( size_t iFirst,
                                   size_t length ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: remove_section_PRE" ) ;
   
   PEL_ASSERT( iFirst+length <= index_limit() ) ;
   PEL_ASSERT( length!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: remove_section_POST( size_t old_index_limit,
                                    size_t old_count,
                                    size_t length,
                                    size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: remove_section_POST" ) ;
   
   PEL_ASSERT( index_limit() <= old_index_limit ) ;
   PEL_ASSERT( count() <= old_count ) ;
   PEL_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: at_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: at_PRE" ) ;
   
   PEL_ASSERT( i < index_limit() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: at_POST( PEL_Object const* resu, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: at_POST" ) ;
   
   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: index_of_PRE( PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: index_of_PRE" ) ;
   
   PEL_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: index_of_POST( size_t result, PEL_Object const* object ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: index_of_POST" ) ;

   PEL_ASSERT( EQUIVALENT( !has(object), result==badIndex ) ) ;
   PEL_ASSERT( result==badIndex || ( result < index_limit() ) ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: destroy_items_and_remove_section_PRE( size_t iFirst,
                                                     size_t length ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: destroy_items_and_remove_section_PRE" ) ;
   
   PEL_ASSERT( iFirst+length<=index_limit() ) ;
   PEL_ASSERT( length!=0 ) ;
   PEL_ASSERT( FORALL( ( size_t i=iFirst ; i<iFirst+length ; ++i ),
                       IMPLIES( at(i)!=0, at(i)->owner()==0 ) ) ) ;
   return( true ) ;
}



//----------------------------------------------------------------------
bool
PEL_Sequence:: clear_POST( size_t old_state_id ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Sequence:: clear_POST" ) ;
   PEL_ASSERT( PEL_Collection::clear_POST( old_state_id ) ) ;
   PEL_ASSERT( count()==0 ) ;
   return( true ) ;
}
