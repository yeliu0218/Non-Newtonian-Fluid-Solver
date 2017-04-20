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

#include <PEL_BalancedBinaryTree.hh>

#include <PEL_BalancedBinaryTreeNode.hh>
#include <PEL_Error.hh>
#include <PEL_assertions.hh>


//-------------------------------------------------------------------------
PEL_BalancedBinaryTree*
PEL_BalancedBinaryTree:: create( PEL_Object* a_owner ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: create" ) ;
   PEL_BalancedBinaryTree* result = new PEL_BalancedBinaryTree( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTree:: PEL_BalancedBinaryTree( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_Set( a_owner ), 
     top( 0 ), 
     nbEntries( 0 )
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: PEL_BalancedBinaryTree" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTree:: ~PEL_BalancedBinaryTree( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: ~PEL_BalancedBinaryTree" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( top!=0 )
   {
      top->destroy() ;
   }
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTree*
PEL_BalancedBinaryTree:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_BalancedBinaryTree* result =  PEL_BalancedBinaryTree::create( a_owner ) ;
   PEL_Iterator* it = create_iterator( 0 ) ;
   for( ; it->is_valid() ; it->go_next() )
   {         
      result->extend( it->item() ) ;
   }
   it->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;

   return( result ) ; 
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTree:: extend( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: extend" ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   PEL_CHECK_PRE( extend_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count() ) ;
   PEL_SAVEOLD( bool, has, has( object ) ) ;

   PEL_Object * ret ;   
   
   if( top==0 )
   {
      top = new PEL_BalancedBinaryTreeNode( object ) ; //?????????????????
      ret = object ;
   }
   else
   {
      ret = top->addChild( object, top ) ;
   }
   if( ret==object )
   {
       nbEntries++ ;
   }
   report_state_change() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
//   return( ret ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_BalancedBinaryTree:: count( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: count" ) ;
   PEL_CHECK_POST( count_POST( nbEntries ) ) ;
   return( nbEntries ) ;
}



//-------------------------------------------------------------------------
PEL_Object*
PEL_BalancedBinaryTree:: item( PEL_Object const* object ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: item" ) ;
   PEL_CHECK_PRE( item_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Object* ret = 0 ;
   if( top!=0 )
   {
      ret = top->search( object ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( ret, object ) ) ;

   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTreeIterator*
PEL_BalancedBinaryTree:: create_iterator( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: create_iterator" ) ;
   PEL_BalancedBinaryTreeIterator* result = 
                     PEL_BalancedBinaryTreeIterator::create( a_owner, this ) ;

   PEL_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTree:: remove( PEL_Object const* object ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: remove" ) ;
   PEL_CHECK_PRE( remove_PRE( object ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, count, count()) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;
   
   PEL_BalancedBinaryTreeNode* node = top->remove( object, top ) ;
   nbEntries-- ;
   node->destroy() ;
   report_state_change() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( remove_POST( OLD( count ), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTree:: destroy_item_and_remove( PEL_Object const* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: destroy_item_and_remove" ) ;
   PEL_Error::object()->raise_not_implemented( this, 
                                               "destroy_item_and_remove" ) ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTree:: clear( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: clear" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, state_id, state_id() ) ;

   if( top != 0 ) top->destroy() ;
   top = 0 ;
   nbEntries = 0 ;
   report_state_change() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( clear_POST( OLD(state_id) ) ) ;   
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTree:: destroy_items_and_clear( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTree:: destroy_items_and_clear" ) ;
   PEL_Error::object()->raise_not_implemented( this, 
                                               "destroy_items_and_clear" ) ;
}



//-------------------------------------------------------------------------
bool
PEL_BalancedBinaryTree:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Collection::invariant() ) ; //??????????????? PEL_Set
   size_t cpt = 0 ;
   
   if( top!=0 )
   {
      cpt = top->getNbElem() ;
   }
   
   PEL_ASSERT( cpt==nbEntries ) ;
   return( true ) ;
}



