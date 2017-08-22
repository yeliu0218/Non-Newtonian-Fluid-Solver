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

#include <PEL_BalancedBinaryTreeIterator.hh>

#include <PEL_BalancedBinaryTree.hh>
#include <PEL_BalancedBinaryTreeNode.hh>
#include <PEL_List.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_BalancedBinaryTreeIterator*
PEL_BalancedBinaryTreeIterator:: create( PEL_Object* a_owner,
                                         PEL_BalancedBinaryTree const* a_tree )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: create" ) ;
   PEL_BalancedBinaryTreeIterator* result = 
                         new PEL_BalancedBinaryTreeIterator( a_owner, a_tree ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( EQUIVALENT( a_tree->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_BalancedBinaryTreeIterator:: PEL_BalancedBinaryTreeIterator(
                                          PEL_Object* a_owner,
                                          PEL_BalancedBinaryTree const* aTree )
//----------------------------------------------------------------------
   : PEL_Iterator( a_owner, aTree ), 
     tree( aTree ),
     valid( false )
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: PEL_BalancedBinaryTreeIterator" ) ;
   pathList = PEL_List::create( this ) ;
   reset() ;
   valid = next() ;
}



//----------------------------------------------------------------------
PEL_BalancedBinaryTreeIterator:: ~PEL_BalancedBinaryTreeIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: ~PEL_BalancedBinaryTreeIterator" ) ;
}



//----------------------------------------------------------------------
bool
PEL_BalancedBinaryTreeIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: is_valid" ) ;
   PEL_CHECK_PRE( is_valid_PRE() ) ;
   return( valid ) ;
}



//----------------------------------------------------------------------
void
PEL_BalancedBinaryTreeIterator:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: start" ) ;
   reset() ;
   valid = next() ;
   record_container_state_id() ;
   PEL_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
PEL_BalancedBinaryTreeIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   valid = next() ;
}



//----------------------------------------------------------------------
bool
PEL_BalancedBinaryTreeIterator:: next( void ) 
//----------------------------------------------------------------------
{
   bool ret = false ;
   if( where )
   {
      PEL_BalancedBinaryTreeNode * less = where->less ;
      PEL_BalancedBinaryTreeNode * more = where->more ;
      
      if( state==LEFT_TO_DO && less!=0 )
      {
         pathList->append( where ) ;
         where = less ;
         ret = next() ;
      }
      else if( ( state==LEFT_TO_DO && less==0 )
               || state == CURRENT_TO_DO )
      {
         ret = true ;
         state = RIGHT_TO_DO ;
      }
      else if( state == RIGHT_TO_DO && more!=0 )
      {
         state = LEFT_TO_DO ;
         pathList->append( where ) ;
         where = more ;
         ret = next() ;
      }
      else  // state == RIGHT_TO_DO && more==0 || state == FINISHED
      {
         if( pathList->count() != 0 )
         {
            PEL_BalancedBinaryTreeNode * up = 
            static_cast<PEL_BalancedBinaryTreeNode * >( pathList->at(pathList->count()-1) ) ;
            pathList->remove_at( pathList->count()-1 ) ;
//???????????????????
//            PEL_BalancedBinaryTreeNode * up =
//               static_cast<PEL_BalancedBinaryTreeNode * >( pathList->remove_at( pathList->count()-1) ) ;
            if( up->more == where )
            {
               where = up ;
               state = FINISHED ;
            }
            else
            {
               PEL_ASSERT( up->less == where ) ;
               where = up ;
               state = CURRENT_TO_DO ;
            }
            ret = next() ;
         }
         else
         {
            where = 0 ;
         }
      }
   }
   return ret ;
}



//----------------------------------------------------------------------
void
PEL_BalancedBinaryTreeIterator:: reset( void )
//----------------------------------------------------------------------
{
   init = true ;
   pathList->clear() ;
   where = tree->top ;
   state = LEFT_TO_DO ;
}


//----------------------------------------------------------------------
PEL_Object*
PEL_BalancedBinaryTreeIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BalancedBinaryTreeIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;

   PEL_Object* result = where->val ;

   PEL_CHECK_POST( item_POST( result ) ) ;
   return( result ) ;
}



