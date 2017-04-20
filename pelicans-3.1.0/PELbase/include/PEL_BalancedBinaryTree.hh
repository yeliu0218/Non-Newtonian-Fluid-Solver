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

#ifndef PEL_BALANCED_BINARY_TREE_HH
#define PEL_BALANCED_BINARY_TREE_HH

#include <PEL_Set.hh>

#include <PEL_BalancedBinaryTreeIterator.hh>

//--------------------------------------------------------------------------
// Sorted sets, implemented as balanced binary search trees, to allow
// fast inserting searching operations
//--------------------------------------------------------------------------
// For each item, three pointers and a size_t are stored.
//--------------------------------------------------------------------------

class PEL_BalancedBinaryTreeNode ;

class PEL_EXPORT PEL_BalancedBinaryTree : public PEL_Set
{
      
   public : //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_BalancedBinaryTree* create( PEL_Object* a_owner ) ;
      
      virtual PEL_BalancedBinaryTree* create_clone( 
                                             PEL_Object* a_owner ) const ;
      
   //-- Measurement   

      virtual size_t count( void ) const ;

   //-- Access

      virtual PEL_Object* item( PEL_Object const* object ) const ;

      virtual PEL_BalancedBinaryTreeIterator* create_iterator(
                                             PEL_Object* a_owner ) const ;
      
   //-- Element change

      virtual void extend( PEL_Object* object ) ;

      
   //-- Removal
      
      virtual void remove( PEL_Object const* object ) ;

      virtual void destroy_item_and_remove( PEL_Object const* object ) ;

      virtual void clear( void ) ;

      virtual void destroy_items_and_clear( void ) ;
     
   protected ://------------------------------------------------------------
      
      virtual ~PEL_BalancedBinaryTree( void ) ;

      PEL_BalancedBinaryTree( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
   
      virtual bool invariant( void ) const ;

   private :  //------------------------------------------------------------

      PEL_BalancedBinaryTree( void ) ;
      PEL_BalancedBinaryTree( PEL_BalancedBinaryTree const& other ) ;
      PEL_BalancedBinaryTree const& operator=( 
                              PEL_BalancedBinaryTree const& other ) ;

      friend class PEL_BalancedBinaryTreeIterator ;
      
   //-- Attributes

      PEL_BalancedBinaryTreeNode* top ;
      size_t nbEntries ;
      
};

#endif
