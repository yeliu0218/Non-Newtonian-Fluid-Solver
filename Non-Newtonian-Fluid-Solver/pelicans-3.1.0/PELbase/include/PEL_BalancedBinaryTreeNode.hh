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

#ifndef PEL_BALANCED_TREE_NODE_HH
#define PEL_BALANCED_TREE_NODE_HH

#include <PEL_Object.hh>

// Balanced binary tree node class.
// Node realizing a binary tree with balance behaviour.
// This class is used only by PEL_BalancedBinaryTree class.

// A tree node is defined by :
// . a pointer to a PEL_Object
// . the left and right branches such that all elements contained in the
//    left sub-tree are lower than the object of the node and thoses in
//    right sub-tree are larger ( in the sense of three_way_comparison method ).
// . the height (0 for a leap) of the node.
//
// The binary tree is said balanced because it is assumed that left height
// and right one don't differ for more that 2.
// The balance allows then fast inserting and searching.

class PEL_EXPORT PEL_BalancedBinaryTreeNode : public PEL_Object
{
   public: //---------------------------------------------------------------
   protected: //-----------------------------------------------------------

      /* Contructs a leap associated with PEL_Object object. */
      //??????????? mettre le create qui va bien
      PEL_BalancedBinaryTreeNode( PEL_Object * object ) ;
      
      /* Delete self. */ //??????? private.....
      virtual ~PEL_BalancedBinaryTreeNode( void ) ;

      /* Searches in the sub-tree for an item matching value */ 
      PEL_Object * search( const PEL_Object * value ) const ;

      /* Inserts a new item in the sub-tree */
      PEL_Object * addChild( PEL_Object * value,
                             PEL_BalancedBinaryTreeNode *& top ) ;

      /* Delete in the sub-tree the item matching value
         top is the pointer referencing self. */
      PEL_BalancedBinaryTreeNode * remove( const PEL_Object * value,
                                           PEL_BalancedBinaryTreeNode *& top ) ;

      /* Returns cardinal */
      size_t getNbElem( void ) const ;
      
      friend class PEL_BalancedBinaryTree ; //????????? trop de friend
      friend class PEL_BalancedBinaryTreeIterator ;
      
      virtual bool invariant( void ) const ;
     
   private: //--------------------------------------------------------
      
      /* Delete sub-tree */
      void clear( void ) ;
      
      /* Balances sub-tree */ 
      void balance( PEL_BalancedBinaryTreeNode *& top ) ;

      /* Copy items of node in sub-tree */
      void copySubTree( PEL_BalancedBinaryTreeNode * node,
                        PEL_BalancedBinaryTreeNode *& top ) ;

      /* Estimates current height */
      void calculHeight( void ) ;
      
      //--------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------

      size_t height ;
      PEL_Object * val ;
      PEL_BalancedBinaryTreeNode * less ;
      PEL_BalancedBinaryTreeNode * more ;
      
};

#endif // PEL_BALANCED_TREE_NODE_HH
