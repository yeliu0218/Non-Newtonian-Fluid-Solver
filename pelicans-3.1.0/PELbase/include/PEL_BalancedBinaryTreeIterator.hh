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

#ifndef PEL_BALANCED_BINARY_TREE_ITERATOR_HH
#define PEL_BALANCED_BINARY_TREE_ITERATOR_HH

#include <PEL_Iterator.hh>

class PEL_BalancedBinaryTree ;
class PEL_BalancedBinaryTreeNode ;
class PEL_List ;

//---------------------------------------------------------------------------
// Iterators on items of PEL_BalancedBinaryTree collections
//---------------------------------------------------------------------------
// Items are accessed from the smallest to the biggest, with respect to
// the total order defined by `three_way_comparison'
//---------------------------------------------------------------------------

class PEL_EXPORT PEL_BalancedBinaryTreeIterator : public PEL_Iterator
{
   public: //---------------------------------------------------------------

   //-- Initialization

      // Create and return an iterator over items of `a_tree'.
      static PEL_BalancedBinaryTreeIterator* create( 
                                       PEL_Object* a_owner,
                                       PEL_BalancedBinaryTree const* a_tree ) ;
       
   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual PEL_Object* item( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;

      virtual void go_next( void ) ;
             
   protected: //------------------------------------------------------------

      virtual ~PEL_BalancedBinaryTreeIterator( void ) ;

      PEL_BalancedBinaryTreeIterator( PEL_Object* a_owner,
                                      PEL_BalancedBinaryTree const* aTree ) ;
      
   private: //--------------------------------------------------------------

      PEL_BalancedBinaryTreeIterator( void ) ;
      PEL_BalancedBinaryTreeIterator( 
                               PEL_BalancedBinaryTreeIterator const& other ) ;
      PEL_BalancedBinaryTreeIterator const& operator=( 
                               PEL_BalancedBinaryTreeIterator const& other ) ;

      bool next( void ) ;
      void reset( void ) ;

      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      PEL_BalancedBinaryTree const* tree ;
      PEL_BalancedBinaryTreeNode* where ;
      enum state_type{ LEFT_TO_DO, CURRENT_TO_DO, RIGHT_TO_DO, FINISHED } ;
      state_type state ;
      PEL_List* pathList ;
      bool init ;
      bool valid ;
} ;

#endif
