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

#ifndef PEL_NUMBERED_DOUBLE_VECTOR_S_HH
#define PEL_NUMBERED_DOUBLE_VECTOR_S_HH

#include <PEL_Object.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>

class PEL_BalancedBinaryTree ;
class PEL_BalancedBinaryTreeIterator ;
class PEL_DoubleComparator ;
class PEL_NumberedDoubleVectorsItem ;

/*
Sets of numbered instances of doubleVector, each of them
having the same size.

PUBLISHED
*/

class PEL_EXPORT PEL_NumberedDoubleVectors : public PEL_Object
{
   public: //------------------------------------------------------------------
   
   //-- Instance delivery and initialization

      // Create, initialize and return an instance devoted to store 
      // doubleVector instances of size `a_size_of_items'.
      static PEL_NumberedDoubleVectors* create(
                            PEL_Object* a_owner,
                            PEL_DoubleComparator const* dbl_comp,
                            size_t a_size_of_items ) ;

   //-- Status

      // number items
      size_t nb_items( void ) const ;

      // number of coordinates of all items
      size_t size_of_items( void ) const ;

   //-- Access with indices

      // Does `self' contains an item comparing equal to `item' ?
      bool has( doubleVector const& item ) const ;

      // index of the point matching `item'
      size_t index( doubleVector const& item ) const ;

   //-- Global access
      
      // array containing all items ordered independently from their 
      // introduction 
      doubleArray2D const& ordered_items( void ) const ;    

      // order of items independent from introduction order
      size_t_vector const& order( void ) const ;    
      
   //-- Element change

      // Ensure that self includes an item equal to `item'.
      void extend( doubleVector const& item ) ;
      
   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      PEL_NumberedDoubleVectors( void ) ;
     ~PEL_NumberedDoubleVectors( void ) ;
      PEL_NumberedDoubleVectors( PEL_NumberedDoubleVectors const& other ) ;
      PEL_NumberedDoubleVectors& operator=( 
                                 PEL_NumberedDoubleVectors const& other ) ;
      
      PEL_NumberedDoubleVectors( PEL_Object* a_owner, 
                                 PEL_DoubleComparator const* dbl_comp,
                                 size_t a_size_of_items ) ;

   //-- Attributes

      friend class PEL_NumberedDoubleVectorsItem ;
      
      PEL_DoubleComparator const* const DBL_COMP ;
      
      size_t const DIM ;
      size_t NB_PTS ;
      
      PEL_BalancedBinaryTree* PT_TREE ;
      mutable PEL_BalancedBinaryTreeIterator* PT_TREE_IT ;
      mutable doubleArray2D ALL_ITEMS ;
      mutable size_t_vector ORDER ;

      mutable PEL_NumberedDoubleVectorsItem* TMP ;
} ;

#endif
