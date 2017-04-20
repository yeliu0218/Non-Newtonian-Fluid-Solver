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

#ifndef PEL_VECTOR_ITERATOR_HH
#define PEL_VECTOR_ITERATOR_HH

#include <PEL_Iterator.hh>

class PEL_Vector ;

/*
Iterators on items of PEL_Vector collections.

The sequence of access is the order of increasing indices of items.
*/

class PEL_EXPORT PEL_VectorIterator : public PEL_Iterator
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an iterator over non NULL items of `vector'.
      static PEL_VectorIterator* create( PEL_Object* a_owner,
                                         PEL_Vector const* vector ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      size_t index_of_item( void ) const ;

      virtual PEL_Object* item( void ) const ;

   //-- Cursor movement

      void go_i_th( size_t i ) ;

      virtual void start( void ) ;

      virtual void go_next( void ) ;
             
   protected: //------------------------------------------------------------

      virtual ~PEL_VectorIterator( void ) ;

      PEL_VectorIterator( PEL_Object* a_owner,
                          PEL_Vector const* vector ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   private: //--------------------------------------------------------------

      PEL_VectorIterator( void ) ;
      PEL_VectorIterator( PEL_VectorIterator const& other ) ;
      PEL_VectorIterator const& operator=( PEL_VectorIterator const& other ) ;

   //-- Attributes

      PEL_Vector const* vec ;
      bool valid ;
      size_t current_i ;
      PEL_Object* current_item ;
} ;

#endif
