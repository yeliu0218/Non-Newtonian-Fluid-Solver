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

#ifndef PEL_ITERATOR_HH
#define PEL_ITERATOR_HH

#include <PEL_Object.hh>

class PEL_Container ;

/*
Iterators to traverse data structures and visit every non NULL items 
exactly once

Items are accessed sequentially, one-way.
The sequence of access depends on the underlying data structure.
*/

class PEL_EXPORT PEL_Iterator : public PEL_Object
{
   public: //---------------------------------------------------------------

   //-- Status report

      // Has the traversed container been modified since last call 
      // to `::start' ?
      virtual bool container_has_been_modified( void ) const ;

      // Is current position valid ?
      virtual bool is_valid( void ) const = 0 ;

   //-- Access

      // item at the current position
      virtual PEL_Object* item( void ) const = 0 ;

   //-- Cursor movement

      // Move to first position.
      virtual void start( void ) = 0 ;

      // Move to next position.
      virtual void go_next( void ) = 0 ;
      
   protected: //------------------------------------------------------------
      
      virtual ~PEL_Iterator( void ) ;

      PEL_Iterator( PEL_Object* a_owner, PEL_Container const* list ) ;
      
      PEL_Iterator( PEL_Object* a_owner ) ;
      
      void base_initialize( PEL_Container const* list ) ;

      void record_container_state_id( void ) ;
      
   //-- Preconditions, Postconditions, Invariant   
   
      virtual bool invariant( void ) const ;

      virtual bool start_POST( void ) const ;
      
      virtual bool item_PRE( void ) const ;
      virtual bool is_valid_PRE( void ) const ;
      virtual bool go_next_PRE( void ) const ;
      virtual bool item_POST( PEL_Object const* result ) const ;
      
   private: //--------------------------------------------------------------

      PEL_Iterator( void ) ;
      PEL_Iterator( PEL_Iterator const& other ) ;
      PEL_Iterator const& operator=( PEL_Iterator const& other ) ;

   //-- Attributes

      PEL_Container const* my_list ;
      size_t list_state ;

} ;

#endif
