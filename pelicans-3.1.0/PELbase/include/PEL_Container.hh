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

#ifndef PEL_CONTAINER_HH
#define PEL_CONTAINER_HH

#include <PEL_Object.hh>

class PEL_Iterator ;

/*
Data structures
  - used to hold zero or more items, where item denotes a possibly NULL
    object ;
  - with a finite item count ;
  - for which there exist a traversal policy that will visit every item
    exactly once.
The items are referred to, accessed and manipulated exclusively through
pointers. The NULL object is thus represented by the 0 pointer.
*/

class PEL_EXPORT PEL_Container : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Container* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Status

      // Do `object1' and `object2' compare equal for search operations ?
      // IMPLEMENTATION : use of is_equal()
      virtual bool matching_items( PEL_Object const* object1,
                                   PEL_Object const* object2 ) const ;

      // identification of the object state (two different states have
      // different identifiers)
      size_t state_id( void ) const ;

   //-- Measurement

      // number of non NULL items
      virtual size_t count( void ) const = 0 ;

   //-- Access

      // Does self include an item matching `object' ?
      bool has( PEL_Object const* object ) const ;

      // first item matching object if any ; O otherwise
      virtual PEL_Object* item( PEL_Object const* object ) const = 0 ;

      // Create and return an iterator on non NULL items.
      virtual PEL_Iterator* create_iterator( PEL_Object* a_owner ) const = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;      
      
   protected: //--------------------------------------------------------

      virtual ~PEL_Container( void ) ;

      PEL_Container( PEL_Object* a_owner ) ;
      
      bool create_clone_POST( PEL_Container* result,
                              PEL_Object* a_owner ) const ;

      virtual bool count_POST( size_t result ) const ;

      virtual bool matching_items_PRE( PEL_Object const* object1,
                                       PEL_Object const* object2 ) const ;

      virtual bool item_PRE( PEL_Object const* object ) const ;
      virtual bool item_POST( PEL_Object const* result, 
                              PEL_Object const* object ) const ;

      virtual bool create_iterator_POST( PEL_Iterator* result,
                                         PEL_Object* a_owner ) const ;

      void report_state_change( void ) ;
      
   private: //----------------------------------------------------------

      PEL_Container( void ) ;
      PEL_Container( PEL_Container const& other ) ;
      PEL_Container const& operator=( PEL_Container const& other ) ;

   //-- Attributes
      
      size_t MY_STATE ;
} ;

#endif
