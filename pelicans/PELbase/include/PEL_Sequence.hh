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

#ifndef PEL_SEQUENCE_HH
#define PEL_SEQUENCE_HH

#include <PEL_Collection.hh>

/*
Collections whose items are in one-to-one correspondance with
integers of a contiguous interval with zero as lower bound.
items can be the NULL object.
*/

class PEL_EXPORT PEL_Sequence : public PEL_Collection
{

   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Sequence* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Comparison

      virtual bool comparable( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

      // IMPLEMENTATION : true if three_way_comparison() returns 0
      virtual bool is_equal( PEL_Object const* other ) const ;

      // IMPLEMENTATION : lexicographic order
      virtual int three_way_comparison( PEL_Object const* other ) const ;

   //-- Measurement

      // exclusive upper limit for indices
      virtual size_t index_limit( void ) const = 0 ;

   //-- Access
      
      /** Class constant returned by index when no element of self match
          is_equal relation. */
      static size_t const badIndex ;

      // item at index `i'
      virtual PEL_Object* at( size_t i ) const = 0 ;

      // index of the first item matching `object' if any ; a number out
      // of the valid index interval if none
      virtual size_t index_of( PEL_Object const* object ) const = 0 ;

   //-- Sorting

      virtual void sort( void ) ;

   //-- Element change

      virtual void extend( PEL_Object* object ) ;

      // Add `object' to end.
      virtual void append( PEL_Object* object ) = 0 ;

      // Add `object' to beginning.
      virtual void prepend( PEL_Object* object ) = 0  ;

      // Assign `object' to the `i'-th entry.
      virtual void set_at( size_t i, PEL_Object* object ) = 0  ;

      // Increase by one the index of all items starting from the `i'-th 
      // position and add `object' to the `i'-th position.
      virtual void insert_at( size_t i, PEL_Object* object ) = 0  ;

   //-- Removal

      virtual void remove( PEL_Object const* object ) ;

      virtual void destroy_item_and_remove( PEL_Object const* object ) ;
      
      // Remove item at index `i'.
      virtual void remove_at( size_t i ) = 0 ;

      // Terminate and remove item at index `i'.
      void destroy_item_and_remove_at( size_t i ) ;

      // Remove `length' contiguous items, starting from the one of index
      // `iFirst'.
      virtual void remove_section( size_t iFirst, size_t length ) = 0 ;

      // Terminate and remove `length' contiguous items, starting from the 
      // one of index `iFirst'.
      virtual void destroy_items_and_remove_section( size_t iFirst,
                                                     size_t length ) = 0 ;

      virtual void destroy_items_and_clear( void ) ;
      
      virtual void clear( void ) = 0 ;
      
   //-- Statics

   protected: //------------------------------------------------------------

      virtual ~PEL_Sequence( void ) ;

      PEL_Sequence( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      bool create_clone_POST( PEL_Sequence* result,
                              PEL_Object* a_owner ) const ;

      virtual bool append_PRE( PEL_Object const* object ) const ;
      virtual bool append_POST( size_t old_index_limit,
                                size_t old_count,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool prepend_PRE( PEL_Object const* object ) const ;
      virtual bool prepend_POST( size_t old_index_limit,
                                 size_t old_count,
                                 PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool set_at_PRE( size_t i, PEL_Object const* object ) const ;
      virtual bool set_at_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t i,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool insert_at_PRE( size_t i, PEL_Object const* object ) const ;
      virtual bool insert_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t i,
                                   PEL_Object const* object,
                                   size_t old_state_id ) const ;

      virtual bool remove_at_PRE( size_t i ) const ;
      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_PRE( size_t iFirst, size_t length ) const ;
      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id ) const ;

      virtual bool at_PRE( size_t i ) const ;
      virtual bool at_POST( PEL_Object const* resu, size_t i ) const ;

      virtual bool index_of_PRE( PEL_Object const* object ) const ;
      virtual bool index_of_POST( size_t result, PEL_Object const* object ) const ;
      
      virtual bool destroy_items_and_remove_section_PRE( size_t iFirst,
                                                         size_t length ) const ;

      virtual bool clear_POST( size_t old_state_id ) const ;

   private: //--------------------------------------------------------------

      PEL_Sequence( void ) ;
      PEL_Sequence( PEL_Sequence const& other ) ;
      PEL_Sequence const& operator=( PEL_Sequence const& other ) ;

} ;

#endif
