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

#ifndef PEL_MAP_HH
#define PEL_MAP_HH

#include <PEL_Collection.hh>

#include <PEL_MapIterator.hh>

//----------------------------------------------------------------------------
// Data structures, used to store non NULL items identified by keys that may
// be hashed into an integer index
//----------------------------------------------------------------------------
// Implemented using a hash table of key/item pairs
//----------------------------------------------------------------------------
// To find an item associated to a given key, that key is first hashed into
// an integer index to identify the bucket in which the key/item is stored.
// Then, the key is searched within that bucket with a comparing criteria 
// based on the is_equal method.
//
// The number of buckets is set at creation. A number of (key/value) 
// pairs greatly exceeding the number of buckets will lead to a
// significant loss of efficiency when retrieving an item associated to a
// given key since the searching algorithm within a bucket is linear.
//----------------------------------------------------------------------------


class PEL_HashTableSet ;

class PEL_EXPORT PEL_Map : public PEL_Container
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return a new instance.
      static PEL_Map* create( PEL_Object* a_owner, size_t size = 20 ) ; 

      virtual PEL_Map* create_clone( PEL_Object* a_owner ) const ;

   //-- Measurement

      // number of buckets of the underlying hash table
      size_t nb_buckets( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      // Is there an item with a key comparing equal to `key' ?
      bool has_key( PEL_Object const* key ) const ;

      virtual PEL_Object* item( PEL_Object const* object ) const ;

      // item of key comparing equal to `key' if any, 0 otherwise 
      PEL_Object* item_at( PEL_Object const* key ) const ;

      virtual PEL_MapIterator* create_iterator( PEL_Object* a_owner ) const ;

   //-- Element change

      // Ensure that the item of key comparing equal to `key' is `a_item'.
      void set_item_at( PEL_Object* key, PEL_Object* a_item ) ;
   
   //-- Removal         

      // Remove the pair key/item of key comparing equal to `key'.
      void remove_at( PEL_Object const* key ) ;

      // Remove all pairs key/item.
      void clear( void ) ;

      
   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------

      friend class PEL_MapIterator ;

      PEL_Map( void ) ;
     ~PEL_Map( void ) ;
      PEL_Map( PEL_Map const& other ) ;
      PEL_Map const& operator=( PEL_Map const& other ) ;

      PEL_Map( PEL_Object* a_owner, size_t aSize ) ;
      PEL_Map( PEL_Object* a_owner, PEL_Map const* other ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_HashTableSet* hList ;
} ;


#endif
