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

#ifndef PEL_HASH_TABLE_SET_HH
#define PEL_HASH_TABLE_SET_HH

#include <PEL_Set.hh>

#include <PEL_HashTableSetIterator.hh>
#include <PEL_List.hh>

//---------------------------------------------------------------------------
//  Sets implemented as hash tables
//---------------------------------------------------------------------------

class PEL_Vector ;

class PEL_EXPORT PEL_HashTableSet : public PEL_Set
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return a new instance.
      static PEL_HashTableSet* create( PEL_Object* a_owner,
                                       size_t size = 20 ) ; 

      virtual PEL_HashTableSet* create_clone( PEL_Object* a_owner ) const ;

   //-- Measurement   

      // number of buckets
      size_t nb_buckets( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      virtual PEL_Object* item( PEL_Object const* object ) const ;

      virtual PEL_HashTableSetIterator* create_iterator( 
                                             PEL_Object* a_owner ) const ;

   //-- Element change

      virtual void extend( PEL_Object* object ) ;

   //-- Removal
      
      virtual void remove( PEL_Object const* object ) ;

      virtual void destroy_item_and_remove( PEL_Object const* object ) ;

      virtual void clear( void ) ;

      virtual void destroy_items_and_clear( void ) ;

   protected: //------------------------------------------------------------

      virtual ~PEL_HashTableSet( void ) ;

      PEL_HashTableSet( PEL_Object* a_owner, size_t size ) ;

   //-- Preconditions, Postconditions, Invariant     
 
      virtual bool invariant( void ) const ;

   private: //--------------------------------------------------------------

      friend class PEL_HashTableSetIterator ;

      PEL_HashTableSet( void ) ;
      PEL_HashTableSet( PEL_HashTableSet const& other ) ;
      PEL_HashTableSet const& operator=( PEL_HashTableSet const& other ) ;
      
      PEL_List * getList( size_t entry ) ;
      PEL_List const* getList( size_t entry ) const ;
      
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      PEL_Vector* hTable ;
      size_t the_size ;
      
} ;


#endif
