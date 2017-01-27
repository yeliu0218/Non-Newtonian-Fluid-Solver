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

#ifndef PEL_VECTOR_HH
#define PEL_VECTOR_HH

#include <PEL_Sequence.hh>

#include <PEL_VectorIterator.hh>

/*
Sequences whose items access from their index if efficient
*/

class PEL_EXPORT PEL_Vector : public PEL_Sequence
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create an return an instance.
      static PEL_Vector* create( PEL_Object* a_owner, size_t size ) ;

      // Reinitialize the internal state, as if the `::create' method
      // was just completed.
      void re_initialize( size_t size ) ;

      virtual PEL_Vector* create_clone( PEL_Object* a_owner ) const ;

      // Reinitialize by copying all the items of `other'.
      void copy( PEL_Vector const* other ) ;
      
   //-- Resizing

      // Change the exclusive upper limit for indices, without losing
      // previously entered items whose index was within the new bound,
      // nor changing their coordinate.
      void resize( size_t size ) ;

  //-- Measurement

      virtual size_t index_limit( void ) const ;

      virtual size_t count( void ) const ;
      
   //-- Access
      
      virtual PEL_Object* item( PEL_Object const* object ) const ;

      virtual PEL_Object* at( size_t i ) const ;

      virtual size_t index_of( PEL_Object const* object ) const ;

      virtual PEL_VectorIterator* create_iterator( PEL_Object* a_owner ) const ;

   //-- Element change

      // Make all items be the 0 pointer.
      void nullify( void ) ;
      
      virtual void append( PEL_Object* object ) ;

      virtual void prepend( PEL_Object* object ) ;

      virtual void set_at( size_t i, PEL_Object* object ) ;

      virtual void insert_at( size_t i, PEL_Object* object ) ; 

   //-- Removal
      
      virtual void remove_at( size_t i ) ;
      
      virtual void remove_section( size_t iFirst, size_t length ) ;

      virtual void destroy_items_and_remove_section( size_t iFirst, 
                                                 size_t length ) ;
      virtual void clear( void ) ;

         
   protected: //--------------------------------------------------------

      virtual ~PEL_Vector( void );

      PEL_Vector( PEL_Object* a_owner, size_t size ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id  ) const ;

  private: //----------------------------------------------------------

      PEL_Vector( void ) ;
      PEL_Vector( PEL_Vector const& other ) ;
      PEL_Vector& operator=( PEL_Vector const& other ) ;
   
   //-- Attributes

      PEL_Object** VECTOR ;
      size_t LENGTH ;
      size_t NB_ENTRIES ;
      size_t CAPACITY ;
} ;

#endif
