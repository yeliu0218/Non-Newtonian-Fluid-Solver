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

#ifndef PEL_LIST_HH
#define PEL_LIST_HH

#include <PEL_Sequence.hh>

#include <PEL_ListIterator.hh>

class PEL_ListItem ;

/*
Sequential, one-way linked lists
*/

class PEL_EXPORT PEL_List : public PEL_Sequence
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static PEL_List* create( PEL_Object* a_owner ) ;

      virtual PEL_List* create_clone( PEL_Object* a_owner ) const ;

      // Reinitialize by copying all the items of `other'.
      void copy( PEL_List const* other ) ;

  //-- Measurement

      virtual size_t index_limit( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      virtual PEL_Object* item( PEL_Object const* object ) const ;

      virtual PEL_Object* at( size_t i ) const ;
         
      virtual size_t index_of( PEL_Object const* object ) const ;

      virtual PEL_ListIterator* create_iterator( PEL_Object* a_owner ) const ;

   //-- Element change

      virtual void append( PEL_Object* object ) ;
      
      virtual void prepend( PEL_Object* object ) ;

      virtual void set_at( size_t i, PEL_Object* object ) ;

      virtual void insert_at( size_t i, PEL_Object* object ) ;

   //-- Removal

      virtual void remove( PEL_Object const* object ) ;
         
      virtual void remove_at( size_t i ) ;
      
      virtual void remove_section( size_t iFirst, size_t length ) ;
      
      virtual void destroy_items_and_clear( void ) ;
      
      virtual void destroy_items_and_remove_section( size_t iFirst, 
                                                 size_t length ) ;

      virtual void clear( void ) ;

   protected: //------------------------------------------------------------

      virtual ~PEL_List( void ) ;

      PEL_List( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      virtual bool set_at_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t i,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;
      virtual bool at_POST( PEL_Object const* result, size_t i ) const ;

      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id ) const ;
      
   private: //--------------------------------------------------------------

      friend class PEL_ListIterator ;

      PEL_List( void ) ;
      PEL_List( PEL_List const& other ) ;
      PEL_List const& operator=( PEL_List const& other ) ;

      PEL_ListItem* theItem( size_t n ) const ;
      
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      PEL_ListItem* theList ;
      PEL_ListItem* theLast ;
      size_t nbElem ;
      
} ;


#endif
