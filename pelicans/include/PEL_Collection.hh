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

#ifndef PEL_COLLECTION_HH
#define PEL_COLLECTION_HH

#include <PEL_Container.hh>

/*
Containers for which basic insertion and removal of items is  meaningfull
*/

class PEL_EXPORT PEL_Collection : public PEL_Container
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Collection* create_clone( PEL_Object* a_owner ) const = 0 ;
           
   //-- Element change 

      // Ensure that self includes an item matching `object'.
      virtual void extend( PEL_Object* object ) = 0 ;

   //-- Removal

      // Remove the first item matching `object'.
      virtual void remove( PEL_Object const* object ) = 0 ;

      // Remove and terminate the first item matching `object'.  */
      virtual void destroy_item_and_remove( PEL_Object const* object ) = 0 ;
      
      // Remove all items.
      virtual void clear( void ) = 0 ;

      // Remove and terminate all items.
      virtual void destroy_items_and_clear( void ) = 0 ;

   protected: //--------------------------------------------------------

      virtual ~PEL_Collection( void ) ;

      PEL_Collection( PEL_Object* a_owner ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool extend_PRE( PEL_Object const* object ) const ;
      virtual bool extend_POST( bool old_has,
                                size_t old_count,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool remove_PRE( PEL_Object const* object ) const ;
      virtual bool remove_POST( size_t old_count,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool clear_POST( size_t old_state_id ) const ;
      
   private: //----------------------------------------------------------

      PEL_Collection( void ) ;
      PEL_Collection( PEL_Collection const& other ) ;
      PEL_Collection const& operator=( PEL_Collection const& other ) ;
      
} ;

#endif
