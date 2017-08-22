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

#ifndef PEL_SET_HH
#define PEL_SET_HH

#include <PEL_Collection.hh>

//---------------------------------------------------------------------------
//   Collections with no NULL items and no matching different items.
//---------------------------------------------------------------------------

class PEL_EXPORT PEL_Set : public PEL_Collection
{
   public : //----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Set* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Element change 

      virtual void extend( PEL_Object* object ) = 0 ;
      
   //-- Removal

      virtual void remove( PEL_Object const* object ) = 0 ;

   protected : //-------------------------------------------------------

      virtual ~PEL_Set( void ) ;

      PEL_Set( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool extend_POST( bool old_has, 
                                size_t old_count, 
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool remove_POST( size_t old_count,
                                PEL_Object const* object,
                                size_t old_state_id ) const ;

   private: //----------------------------------------------------------

      PEL_Set( void ) ;
      PEL_Set( PEL_Set const& other ) ;
      PEL_Set const& operator=( PEL_Set const& other ) ;
} ;

#endif
