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

#ifndef PEL_KEY_ITEM_PAIR_HH
#define PEL_KEY_ITEM_PAIR_HH

#include <PEL_Object.hh>

class PEL_EXPORT PEL_KeyItemPair : public PEL_Object
{
      
   public: //---------------------------------------------------------------


   //-- Initialization

      static PEL_KeyItemPair* create( PEL_Object* a_owner,
                                      PEL_Object* a_key,
                                      PEL_Object* a_item ) ; 

   //-- Access

      virtual size_t hash_code( void ) const ;

      PEL_Object* key( void ) const ;

      PEL_Object* item( void ) const ;

   //-- Element change

      void set_item( PEL_Object* object ) ;

   //-- Comparison

      virtual bool comparable( PEL_Object const* other ) const ;

      virtual bool is_equal( PEL_Object const* other ) const ;

      virtual int three_way_comparison( PEL_Object const* other ) const ;

      
   protected: //----------------------------------------------------

      virtual bool invariant() const ;

   private: //------------------------------------------------------

      PEL_KeyItemPair( void ) ;
     ~PEL_KeyItemPair( void ) ;
      PEL_KeyItemPair( PEL_KeyItemPair const& other ) ;
      PEL_KeyItemPair const& operator=( PEL_KeyItemPair const& other ) ;
                      
      PEL_KeyItemPair( PEL_Object* a_owner,
                       PEL_Object* a_key,
                       PEL_Object* a_item ) ;
      
      //------------------------------------------------------------
      //   ATTRIBUTES
      //------------------------------------------------------------
      PEL_Object* the_key ;
      PEL_Object* the_item ;

} ;


#endif
