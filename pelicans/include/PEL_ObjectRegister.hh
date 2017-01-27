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

#ifndef PEL_OBJECT_REGISTER_HH
#define PEL_OBJECT_REGISTER_HH

#include <PEL_Object.hh>

#include <string>

/*
Registers of `PEL_Object::' instances, eg usable in pluggable factories.

PUBLISHED   
*/

class PEL_Iterator ;
class PEL_Map ;

class PEL_EXPORT PEL_ObjectRegister : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PEL_ObjectRegister* create( PEL_Object* a_owner,
                                         std::string const& a_register_name ) ;
      
   //-- Characteristics

      // name of `self'
      std::string const& register_name( void ) const ;
      
   //-- Access

      // Is there an object registered under the name `a_name' ?
      bool has( std::string const& a_name ) const ;

      // object registered under the name `a_name'
      // (fatal error raised if none)
      PEL_Object* item( std::string const& a_name ) const ;

      // Create and return an iterator on registered objects.
      PEL_Iterator* create_iterator( PEL_Object* a_owner ) const ;
      
   //-- Registration
      
      // Register `an_item' under the name `a_name'.
      // (fatal error raised if such a registration name already exists).
      void register_item( std::string const& a_name, PEL_Object* an_item ) ;

      // Suppress the object that was registered under the name `a_name'
      // (fatal error raised if none).
      void unregister_item( std::string const& a_name ) ;

      // Suppress the registered object `an_item'
      // (fatal error raised if none).
      void unregister_item( PEL_Object* an_item ) ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PEL_ObjectRegister( void ) ;
     ~PEL_ObjectRegister( void ) ;
      PEL_ObjectRegister( PEL_ObjectRegister const& other ) ;
      PEL_ObjectRegister& operator=( PEL_ObjectRegister const& other ) ;

      PEL_ObjectRegister( PEL_Object* a_owner,
                          std::string const& a_register_name ) ;
      
   //-- Attributes
      
      PEL_Map* const REGISTER ;
      std::string const NAME ;
} ;


#endif

