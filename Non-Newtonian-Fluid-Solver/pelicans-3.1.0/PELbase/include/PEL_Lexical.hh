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

#ifndef PEL_LEXICAL_H
#define PEL_LEXICAL_H

#include <PEL_Object.hh>

#include <string>

/* Provides a common interface to Pelicans data structures for reading. */

class PEL_KeywordDataPair ;
class PEL_Module ;
class PEL_List ;
class PEL_Data ;
class PEL_String ;

class PEL_EXPORT PEL_Lexical : public PEL_Object
{
   public: //-------------------------------------------------------
      
      // Initialization
      static PEL_Lexical* create( PEL_Module* aModule ) ;
      static PEL_Lexical* create( PEL_Data * aSimpleValue ) ;
      static PEL_Lexical* create( PEL_List * aList ) ;
      
   //-- Access

      // Is self a module ?
      bool is_module( void ) const;
      
      // Is self a simple value ?
      bool is_data( void ) const ;
      
      // Is self a list ?
      bool is_list( void ) const ;
      
   //-- Retrieving
      
      // Retrieves self as a module.
      PEL_Module * to_module( void ) ;
      
       // Retrieves self as a simple value.
      PEL_Data* to_data( void ) ;
      
       // Retrieves self as a list.
      PEL_List * to_list( void ) ;
      
  //-- Removing   
      static void remove_all_lexical( void ) ;
      
   protected: //-------------------------------------------------------
      
      // Destructor
      virtual ~PEL_Lexical( void ) ;
      
   private: //-------------------------------------------------------
      
      PEL_Lexical( PEL_Lexical const& other  ) ;
      PEL_Lexical const& operator=( PEL_Lexical const& other  ) ;     

      PEL_Lexical( void ) ;
      PEL_Lexical( PEL_Object* a_owner, PEL_Module* aModule ) ;
      PEL_Lexical( PEL_Object* a_owner, PEL_KeywordDataPair * anAssignment ) ;
      PEL_Lexical( PEL_Object* a_owner, PEL_Data * aSimpleValue ) ;
      PEL_Lexical( PEL_Object* a_owner, PEL_List * aSimpleValue ) ;
      
      bool invariant( void ) const ;
      static PEL_Lexical* common_owner( void ) ;
      static PEL_Lexical* the_common_owner ;
      
      //-------------------------------------------------------------
      //   ATTRIBUTES
      //-------------------------------------------------------------
      PEL_Module* myModule ;
      PEL_Data* myData ;
      PEL_List* myList ;
};

typedef PEL_Lexical* PEL_LexicalPtr ;

#endif 
