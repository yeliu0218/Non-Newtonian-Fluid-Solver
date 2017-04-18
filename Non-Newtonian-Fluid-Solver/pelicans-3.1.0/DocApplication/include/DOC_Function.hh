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

#ifndef DOC_Function_HH
#define DOC_Function_HH

#include <DOC_Symbol.hh>
#include <string>

class DOC_Class ;
class DOC_Method ;
class PEL_List ;
class DOC_Writer ;

// DOC_Symbolic expressions.

class PEL_EXPORT DOC_Function : public DOC_Symbol
{
   public://----------------------------------------------------------
      
  //-- Assignment attempt
      
      virtual DOC_Function* to_function( void ) ;

   //-- Creation
      static DOC_Function* create( PEL_Object* a_owner,
                               std::string const& theDOC_Text ) ;
      static DOC_Function* create( PEL_Object* a_owner,
                               std::string const& theDOC_Text, DOC_Symbol * arg ) ;
      static DOC_Function* create( PEL_Object* a_owner,
                               std::string const& theDOC_Text, DOC_Symbol * arg1, DOC_Symbol * arg2, bool est_infixe = true ) ;
      static DOC_Function* create( PEL_Object* a_owner,
                               std::string const& theDOC_Text, DOC_Symbol * arg1, DOC_Symbol * arg2, DOC_Symbol * arg3 ) ;
      static DOC_Function* create( PEL_Object* a_owner,
                               DOC_Function * fonc,
                               PEL_List* lst ) ;

   //-- Status
      std::string name( void ) const ;
      // Is `self' representing some class item of `a_class' ?
      DOC_ClassItem const* is_class_item( DOC_Class const* a_class ) const ;
      
   //-- Input - Output
      
      virtual std::string text( void ) const ;
      
      virtual std::string referenced_string( DOC_Class const* classe,
                                             DOC_Method const* method,
                                             DOC_Writer const* red ) const ;

   protected://-------------------------------------------------------
   private://---------------------------------------------------------
      
      DOC_Symbol* arg( size_t n ) const ;
      DOC_Function( void ) ;
      DOC_Function( PEL_Object* a_owner,
                std::string const& theDOC_Text ) ;
      DOC_Function( PEL_Object* a_owner,
                std::string const& theDOC_Text, DOC_Symbol * arg ) ;
      DOC_Function( PEL_Object* a_owner,
                std::string const& theDOC_Text, DOC_Symbol * arg1, DOC_Symbol * arg2, bool est_infixe = true ) ;
      DOC_Function( PEL_Object* a_owner,
                std::string const& theDOC_Text, DOC_Symbol * arg1, DOC_Symbol * arg2, DOC_Symbol * arg3 ) ;
      DOC_Function( PEL_Object* a_owner,
                DOC_Function * fonc, PEL_List* lst ) ;

      ~DOC_Function( void ) ;
      DOC_Function& operator=( DOC_Function const& f ) ;
      DOC_Function( DOC_Function const& f ) ;
      std::string myDOC_Symbol ;
      PEL_List* mesDOC_Arguments ;
      bool infixe ;
      DOC_Function* symb ;
      bool unaire ;
      bool identificateur ;
      static int nb_blancs ;
      
} ;
      

#endif // DOC_Function_HH
