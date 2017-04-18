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

#ifndef DOC_Symbol_HH
#define DOC_Symbol_HH

#include <PEL_Object.hh>
#include <string>

class DOC_Argument ;
class DOC_Attribute ;
class DOC_Text ;
class DOC_Class ;
class DOC_ClassItem ;
class DOC_Function ;
class DOC_Sequence ;
class DOC_Method ;
class DOC_Writer ;
class DOC_Type ;
class DOC_Typedef ;

// Lexical symbols produced throught lexical analysis. 

class PEL_EXPORT DOC_Symbol : public PEL_Object
{
   public: //---------------------------------------------------------

  //-- Assignment attempt
      
      // assign `self' as an instance of `DOC_Text::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Text* to_text( void ) ;
      
      // assign `self' as an instance of `DOC_Class::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Class* to_class( void ) ;
      
      // assign `self' as an instance of `DOC_Type::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Type* to_type( void ) ;
      
      // assign `self' as an instance of `DOC_Sequence::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Sequence* to_sequence( void ) ;
      
      // assign `self' as an instance of `DOC_Argument::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Argument* to_argument( void ) ;
      
      // assign `self' as an instance of `DOC_Function::'
      // IMPLEMENTATION : return NULL
      virtual DOC_Function* to_function( void ) ;
      
      // assign `self' as an instance of `DOC_ClassItem::'
      // IMPLEMENTATION : return NULL
      virtual DOC_ClassItem* to_class_item( void ) ;

   //-- Creation
   //-- Input - Output
      
      virtual std::string text( void ) const ;
      
      // display `self' using `red' in volatile context defined by
      // `method' and `red' (they can be null) to express dynamic reference
      virtual std::string referenced_string( DOC_Class const* classe,
                                             DOC_Method const* method,
                                             DOC_Writer const* red ) const ;
      
   protected: //------------------------------------------------------
      
      DOC_Symbol( PEL_Object* a_owner ) ;
      virtual ~DOC_Symbol( void ) ;
     
   private: //--------------------------------------------------------
      
      DOC_Symbol& operator=( DOC_Symbol const& s) ;
      DOC_Symbol( DOC_Symbol const& s ) ;
      DOC_Symbol( void ) ;
      void error( void ) ;
      
} ;

      
#endif // DOC_Symbol_HH
