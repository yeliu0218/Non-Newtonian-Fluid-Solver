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

#ifndef ARGUMENT_HH
#define ARGUMENT_HH

#include <iostream>
#include <string>
#include <DOC_Symbol.hh>

class DOC_Type ;

// Class argument.

class PEL_EXPORT DOC_Argument : public DOC_Symbol
{
   public: //--------------------------------------------------

  //-- Assignment attempt
      
      virtual DOC_Argument* to_argument( void ) ;
      
      static DOC_Argument* create( PEL_Object* a_owner,
                                   DOC_Type * typ,
                                   std::string const& name ) ;
      
      // Set initial value to `init'.
      void initialize( std::string const& init ) ;

   //-- Input - Output
      
      DOC_Type const* type( void ) const ;
      std::string const& silent_var( void ) const ;
      virtual std::string text( void ) const ;

      //-- Comparison
      bool is_equal( PEL_Object const* other ) const ;
      
   protected: //-----------------------------------------------

      DOC_Argument( PEL_Object* a_owner, DOC_Type * typ, std::string const& name ) ;
      virtual ~DOC_Argument( void ) ;
      
   private: //-------------------------------------------------

      DOC_Argument& operator=( DOC_Argument const& a ) ;
      DOC_Argument( DOC_Argument const& a ) ;
      DOC_Argument( void ) ;
      
      std::string nameMuet ;
      std::string initialisation ;
      DOC_Type * myDOC_Type ;
} ;

#endif

      
