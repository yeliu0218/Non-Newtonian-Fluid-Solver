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

#ifndef DOC_ClassItem_HH
#define DOC_ClassItem_HH

#include <iostream>
#include <string>

#include <DOC_Symbol.hh>

class DOC_Attribute ;
class DOC_Type ;
class DOC_Category ;
class DOC_Class ;
class DOC_Method ;
class DOC_Writer ;

// Class item.

class PEL_EXPORT DOC_ClassItem : public DOC_Symbol
{
   public: //-------------------------------------------------

   //-- Assignment attempt
      
      virtual DOC_ClassItem* to_class_item( void ) ;      
      void attach( DOC_Class* classe ) ;
      
      enum Protection { Public, Protected, Private }  ;

      virtual DOC_Method* method( void ) ;
      
      virtual DOC_Attribute* attribute( void ) ;

   //-- Status

      size_t index( void ) const ;
      
      DOC_Class const* def_class( void ) const ;
      
      bool visible( Protection level ) const ;
      
      Protection protection( void ) const ;
      
      virtual bool has_to_document( DOC_Class const* cl ) const ;

      DOC_Category const* category( void ) const ;

      std::string category_displaye( void ) const ;
      
      virtual std::string const& name( void ) const = 0 ;
      
      virtual bool is_inheritable( void ) const ;
      
      virtual std::string comment( void ) const ;

      virtual bool declare( std::string const& a__name ) const ;
      
   //-- Input - Output
      
      virtual std::string signature( void ) const = 0 ;
      
      virtual std::string prototype( DOC_Writer& sullizer ) const = 0 ;
      
      virtual void display( std::ostream& out ) const = 0 ;

      std::string const& bookmark( void ) const ;
      
      void set_bookmark( std::string const& a_name ) ;
      
    //-- Comparison

      virtual int three_way_comparison( PEL_Object const* other ) const ;
     
   protected: //-----------------------------------------------

      DOC_ClassItem( Protection maprotection,
                 DOC_Category const* macategory,
                 DOC_Text const* mycomment ) ;
      virtual ~DOC_ClassItem( void ) ;
      
      void inherit_category( DOC_ClassItem const* other) ;
      
   private: //-------------------------------------------------

      static size_t global_index ;
      DOC_ClassItem( void ) ;
      DOC_ClassItem( DOC_ClassItem const& other ) ;
      DOC_ClassItem& operator=( DOC_ClassItem const& other ) ;
      
      Protection niveauProtection ;
      DOC_Category const* maDOC_Category ;
      std::string myComment ;
      DOC_Class const* maDOC_Class ;
      size_t my_index ;
      std::string alias ;
} ;


#endif
