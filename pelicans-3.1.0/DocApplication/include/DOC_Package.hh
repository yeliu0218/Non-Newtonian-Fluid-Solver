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

#ifndef DOC_PACKAGE_HH
#define DOC_PACKAGE_HH

#include <PEL_List.hh>
#include <string>
#include <stringVector.hh>

class DOC_Class ;
class DOC_Writer ;

// DOC_Class packaging.
// Only `Master' package doesn't own father.

class PEL_EXPORT DOC_Package : public PEL_Object
{
   public: //----------------------------------------------------------------
   //-- Input - Output
      // Parse configuration file.
      static void read( std::string const& file ) ;

   //-- Status
      std::string const& name( void ) const ;
      DOC_Package* father( void ) const ;
      DOC_Package* find_sub_package( std::string const& a_name ) const ;      
      std::string const& comment( void ) const ;
      
      // package associated to `classe'.
      static DOC_Package const* attached_to(  DOC_Class const* classe,
                                              std::string& comment_court ) ;

      // Search for a package given its name `a_name'.
      static DOC_Package* find( std::string const& a_name ) ;
      
      // Higher package.
      static DOC_Package* master( void )  ;
      
     PEL_List const* sub_packages( void ) const ;

      // Is `cl' a class contained in `self' or a sub-package of `self' ?
      bool own_class( DOC_Class const* cl ) const ;

      bool own_non_external_classes( void ) const ;

      // top package under master one
      DOC_Package const* top( void ) const ;
      
      
       static std::string const& application( void ) ;
   //-- Input - Output
      void print( std::ostream& os, size_t ident_witdth ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      DOC_Package( void ) ;
     ~DOC_Package( void ) ;
      DOC_Package( DOC_Package const& other ) ;
      DOC_Package& operator=( DOC_Package const& other ) ;

      DOC_Package( PEL_Object* a_owner,
                   std::string const& a_name,
                   DOC_Package* a_father,
                   std::string const& com ) ;

      DOC_Package( PEL_Object* a_owner ) ;

      // list of all packages
      static PEL_List* packages( void ) ;

      //-- Creation
      static DOC_Package* create( std::string const& a_name,
                                  DOC_Package* a_father,
                                  std::string const& comment  ) ;
       // Hidden package.
      static DOC_Package* miscellanous( void )  ;
      
   //-- Class attributes

      static DOC_Package* Master ;
      static std::string application_name ;
      
   //-- Attributes

      std::string my_name ;
      std::string my_comment ;
      DOC_Package* my_father ;
      PEL_List* subpackages ;
      stringVector elements ;
      stringVector element_comments ;
} ;

#endif
