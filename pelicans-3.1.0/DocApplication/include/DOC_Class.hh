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

#ifndef DOC_Class_HH
#define DOC_Class_HH

#include <PEL_List.hh>
#include <string>

#include <DOC_Tools.hh>
#include <DOC_Writer.hh>
#include <DOC_Symbol.hh>

// DOC_Class representation.

class DOC_Attribute ;
class DOC_ClassItem ;
class DOC_Method ;
class DOC_Package ;
class DOC_Writer ;

class PEL_EXPORT DOC_Class : public DOC_Symbol
{
   public: //------------------------------------------------------

  //-- Assignment attempt(1)

      virtual DOC_Class* to_class( void ) ;
      
  //-- Creation(2)
      
      static DOC_Class* create( std::string const& a_name, 
                                DOC_Text* a_comment,
                                std::string const& a_file,
                                int a_line_nb,
                                DOC_Sequence* arg,
                                DOC_Class * inherit_from = 0 ) ;
      static DOC_Class* search( std::string const& a_name ) ;
      
  //-- General setting(3)
      
      // declare `class_name' as the name of a upper class in a
      // hierarchical sense.
      static void declare_upper( stringVector const& class_name ) ;
      
      void inherit( PEL_List* items  ) const ;

  //-- Status(4)
      
      std::string const& name( void ) const ;

      static void display_all( void ) ;

      DOC_ClassItem * find( std::string const& expr ) const ;
      
      DOC_Method * find_method( DOC_Method const* meth ) const ;
      
      static DOC_ClassItem const* search_and_find( std::string const& expr ) ;
      
      bool is_plug_point( void ) const ;
      
      bool is_external( void ) const ;
      
      bool is_upper_class( void ) const ;
      
      bool is_component( void ) const ;
      
      DOC_Class const* mother( void ) const ;
      
      std::string const& comment( void ) const ;
      
      PEL_List const* inherited_classes( void ) const ;
      
      PEL_List const* owned_components( void ) const ;
      
      bool is_virtual( void ) const ;
      
      DOC_Package const* package( void  ) const ;
      
      std::string short_description( void ) const ;

      // All classes list.
      static PEL_List* list( void ) ;
      
      // all items contained by `self'
      PEL_List const* all_elements( void ) const ;
      
  //-- Documentation processing(5)
      
      void process_documentation( DOC_Writer& sullizer ) ;
      
      bool has_to_document( DOC_ClassItem const* item ) const ;      
      
      bool publish_source( void ) const ;

  //-- Diagnostic(6)
      
      void verify( void  ) const ;

      void warning( std::string const& mess ) const ;    
      
  //-- Comparison(7)

      virtual int three_way_comparison( PEL_Object const* other ) const ;
     
   protected: //---------------------------------------------------

   private: //-----------------------------------------------------
      
      DOC_Class( void ) ;
      DOC_Class( PEL_Object* a_owner,
                 std::string const& a_name, 
                 DOC_Text* a_comment,
                 std::string const& a_file,
                 int a_line_nb,
                 DOC_Sequence* arg,
                 DOC_Class * inherit_from ) ;
      ~DOC_Class( void ) ;
     
      DOC_Class& operator=( DOC_Class const& c ) ;
      DOC_Class( DOC_Class const& c ) ;
      
      void add_descendant( DOC_Class const* element ) ;
      void add_element(  DOC_ClassItem* element ) ;
      
      void verify_default_implementation( void ) const ;
      void verify_protected_area( void ) const ;
      
      std::string myName ;
      DOC_Class * mother_class ;
      std::string strcomment ;
      bool virtuel ;
      bool abstract ;

      std::string file ;
      int line_number ;
      bool modele ;
      bool upper_DOC_Class ;
      
      static DOC_Class * courant ;
      static stringVector upper_DOC_Classs ;
      PEL_List* elements ;
      PEL_List* inherited_list ;
      bool pub_src ;
      mutable int COMPONENT ;
} ;



#endif


