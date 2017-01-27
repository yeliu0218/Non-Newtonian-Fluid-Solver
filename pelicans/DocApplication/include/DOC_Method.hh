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

#ifndef DOC_Method_HH
#define DOC_Method_HH
#include <PEL_List.hh>

#include <DOC_Attribute.hh>

// DOC_Method representation.

class DOC_Category ;
class DOC_Function ;
class DOC_Type ;

class PEL_EXPORT DOC_Method : public DOC_Attribute
{
   public: //-------------------------------------------------

  //-- Assignment attempt

      virtual DOC_Method* method( void ) ;

   //-- Input - Output
      
      virtual std::string prototype( DOC_Writer& sullizer ) const ;
      virtual void display( std::ostream& out ) const ;
      virtual std::string signature( void ) const ;
      
  //-- Status
      virtual std::string comment( void ) const ;
      virtual bool is_inheritable( void ) const ;

  //-- Comparison
      virtual bool is_equal( PEL_Object const* other ) const ;
      
  //-- Creation
      
      // Depending on the fact that current parsed file has implemetation
      // kind or not, a new method is built or a method whose
      // prototype is already defined is completed.
      static DOC_Method* create(  std::string const& a_name,
                              DOC_Type * a_type,
                              Protection protection,
                              DOC_Category const* category,
                              DOC_Text * comment,
                              DOC_Sequence* arguments,
                              DOC_Sequence* modifiers,
                              int line = -1,
                              std::string const& file = "" ) ;
      
      void set_constant( void ) ;
      void set_virtual( void ) ;
      void overide( DOC_Method const* precedent ) ;
      void set_label( std::string const& a_label ) ;
      std::string const& label( void ) const ;
      enum Condition { pre, post, assert } ;
      void recover_conditions( DOC_Class * aclass,
                               Condition preCondition ) ;
      PEL_List const* conditions_list( Condition preCondition ) const ;
      void add_condition( DOC_Function * test, Condition cond ) ;
      
   //-- Status
 
      DOC_Method* has_method( DOC_Class const* a_classe,
         std::string const& method_name, std::string const& ext ) const ;
      
      PEL_List const* arguments( void ) const ;
      bool is_virtual( void ) const ;
      bool is_abstract( void ) const ;
      bool is_constructor( void ) const ;
      bool is_destructor( void ) const ;
      // name whithout suffix
      std::string short_name( void ) const ;
      std::string full_name( void ) const ;
      bool has_argument( std::string const& argstr ) const ;
      bool is_postcondition( void ) const ;
      bool is_precondition( void ) const ;
      bool is_invariant( void ) const ;
      DOC_Method const* overiden_method( void ) const ;
      PEL_List const* pre_conditions( void ) const ;
      PEL_List const* post_conditions( void ) const ;
      DOC_Method const* referring_method( void ) const ;
      
   //-- Diagnostic
      
      virtual void verify( void  ) const ;
      
      enum File { Implementation, Definition } ;
      
      void warning( File fich, std::string mess ) const ;

      std::string const& filename( File fich ) const ;

      // Format comment `com' using `red' to underline referenced
      //  item between `debut' and `fin' (class, class item)
      //  with `marqueur' level of output.
      static std::string format_comment( DOC_Writer const& red,
                                         size_t marqueur,
                                         std::string const& com,
                                         char debut,
                                         char fin,
                                         DOC_Class const* def ) ;
      
   protected: //-----------------------------------------------

   private: //-------------------------------------------------

      DOC_Method( void ) ;
      DOC_Method( DOC_Method const& other ) ;
      DOC_Method& operator=( DOC_Method const& other ) ;
      
      DOC_Method( std::string const& a_name,
              DOC_Type * a_type,
              Protection protection,
              DOC_Category const* category,
              DOC_Text * comment,
              DOC_Sequence* arguments,
              DOC_Sequence* modifiers ) ;
      
      ~DOC_Method( void ) ;


      void check_header( std::string const& header ) const ;

      void check_query_header( std::string const& header ) const ;

      void check_command_header( std::string const& header ) const ;

      bool has_compatible_arguments( DOC_Method const* called ) const ;
      
      DOC_Category const* laDOC_Category ;
      int index_cat ;
      PEL_List* args ;
      bool estVirtuel ;
      bool estVirtuel_pur ;
      bool is_constant ;
      DOC_Method const* surdefinition ;
      bool fait_PRE ;
      bool fait_POST ;
      bool constructeur ;
      bool destructeur ;
      DOC_Method const* implementation ;
      int def_line ;
      std::string def_fich ;
      int imp_line ;
      std::string imp_fich ;
      PEL_List* preConditions ;
      PEL_List* postConditions ;
      PEL_List* proceded_preConditions ;
      PEL_List* proceded_postConditions ;
      std::string my_label ;
      mutable DOC_ClassItem const* referring_item ; // For post and precondition only.
} ;


#endif
