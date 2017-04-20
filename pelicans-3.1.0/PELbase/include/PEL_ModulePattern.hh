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

#ifndef PEL_MODULE_PATTERN_HH
#define PEL_MODULE_PATTERN_HH

#include <PEL_Object.hh>
#include <PEL_Data.hh>

#include <string>

class PEL_Module ;
class stringVector ;

/*
  Provides functionalities to analyse and/or verify data deck.
*/

class PEL_EXPORT PEL_ModulePattern : public PEL_Object
{
   public: //----------------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance attached to base pattern `a_name'.
      static PEL_ModulePattern* create( PEL_Object* a_owner,
                                        PEL_Module const* way,
                                        std::string const& a_name ) ;
      // Duplicate `self'.
      virtual PEL_ModulePattern* create_clone( PEL_Object* a_owner ) const ;

      // Search for description of class `class_name'.
      // Return NULL if no description is found.
      static PEL_ModulePattern* create_pattern_description(
                                        PEL_Object* a_owner,
                                        std::string const& class_name ) ;
      
  //-- Status report

      // Does current pattern is an indirection on sub-pattern ?
      static bool is_mutable( PEL_Module const* module ) ;
      
      
      // refering pattern description ( follow or the class name )
      static std::string const& indirection_to( PEL_Module const* module ) ;
      
      // does `module' refer to a module description ?
      static bool is_module_description( PEL_Module const* module ) ;
      
      // does `module' refer to an entry description ?
      static bool is_entry_description( PEL_Module const* module ) ;
      
      // does `module' refer to a variable description ?
      static bool is_variable_description( PEL_Module const* module ) ;
      
      // does `module' refer to a condition description ?
      static bool is_condition_description( PEL_Module const* module ) ;
      
      // valid choices for indirection
      static stringVector const& valid_indirections( PEL_Module const* module ) ;
      // needed module for current pattern
      stringVector const& mandatory_modules( PEL_Module const* way, bool first=true ) const ;
      
      // needed entries for current pattern
      stringVector const& mandatory_entries( PEL_Module const* way, bool first=true ) const ;
      
      // type of data
      std::string const& type_of_entry( std::string const& a_name ) const ;
      
      // name of `self'
      std::string const& name( void ) const ;
      
      // validity of `checked' with current pattern
      PEL_Module* validity( PEL_Object* a_owner,
                            PEL_Module const* checked,
                            bool recurse,
                            PEL_Module* result=0 ) const ;
      
      stringVector const&  provided_variables( void ) const ;
      
      static std::string const& variable_access( PEL_Module const* module,
                                                 PEL_Module const* way )  ;
      
   //-- Modifier
      
      enum Access { mandatory, optional, generic, unspecified } ;

      // Add module `path_and_name' to current pattern with accessibility
      //  `access' following `way'.
      void add_pattern( std::string const& path_and_name,
                        PEL_Module const* way,
                        Access access ) const ;
      
      // Add generic entry in current pattern following `way'.
      void add_generic_keyword( std::string const& a_keyword,
                                PEL_Module const* way,
                                PEL_Data::Type type ) ;
      
      // Add entry `a_keyword' in current pattern following `way'.
      void add_entry( std::string const& a_keyword,
                      PEL_Module const* way,
                      PEL_Data::Type type,
                      Access acc ) ;
      
      // Add provided entry `a_keyword' in current pattern following `way'.
      void add_provided_variable( std::string const& a_keyword,
                                  PEL_Module const* way,
                                  PEL_Data::Type type ) ;
      
      // Create subpattern of `self' following `way' at `path_and_name' position.
      PEL_ModulePattern* create_subpattern(
                                  PEL_Object* a_owner,
                                  PEL_Module const* way,
                                  std::string const& path_and_name ) const ;
      
      // pattern on generic module of current pattern following `way'
      PEL_ModulePattern* generic_pattern( PEL_Object* a_owner,
                                          PEL_Module const* way ) const ;
      
      static void expand_then_simplify( void ) ;
      
      void attach_help_data( std::string const& keyword,
                             std::string const& help ) ;
      
      void attach_default_data( std::string const& keyword,
                                std::string const& value ) ;
      
      void attach_verify_data( std::string const& keyword,
                               std::string const& expression ) ;
      
      void attach_list_of_valid_choices( std::string const& keyword,
                                         stringVector const& valid_choices ) ;
      
      void attach_list_of_dynamic_choices( std::string const& keyword,
                                           std::string const& regexp,
                                           std::string const& where ) ;
      
      void attach_file_extension( std::string const& filename,
                                  std::string const& mode ) ;
      
      // Is `a_name' an allowed module of current pattern ?
      PEL_Module* allowed_module( std::string const& a_name,
                                  PEL_Module const* way ) const ;
      
      PEL_Module* generic_module( PEL_Module const* way ) const ;
      
      PEL_Module* generic_entry( PEL_Module const* way ) const ;
      
      // Is `a_name' an allowed entry of current pattern ?
      PEL_Module const*  allowed_entry( std::string const& a_name,
                                        PEL_Module const* way  ) const ;
      
   //-- Module pattern base management

      // Attach module pattern base to file `filename' and start
      //  pattern recognition process.
      // If the file exists, it is read and inserted in base.
      static void build_pattern_base( std::string const& filename ) ;
      
      // Attach module pattern base to file `filename'.
      static void open_pattern_base( std::string const& filename ) ;
      
      // Close module pattern base.
      static void close_pattern_base( void ) ;
      
      // Save base to file.
      static void save_pattern( void ) ;
      
      // Has pattern to be built ?
      static bool build_pattern( void ) ;
      
      // Has pattern base ?
      static bool has_pattern( void ) ;
     
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------
      
      PEL_ModulePattern( void ) ;
     ~PEL_ModulePattern( void ) ;
      PEL_ModulePattern( PEL_ModulePattern const& other ) ;
      PEL_ModulePattern& operator=( PEL_ModulePattern const& other ) ;

      PEL_ModulePattern( PEL_Object* a_owner, PEL_Module* mm, PEL_Module const* a_father=0 ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Internal
      
      static void set_property( PEL_Module* mod,
                                std::string const& name,
                                std::string const& property ) ;
      
      static std::string const& property( PEL_Module const* mod,
                                          std::string const& name ) ;
      
      static bool is_module_matching( std::string const& module_name,
                                      std::string const& pattern_name ) ;
      
      void add_pattern_simple( std::string const& a_name,
                               PEL_Module const* way,
                               Access access ) const ;
      
      void add_entry_simple( std::string const& a_keyword,
                             PEL_Data::Type type,
                             Access access  ) ;
      
      static PEL_List* list_of_selected_conditional_modules(
                             PEL_Object * a_owner,
                             PEL_Module const* module,
                             PEL_Module const* way ) ;
      
      static PEL_List* list_of_conditional_modules(
                             PEL_Object * a_owner,
                             PEL_Module const* module ) ;
      
      static void check_name_validity( std::string const& a_name ) ;
      
      static PEL_Module*  module_of_pattern(
                             std::string const& class_name,
                             std::string const& instance_name,
                             bool build_if_not_exist ) ;
      
      static stringVector const& access_name( void ) ;
      
      PEL_Module* find_subpattern(
                             PEL_Module const* way,
                             std::string const& a_name,
                             PEL_Module *& a_father ) const ;
      
      PEL_Module* find_subpattern_simple(
                             PEL_Module const* way,
                             std::string const& a_name,
                             PEL_Module *& a_father ) const ;
      
      static void split( std::string const& a_name,
                         std::string& first_dir,
                         std::string& last ) ;
      
      static std::string class_part( std::string const& class_name ) ;
      
      void expand_then_simplify_one( void ) ;
      
      // Is value of `a_name' as an entry of `way' a valid
      //  with respect to `model' ?
      void check_allowed_value_for_entry(
                             PEL_Module* result,
                             PEL_Module const* model,
                             std::string const& a_name,
                             PEL_Module const* way ) const ;
      
      std::string inferred_generic_name( std::string const& a_name ) const ;
      
      static PEL_Module* child( PEL_Module const* module,
                                std::string const& name ) ;
      
      static void extend_entries( PEL_Module * module,
                                  std::string const entrie_name,
                                  std::string const& value ) ;
      
      static void complete_indirection_choices( PEL_Module * root ) ;

      static void check_evaluable( bool pattern_mod,
                                   PEL_Module const* mod,
                                   std::string const& keyword,
                                   PEL_Data::Type data_type ) ;
      
      static void check_in_spec( PEL_Module* result,
                                 PEL_Module const* model,
                                 std::string const& a_name,
                                 PEL_Module const* way ) ;
      
      static void check_vector_in_spec( PEL_Module* result,
                                        PEL_Module const* model,
                                        std::string const& a_name,
                                        PEL_Module const* way ) ;
      static void check_unique_spec( PEL_Module* result,
                                     PEL_Module const* model,
                                     std::string const& a_name,
                                     PEL_Module const* way ) ;
      static void merge( PEL_Module* current, PEL_Module const* to_add ) ;
      
   //-- Static attribute
      
      static PEL_Module* MP_File ;
      static PEL_Module* MP_Base ;
      static PEL_Module* MP_Pattern ;
      static std::string pattern_filename ;
      static bool HAS ;
      static bool BUILD ;

   //-- Attribute
      
      PEL_Module* const mod ;
      PEL_Module const* father ;
} ;

#endif
