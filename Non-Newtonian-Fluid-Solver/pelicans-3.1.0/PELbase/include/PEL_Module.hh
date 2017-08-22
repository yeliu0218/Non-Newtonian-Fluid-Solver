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

#ifndef PEL_MODULE_HH
#define PEL_MODULE_HH

#include <PEL_Object.hh>

#include <string>
#include <iosfwd>

class PEL_KeywordDataPair ;
class PEL_KeywordDataIterator ;
class PEL_List ;
class PEL_ModuleIterator ;
class PEL_ModuleExplorer ;
class PEL_Context ;
class PEL_ContextPair ;
class PEL_ContextSimple ;
class PEL_Data ;
class PEL_DataWithContext ;
class PEL_String ;
class PEL_Variable ;

/* 
Nodes of the PELICANS Hierarchical Data System.

The PELICANS Hierarchical Data System is described in `PEL_ModuleExplorer::'.
Note on context : all childs module share context from their father with
 their own context.
*/

class PEL_EXPORT PEL_Module : public PEL_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Creates and return an instance containing neither modules nor entries.
      static PEL_Module* create( PEL_Object* a_owner,
                                 std::string const& a_name ) ;
      
      // Create and return an instance being the root a the module hierarchy
      // stored in `file_name'.
      static PEL_Module* create( PEL_Object* a_owner,
                                 std::string const& a_name,
                                 std::string const& file_name,
                                 PEL_Context const* ct = 0 ) ;
      
      static PEL_Module* create( PEL_Object* a_owner,
                                 std::string const& a_name,
                                 std::istream& input_stream,
                                 PEL_Context const* ct = 0 ) ;

      // Create and return a clone of `self': modules and entries are clones,
      // all variables of `::father'() context needed is retrieved.
      virtual PEL_Module* create_clone( PEL_Object* a_owner ) const ;

      // Create and return a copy of `self': modules and entries are clones,
      // `::father'() context is lost.
      PEL_Module* create_copy( PEL_Object* a_owner ) const ;
      
      static PEL_Module* create_as_difference( PEL_Object* a_owner,
                                               std::string const& a_name,
                                               PEL_Module const* m1,
                                               PEL_Module const* m2,
                                               PEL_ModuleExplorer const* exp ) ;

   //-- Modifiers

      void modify_module_name( std::string const& a_name ) ;
      
      // Make `self' contain the keyword-data pair defined by `keyword'
      // and `data' as an entry. If `self' already contains an entry
      // with the same keyword, a fatal error is raised.
      void add_entry( std::string const& keyword, PEL_Data const* data ) ;

      // Make `self' contain `a_module'. 
      // If `self' already contains a module with the same name as `a_module', 
      // the two modules are merged (if it happens that two leaves have the 
      // same keyword, a fatal error is raised).
      void add_module( PEL_Module* a_module ) ;
      
      // Add all the modules and entries of `a_module' to `self'.
      // If it happens that two leaves have the same keyword, the value of
      // `a_module' is taken).
      void merge_module( PEL_Module* a_module ) ;

      // Remove the module identified by the `path_and_name' in the 
      // module hierarchy of root `self'.
      void remove_module( std::string const& path_and_name ) ;

      // Remove the leaf identified by the `path_and_name' in the 
      // module hierarchy of root `self'.
      void remove_entry( std::string const& path_and_name ) ;

   //-- Access
      
      // name of `self'
      std::string const& name( void ) const ;

      // absolute path-name of `self'
      std::string const& absolute_path_name( void ) const ;

      // Is `self' empty ?
      bool is_empty( void ) const ;

      // Is there a module in the 
      // module hierarchy of root `self' ? 
      bool has_module( void ) const ;

      // Is there a module identified by the path-name `path_and_name' in the 
      // module hierarchy of root `self' ? 
      bool has_module( std::string const& path_and_name ) const ;
      
      // module containing this if any (NULL otherwise)
      PEL_Module const* father( void ) const ;
      
      // first module whose path-name is `path_and_name' in the module
      // hierarchy of root `self' (if none, a fatal error is raised)
      PEL_Module* module( std::string const& path_and_name ) const ;
      
      // Is there an entry in the 
      // module hierarchy of root `self' ? 
      bool has_entry( void ) const ;

      // Is there an entry identified by the path-name `path_and_name' in the 
      // module hierarchy of root `self' ? 
      bool has_entry( std::string const& path_and_name ) const ;
      
      // data of the first entry identified by the path-name `path_and_name' 
      // in the module hierarchy of root `self' 
      // (if none, a fatal error is raised)
      PEL_Data const* data_of_entry( std::string const& path_and_name ) const ;

      // Make `new_data' be the data of the first  entry identified by the 
      // path-name `path_and_name' in the module hierarchy of root `self'
      // (if none, a fatal error is raised).
      void replace_data_of_entry( std::string const& path_and_name,
                                  PEL_Data const* new_data ) const ;

      // Create and return a list of entries such that the regular expression
      // `regexp' matches their path name relative to `self'.
      PEL_List* create_data_selection( PEL_Object* a_owner,
                                       std::string const& regexp,
                                       PEL_List* result,
                                       std::string const& where ) const ;
   //-- Context management
      
      // variable context associated to `self'
      PEL_Context const* context( void ) const ;
      
      // Add to `::context()' variable `variable' with default
      // `value'.
      void add_variable( PEL_Variable const* variable,
                         PEL_Data* value ) ;
      
      // Modify `::context()' variable `variable' with 
      // `value'.
      void modify_variable( PEL_Variable const* variable,
                            PEL_Data* value ) ;
      
   //-- Input - Output
      
      // Append the Module Hierarchy of root `self' to file whose name is
      // `file' with either "text" or "hybrid" `format'.
      // In hybrid format, whenever possible, each entry is replaced
      // by reference to an associated binary file named `file'.bin
      // The hybrid format is used :
      //   - to reduce file size ;
      //   - to ensure exact float recovering.
      // (binary files generated with hybrid format might not be readable on 
      // system different from these used to produce them).
      void write( std::string const& file,
                  std::string const& format ) const ;
      
      // IMPLEMENTATION : write the Module Hierarchy of root `self' with the
      // format of a PELICANS data file.
      virtual void print( std::ostream& os, size_t indent_width ) const ;

      virtual void display_info( std::ostream& os, size_t indent_width ) const ;

      // name of module beeing parsed: usefull in error treatment
      static std::string const& current_parsed_module_path_name( void ) ;

      // string representation of owned entry given by its path
      std::string data_as_string( std::string const& path_and_name,
                                  PEL_Context const* ct,
                                  bool & failed ) const ;
      
   //-- Iterators
      
      // Create and return an iterator on modules contained in `self'.
      PEL_ModuleIterator* create_module_iterator( PEL_Object* a_owner ) const ;

      // Create and return an iterator on entries contained in `self'.
      PEL_KeywordDataIterator* create_entry_iterator( PEL_Object* a_owner ) const ;

   //-- Path utilities

      // name of the module or entry defined by the path-name `path_and_name'
      // (similar to the unix command basename associated to file names)
      static std::string basename( std::string const& path_and_name ) ;
      
      // name of the module containing the module or entry defined by the
      // path-name `path_and_name'
      // (similar to the unix command dirname associated to file names)
      static std::string dirname( std::string const& path_and_name ) ;

      bool substitute_variables( std::string& replaced,
                                 PEL_Context const* ct ) const ;

      PEL_DataWithContext const* create_evaluation(
                                   PEL_Object * a_owner,
                                   std::string const& expression,
                                   PEL_Context const* ct ) const ;
   
   //-- Comparison

      // IMPLEMENTATION : `other' must be a `PEL_Module::' object
      // or a `PEL_String::' object
      virtual bool comparable( PEL_Object const* other ) const ;
      
      // IMPLEMENTATION : if `other' is a `PEL_Module::' object, it is equal to
      // to `self' if it has the same name ; if `other' is a `PEL_String::' 
      // object, it is equal to `self' if it represents its name
      virtual bool is_equal( PEL_Object const* other ) const ;
      
      // IMPLEMENTATION : compare `self' name to `other' name
      virtual int three_way_comparison( PEL_Object const* other ) const ;

      // IMPLEMENTATION : hash code of a `PEL_String::' object defined from
      // `self' name
      virtual size_t hash_code( void ) const ;

   protected: //-------------------------------------------------------
    	    	
   private: //-------------------------------------------------------

      PEL_Module( void ) ;
     ~PEL_Module( void ) ;
      PEL_Module( PEL_Module const& other ) ;
      PEL_Module const& operator=( PEL_Module const& other ) ;
  
      PEL_Module( PEL_Object* a_owner,
                  std::string const& a_name,
                  PEL_Context const* ct ) ;

      void complete_context( PEL_Module* root,
                             PEL_Module* dup,
                             PEL_Context const* ref_ctx ) const ;
      
      void complete_context( PEL_Module* root,
                             PEL_Data const* dat,
                             PEL_Context const* ref_ctx,
                             bool& has_modif ) const ;

      bool find( std::string const& nom,
                 PEL_Module*& theModule,
                 PEL_KeywordDataPair*& theAssignment ) const ;
      static void split( std::string const& nom,
                         std::string& token,
                         std::string& otherToken ) ;
      
      // Print entire tree and inserts tabulation for each level
      std::ostream& recursive_print( std::ostream& s ,
                                     int n,
                                     bool hybrid,
                                     std::string const& bin_file ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Attributes

      PEL_String* const NAME ;

      // Children structure:
      PEL_Module const* FATHER ;
      PEL_List* const MODS ;     // List of PEL_Module*
      PEL_List* const ENTRIES ;  // List of PEL_KeywordDataPair*

      // Context:
      PEL_ContextSimple* const CTX ;
      PEL_ContextPair* const TMP_CTX ;

};

#endif 
