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

#ifndef PEL_MODULE_EXPLORER_HH
#define PEL_MODULE_EXPLORER_HH

#include <PEL_Object.hh>
#include <PEL_Data.hh>

#include <string>

class PEL_Application ;
class PEL_Container ;
class PEL_Context ;
class PEL_ContextPair ;
class PEL_DataWithContext ;
class PEL_List ;
class PEL_Module ;
class PEL_ModuleIterator ;
class PEL_ModulePattern ;
class PEL_KeywordDataIterator ;
class boolVector ;
class doubleVector ;
class doubleArray2D ;
class boolArray2D ;
class stringArray2D ;
class doubleArray3D ;
class intArray2D ;
class intArray3D ;
class intVector ;
class stringVector ;

/*
Navigators to interrogate a database of the PELICANS Hierarchical Data System.

Each `PEL_ModuleExplorer::' object is attached to a `PEL_Module::' object
and offers ad hoc navigational and query facilities for accessing the data
of the associated Module Hierarchy.
Moreover, PEL_ModuleExplorer objects can be linked with a pattern recognition
mecanism. Two modes are then available :
build mode : pattern describing current HDS is build ;
verify mode : current HDS is compared to match a pattern and `validity' feature
 become available.
Since each `PEL_Module::' object defines a Module Hierarchy, each 
`PEL_ModuleExplorer::' object is said to be attached to a Module Hierarchy.

DESCRIPTION OF THE PELICANS HIERARCHICAL DATA SYSTEM
----------------------------------------------------

Data are organized into a hierarchy. The nodes of the hierarchy are called 
modules (hence the name: Module Hierarchy) while the leaves are the 
data themselves.

Each module may contain :
   - other modules, and 
   - entries represented by a keyword-data pair (thus, within a module,
     data are uniquely identified by a keyword).

Each module defines a Module Hierarchy (a module contains other modules
themselves containing other modules and so on).
The upper module of a Module Hierarchy is denoted the root module.
 
Given the root module of a Module Hierarchy, modules and entries can be 
retrieved knowing a path-name. A path-name is the name of a module or 
of the keyword of an entry preceeded by a concatenation of 
one or several module names separated with "/" . 
The starting node of a path is the root module if the first character is "/" 
(absolute path-name), otherwise (relative path-name) that
starting node is unspecified  (thus path-name looks like a unix file name). 

Example : in the Module Hierarchy defined by :
   MODULE mod1
     MODULE mod2
       k = "v"
     END MODULE mod2
   END MODULE mod1 
the keyword-data pair  {k,"v"} can be retrieved in module mod1 with either of 
the following path-names : 
        "/mod1/mod2/k" (absolute path-name)
        "mod2/k"       (relative path-name)
        "k"            (relative path-name)

Note that modules and entries might not be uniquely defined by
relative path-names. They are uniquely defined by absolute path-names.

Each `PEL_ModuleExplorer::' object is attached to a `PEL_Module::' object.
*/

class PEL_EXPORT PEL_ModuleExplorer : public PEL_Object
{
   public: //----------------------------------------------------------------

   //-- Module pattern

      enum PatternStatus { ignore, build, verify } ;
      
      // status of `self' relative to pattern recognition
      PatternStatus pattern_status( void ) const ;
      
      // pattern module used to describe pattern ( can be null if no pattern 
      // exist for associated module)
      PEL_ModulePattern* pattern( void ) const ;
      
      // compare `self' to underlying pattern and return diagnostic message
      //  if pattern is not respected
      // `recurse' is used to checked recursively all sub-modules
      PEL_ModuleExplorer* validity( PEL_Object* a_owner,
                                    bool recurse=true ) const ;
      
   //-- Instance delivery and initialization

      // Create and return an instance attached to `mm'.
      static PEL_ModuleExplorer* create( PEL_Object* a_owner,
                                         PEL_Module const* mm,
                                         PatternStatus status = ignore ) ;

      // Create and return an instance attached to the first module identified
      // by the path-name `path_and_name' in the attached Module Hierarchy.
      PEL_ModuleExplorer* create_subexplorer( 
                                 PEL_Object* a_owner,
                                 std::string const& path_and_name ) const ;

      virtual PEL_ModuleExplorer* create_clone( PEL_Object* a_owner ) const ;

   //-- Status report

      // name of the attached module
      std::string const& name( void ) const ;

      // absolute path-name of `self'
      std::string const& absolute_path_name( void ) const ;

      // owner of the attached module
      PEL_Object const* owner_of_attached_module( void ) const ;

      // Is the attached module empty ?
      bool is_empty( void ) const ;
      
      // Does the attached module contains a module of name `a_path_and_name' ?
      bool has_module( std::string const& a_path_and_name ) const ;

      // Does the attached module contains an entry whose keyword is 
      // `a_path_and_keyword' ?
      bool has_entry( std::string const& a_path_and_keyword ) const ;

  //-- Iterator on modules contained in the attached module
      
      // Move module iterator to first position.
      void start_module_iterator( void ) ;

      // Is module iterator position valid ?
      bool is_valid_module( void ) const ;

      // Move module iterator one position within the modules contained
      // in the attached module.
      void go_next_module( void ) ;

      // Create and return an instance attached to the module at
      // module iterator position.      
      PEL_ModuleExplorer* create_subexplorer( PEL_Object* a_owner ) const ;

  //-- Iterator on entries contained in the attached module
      
      // Move entry iterator to first position.
      void start_entry_iterator( void ) ;

      // Is entry iterator position valid ?
      bool is_valid_entry( void ) const ;

      // Move entry iterator one position within the entries contained in the
      // attached module.
      void go_next_entry( void ) ;

      // keyword of the entry at current entry iterator position
      std::string const& keyword( void ) const ;
      
      // data of the entry at current entry iterator position
      PEL_DataWithContext* data( PEL_Object* a_owner,
                                 PEL_Context const* ct = 0 ) const ;

  //-- Access to data contained in the attached module
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none)
      PEL_DataWithContext*
      abstract_data( PEL_Object* a_owner,
                     std::string const& a_path_and_keyword,
                     PEL_Context const* ct = 0 ) const  ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      bool bool_data( std::string const& a_path_and_keyword,
                      PEL_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      int int_data( std::string const& a_path_and_keyword,
                    PEL_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      double double_data( std::string const& a_path_and_keyword,
                          PEL_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      std::string const& string_data( std::string const& a_path_and_keyword,
                                      PEL_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intVector const&  intVector_data( std::string const& a_path_and_keyword,
                                        PEL_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleVector const& doubleVector_data( std::string const& a_path_and_keyword,
                                             PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleArray2D const& doubleArray2D_data( std::string const& a_path_and_keyword,
                                               PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      stringArray2D const& stringArray2D_data( std::string const& a_path_and_keyword,
                                               PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      boolArray2D const& boolArray2D_data( std::string const& a_path_and_keyword,
                                               PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleArray3D const& doubleArray3D_data( std::string const& a_path_and_keyword,
                                               PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intArray2D const& intArray2D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intArray3D const& intArray3D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      boolVector const& boolVector_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      stringVector const& stringVector_data( std::string const& a_path_and_keyword,
                                             PEL_Context const* ct = 0  ) const  ;

      PEL_Module* create_clone_of_attached_module( PEL_Object* a_owner ) const ;
      
   //-- Meta data information
      
      void set_help( std::string const& a_keyword,
                     std::string const& expression ) const ;
      void set_default( std::string const& a_keyword,
                        std::string const& expression ) const ;
      void test_data( std::string const& a_keyword,
                      std::string const& expression,
                      PEL_Context const* ct = 0 ) const ;
      void test_data_in( std::string const& a_keyword,
                         stringVector const& choices ) const ;
      void test_data_as( std::string const& a_keyword,
                         std::string const& regexp,
                         std::string const& where="true" ) const ;
      void test_file( std::string const& filename,
                      std::string const& mode ) const ;
      
   //-- Input - Output

      // IMPLEMENTATION : write the Module Hierarchy of root `self' with the
      // format of a PELICANS data file.
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      void write( std::string const& file,
                  std::string const& format ) const ;

      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Module pattern
      
      bool build_pattern( void ) const ;

      void declare_data( std::string const& a_path_and_keyword,
                         PEL_Data::Type type ) const ;
      
   //-- Access to context
      
      PEL_Context const* context( std::string const& a_path_and_name,
                                  PEL_Context const* ct ) const ;
      
      PEL_ModuleExplorer( void ) ;
     ~PEL_ModuleExplorer( void ) ;
      PEL_ModuleExplorer( PEL_ModuleExplorer const& other ) ;
      PEL_ModuleExplorer const& operator=(
          	          PEL_ModuleExplorer const& other ) ;

      PEL_ModuleExplorer( PEL_Object* a_owner,
                          PEL_Module const* mm,
                          PatternStatus status ,
                          PEL_ModulePattern* pat ) ;

      //---------------------------------------------------------------------
      //   ATTRIBUTES
      //---------------------------------------------------------------------
      PEL_Module const* mod ;
      PEL_ModuleIterator * submodule_iterator ;
      PEL_KeywordDataIterator * keyword_iterator ;
      mutable PEL_ContextPair* tmp_context ;
      
      PEL_ModulePattern* MP ;
      bool keyword_iterator_started ;
      PatternStatus my_status ;
} ;

#endif
