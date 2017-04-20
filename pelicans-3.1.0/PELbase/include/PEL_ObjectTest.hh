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

#ifndef PEL_OBJECT_TEST_HH
#define PEL_OBJECT_TEST_HH

#include <PEL_Object.hh>

#include <iosfwd>
#include <string>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Unit tests of a class.

Each instance is designed to perform a set of elementary tests

FRAMEWORK INSTANTIATION

   1. Derive a concrete subclass, say MyClassTest, devoted to perform
      unit tests of a given class, say "MyClass".
   2. Choose a name for MyClassTest, say "my_test".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the unique instance of MyClassTest to be registered :
      5.1 Implement a default constructor that initializes the
          `PEL_ObjectTest::' subobject by calling
             `PEL_ObjectTest( std::string const&, std::string const& )'
          with "MyClass" and "my_test" as arguments.
          Example of pseudo-code :
          | MyClassTest:: MyClassTest( void )
          |     : PEL_ObjectTest( "MyClass", "my_test" ) {}
      5.2 Define and initialize a static unique instance by calling the default
          constructor.
             declaration (in the header file, eg MyClassTest.hh) :
             | static MyClassTest const* UNIQUE_INSTANCE ;
             definition (in the implementation file, eg MyClassTest.cc) :
             | MyClassTest const* 
             | MyClassTest::UNIQUE_INSTANCE = new MyClassTest() ;
   6. If MyClassTest is to be used in conjonction with a data deck,
      implement `process_one_test' and possibly `::reset_all_tests'
      If not, implement `::process_all_tests'
      In any case, the result of each elementary test performed should
      be communicated by calling `::notify_one_test_result'

PUBLISHED */

class PEL_EXPORT PEL_ObjectTest : public PEL_Object
{
   public: //--------------------------------------------------------------

   //-- Instance delivery and initialization

      // an instance registered as `a_name'
      static PEL_ObjectTest* object( std::string const& a_name ) ;

   //-- Identification

      // registration name of `self'
      std::string const& registration_name( void ) const ;

   //-- Possible data deck

      // Indicate that the Hierarchical Data Structure accessible through
      // `exp' is a data deck for `self'.
      void set_data_deck_explorer( PEL_ModuleExplorer const* exp ) ;

      // Is there a data deck for `self' and an associated explorer ?
      bool has_data_deck_explorer( void ) const ;

      // explorer on the data deck of `self'
      PEL_ModuleExplorer const* data_deck_explorer( void ) const ;

   //-- Elementary tests execution

      /* Run all elementary tests. calling in sequence
            1. `::reset_all_tests'
            2. `::process_all_tests'. */
      void run_all_tests( void ) ;

      // Does the previous calls to `::run_all_tests()' on behalf of
      // any instance lead to successful tests ?
      static bool tests_of_all_instances_are_successful( void ) ;
            
   protected: //-----------------------------------------------------------

   //-- Plug-in

      virtual ~PEL_ObjectTest( void ) ;

      // Registration of an instance as `my_name' devoted to
      // testing the class called `tested_class_name'.
      PEL_ObjectTest( std::string const& tested_class_name,
                      std::string const& my_name ) ;
      
   //-- Elementary tests management
      
      // Perform the first step of `::run_all_tests'.
      // IMPLEMENTATION : do nothing.
      virtual void reset_all_tests( void ) ;
      
      // Perform the second step of `::run_all_tests'.
      // IMPLEMENTATION : call in sequence `::process_one_test' with an 
      // argument refering to all subexplorers accessible 
      // by `::data_deck_explorer()'.
      virtual void process_all_tests( void ) ;
      
      // Perform the steps of `::process_all_tests' where `exp'
      // refers  in sequence to all the subexplorers accessible
      // by `::data_deck_explorer()'.
      // IMPLEMENTATION : do nothing.
      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;
      
      // Notify that some test, denoted by `displayed_name', was
      // successful if `success' is true or failed if `success' is false.
      void notify_one_test_result( std::string const& displayed_name, 
                                   bool success ) ;

      void print_time_result( std::string const& name, double tt ) const ;

      void print_memory_result( std::string const& name, 
                                size_t memory_size ) const ;

      // stream used for outputs
      static std::ostream& out( void ) ;

   private: //-------------------------------------------------------------

      PEL_ObjectTest( void ) ;
      PEL_ObjectTest( PEL_ObjectTest const& other ) ;
      PEL_ObjectTest& operator=( PEL_ObjectTest const& other ) ;
      
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Class attributes

      static bool FAILURE ;

   //-- Attributes

      std::string const NAME ;
      std::string const TESTED_CLASS ;
      
      PEL_ModuleExplorer* EXP ;
      size_t NB_ELEMENTARY_TESTS ;
      size_t NB_ELEMENTARY_TESTS_OK ;
} ;

#endif 
