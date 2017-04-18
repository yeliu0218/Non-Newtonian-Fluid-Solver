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

#ifndef PEL_APPLICATION_HH
#define PEL_APPLICATION_HH

#include <PEL_Object.hh>
#include <string>
#include <stringVector.hh>

class PEL_ListIdentity ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Applications, performing their specific tasks.

The program execution consists of five stages :
  1. Initial stage (Big-Bang time) : all statics are initialized and
     the only instance of PEL_Root is created.
  2. The data deck is read and stored in memory.
  3. An instance of a concrete subclass of PEL_Application is created.
  4. Program core execution : the program execution proceeds by performing
     its specific tasks.
  5. Final stage : termination of the only instance of PEL_Root, leading 
     to the termination of all objects belonging to a ownership tree whose
     root node is not the NULL object.
The PEL_Application class provides an interface for executing the 
specific tasks of the above fouth point.

FRAMEWORK INSTANTIATION

   CASE 1 : derivation of a concrete subclass

   1. Derive a concrete subclass, say MyAppli.
   2. Choose a name for MyAppli, say "my_appli".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `PEL_Application::' subobject by calling
               `PEL_Application( std::string const& )'
          with "my_appli" as argument.
          Example of pseudo-code :
          | MyAppli:: MyAppli( void ) : PEL_Application( "my_appli" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyAppli.hh) :
             | static MyAppli const* PROTOTYPE ;
             definition (in the implementation file, eg MyAppli.cc) :
             | MyAppli const* MyAppli::PROTOTYPE = new MyAppli() ;'
   6. Implement a private constructor that initializes the 
      `PEL_Application::' subobject by calling
                 `PEL_Application( PEL_Object*, PEL_ModuleExplorer const* )'
      or
                 `PEL_Application( PEL_Object*, stringVector& )'
      Example of pseudo-code :
      | MyAppli:: MyAppli( PEL_Object* a_owner,
      |                    PEL_ModuleExplorer const* exp )
      |    : PEL_Application( a_owner, exp ), ...
      | { ... }
   7. Implement the `::create_replica' method that allocates an object
      of type `MyAppli' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | MyAppli* MyAppli::create_replica( PEL_Object* a_owner,
      |                                   PEL_ModuleExplorer const* exp ) const
      | {
      |    PEL_LABEL( "MyAppli::create_replica" ) ;
      |    PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
      |    MyAppli* result = new MyAppli( a_owner, exp ) ;
      |    PEL_CHECK( create_replica_POST( result, a_owner, exp ) ;
      |    return result ;
      | }
   8. Implement the `::run' method

   CASE 2 : derivation of an abstract subclass

   1. Derive an abstract subclass, say MyAppli.
   2. Implement a protected virtual destructor.
   3. Implement a protected constructor that initializes the
      `PEL_Application::' subobject by calling
               `PEL_Application( std::string const& )'
      Example of pseudo-code :
      | MyAppli:: MyAppli( std::string const& name ) 
      |    : PEL_Application( name ) {}
      This constructor is devoted to be used by the concrete subclasses 
      of MyAppli for the registration of their prototype.
   4. Implement a protected constructor that initializes the 
      `PEL_Application::' subobject by calling
                 `PEL_Application( PEL_Object*, PEL_ModuleExplorer const* )'.
      Example of pseudo-code :
      | MyAppli:: MyAppli( PEL_Object* a_owner,
      |                    PEL_ModuleExplorer const* exp )
      |    : PEL_Application( a_owner, exp ), ...
      | { ... }
      This constructor is devoted to be used to initialize the MyAppli
      base class subobject when creating objects of concrete subclasses
      of MyAppli (such creations are performed in the `create_replica::'
      method whose implementation is deferred into those concrete subclasses).
*/

class PEL_EXPORT PEL_Application : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance of `PEL_Application::' according
      // to the data attainable by `exp'.
      static PEL_Application* make( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) ;

      // Create and return an instance according to the command-line
      // arguments gathered in `args'. 
      static PEL_Application* make( PEL_Object* a_owner,
                                    stringVector& args ) ;

   //-- Program core execution(0.2)

      // Perform the specific tasks of the application (Called by main()).
      virtual void run( void ) = 0 ;

   //-- Persistence

      void register_storable_objects( void ) ;

      void write_storable_objects( void ) const ;

      void restore_registered_objects( PEL_ObjectReader* ret ) const ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PEL_Application( void ) ;

      // Registration of an instance as `name'.
      PEL_Application( std::string const& name ) ;

      // In the constructor called by `::create_replica' or 
      // `::create_replica_from_args' : initialization the base class subobject
      // (`exp' can possibly be 0 ).
      PEL_Application( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp ) ;

      virtual PEL_Application* create_replica( 
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const = 0 ;

      // IMPLEMENTATION : raise a fatal error.
      virtual PEL_Application* create_replica_from_args( 
                                   PEL_Object* a_owner,
                                   stringVector& args ) const ;

      bool is_a_prototype( void ) const ;

   //-- Command line(1010.0)

      void notify_error_in_arguments( void ) const ;

      virtual void print_usage( void ) const ;
      virtual void print_options( void ) const ;
      virtual void print_operands( void ) const ;
      virtual void print_exit_status( void ) const ;

      std::string usage_title( std::string const& name ) const ;
      std::string options_title( void ) const ;
      std::string operands_title( void ) const ;
      std::string exit_status_title( void ) const ;

   //-- Preconditions, Postconditions, Invariant    

      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( PEL_Application const* result,
				PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) const ;

      bool create_replica_from_args_POST( PEL_Application const* result,
		         		  PEL_Object* a_owner,
				          stringVector& args ) const ;

      virtual bool invariant( void ) const ;

   //-- Persistence
      
      // Extend `list' (with the `PEL_ListIdentity::extend' method) so that it
      // contains all objects required by the storage and retrieval mechanisms.
      // IMPLEMENTATION : do nothing, i.e. leave `list' unchanged.
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;

      // name of the module containing the data related to
      // the storage mechanism (of persistence)
      // IMPLEMENTATION : "PEL_ObjectWriter"
      virtual std::string const& object_writer_module_name( void ) const ;

   private: //----------------------------------------------------------

      PEL_Application( void ) ;
      PEL_Application( PEL_Application const& other ) ;
      PEL_Application& operator=( PEL_Application const& other ) ;

      void print_usage_then_exit( int exit_status = 0 ) const ;

      void initialize_objects_storage( PEL_ModuleExplorer const* exp ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const IS_PROTO ;

      PEL_ObjectWriter* SAVER ;

      // List of the persistent objects :
      PEL_ListIdentity* persistent_objects ;
      mutable PEL_ListIterator* persistent_objects_it ;
} ;

#endif



