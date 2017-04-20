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

// LEVEL macro is undefined in PEL_assertions.hh
int const a_compilation_level = LEVEL ;

#include <PEL.hh>
#include <PEL_Exec.hh>

#include <PEL_assertions.hh>
#include <PEL_Application.hh>
#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Error.hh>
#include <PEL_Exceptions.hh>
#include <PEL_ExternalAPI.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_Object.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Timer.hh>
#include <PEL_Variable.hh>

#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>
using std::endl ;

//------------------------------------------------------------------------
bool PEL_Exec::VERBOSE = false ;
bool PEL_Exec::ONLY_HELP = false ;
bool PEL_Exec::SIGNAL_HANDLING = true ;
PEL_Timer* PEL_Exec::APPLI_TIMER = 0 ;
std::string PEL_Exec::EXE_FILE = "" ;
int PEL_Exec::EXIT_STATUS = 0 ;
PEL_ModuleExplorer::PatternStatus PEL_Exec::PATTERN_STATUS =
                                              PEL_ModuleExplorer::ignore ;
PEL_ContextSimple* PEL_Exec::EXEC_CONTEXT = 0 ;
std::istream* PEL_Exec::PEL_IN = 0 ;
std::ostream* PEL_Exec::PEL_OUT = 0 ;
std::ostream* PEL_Exec::PEL_ERR = 0 ;
bool PEL_Exec::EXTERNAL_API = true ;
//------------------------------------------------------------------------

//------------------------------------------------------------------------
size_t
PEL_Exec:: compilation_level( void )
//-------------------------------------------------------------------------
{
   return( (size_t) a_compilation_level ) ;
}

//------------------------------------------------------------------------
std::string const&
PEL_Exec:: date_of_compilation( void )
//-------------------------------------------------------------------------
{
   static std::string const result = __DATE__ ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
PEL_Exec:: name_of_exe( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: name_of_exe" ) ;
   return( EXE_FILE );
}

//------------------------------------------------------------------------
std::string
PEL_Exec:: dynamic_check_list( void )
//-------------------------------------------------------------------------
{
   std::string result ;
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Precondition ) )
   {
      result += "preconditions " ;
   }
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Postcondition ) )
   {
      result += "postconditions " ;
   }
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Invariant ) )
   {
      result += "invariants " ;
   }
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Check ) )
   {
      result += "checks " ;
   }
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      result += "objects " ;
   }
   
   return( result ) ;
}

//------------------------------------------------------------------------
int
PEL_Exec::exit_code( void )
//------------------------------------------------------------------------
{
   return( EXIT_STATUS ) ;
}

//------------------------------------------------------------------------
void
PEL_Exec::set_exit_code( int status )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec::set_exit_code" ) ;
   EXIT_STATUS = status ;
   PEL_CHECK_POST( exit_code() == status ) ;
}

//------------------------------------------------------------------------
PEL_Module const* 
PEL_Exec:: create_module_with_data( PEL_Object* a_owner,
                                    stringVector& args )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: create_module_with_data" ) ;
   PEL_Module* result = 0 ;
   try 
   { 
      PEL_CHECK_PRE( a_owner!=0 ) ; 

      if( !ONLY_HELP )
      {
         std::string dataFile = "" ;
         PEL_Module* main_module = 0 ;

         if( args.size()!=0 && args(0)=="-PEL_Application" )
         {
            // la structure de modules est écrite en ligne....
            PEL_ASSERT( args.size() > 1 ) ;
            args.remove_at(0) ;
            std::string api = args(0) ;
            args.remove_at(0) ;
            main_module = create_module_from_args( 0, api, args ) ;
            dataFile = "" ;
         }
         else if( args.size() == 1 )
         {
            dataFile = args(0) ;
            args.remove_at(0) ;
         }
      
         // Disiplay PELICANS header.
         if( VERBOSE ) print( dataFile, PEL_Exec::out() ) ;
   
         if( ( main_module==0 ) && ( !dataFile.empty() ) )
         {
            // data-deck reading and storage in memory
            main_module = PEL_Module::create( 0, "MAIN", dataFile,
                                              execution_context() ) ;
         }

         if( main_module!=0 ) 
         {
            PEL_ModuleIterator* it = main_module->create_module_iterator( 0 ) ;
            it->start() ;
            if( !it->is_valid() ) 
            {
               PEL_Error::object()->raise_plain(
                     "Empty file " + dataFile +"\nNo root module found" ) ;
            }
            
            result = it->item() ; it->go_next() ;
            if( it->is_valid() ) 
            {
               PEL_Error::object()->raise_plain(
                     "No only one root module found in "+dataFile ) ;
            }
            it->destroy() ;
            main_module->remove_module( result->name() ) ;
            main_module->change_owner( a_owner, result ) ;
            main_module->destroy() ; main_module = 0 ;
         }
      }

      PEL_CHECK_POST( IMPLIES( result!=0, result->owner() == a_owner ) ) ;
   }
   catch( PEL_Exceptions::Error )
   {
      result = 0 ;
   }
   return( result ) ;  
}

//------------------------------------------------------------------------
PEL_Application* 
PEL_Exec:: create_application( PEL_Object* a_owner,
                               PEL_Module const* appli_mod,
                               stringVector& args )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: create_application" ) ;
   
   PEL_Application* result = 0 ;

   if( !ONLY_HELP && appli_mod!=0 )
   {
      // creation of an instance of a concrete subclass of PEL_Application
      try
      {
         PEL_ModuleExplorer* exp =
            PEL_ModuleExplorer::create( a_owner,
                                        appli_mod,
                                        PATTERN_STATUS ) ;
         if( PATTERN_STATUS==PEL_ModuleExplorer::verify )
         {
            PEL_ModuleExplorer* validity = exp->validity(0) ;
            if( !validity->is_empty() )
            {
               PEL_Error::object()->display_data_checking(
                  validity ) ;
               PEL_Error::object()->raise_plain( "Pattern check aborts" ) ;
            }
            validity->destroy() ; validity = 0 ;
         }
         result = PEL_Application::make( a_owner, exp ) ;
      }
      catch( PEL_Exceptions::Error )
      {
         result = 0 ;
      }
   }
   else if( !ONLY_HELP && args.size()!=0 )
   {
      try
      {
         result = PEL_Application::make( a_owner, args ) ;
      }
      catch( PEL_Exceptions::Error )
      {
         result = 0 ;
      }
   }

   try
   {
      PEL_CHECK_POST( IMPLIES( result!=0, result->owner() == a_owner ) ) ;
   }
   catch( PEL_Exceptions::Error )
   {
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_Communicator const*
PEL_Exec:: communicator( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: communicator" ) ;
   static PEL_Communicator const* result = 0 ;
   if( result==0 )
   {
      std::string com_name = "PEL_SequentialCommunicator" ;
      PEL_Variable const* var = PEL_Variable::object( "BS_with_MPI" ) ;
      if( execution_context()->has_variable( var ) )
      {
         PEL_ASSERT( execution_context()->value( var )->to_bool() ) ;
         com_name = "EXT_MPIcommunicator" ;
      }
      result = PEL_Communicator::object( com_name ) ;
   }
   PEL_CHECK_POST( result != 0 ) ;
   
   // This post-condition verifies that no static objet needs communicator
   // before MPI one has been created.
   PEL_CHECK_POST( IMPLIES(
                      execution_context()->has_variable(PEL_Variable::object("BS_with_MPI") ),
                      result->name()=="EXT_MPIcommunicator" ) ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_Context const*
PEL_Exec:: execution_context( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: execution_context" ) ;
   
   static bool first = true ;
   if( first )
   {
      EXEC_CONTEXT = PEL_ContextSimple::create( PEL_Root::object() ) ;
      first = false ;
   }
   PEL_Context const* result = EXEC_CONTEXT ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==PEL_Root::object() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PEL_Exec:: add_variable_to_execution_context( PEL_Variable const* a_variable,
                                              PEL_Data* a_value )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: add_variable_to_execution_context" ) ;
   PEL_CHECK_PRE( a_variable!=0 ) ;
   PEL_CHECK_PRE( a_value!=0 && a_value->owner()==0 ) ;
   PEL_CHECK_PRE( !PEL_Exec::execution_context()->has_variable( a_variable ) ) ;

   // Creation of the context if needed :
   PEL_Context const* ct = PEL_Exec::execution_context() ;
   PEL_ASSERT( ct != 0 ) ;

   a_value->set_owner( EXEC_CONTEXT ) ;
   EXEC_CONTEXT->extend( a_variable, a_value ) ;

   PEL_CHECK_POST( PEL_Exec::execution_context()->has_variable( a_variable ) ) ;
   PEL_CHECK_POST( PEL_Exec::execution_context()->value( a_variable )==a_value ) ;  
}

//------------------------------------------------------------------------
std::ostream&
PEL_Exec::out( void )
//------------------------------------------------------------------------
{
   return( PEL_OUT==0 ? std::cout : *PEL_OUT ) ;
}

//------------------------------------------------------------------------
std::ostream&
PEL_Exec::err( void )
//------------------------------------------------------------------------
{
   return( PEL_ERR==0 ? std::cerr : *PEL_ERR ) ;
}

//------------------------------------------------------------------------
std::istream&
PEL_Exec::in( void )
//------------------------------------------------------------------------
{
   return( PEL_IN==0 ? std::cin : *PEL_IN ) ;
}

//------------------------------------------------------------------------
int
PEL_Exec:: initialize( int argc, char* argv[], stringVector& args )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: initialize" ) ;
   PEL_CHECK( args.size()==0 ) ;

   int result = 0 ;

   try
   {
      ONLY_HELP = false ;

      for( size_t i=0 ; i<(size_t)argc ; i++ )
      {
         std::string arg =  argv[i] ;
         
         if( arg == "-no_external_API" )
            EXTERNAL_API = false ;
      }
      if( EXTERNAL_API)  PEL_ExternalAPI::initialize_all_APIs( argc, argv ) ;

      args.re_initialize( argc ) ;
      for( size_t i=0 ; i<(size_t)argc ; i++ )
      {
         args(i) = argv[i] ;
      }

      // Parse command line arguments specific to PELICANS.
      parse_arguments( args ) ;

      if( SIGNAL_HANDLING ) PEL_System::exception_trapping() ;
      
      if( VERBOSE )
      {
         APPLI_TIMER = PEL_Timer::create( PEL_Root::object() ) ;
         APPLI_TIMER->start() ;
      }
      if( communicator()->nb_ranks()>1 )
      {
         std::cout << "Process #" << communicator()->rank()
                   << " : launched on " << PEL_System::host_name()
                   << " (" << PEL_System::process_id() << ")." <<  std::endl ;
      }
   }
   catch( PEL_Exceptions::Error )
   {
      result = 1 ;
   }

   return( result ) ;
}

//------------------------------------------------------------------------
void
PEL_Exec:: run_application( PEL_Application* appli  )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec::run_application" ) ;

   if( ONLY_HELP )
   {
      PEL_Exec::print_usage( PEL_Exec::err() ) ;
   }
   else
   {
      try
      {
         if( appli!=0 )
         {
            // program core exection : the tasks specific to the concrete 
            // application are performed
            appli->run() ;
         }
      }
      catch( PEL_Exceptions::Error )
      {
      }
   }
}

//------------------------------------------------------------------------
void
PEL_Exec:: terminate( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec::terminate" ) ;

   try
   {
      size_t const rank = communicator()->rank() ;
      size_t const nb_ranks = communicator()->nb_ranks() ;
   
      if( APPLI_TIMER !=0 )
      {
         APPLI_TIMER->stop() ;
         PEL_Exec::out() << "*** Elapsed time in second: " << std::endl ;
         PEL_Exec::out() << "user " << APPLI_TIMER->time() ;
         PEL_Exec::out() << std::endl << std::endl ;
      }
      if( PATTERN_STATUS == PEL_ModuleExplorer::build &&
          PEL_ModulePattern::build_pattern() )
      {
         PEL_ModulePattern::save_pattern() ;
      }
      // termination of all objects belonging to a ownership tree whose
      // root node is not the NULL object
      PEL_Root::cleanup() ;

      if(EXTERNAL_API) PEL_ExternalAPI::terminate_all_APIs() ;
      PEL_Exec::check_for_remaining_objects() ;

      if( nb_ranks>1 ) // communicator() is destroyed...
      {
         std::cout << "Process #" << rank << " : terminated." << std::endl ;
      }
   }
   catch( PEL_Exceptions::Error )
   {
      if( VERBOSE )
      {
         PEL_Exec::out() << "Unable to achieve PELICANS termination !!!"
                         << std::endl ;
      }   
   }

   out().flush() ;
   err().flush() ;
   if( PEL_OUT!=0 )
   {
      delete PEL_OUT ; PEL_OUT = 0 ;
   }
   if( PEL_ERR!=0 )
   {
      delete PEL_ERR ; PEL_ERR = 0 ;
   }
   if( PEL_IN!=0 )
   {
       delete PEL_IN ; PEL_IN = 0 ;
   }
}

//------------------------------------------------------------------------
void
PEL_Exec:: check_for_remaining_objects( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec::check_for_remaining_objects" ) ;
   if( PEL_Object::GetNumberOf_PEL_objects() > 0 )
   {
      PEL_Exec::out() << std::endl << std::endl ;
      PEL_Exec::out() << "after cleaning up, number of remaining P-objects : "
                 << PEL_Object::GetNumberOf_PEL_objects()
                 << std::endl ;
      PEL_Object::TraceRemainingObjects( PEL_Exec::out() ) ;
      set_exit_code( 2 ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Exec::print( std::string const& data, std::ostream& os ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: print" ) ;
   
   os << endl ;
   os << "*** Operating system: " <<
      PEL_System::sysname() << endl  ;
   os << endl ;
   os << "*** Executable: " << name_of_exe() << endl << endl ;

   if( !dynamic_check_list().empty() )
   {
      os << "*** Built-in tests that are dynamically enabled:" << endl ; 
      os << "       " << dynamic_check_list() << endl << endl ;
   }

   os << "*** Data file: " ;
   if( data== "-stdin" )
   {
      os << "read from standard input" ;
   }
   else
   {
      os << data ;
   }
   os << endl << endl;
   
   os << "*** PELICANS library" << endl ;
   os << "       compiler          : " << PEL_System::compiler_name() << endl ;
   os << "       compilation date  : " << date_of_compilation() << endl ;
   os << "       compilation level : opt" << compilation_level() << endl ;
   os << endl ;
}

//-------------------------------------------------------------------------
void
PEL_Exec:: parse_arguments( stringVector& args )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: parse_arguments" ) ;

   EXE_FILE = args(0) ;
   args.remove_at( 0 ) ;
   size_t iArgs=0 ;
   while( iArgs<args.size() )
   {
      std::string const& argument = args(iArgs) ;
      bool done = true ;
      // Here take place loop for argument scanning
      if( argument == "-Cpost" )
      {
         PEL_Assertion::add_handled_check( PEL_Assertion::Postcondition ) ;
      }
      else if( argument == "-Cobjects" )
      {
         PEL_Assertion::add_handled_check( PEL_Assertion::Objects ) ;
      }
      else if( argument == "-catch" )
      {
         PEL_Assertion::add_handled_check( PEL_Assertion::Objects ) ;
         args.remove_at( iArgs ) ;
         std::istringstream is( args(iArgs) ) ;
         size_t r ;
         is >> r ;
         
         PEL_Object::catch_object_by_rank( r ) ;
      }
      else if( argument == "-build_pattern" )
      {
         args.remove_at( iArgs ) ;
         PEL_ModulePattern::build_pattern_base( args(iArgs).c_str() ) ;
         PATTERN_STATUS = PEL_ModuleExplorer::build ;
      }
      else if( argument == "-check_pattern" )
      {
         args.remove_at( iArgs ) ;
         PEL_ModulePattern::build_pattern_base( args(iArgs).c_str() ) ;
         PATTERN_STATUS = PEL_ModuleExplorer::verify ;
      }
      else if( argument == "-Call" )
      {
         PEL_Assertion::add_handled_check( PEL_Assertion::Postcondition ) ;
         PEL_Assertion::add_handled_check( PEL_Assertion::Invariant ) ;
         PEL_Assertion::add_handled_check( PEL_Assertion::Check ) ;         
      }
      else if( argument == "-v" )
      {
         VERBOSE = true ;
      }
      else if( argument == "-o" )
      {
         args.remove_at( iArgs ) ;
         std::string outfile = args(iArgs) ;
         if( communicator()->nb_ranks()>1 )
         {
            std::ostringstream os ;
            os << "." <<  (int)communicator()->rank() ;
            outfile += os.str() ;
         }
         PEL_OUT = new std::ofstream( outfile.c_str() ) ;
         if( !(*PEL_OUT) )
         {
            PEL_Error::object()->raise_plain(
               "Unable to open standard output file "+outfile ) ;
         }
      }
      else if( argument == "-H" )
      {
         ONLY_HELP = true ;
      }
      else if( argument == "-no_signal_handling" )
      {
         SIGNAL_HANDLING = false ;
      }
      else if( argument == "-no_external_API" )
      {
         // done
      }
      else if( argument == "-notify_parallel" )
      {  
         PEL_Context const* ct = PEL_Exec::execution_context() ;
         if( !ct->has_variable( PEL_Variable::object( "BS_with_MPI" ) ) )
         {
            PEL_Error::object()->raise_plain(
               "Parallel execution not allowed: no MPI tools defined" ) ;
         }
      }
      else
      {
         done = false ;
      }
      if( done )
      {
         args.remove_at( iArgs ) ;
      }
      else
      {
         iArgs++ ;
      }
   }
   if( args.size() == 0 )
   {
      ONLY_HELP = true ;
   }
}

//------------------------------------------------------------------------
PEL_Module* 
PEL_Exec:: create_module_from_args( PEL_Object* a_owner,
                                    std::string const& application_name,
                                    stringVector& args )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Exec:: create_module_from_args" ) ;
   
   std::stringstream str ;
   str << "MODULE PEL_Application" << std::endl ;
   str << "concrete_name=\"" << application_name << "\"" << std::endl ;
   size_t i=0 ;
   
   while( i<args.size() )
   {
      std::string const& arg = args(i) ;
      size_t idx = arg.find_first_of( "=" ) ;
      if( idx<arg.length() )
      {
         str << arg << std::endl;
         args.remove_at(i) ;
      }
      else
      {
         i++ ;
      }
   }
   str << "END MODULE PEL_Application" << std::endl ;
   
   PEL_Module* result = PEL_Module::create( a_owner, "MAIN", str ) ;
   
   return( result ) ;
}

//------------------------------------------------------------------------
void
PEL_Exec:: print_usage( std::ostream& os )
//------------------------------------------------------------------------
{
   os << "Usage: " << endl ;
   os << "   [executable] [-H] [-V] [-Cobjects] [-catch <obj_nb>] "
      << "[-Cpost|-Call] [ -build_pattern <file> ] ( data_file | -stdin )" 
      << endl << endl ;
   os << "Options:" << endl ;
   os << "   [-h] : Display help and exit" << endl ;
   os << "   [-v] : verbosity"  << endl ;
   os << "   [-Cobjects]    : orphean objects tracking" << endl ;
   os << "   [-catch <obj_nb>]: object creation catching" << endl ;
   os << "   [-Cpost|-Call] : assertion level modification" << endl ;
   os << "   [ -build_pattern <file> ] : extend given pattern" << endl ;
   os << "   [ -check_pattern <file> ] : check from given pattern" << endl ;
   os << "   [ -no_signal_handling ] : doesn't provide signal handling facility" << endl ;
   os << "   [ -no_external_API ] : doesn't provide extra libraries management" << endl ;
   os << endl << "Arguments:" << endl ;
   os << "   ( data_file | -stdin ) : data file or standard input"<< endl ;
}
