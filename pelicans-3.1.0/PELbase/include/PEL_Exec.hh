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

#ifndef PEL_EXEC_HH
#define PEL_EXEC_HH

#include <cstddef>
#include <iosfwd>

#include <PEL_ModuleExplorer.hh>

class stringVector ;

class PEL_Application ;
class PEL_Communicator ;
class PEL_Context ;
class PEL_ContextSimple ;
class PEL_Module ;
class PEL_Object ;
class PEL_Timer ;

class PEL_EXPORT PEL_Exec
{
   public: //-----------------------------------------------------------------
      
   //-- PELICANS library self-description
      
      // Compilation level of PELICANS library.
      static size_t compilation_level( void ) ;

      // Date when was compiled PELICANS library currently used.
      static std::string const& date_of_compilation( void ) ;
      
      // Name of executable.
      static std::string const& name_of_exe( void ) ;

      // List of dynamic checks activated during execution.
      static std::string dynamic_check_list( void ) ;

   //-- PELICANS core execution

      // Initialize platform from argument list `args', 
      // and return 0 if successful.
      // Arguments not recognized by PELICANS are set in args.
      static int initialize( int argc, char * argv[],
                              stringVector& args ) ;

      // Terminate PELICANS session.
      static void terminate( void ) ;

      // Return to system exit code.
      // >=1 : exit on user-defined user error
      // 1 : exit on default user error
      // 0 : normal exit
      // -1 : exit on internal error
      // -2 : exit on internal error du to remaining objects after termination.
      static int exit_code( void ) ;
      
      // Modify exit code.
      static void set_exit_code( int status ) ;
      
   //-- `::PEL_Application' factory
      
      // Create application module from file whose name is given in `args'
      //  remaining argument list.
      static PEL_Module const* create_module_with_data( PEL_Object* a_owner,
                                                        stringVector& args ) ;

      // Create application from module `mod'
      static PEL_Application* create_application( PEL_Object* a_owner,
                                                  PEL_Module const* appli_mod,
                                                  stringVector& args ) ;
      
      // Launch application `appli'.
      static void run_application( PEL_Application* appli  ) ;

   //-- Parallel tools

      static PEL_Communicator const* communicator( void ) ;

   //-- Context of execution

      static PEL_Context const* execution_context( void ) ;
      static void add_variable_to_execution_context(
                          PEL_Variable const* a_variable, PEL_Data* a_value ) ;
      
   //-- Input - Output
      
      // Pelicans standard output stream.
      static std::ostream& out( void ) ;
      
      // Pelicans standard error stream.
      static std::ostream& err( void ) ;
      
      // Pelicans standard input stream.
      static std::istream& in( void ) ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      PEL_Exec( void ) ;
     ~PEL_Exec( void ) ;
      PEL_Exec( PEL_Exec const& other ) ;
      PEL_Exec& operator=( PEL_Exec const& other ) ;

      // Check for remaining objects.
      static void check_for_remaining_objects( void ) ;

      static void print( std::string const& data,
                         std::ostream& os ) ;
      
      static void parse_arguments( stringVector& args ) ;
      static PEL_Module* create_module_from_args(
                                     PEL_Object* a_owner,
                                     std::string const& application_name,
                                     stringVector& args ) ;
      
      static void print_usage( std::ostream& os ) ;
      
      static bool VERBOSE ;
      static bool ONLY_HELP ;
      static bool SIGNAL_HANDLING ;
      static bool EXTERNAL_API ;
      
      static PEL_Timer* APPLI_TIMER ;
      
      static std::string EXE_FILE ;

      static int EXIT_STATUS ;
      
      static PEL_ModuleExplorer::PatternStatus PATTERN_STATUS ;

      static PEL_ContextSimple* EXEC_CONTEXT ;

      static std::istream* PEL_IN ;
      static std::ostream* PEL_OUT ;
      static std::ostream* PEL_ERR ;
      
} ;

#endif
