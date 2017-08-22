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

#include <PEL_RunTest.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Comparator.hh>
#include <PEL_Context.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_System.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <PEL_Variable.hh>
#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

using std::endl ;
using std::istringstream ;
using std::ostringstream ;
using std::string ;

PEL_RunTest const* PEL_RunTest::PROTOTYPE = new PEL_RunTest() ;

struct PEL_RunTest_ERROR
{
   static void n0( std::string const& f_name ) ;
   static void n1( std::string const& config_file ) ;
   static void n2( std::string const& config_file ) ;
   static void n3( std::string const& config_file, 
                   std::string const& nb_procs ) ;
   static void n4( std::string const& config_file, size_t nb_procs ) ;
} ;

//----------------------------------------------------------------------------
PEL_RunTest:: PEL_RunTest( void )
//----------------------------------------------------------------------------
   : PEL_Application( "peltest" )
   , VERBOSE( false )
   , TESTED_EXE( "" )
   , PELCMP_EXE( "" )
   , DATA_FILE( "" )
   , CONFIG_FILE( "" )
   , COMMON_ROOT( "" )
   , LIST( 0 )
   , DIR_IGNORED( 0 )
   , NAME( 0 )
   , DIR( 0 )
   , OPTIONS( 0 )
   , ACTION( ignore )
   , VERIFY( false )
   , PATTERN_FILE( "" )
   , TEST_RESULT( test_successful )
   , NB_SUCCESSFUL( PEL::bad_index() )
   , NB_AMBIGUOUS( PEL::bad_index() )
   , NB_FAILED( PEL::bad_index() )
   , GLOB_IGNORED( 0 )
   , RESU( "" )
   , MPI_RUN( "" )
   , DBL_EXACT( false )
   , MY_DBL_EPS( PEL::bad_double() )
   , MY_DBL_MIN( PEL::bad_double() )
{
}

//----------------------------------------------------------------------------
PEL_RunTest*
PEL_RunTest:: create_replica( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   PEL_RunTest* result = new PEL_RunTest( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_RunTest:: PEL_RunTest( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , VERBOSE( false )
   , TESTED_EXE( PEL_Exec::name_of_exe() )
   , PELCMP_EXE( PEL_Exec::name_of_exe() )
   , DATA_FILE( "data.pel" )
   , CONFIG_FILE( "config.pel" )
   , COMMON_ROOT( "" )
   , LIST( exp->stringVector_data( "test_directories" ) )
   , DIR_IGNORED( exp->has_entry( "ignored_directories" ) ?
         exp->stringVector_data( "ignored_directories" ) : 0 ) 
   , NAME( 0 )
   , DIR( 0 )
   , OPTIONS( 0 )
   , ACTION( ignore )
   , VERIFY( false )
   , PATTERN_FILE( "" )
   , TEST_RESULT( test_successful )
   , NB_SUCCESSFUL( PEL::bad_index() )
   , NB_AMBIGUOUS( PEL::bad_index() )
   , NB_FAILED( PEL::bad_index() )
   , GLOB_IGNORED( 0 )
   , RESU( "resu" )
   , MPI_RUN( "" )
   , DBL_EXACT( false )
   , MY_DBL_EPS( PEL::bad_double() )
   , MY_DBL_MIN( PEL::bad_double() )
   , MODE( PEL_String::create( 0, "" ) )
{
   PEL_LABEL( "PEL_RunTest:: PEL_RunTest" ) ;

   set_mpi_run() ;
   
   if( exp->has_entry( "executable" ) )
   {
      TESTED_EXE = exp->string_data( "executable" ) ;
      exp->test_file( "executable", "read" ) ;
   }
    
   if( exp->has_entry( "mpirun" ) )
   {
      MPI_RUN = exp->string_data( "mpirun" ) ;
      exp->test_file( "mpirun", "read" ) ;
   }
   
   if( exp->has_entry( "data_filename" ) )
   {
      DATA_FILE = exp->string_data( "data_filename" ) ;
   }
   
   if( exp->has_entry( "configuration_filename" ) )
   {
      CONFIG_FILE = exp->string_data( "configuration_filename" ) ;
   }
   
   if( exp->has_entry( "verbose" ) )
   {
      VERBOSE = exp->bool_data( "verbose" ) ;
   }

   if( exp->has_entry( "run_options" ) )
   {
      OPTIONS = exp->stringVector_data( "run_options" ) ;
   }
   
   if( exp->has_module( "pattern" ) ) 
   {
      PEL_ModuleExplorer const* sub_exp =
         exp->create_subexplorer( this, "pattern" ) ;
      PATTERN_FILE = PEL_System::absolute_path(
                            sub_exp->string_data( "pattern_filename" ) ) ;
      std::string const& type = sub_exp->string_data( "type" ) ;
      sub_exp->test_data_in( "type", "build,verify,build_then_verify" ) ;
      if( type == "build" )
      {
         sub_exp->test_file( "pattern_filename", "write" ) ;
         ACTION = build ;
      }
      else if( type == "verify" )
      {
         sub_exp->test_file( "pattern_filename", "read" ) ;
         ACTION = verify ;
      }
      else if( type == "build_then_verify" )
      {
         sub_exp->test_file( "pattern_filename", "write" ) ;
         ACTION = build_then_verify ;
      }
   }
   
   if( exp->has_module( "pel_compare" ) )
   {
      PEL_ModuleExplorer const* se = 
                         exp->create_subexplorer( 0, "pel_compare" ) ;
      PEL_Comparator::set_preferred_motifs_formats( se ) ;
      se->destroy() ;
   }   
   if( exp->has_entry( "output_directories" ) )
   {
      DIR = exp->stringVector_data( "output_directories" ) ;
      set_resu_dirs( LIST, COMMON_ROOT, NAME ) ;
      if( DIR.size() != LIST.size() )
      {
         PEL_Error::object()->raise_data_error(
            exp,
            "output_directories",
            "should have the same size than \"test_directories\" table" ) ;
      }     
   }
   else
   {
      std::string const default_root = "PELICANS_TESTS" ;
      set_resu_dirs( LIST, COMMON_ROOT, NAME ) ;
      DIR = NAME ;
      for( size_t i=0 ; i<DIR.size() ; ++i )
      {
         std::string d = DIR(i) ;
         DIR(i) = default_root+PEL_System::path_name_separator()+d ;
      }
      PEL_ASSERT( LIST.size() == DIR.size() ) ;
   }
   
   if( exp->has_entry( "files_to_ignore" ) )
   {
      GLOB_IGNORED = exp->stringVector_data( "files_to_ignore" ) ;
   }
   
   if( exp->has_entry( "output_file" ) )
   {
      RESU = exp->string_data( "output_file" ) ;
   }
   
   if( exp->has_module( "double_comparison" ) )
   {
      PEL_ModuleExplorer const* ee = 
                         exp->create_subexplorer( this ,"double_comparison" ) ;
      std::string const& tt = ee->string_data( "type" ) ;
      ee->test_data_in( "type", "PEL_double_equality,exact" ) ;
      if( tt == "PEL_double_equality" )
      {
         DBL_EXACT = false ;
         MY_DBL_EPS = ee->double_data( "dbl_eps" ) ;
         ee->test_data( "dbl_eps", "dbl_eps>=0.0" ) ;
         MY_DBL_MIN = ee->double_data( "dbl_min" ) ;
         ee->test_data( "dbl_min", "dbl_min>=0.0" ) ;
      }
      else if( tt == "exact" )
      {
         DBL_EXACT = true ;
      }
   }   
   PEL_Exec::add_variable_to_execution_context(
                      PEL_Variable::object( "SS_peltest_MODE" ),
                      MODE ) ;
}

//----------------------------------------------------------------------------
PEL_RunTest*
PEL_RunTest:: create_replica_from_args( PEL_Object* a_owner,
                                        stringVector& args ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: create_replica_from_args" ) ;
   
   PEL_RunTest* result = new PEL_RunTest( a_owner, args ) ;

   PEL_CHECK( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_RunTest:: PEL_RunTest( PEL_Object* a_owner, stringVector& args )
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , VERBOSE( false )
   , TESTED_EXE( PEL_Exec::name_of_exe() )
   , PELCMP_EXE( PEL_Exec::name_of_exe() )
   , DATA_FILE( "data.pel" )
   , CONFIG_FILE( "config.pel" )
   , COMMON_ROOT( "" )
   , LIST( 0 )
   , DIR_IGNORED( 0 )
   , NAME( 0 )
   , DIR( 0 )
   , OPTIONS( 0 )
   , ACTION( ignore )
   , VERIFY( false )
   , PATTERN_FILE( "" )
   , TEST_RESULT( test_successful )
   , NB_SUCCESSFUL( PEL::bad_index() )
   , NB_AMBIGUOUS( PEL::bad_index() )
   , NB_FAILED( PEL::bad_index() )
   , GLOB_IGNORED( 0 )
   , RESU( "resu" )
   , MPI_RUN( "" )
   , DBL_EXACT( false )
   , MY_DBL_EPS( PEL::bad_double() )
   , MY_DBL_MIN( PEL::bad_double() )
   , MODE( PEL_String::create( 0, "" ) )
{
   set_mpi_run() ;
   
   while( args.size() != 0 )
   {
      std::string arg = args( 0 ) ;
      args.remove_at( 0 ) ;
      if( arg == "-pattern_build" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         PATTERN_FILE = args( 0 ) ;
         args.remove_at( 0 ) ;
         ACTION = build ;
      }
      else if( arg == "-pattern_verify" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         PATTERN_FILE = args( 0 ) ;
         args.remove_at( 0 ) ;
         ACTION = verify ;
      }
      else if( arg == "-pattern_build_then_verify" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         PATTERN_FILE = args( 0 ) ;
         args.remove_at( 0 ) ;
         ACTION = build_then_verify ;
      }
      else if( arg == "-post" )
      {
         OPTIONS.append( "-Cpost" ) ;
      }
      else if( arg == "-verbose" )
      {
         VERBOSE = true ;
      }
      else if( arg == "-all" )
      {
         OPTIONS.append( "-Call" ) ;
      }
      else if( arg == "-test_directories" )
      {
         while( args.size() != 0 && args(0)[0] != '-' )
         {
            LIST.append( args(0) ) ;
            args.remove_at( 0 ) ;
         }
      }
      else if( arg == "-files_to_ignore" )
      {
         while( args.size() != 0 && args(0)[0] != '-' )
         {
            GLOB_IGNORED.append( args(0) ) ;
            args.remove_at( 0 ) ;
         }
      }
      else if( arg == "-output_file" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         RESU = args( 0 ) ;
         args.remove_at( 0 ) ;
      }
      else if( arg == "-output_directories" )
      {
         while( args.size() != 0 && args(0)[0] != '-' )
         {
            DIR.append( args(0) ) ;
            args.remove_at( 0 ) ;
         }
         set_resu_dirs( LIST, COMMON_ROOT, NAME ) ;
      }
      else if( arg == "-executable" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         TESTED_EXE = args( 0 ) ;
         args.remove_at( 0 ) ;
      }
      else if( arg == "-exact" )
      {
         DBL_EXACT = true ;
      }
      else if( arg == "-dbl_eps" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         std::istringstream is( args( 0 ) ) ;
         is >> MY_DBL_EPS ;         
         args.remove_at( 0 ) ;
         if( MY_DBL_EPS<0.0 ) 
         {
            notify_error_in_arguments() ;
         }
      }
      else if( arg == "-dbl_min" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         std::istringstream is( args( 0 ) ) ;
         is >> MY_DBL_MIN ;         
         args.remove_at( 0 ) ;
         if( MY_DBL_MIN<0.0 ) 
         {
            notify_error_in_arguments() ;
         }
      }
      else if( arg == "-mpirun" )
      {
         if( args.size()==0 ) notify_error_in_arguments() ;
         MPI_RUN = args( 0 ) ;
         args.remove_at( 0 ) ;
      }
      else
      {
         notify_error_in_arguments() ;
      }
   }
   
   if( LIST.size() == 0 )
   {
      notify_error_in_arguments() ;
   }
   if( DIR.size() == 0 )
   {
      std::string const default_root = "PELICANS_TESTS" ;
      set_resu_dirs( LIST, COMMON_ROOT, NAME ) ;
      DIR = NAME ;
      for( size_t i=0 ; i<DIR.size() ; ++i )
      {
         std::string d = DIR(i) ;
         DIR(i) = default_root+PEL_System::path_name_separator()+d ;
      }
   }
   
   if( ! ( ( MY_DBL_EPS==PEL::bad_double() && MY_DBL_MIN==PEL::bad_double() ) ||
           ( MY_DBL_EPS!=PEL::bad_double() && MY_DBL_MIN!=PEL::bad_double() ) ) )
      notify_error_in_arguments() ;
   if( DBL_EXACT && ( MY_DBL_EPS!=PEL::bad_double() || MY_DBL_MIN!=PEL::bad_double() ) )
      notify_error_in_arguments() ;      
   if( DIR.size() != LIST.size() )
   {
      notify_error_in_arguments() ;
   }

   PEL_Exec::add_variable_to_execution_context(
                      PEL_Variable::object( "SS_peltest_MODE" ),
                      MODE ) ;
}

//----------------------------------------------------------------------------
PEL_RunTest:: ~PEL_RunTest( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PEL_RunTest:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest::run" ) ;   

   out() << "executable : " << TESTED_EXE << std::endl ;
   out() << "tests      : " << COMMON_ROOT << std::endl ;

   NB_SUCCESSFUL = 0 ;
   NB_AMBIGUOUS = 0 ;
   NB_FAILED = 0 ;
   
   VERIFY = ( ACTION == verify ) ;
   do_tests() ;
   
   if( ACTION == build_then_verify )
   {
      VERIFY = true ;
      do_tests() ;
   }
   
   out() << std::endl ;
   std::string const stars( 50, '*' ) ;
   out() << stars << std::endl ;
   out() << "* number of successful tests: " << NB_SUCCESSFUL << std::endl ;
   out() << "* number of  ambiguous tests: " << NB_AMBIGUOUS << std::endl ;
   out() << "* number of    failing tests: " << NB_FAILED << std::endl ;
   out() << stars << std::endl ;
   out() << std::endl ;
}

//----------------------------------------------------------------------------
void
PEL_RunTest:: do_tests( void ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest::do_tests" ) ;
   if( VERIFY )
   {
      MODE->set( "verify" ) ;
   }
   else
   {
      MODE->set( "execution" ) ;
   }
   
   for( size_t i=0 ; i<LIST.size() ; i++ )
   {
      stringVector datas( 0 ) ;
      PEL_System::find( LIST( i ), DATA_FILE, datas, true ) ;
      if( VERBOSE )
         PEL::out() << "List of file "<<DATA_FILE<<" in "<<LIST(i)
                    <<std::endl<<datas<<std::endl<<std::endl  ;
      
      datas.sort() ;
      for( size_t j=0 ; j<datas.size() ; j++ )
      {
         std::string const& data_file = datas( j ) ;
         bool ignored = false ;
         for( size_t kk=0 ; kk<DIR_IGNORED.size() ; ++kk )
         {
            if( data_file.find( DIR_IGNORED( kk ) ) < data_file.length() )
            {
               ignored = true ;
               break ;
            }
         }
         
         if( !ignored ) do_test( LIST( i ), NAME( i ), DIR( i ), datas( j ) ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
PEL_RunTest:: do_test( std::string const& base,
                       std::string const& print_name,
                       std::string const& test_directory,
                       std::string const& filename ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest::do_test" ) ;
   PEL_CHECK( filename.find( base )==0 ) ;
   
   TEST_RESULT = test_successful ;
 
   std::string const tirets( 50, '-' ) ;
   std::string const rel = filename.substr( base.length()+1,
                                            filename.length()-base.length()-1 ) ;

   // Banner: begin test:
   out() << tirets << std::endl << "| " ;
   std::string nn = print_name ;
   if( PEL_System::dirname( rel ) != "." )
   {
      nn += PEL_System::path_name_separator() + PEL_System::dirname( rel ) ;
   }
   if( VERIFY )
   {
      out() << "verifying " << PEL_System::basename( filename )
            << " of " << nn << std::endl ;
   }
   else
   {
      out() << nn << std::endl ;
   }
   out() << tirets << std::endl ;

   // Configuration data base:
   PEL_ModuleExplorer const* conf =
                         create_config( 0, PEL_System::dirname( filename ) ) ;
   if( conf->has_entry( "condition_for_execution" ) )
   {
      bool do_it = conf->bool_data( "condition_for_execution",
                                    PEL_Exec::execution_context() ) ;
      if( !do_it )
      {
         out() << "| Test not executed" << std::endl << tirets << std::endl ;
         conf->destroy() ; conf = 0 ;
         return ; // <----------
      }
   }
   
   // Build command line:
   stringVector cmd( 0 ) ;
   std::string dir ;
   std::string output_file ;
      
   bool const parallel = conf->has_entry( "mpi_options" ) ;
   // mpi:
   if( parallel )
   {
      // Check mpi:
      if( MPI_RUN.empty() )
      {
         notify_test_failure( "Parallel test not allowed: no MPI tools defined" ) ;
      }
      cmd.append( MPI_RUN ) ;
      stringVector const& opts = conf->stringVector_data( "mpi_options" ) ;

      if( conf->has_entry( "mpi_machinefile" ) )
      {
         stringVector const& mm = 
            conf->stringVector_data( "mpi_machinefile" ) ;
         check_mpi_opts( mm, opts ) ;
         std::ofstream ff( "machines.tmp" ) ;
         if( !ff ) PEL_Error::object()->raise_file_handling( "machines.tmp",
                                                             "open" ) ;
         for( size_t i=0 ; i<mm.size() ; ++i )
            ff << mm(i) << std::endl ;
         cmd.append( "-machinefile" ) ;
         cmd.append( PEL_System::absolute_path( "machines.tmp" ) ) ;
      }

      for( size_t i=0 ; i<opts.size() ; ++i ) cmd.append( opts(i) ) ;
   }
   if( VERIFY ) 
   {
      cmd.append( PEL_System::absolute_path( PELCMP_EXE ) ) ;
      cmd.append( "-A" ) ;
      cmd.append( "check" ) ;
      cmd.append( "-s" ) ;
      cmd.append( PATTERN_FILE ) ;
      cmd.append( PEL_System::absolute_path( filename ) ) ;
      dir = PEL_System::dirname( PEL_System::absolute_path( filename ) ) ;
      output_file = "" ;
   }
   else
   {
      output_file = RESU ;
      dir = test_directory + PEL_System::path_name_separator()
                           + PEL_System::dirname( rel ) ;
      PEL_System::mkdir( dir ) ;
      stringVector oldfiles( 0 ) ;
      PEL_System::find( dir, "*", oldfiles, false ) ;
      for( size_t i=0 ; i<oldfiles.size() ; ++i )
      {
         PEL_System::erase( oldfiles( i ) ) ;
      }
      
      // Executable:
      {
         cmd.append( PEL_System::absolute_path( TESTED_EXE ) ) ;
      }
      // Data file:
      {
         cmd.append( PEL_System::absolute_path( filename ) ) ;
      }
      // Output file:
      if( parallel )
      {
         cmd.append( "-o" ) ;
         cmd.append( output_file ) ;
         output_file = "" ;
      }
      // Verbose flag:
      {
         cmd.append( "-v" ) ;
      }
      // Pattern built:
      if( !PATTERN_FILE.empty() )
      {
         cmd.append( "-build_pattern" ) ;
         cmd.append( PATTERN_FILE ) ;
      }
      // Additional options:
      for( size_t i=0 ; i<OPTIONS.size() ; ++i )
      {
         cmd.append( OPTIONS(i) ) ;
      }
      if( conf->has_entry( "run_options" ) ) 
      {
         stringVector const& extra = conf->stringVector_data( "run_options" ) ;
         for( size_t i=0 ; i<extra.size() ; ++i )
         {
            cmd.append( extra(i) ) ;
         }
      }
   }
   
   if( VERBOSE ) out() << "| Running : " << cmd << std::endl ;
   
   if( TEST_RESULT == test_successful )
   {
      int res = PEL_System::run( cmd, dir, output_file ) ;

      do_report( conf,
                 res,
                 PEL_System::dirname( PEL_System::absolute_path( filename ) ),
                 dir,
                 PEL_System::basename( filename ) ) ;
   }

   // Banner: test diagnostic
   if( TEST_RESULT == test_successful )
   {
      out() << "| Test is successful" << std::endl ;
      ++NB_SUCCESSFUL ;
   }
   else if( TEST_RESULT == test_failed )
   {
      
      out() << "| Test failed" << std::endl ;
      PEL_Exec::set_exit_code( 1 ) ;
      ++NB_FAILED ;
   }
   else  if( TEST_RESULT == test_ambiguous )
   {
      out() << "| Strict comparison failed" << std::endl ;
      out() << "| Test success ?... to be analyzed..." << std::endl ;
      if( PEL_Exec::exit_code() == 0 )
      {
         PEL_Exec::set_exit_code( 2 ) ;
      }
      ++NB_AMBIGUOUS ;
   }    
   out() << tirets << std::endl ;

   PEL_System::erase( "machines.tmp" ) ;
   conf->destroy() ; conf = 0 ;
}

//----------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_RunTest:: create_config( PEL_Object* a_owner,
                             std::string const& base ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest::create_config" ) ;
   PEL_CHECK( !base.empty() ) ;

   PEL_ModuleExplorer* result ;
   PEL_Module* mod_result ;
   
   std::string const conf_file =
               base + PEL_System::path_name_separator() + CONFIG_FILE ;
   
   std::ifstream file( conf_file.c_str() ) ;
   if( file.good() )
   {
      file.close() ;
      if( VERBOSE )
      {
         out() << "| reading configuration file: " <<  conf_file << std::endl ;
      }
      PEL_Module* mod = PEL_Module::create( 0, "main", conf_file ) ;
      PEL_ModuleIterator * it = mod->create_module_iterator(0) ;
      it->start() ;
      mod_result = it->item() ;
      result = PEL_ModuleExplorer::create( a_owner, mod_result ) ;
      mod->set_owner( result ) ;
      it->destroy() ; it = 0 ;
   }
   else
   {
      mod_result = PEL_Module::create( 0, "run_test" ) ;
      std::string const old_conf_file =
                 base + PEL_System::path_name_separator() + "pel_test.pl" ;
      std::ifstream old_file( old_conf_file.c_str() ) ;

      // Try to read obsolete perl configuration file:
      if( old_file.good() )
      {
         if( VERBOSE )
         {
            out() << "| reading configuration file: " <<  old_conf_file << std::endl ;
         }
         out() << "| WARNING: configuration file \"pel_test.pl\" is obsolete" << std::endl ;
         std::string line ;
         while( !old_file.eof() )
         {
            getline( old_file, line ) ;
            if( !old_file.eof() )
            {
               if( !old_file.good() )
               {
                  PEL_Error::object()->raise_plain(
                     "*** syntax error file: \""+old_conf_file+"\"" ) ;
               }
               if( line.find( "$pel_run_options" )==0 )
               {
                  size_t idx = line.find( "\"" ) ;
                  size_t jdx = line.rfind( "\"" ) ;
                  if( !( idx < jdx && jdx < line.length() ) )
                  {
                     PEL_RunTest_ERROR:: n0( old_conf_file ) ;
                  }
                  stringVector extra(1) ;
                  extra(0) = line.substr( idx+1, jdx-idx-1 ) ;
                  mod_result->add_entry(
                     "run_options",
                     PEL_StringVector::create( mod_result, extra ) ) ;
               }
               else if( line.find( "@files_to_ignore" )==0 )
               {
                  size_t idx = line.find( "(" ) ;
                  size_t jdx = line.rfind( ")" ) ;
                  if( !( idx < jdx && jdx < line.length() ) )
                  {
                     PEL_RunTest_ERROR:: n0( old_conf_file ) ;
                  }
                  std::string list = line.substr( idx+1, jdx-idx-1 ) ;
                  size_t dx = 0 ;
                  size_t odx = 0 ;
                  stringVector list_of_files(0) ;
                  while( ( dx = list.find( ",", odx ) ) < list.length() )
                  {
                     std::string afile = list.substr( odx, dx-odx ) ;
                     odx = dx+1 ;
                     size_t iidx = afile.find( "\"" ) ;
                     size_t jjdx = afile.rfind( "\"" ) ;
                     if( !( iidx < jjdx && jjdx < afile.length() ) )
                     {
                        PEL_RunTest_ERROR:: n0( old_conf_file ) ;
                     }
                     list_of_files.extend( afile.substr( iidx+1, jjdx-iidx-1 ) ) ;
                  }
                  std::string afile = list.substr( odx, list.length()-odx ) ;
                  size_t iidx = afile.find( "\"" ) ;
                  size_t jjdx = afile.rfind( "\"" ) ;
                  if( !( iidx < jjdx && jjdx < afile.length() ) )
                  {
                     PEL_RunTest_ERROR:: n0( old_conf_file ) ;
                  }
                  list_of_files.extend( afile.substr( iidx+1, jjdx-iidx-1 ) ) ;
                  if( !( list_of_files.size() > 0 ) )
                  {
                     PEL_RunTest_ERROR:: n0( old_conf_file ) ;
                  }
                  mod_result->add_entry(
                     "files_to_ignore",
                     PEL_StringVector::create( mod_result, list_of_files ) ) ;
               }
            }
         }
      }
      stringVector error(0) ;
      PEL_System::find( base, "expected.err*", error, false ) ;
      if( error.size() != 0 )
      {
         mod_result->add_entry( "failure_expected",
                                PEL_Bool::create( mod_result, true ) ) ;
      }
      result = PEL_ModuleExplorer::create( a_owner, mod_result ) ;
      mod_result->set_owner( result ) ;
   }

   PEL_CHECK_POST( result !=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PEL_RunTest:: do_report( PEL_ModuleExplorer const* conf,
                         int err_code,
                         std::string const& base_dir,
                         std::string const& test_dir,
                         std::string const& filename ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: do_report" ) ;
   bool failure = !VERIFY && ( conf->has_entry( "failure_expected" ) &&
                               conf->bool_data( "failure_expected" ) ) ;
   
   if( err_code!=0 && !failure )
   {
      std::ostringstream is ;
      is << "Test execution failed : (exit code " ;
      is << err_code << ")" ;
      notify_test_failure( is.str() ) ;
   }
   else if( err_code==0 && failure )
   {
      notify_test_failure( "Test execution must fail, and it does not !" ) ;
   }
   else if( !VERIFY ) 
   {
      stringVector generated(0) ;
      PEL_System::find( test_dir, "*", generated, false ) ;
      stringVector original(0) ;
      PEL_System::find( base_dir, "*", original, false ) ;
      stringVector ignored(0) ;
      if( conf->has_entry( "files_to_ignore" ) )
         ignored  = conf->stringVector_data( "files_to_ignore" ) ;
      bool resu_build = false ;
      for( size_t j=0 ; j<generated.size() ; ++j )
      {
         std::string fg = PEL_System::basename( generated(j) ) ;
         if( VERBOSE ) out() << "| checking file " << fg << std::endl ;
         if( ignored.has( fg ) )
         {
            ignored.remove_at( ignored.index_of( fg ) ) ;
            if( VERBOSE ) out() << "|    ignored" << std::endl ;
         }
         else if( GLOB_IGNORED.has( fg ) )
         {
            if( VERBOSE ) out() << "|    ignored" << std::endl ;
         }
         else if( fg == RESU || fg.find( RESU + "." )==0 )
         {
            resu_build = true ;
         }
         else
         {
            bool found = false ;
            for( size_t i=0 ; !found && i<original.size() ; i++ )
            {
               std::string fn = PEL_System::basename( original(i) ) ;
               found = ( fn == fg ) ;
               if( found )
               {
                  PEL_ModuleExplorer const* ee = specific_option( conf, fg ) ;
                  std::string format = "UNKNOWN" ;

                  PEL_Comparator::detect_file_format( ee, fg, format ) ;
                  do_compare( original(i), generated(j), ee, format ) ;
                  
               }
            }
            if( !found ) notify_test_failure( "Unknown generated file: "+fg ) ;
         }
      }
      if( !resu_build )
      {
         notify_test_failure( "No output file " + RESU + " generated" ) ;
      }
      for( size_t j=0 ; j<ignored.size() ; ++j )
      {
         notify_test_failure( "File not generated: "+ignored(j) ) ;
      }
   }
}

//----------------------------------------------------------------------------
PEL_ModuleExplorer const*
PEL_RunTest:: specific_option( PEL_ModuleExplorer const* conf,
                               std::string const& filename ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: specific_option" ) ;

   PEL_ModuleExplorer const* result = 0 ;

   if( conf->has_module( "PEL_Comparator" ) ) 
   {
      PEL_ModuleExplorer* e = conf->create_subexplorer( 0, "PEL_Comparator" ) ;
      e->start_module_iterator() ;
      for( ; e->is_valid_module() ; e->go_next_module() ) 
      {
         PEL_ModuleExplorer const* ee = e->create_subexplorer( e ) ;
         if( ee->string_data( "filename" ) == filename ) 
         {
            result = ee->create_clone( this ) ;
            break ;
         }
      }
      e->destroy() ;
   }

   PEL_CHECK( result==0 ||
              ( result->owner()==this &&
                result->string_data( "filename" )==filename ) ) ;
   return result ;
}

//----------------------------------------------------------------------------
void
PEL_RunTest:: do_compare( std::string const& first,
                          std::string const& second,
                          PEL_ModuleExplorer const* options,
                          std::string const& format ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: do_compare" ) ;
   stringVector cmd(5) ;
   cmd(0) = PEL_System::absolute_path( PELCMP_EXE ) ;
   cmd(1) = "-A";
   cmd(2) = "pelcmp"  ;
   cmd(3) = PEL_System::absolute_path( first ) ;
   cmd(4) = PEL_System::absolute_path( second ) ;
   
   if( VERBOSE ) out() << "|    format: " << format << endl ;

   if( MY_DBL_EPS != PEL::bad_double() ) 
   {
      {
         std::ostringstream os ;
         os << MY_DBL_EPS ;
         cmd.append( "-dbl_eps" ) ;
         cmd.append( os.str() ) ;
      }
      {
         PEL_ASSERT( MY_DBL_MIN != PEL::bad_double() ) ;
         std::ostringstream os ;
         os << MY_DBL_MIN ;
         cmd.append( "-dbl_min" ) ;
         cmd.append( os.str() ) ;
      }
   }
   if( VERBOSE )
   {
      cmd.append( "-verbose" ) ;
   }
   if( options != 0 ) 
   {
      if( options->has_entry( "ignore_data" ) ) 
      {
         stringVector const& id = options->stringVector_data( "ignore_data" ) ;
         for( size_t i=0 ; i<id.size() ; i++ )
         {
            cmd.append( "-ignore_data" ) ;
            cmd.append( id(i) ) ;
         }
      }
      if( !DBL_EXACT && MY_DBL_EPS == PEL::bad_double() && 
          options->has_module( "double_comparison" ) )
      {
         PEL_ModuleExplorer const* e = 
                   options->create_subexplorer( 0, "double_comparison" ) ;
         {
            std::ostringstream os ;
            os << e->double_data( "dbl_eps" ) ;
            cmd.append( "-dbl_eps" ) ;
            cmd.append( os.str() ) ;
         }
         {
            std::ostringstream os ;
            os << e->double_data( "dbl_min" ) ;
            cmd.append( "-dbl_min" ) ;
            cmd.append( os.str() ) ;
         }
         e->destroy() ;
      }
   }
   cmd.append( "-format" ) ;
   cmd.append( format ) ;
   if( VERBOSE )  out() << "|    do_compare : " << cmd << std::endl ;
   
   int res = PEL_System::run( cmd, ".", "cmp.out" ) ;
   if( res!=0 )
   {
      notify_test_ambiguous( PEL_System::basename(first) + 
                             " differs from the original" ) ;
      std::ifstream ares( "cmp.out" ) ;
      while( !ares.eof() )
      {
         std::string line ;
         getline( ares, line ) ;
         notify_test_ambiguous( line ) ;
      }
   }
   PEL_System::erase( "cmp.out" ) ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: notify_test_failure( std::string const& msg ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: notify_test_failure" ) ;
   
   out() << "| " << msg << std::endl ;
   TEST_RESULT = test_failed ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: notify_test_ambiguous( std::string const& msg ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: notify_test_ambiguous" ) ;
   
   out() << "| " << msg << std::endl ;
   TEST_RESULT = test_ambiguous ;
   
}

//---------------------------------------------------------------------------
std::ostream&
PEL_RunTest:: out( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: out" ) ;
   
   return PEL::out() ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << usage_title( "peltest" )  ;
   PEL::out() << "[options...]  -test_directories dir1 ... dirn" << endl ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << "     -executable <exe>" << endl ;
   PEL::out() << "          path name of the executable to be tested" 
              << endl << endl ;
   PEL::out() << "     -files_to_ignore file_1 ... file_n" << endl ;
   PEL::out() << "          generated files to ignore" 
              << endl << endl ;
   PEL::out() << "     -output_file resu_file [default \"resu\"]" << endl ;
   PEL::out() << "          name of the file from standard output of each "
              << "test running" << endl << endl ;
   PEL::out() << "     -output_directories dir1 ... dirn" << endl ;
   PEL::out() << "          path names of the directories containing the "
              << "result of the test cases" << endl << endl ;
   PEL::out() << "     -pattern_build pattern_file" << endl ;
   PEL::out() << "          build pattern file corresponding to all data files"
              << endl << endl ;
   PEL::out() << "     -pattern_verify pattern_file" << endl ;
   PEL::out() << "          verify all data files" 
              << endl << endl ;
   PEL::out() << "     -pattern_build_then_verify pattern_file" << endl ;
   PEL::out() << "          build pattern file corresponding to all data files"
              << endl << "          and then verify all data files" 
              << endl << endl ;
   PEL::out() << "     -all and -post" << endl ;
   PEL::out() << "          test cases additional option" 
              << endl << endl ;
   PEL::out() << "     -dbl_eps <value>" << endl
	<< "          argument a_dbl_eps in calls to PEL::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   PEL::out() << "     -dbl_min <value>" << endl
	<< "          argument a_dbl_min in calls to PEL::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   PEL::out() << "     -exact" << endl 
        << "          do not use any tolerance when comparing values of " 
        << "type double" << endl << endl ;
   PEL::out() << operands_title() ;
   PEL::out() << "     -test_directories dir1 ... dirn" << endl ;
   PEL::out() << "          path names of the directories containing the "
              << "test cases" << endl << endl ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: print_exit_status( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << exit_status_title() ;
   PEL::out() << "     0    All the tests where successful" << endl ;
   PEL::out() << "    >0    Some tests failed" << endl ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: set_resu_dirs( stringVector const& test_dirs,
                             std::string& common_root,
                             stringVector& resu_dirs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: set_resu_dirs" ) ;
   PEL_CHECK( test_dirs.size() != 0 ) ;

   common_root = "" ;
   resu_dirs = test_dirs ;
   for( size_t i=0 ; i<resu_dirs.size() ; ++i )
   {
      std::string d = PEL_System::absolute_path( resu_dirs(i) ) ;
      if( d[d.length()-1] == PEL_System::path_name_separator() )
      {
         d = d.substr( 0, d.length()-1 ) ;
      }
      resu_dirs(i)= d ;
   }

   if( resu_dirs.size() == 1 )
   {
      common_root = PEL_System::dirname( resu_dirs(0) ) ;
      resu_dirs(0) = PEL_System::basename( resu_dirs(0) ) ;
   }
   else
   {
      bool cont = true ;
      std::string root_dir1 ;
      std::string base_dir1 ;
      std::string root_dir2 ;
      std::string base_dir2 ;
      while( cont )
      {
         root_directory( resu_dirs(0), root_dir1, base_dir1 ) ;
         for( size_t i=1 ; cont && i<resu_dirs.size() ; ++i )
         {
            root_directory( resu_dirs(i), root_dir2, base_dir2 ) ;
            cont = ( root_dir2 == root_dir1 ) ;
         }
         if( cont )
         {
            common_root += PEL_System::path_name_separator() + root_dir1 ;
         }
         for( size_t i=0 ; i<resu_dirs.size() ; ++i )
         {
            root_directory( resu_dirs(i), root_dir2, base_dir2 ) ;
            resu_dirs(i) = "" ;
            if( !cont && !root_dir2.empty() )
            {
               resu_dirs(i) = root_dir2 ;
            }
            resu_dirs(i) += base_dir2 ;
         }
      }
   }
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: set_mpi_run( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: set_mpi_run" ) ;
   
   PEL_Context const* exe_ctx = PEL_Exec::execution_context() ;
   PEL_Variable const* with_mpi =  PEL_Variable::object( "BS_with_MPI" ) ;
   if( exe_ctx->has_variable( with_mpi ) &&
       exe_ctx->value( with_mpi )->to_bool() )
   {
      PEL_Variable const* mpi_run = PEL_Variable::object( "SS_MPI_RUN" ) ;
      if( !exe_ctx->has_variable( mpi_run ) )
      {
         PEL_Error::object()->raise_internal(
            "*** PEL_RunTest error:\n"
            "    variable SS_MPI_RUN not set" ) ;
      }
      MPI_RUN = exe_ctx->value( mpi_run )->to_string() ;
   }
}
   
//---------------------------------------------------------------------------
void
PEL_RunTest:: check_mpi_opts( stringVector const& machines,
                              stringVector const& mpi_options )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: check_mpi_opts" ) ;

   size_t idx = mpi_options.index_of( "-np" ) ;
   if( idx >= mpi_options.size() ) PEL_RunTest_ERROR::n1( CONFIG_FILE ) ;
   if( idx == mpi_options.size()-1 ) PEL_RunTest_ERROR::n2( CONFIG_FILE ) ;
   istringstream istr( mpi_options( idx+1 ) ) ;
   size_t nb_procs ;
   istr >> nb_procs ;
   if( !istr || !istr.eof() ) 
      PEL_RunTest_ERROR::n3( CONFIG_FILE, mpi_options( idx+1 ) ) ;
   if( machines.size() < nb_procs ) 
      PEL_RunTest_ERROR::n4( CONFIG_FILE, nb_procs ) ;
}

//---------------------------------------------------------------------------
void
PEL_RunTest:: root_directory( std::string const& dir,
                              std::string& root_dir,
                              std::string& base_dir )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest:: root_directory" ) ;
   PEL_CHECK( dir.length() > 0 ) ;
   
   size_t idx = dir.find_first_of( PEL_System::path_name_separator(), 1 ) ;
   if( idx < dir.length() )
   {
      root_dir = dir.substr( 1, idx-1 ) ;
      base_dir = dir.substr( idx, dir.length()-idx ) ;
   }
   else
   {
      root_dir = dir.substr( 1, dir.length()-1 ) ;
      base_dir = "" ;
   }
}

//internal------------------------------------------------------------------
void 
PEL_RunTest_ERROR:: n0( std::string const& f_name )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PEL_RunTest error:" << endl ;
   mesg << "    syntax error in the test configuration file" << endl ;
   mesg << "       \"" << f_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal------------------------------------------------------------------
void 
PEL_RunTest_ERROR:: n1( std::string const& config_file )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PEL_RunTest error:" << endl ;
   mesg << "    in the file " << config_file << endl ;
   mesg << "    the entry of keyword \"mpi_options\" should" << endl ;
   mesg << "    contain the \"-np\" item" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal------------------------------------------------------------------
void 
PEL_RunTest_ERROR:: n2( std::string const& config_file )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PEL_RunTest error:" << endl ;
   mesg << "    in the file " << config_file << endl ;
   mesg << "    and in the entry of keyword \"mpi_options\"" << endl ;
   mesg << "    the \"-np\" item should be followed by a string" << endl ;
   mesg << "       giving the number of processors" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal------------------------------------------------------------------
void 
PEL_RunTest_ERROR:: n3( std::string const& config_file,
                        std::string const& nb_procs )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PEL_RunTest error:" << endl ;
   mesg << "    in the file " << config_file << endl ;
   mesg << "    and in the entry of keyword \"mpi_options\" the" << endl ;
   mesg << "    \"-np\" item is followed the invalid item:" << endl ;
   mesg << "       \"" << nb_procs << "\"" << endl ;
   mesg << "    it should be followed by a string giving the" << endl ;
   mesg << "       number of processors" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal------------------------------------------------------------------
void 
PEL_RunTest_ERROR:: n4( std::string const& config_file,
                        size_t nb_procs )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PEL_RunTest error:" << endl ;
   mesg << "    in the file " << config_file << endl ;
   mesg << "    the number of processors requested by" << endl ;
   mesg << "    the entry of keyword \"mpi_options\" is " 
        << nb_procs << "," << endl ;
   mesg << "    hence the entry of keyword \"mpi_machinefile\"" << endl ;
   mesg << "    should contain a list of at least " << endl ;
   mesg << "    " << nb_procs << " machines (possibly duplicated)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

