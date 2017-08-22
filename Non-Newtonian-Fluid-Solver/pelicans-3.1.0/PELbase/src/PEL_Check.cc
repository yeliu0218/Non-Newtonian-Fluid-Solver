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

#include <PEL_Check.hh>

#include <PEL.hh>
#include <PEL_Exec.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_System.hh>
#include <PEL_assertions.hh>

#include <fstream>
#include <iostream>

PEL_Check const* PEL_Check::PROTOTYPE = new PEL_Check( "check" ) ;

/*
  Check validity of several data files from a pattern model.
  
  Two modes are available :
   - in interactive mode, application waits for a command line on standard
     input stream of the form : <filename> ;
   - else files to be checked are given on application command line.

   COMMAND LINE :
   Syntaxes recognized in command line mode are :
   check [-s] -i pattern.pel                        interactive mode
   check [-s] pattern.pel file1.pel ... filen.pel   non interactive mode

   DATA DECK :
   Data files recognized are :
   
   MODULE PEL_Application
      concrete_name = "check"
      pattern_filename = "pattern.pel"
      checked_files = < "file1.pel" file2.pel" ... > // Optionnal
      output_file = "expected.err"                   // Optionnal
   END MODULE PEL_Application
   interactive mode is actived when checked_files is not provided.
   FRAMEWORK INSTANTIATION
   Class wishing inherit from check must implement following methods :
    `do_check(a_owner, exp)' : method that process checking ;
                               
     Plug in methods : create_replica, create_replica_from_args, constructors.
*/

//----------------------------------------------------------------------
PEL_Check:: PEL_Check( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Application( a_name )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE()
   , SILENT( false )
   , NAME( a_name )
{
}

//----------------------------------------------------------------------
PEL_Check:: PEL_Check( PEL_Object* a_owner,
                       std::string const& a_name,
                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE( "" )
   , SILENT( false )
   , INTERACTIVE( false )
   , NAME( a_name )
{
   PEL_LABEL( "PEL_Check:: PEL_Check" ) ;
   PATTERN = exp->string_data( "pattern_filename" ) ;
   if(exp->has_entry("checked_files"))
   {
      MY_ARGS = exp->stringVector_data( "checked_files" ) ;
   }
   else
   {
      INTERACTIVE = true ;
   }
   if( exp->has_entry( "silent_mode" ) )
   {
      SILENT = exp->bool_data( "silent_mode" ) ;
   }
   
   if( exp->has_entry("output_file") )
   {
      FILE=exp->string_data( "output_file" ) ;
   }
}

//----------------------------------------------------------------------
PEL_Check:: PEL_Check( PEL_Object* a_owner,
                       std::string const& a_name,
                       stringVector& args )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE( "" )
   , SILENT( false )
   , INTERACTIVE( false )
   , NAME( a_name )
{
   PEL_LABEL( "PEL_Check:: PEL_Check" ) ;

   for( size_t i=0 ; i<args.size() ; i++ )
   {
      if(args(i)=="-s" )
      {
         SILENT = true ;
         args.remove_at(i) ;
         break ;
      }
      else if(args(i)=="-i" )
      {
         INTERACTIVE = true ;
         args.remove_at(i) ;
         break ;
      }
   }
   
   if( !INTERACTIVE && args.size() < 2 )
   {
      PEL_Error::object()->raise_plain(
         "usage : " + NAME + " pattern.pel data1.pel [datai.pel...]" ) ;
   }
   else if( INTERACTIVE && args.size() != 1 )
   {
      PEL_Error::object()->raise_plain(
         "usage : " + NAME + " -i pattern.pel " ) ;
   }
   PATTERN = args(0) ;
   args.remove_at(0) ;
   MY_ARGS = args ;
   args.re_initialize(0) ;
}

//----------------------------------------------------------------------
PEL_Check:: ~PEL_Check( void )
//----------------------------------------------------------------------
{
   if( this == PROTOTYPE )
   {
      PROTOTYPE = 0 ;  
   }
}

//----------------------------------------------------------------------
PEL_Check* 
PEL_Check:: create_replica( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Check:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_Check* result = new PEL_Check( a_owner, NAME, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Check* 
PEL_Check:: create_replica_from_args( PEL_Object* a_owner,
                                      stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Check:: create_replica_from_args" ) ;

   PEL_Check* result = new PEL_Check( a_owner, NAME, args ) ;

   PEL_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_Check:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Check:: run" ) ;
   if( PEL_ModulePattern::build_pattern() )
   {
      PEL_ModulePattern::save_pattern() ;
   }
   
   PEL_ModulePattern::open_pattern_base( PATTERN ) ;
   
   if( !SILENT )
      PEL_Exec::out() << "****** " << NAME << " output" << std::endl << std::endl ;

   if( !INTERACTIVE )
   {
      if( !SILENT )
         PEL_Exec::out() << NAME << " : interactive mode is off"
                         << std::endl << std::endl ;
      for( size_t i=0 ; i<MY_ARGS.size() ; i++ )
      {
         process( MY_ARGS(i) ) ;         
      }
   }
   else 
   {
      bool end = false ;
      PEL::out() << NAME << " : interactive mode is on" << std::endl ;
      PEL::out() << " Enter name of file to check (QUIT to quit)" << std::endl ;
      while( !end ) 
      {
         std::string a_file("QUIT") ;
      
         std::cin >> a_file ;
         if( a_file == "QUIT" ) 
         {
            end = true ;
         }
         else
         {
            process( a_file ) ;
            PEL::out() << "That's all folks" << std::endl ;
            PEL::out().flush() ;
         }
      }
   }

   PEL_ModulePattern::close_pattern_base() ;
}

//----------------------------------------------------------------------
void
PEL_Check:: process( std::string const& file_to_parse )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Check:: process" ) ;
   
   if( !SILENT )
      PEL_Exec::out() << "*** Checking " << file_to_parse << std::endl ;
      
   stringVector args(1) ;
   args(0) = file_to_parse ;
   PEL_ModuleExplorer* validity = 0 ;
   PEL_Module const* mod = PEL_Exec::create_module_with_data( this, args ) ;
   if( mod!=0 ) // reading successful
   {
      PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( 0, mod ) ;
      validity = do_check( 0, exp ) ;
      exp->destroy() ; exp=0 ;
   }
   else
   {
      std::string const& parsed =
         PEL_Module::current_parsed_module_path_name() ;
         
      if( !parsed.empty() )
      {
         PEL_Module* res = PEL_Module::create( this, "ParsingError" ) ;
         PEL_Error::notify( res, "Parse error", parsed ) ;
         validity = PEL_ModuleExplorer::create(0, res ) ;
      }
   }
   if( !validity->is_empty() )
   {
      if( FILE.empty() )
      {
         PEL_Exec::out() << std::endl ;
         PEL_Error::object()->display_data_checking( validity ) ;
      }
      else
      {
         std::ofstream out( FILE.c_str() ) ;
         if( !out )
            PEL_Error::object()->raise_plain( "Unable to open "+FILE ) ;
            
         validity->print( out, 0 ) ;
         out.close() ;         
      }
      PEL_Exec::set_exit_code( 1 ) ;
   }
   else
   {
      if(FILE.empty())
      {
         if( !SILENT )
         PEL_Exec::out() << std::endl << "Pattern checking is successful"
                         << std::endl ;
      } else 
      {
         PEL_System::erase( FILE ) ;
      }
      
   }
   
   validity->destroy() ;
      
   if( mod!=0 ) destroy_possession( mod ) ;
   mod = 0 ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_Check:: do_check( PEL_Object* a_owner,
                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Check:: do_check" ) ;
   PEL_CHECK_PRE( do_check_PRE( a_owner, exp ) ) ;

   PEL_Module* a_mod = exp->create_clone_of_attached_module(0) ;
   
   PEL_ModuleExplorer* ctrl =
      PEL_ModuleExplorer::create( 0,
                                  a_mod,
                                  PEL_ModuleExplorer::verify ) ;
   
   PEL_ModuleExplorer* result = ctrl->validity(a_owner) ;

   ctrl->destroy() ;
   a_mod->destroy() ;
   
   PEL_CHECK_POST( do_check_POST( a_owner, result ) ) ;
   return result ;
}


//----------------------------------------------------------------------
bool
PEL_Check:: do_check_PRE( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( exp!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Check:: do_check_POST( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->owner()==a_owner ) ;
   return true ;
}

