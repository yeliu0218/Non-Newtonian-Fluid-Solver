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

#include <FE_Launcher.hh>

#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::string ; using std::ostringstream ; using std::ifstream ;

struct FE_Launcher_ERROR {
   static void n0( string const& filename ) ;
   static void n1( PEL_ModuleExplorer const* exp, size_t nbv ) ;
   static void n2( void ) ;
   static void n3( std::string const& name,
                   std::string const& type ) ;
} ;

FE_Launcher const* FE_Launcher::PROTOTYPE = new FE_Launcher( "FE_Launcher" ) ;

//----------------------------------------------------------------------------
FE_Launcher:: FE_Launcher( std::string const& name )
//----------------------------------------------------------------------------
   : PEL_Application( name )
   , PATTERN( PEL_ModuleExplorer::ignore )
   , CONTROLER_FILE()
   , DATA_FILE()
   , NB_CALC( 0 )
{
}

//----------------------------------------------------------------------------
FE_Launcher*
FE_Launcher:: create_replica( PEL_Object* a_owner, 
                              PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Launcher:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   std::string controler_file = "" ;
   if( exp->has_entry( "controler_file" ) )
   {
      controler_file = exp->string_data( "controler_file" ) ;
      exp->test_file( "controler_file", "read" ) ;
      if( exp->pattern_status() != PEL_ModuleExplorer::ignore )
      {
         controler_file = "" ;
      }
   }
   
   FE_Launcher* result = new FE_Launcher( a_owner, exp, controler_file ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
FE_Launcher:: FE_Launcher( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp,
                           std::string const& controler_file  )
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , PATTERN( exp->pattern_status() )
   , CONTROLER_FILE( controler_file )
   , DATA_FILE( exp->string_data( "data_file" ) )
   , NB_CALC( exp->int_data( "nb_calculations" ) )
{
   PEL_LABEL( "FE_Launcher:: FE_Launcher" ) ;
   PEL_CHECK_PRE(
      IMPLIES( exp->pattern_status() != PEL_ModuleExplorer::ignore,
               controler_file.empty() ) ) ;

   ifstream ff( DATA_FILE.c_str(), std::ios::in  ) ;
   if( !ff ) FE_Launcher_ERROR::n0( DATA_FILE ) ;
   ff.close() ;

   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "variables" ) ;
   size_t iv = 0 ;
   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      ++iv ;
   }
 
   TYPE.resize( iv, PEL_Data::Undefined ) ;
   BOOL_VALUES.resize( iv, 
                       boolVector( exp->int_data( "nb_calculations" ) ) ) ;
   BOOL_DATA.resize( iv, (PEL_Bool*)0 ) ;
   INT_VALUES.resize( iv, 
                      intVector( exp->int_data( "nb_calculations" ) ) ) ;
   INT_DATA.resize( iv, (PEL_Int*)0 ) ;
   DBL_VALUES.resize( iv, 
                      doubleVector( exp->int_data( "nb_calculations" ) ) ) ;
   DBL_DATA.resize( iv, (PEL_Double*)0 ) ;
   STRING_VALUES.resize( iv, 
                         stringVector( exp->int_data( "nb_calculations" ) ) ) ;
   STRING_DATA.resize( iv, (PEL_String*)0 ) ;

   iv = 0 ;
   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      PEL_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ;
      string const& name = se->string_data( "name" ) ;
      string const& type = se->string_data( "type" ) ;
      if( name.substr(0,2) == "BS" )
      {
         if( type != "Bool" ) 
            FE_Launcher_ERROR::n3( name, "Bool" ) ;
         TYPE[iv]= PEL_Data::Bool ;
         boolVector const& values = se->boolVector_data( "values" ) ;
         if( values.size() != NB_CALC ) FE_Launcher_ERROR::n1( se, NB_CALC ) ; 
         BOOL_VALUES[iv]= values ;
         PEL_Bool* data = PEL_Bool::create( 0, false ) ;
         BOOL_DATA[iv]= data ;
         PEL_Exec::add_variable_to_execution_context( 
                            PEL_Variable::object( name ), data ) ;
      }
      else if( name.substr(0,2) == "IS" )
      {
         if( type != "Int" ) 
            FE_Launcher_ERROR::n3( name, "Int" ) ;
         TYPE[iv]= PEL_Data::Int ;
         intVector const& values = se->intVector_data( "values" ) ;
         if( values.size() != NB_CALC ) FE_Launcher_ERROR::n1( se, NB_CALC ) ; 
         INT_VALUES[iv]= values ;
         PEL_Int* data = PEL_Int::create( 0, 0 ) ;
         INT_DATA[iv]= data ;
         PEL_Exec::add_variable_to_execution_context( 
                            PEL_Variable::object( name ), data ) ;
      }
      else if( name.substr(0,2) == "DS" )
      {
         if( type != "Double" ) 
            FE_Launcher_ERROR::n3( name, "Double" ) ;
         TYPE[iv]= PEL_Data::Double ;
         doubleVector const& values = se->doubleVector_data( "values" ) ;
         if( values.size() != NB_CALC ) FE_Launcher_ERROR::n1( se, NB_CALC ) ; 
         DBL_VALUES[iv]= values ;
         PEL_Double* data = PEL_Double::create( 0, 0.0 ) ;
         DBL_DATA[iv]= data ;
         PEL_Exec::add_variable_to_execution_context( 
                            PEL_Variable::object( name ), data ) ;
      }
      else if( name.substr(0,2) == "SS" )
      {
         if( type != "String" ) 
            FE_Launcher_ERROR::n3( name, "String" ) ;
         TYPE[iv]= PEL_Data::String ;
         stringVector const& values = se->stringVector_data( "values" ) ;
         if( values.size() != NB_CALC ) FE_Launcher_ERROR::n1( se, NB_CALC ) ; 
         STRING_VALUES[iv]= values ;
         PEL_String* data = PEL_String::create( 0, "" ) ;
         STRING_DATA[iv]= data ;
         PEL_Exec::add_variable_to_execution_context( 
                            PEL_Variable::object( name ), data ) ;
      }
      else FE_Launcher_ERROR::n2() ;

      se->destroy() ; se=0 ;
      ++iv ;
   }
   ee->destroy() ;

}

//----------------------------------------------------------------------------
FE_Launcher:: ~FE_Launcher( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
FE_Launcher:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Launcher:: run" ) ;
   
   for( size_t i=0 ; i<NB_CALC ; ++i )
   {
      for( size_t iv=0 ; iv<TYPE.size() ; ++iv )
      {
         if( TYPE[iv] == PEL_Data::Bool )
         {
            BOOL_DATA[iv]->set( BOOL_VALUES[iv](i) ) ;
         } 
         else if( TYPE[iv] == PEL_Data::Int )
         {
            INT_DATA[iv]->set( INT_VALUES[iv](i) ) ;
         } 
         else if( TYPE[iv] == PEL_Data::Double )
         {
            DBL_DATA[iv]->set( DBL_VALUES[iv](i) ) ;
         }
         else if( TYPE[iv] == PEL_Data::String )
         {
            STRING_DATA[iv]->set( STRING_VALUES[iv](i) ) ;
         }
      }

      PEL_Module* m = PEL_Module::create( PEL_Root::object(),
                                          "TOP",
                                          DATA_FILE,
                                          PEL_Exec::execution_context() ) ;
      
      PEL_Module* m_appli = m->module( "PEL_Application" ) ;

      // Particular case of the saving files:
      PEL_Module* m_saver = 
                  m_appli->module( "PDE_DomainAndFields/PDE_ResultSaver" ) ;
      std::string ini_nn = 
                   m_saver->data_of_entry( "files_basename" )->to_string() ;
      ostringstream nn ;
      nn << ini_nn << "_" << i ;
      m_saver->replace_data_of_entry(
               "files_basename", PEL_String::create( m_saver, nn.str() ) ) ;
      
      if( !CONTROLER_FILE.empty() && PATTERN != PEL_ModuleExplorer::build )
      {
         std::ostringstream tmp ;
         tmp << "temp_file.pel" ;
         if( PEL_Exec::communicator()->nb_ranks()>1 )
         {
            tmp << "." << PEL_Exec::communicator()->rank() ;
         }
         std::string const tmp_file = tmp.str() ;
         std::ofstream out( tmp_file.c_str() ) ;
         if( !out ) FE_Launcher_ERROR:: n0( tmp_file ) ;
         out.close() ;    
         m_appli->write( tmp_file, "text" ) ;
         stringVector args(0) ;
         args.append( "-A" ) ;
         args.append( "check" ) ;
         args.append( CONTROLER_FILE ) ;
         args.append( tmp_file ) ;
         PEL_Application* controler = PEL_Application::make( 0, args ) ;
         controler->run() ;
         controler->destroy() ; controler = 0 ;
         PEL_System::erase( tmp_file ) ;
         if( PEL_Exec:: exit_code() != 0 )
         {
            PEL_Error::object()->raise_plain( "Launcher generated data deck checking failed" ) ;
         }
         PEL::out() << std::endl ;
      }
      PEL_ModuleExplorer const* ee =
                           PEL_ModuleExplorer::create( 0, m_appli, PATTERN ) ;
      PEL_Application* aa = PEL_Application::make( 0, ee ) ;
      aa->run() ;
      aa->destroy() ;
      ee->destroy() ;
   }
}

//internal-------------------------------------------------------------------
void
FE_Launcher_ERROR:: n0( std::string const& filename )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_Launcher :" << endl ;
   mesg << "   unable to open file : \"" << filename << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_Launcher_ERROR:: n1( PEL_ModuleExplorer const* exp, size_t nbv )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MODULE " << exp->name() << endl ;
   mesg << "   the data of keyword \"values\" should be" << endl ;
   mesg << "   a vector of length " << nbv ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_Launcher_ERROR:: n2( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_Launcher :" << endl ;
   mesg << "   the only recognized variables are of type" << endl ;
   mesg << "      - PEL_Data::Bool   (name starting with \"BS\")" << endl ;
   mesg << "      - PEL_Data::Int    (name starting with \"IS\")" << endl ;
   mesg << "      - PEL_Data::Double (name starting with \"DS\")" << endl ;
   mesg << "      - PEL_Data::String (name starting with \"SS\")" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_Launcher_ERROR:: n3( std::string const& name, std::string const& type )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_Launcher :" << endl ;
   mesg << "   a variable of name : \"" << name << "\"" << endl ;
   mesg << "   should be of type  : " << type ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
