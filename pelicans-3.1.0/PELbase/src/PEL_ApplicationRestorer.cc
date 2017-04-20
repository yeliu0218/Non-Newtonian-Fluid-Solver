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

#include <PEL_ApplicationRestorer.hh>

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

using std::string ;

PEL_ApplicationRestorer const* 
PEL_ApplicationRestorer:: PROTOTYPE = new PEL_ApplicationRestorer() ;

//-------------------------------------------------------------------------
PEL_ApplicationRestorer:: PEL_ApplicationRestorer( void )
//-------------------------------------------------------------------------
   : PEL_Application( "PEL_ApplicationRestorer" )
   , READER( 0 )
   , APPLI( 0 )
{
}
   
//---------------------------------------------------------------------------
PEL_ApplicationRestorer*
PEL_ApplicationRestorer:: create_replica( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ApplicationRestorer:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   PEL_ApplicationRestorer* result = new PEL_ApplicationRestorer( a_owner,
                                                                  exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ApplicationRestorer:: PEL_ApplicationRestorer( 
                                             PEL_Object* a_owner,
                                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , READER( 0 )
   , APPLI(0)
{
   PEL_ModuleExplorer* ee = 0 ;

   ee = exp->create_subexplorer( 0, "PEL_ObjectReader" ) ;
   READER = PEL_ObjectReader::create( this, ee ) ;
   ee->destroy() ;

   string appendum_file ;
   if( exp->has_entry( "appendum_file_name" ) )
   {
      appendum_file = exp->string_data( "appendum_file_name" ) ;
   }

   size_t cycle_number = 0 ;
   if( exp->has_entry( "cycle_number" ) )
   {
      cycle_number = exp->int_data( "cycle_number" ) ;
   }
      
   PEL_Module* mod = create_modified_data_deck_module( appendum_file ) ;
   ee = PEL_ModuleExplorer::create( 0, mod ) ;
   APPLI = PEL_Application::make( this, ee ) ;
   ee->destroy() ;

   READER->seek_cycle( cycle_number ) ;
   APPLI->restore_registered_objects( READER ) ;
   READER->close_cycle() ;

   PEL_CHECK_POST( APPLI!=0 ) ;
}

//---------------------------------------------------------------------------
PEL_ApplicationRestorer:: ~PEL_ApplicationRestorer( void )
//---------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_ApplicationRestorer:: run( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ApplicationRestorer:: run" ) ;

   APPLI->run() ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ApplicationRestorer:: create_modified_data_deck_module( 
                                            std::string const& appendum_file )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ApplicationRestorer:: create_modified_data_deck_module" ) ;

   PEL_Module* m = 0 ;

   PEL_Module* header = READER->header_module() ;
   if( header->has_module( "PEL_Application" ) )
   {
      m = header->module( "PEL_Application" ) ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "invalid restart file" ) ; 
   }
   PEL_ASSERT( m != 0 ) ;
   PEL_Module* result = m->create_clone( this ) ;

   result->remove_module( "PEL_ObjectWriter" ) ;
   change_owner( PEL_Root::object(), result ) ;

   if( !appendum_file.empty() )
   {
      PEL_Module* appendum =
         PEL_Module::create( 0, "ROOT", appendum_file,
                             PEL_Exec::execution_context() ) ;
      PEL_Module* o = appendum->module( "PEL_Application" ) ;
      result->merge_module( o ) ;
      appendum->destroy() ; appendum = 0 ;
   }
   
   PEL_CHECK( result != 0 ) ;
   PEL_CHECK( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}
