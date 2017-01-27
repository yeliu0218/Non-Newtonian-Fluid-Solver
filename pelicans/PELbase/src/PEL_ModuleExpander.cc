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

#include <PEL_ModuleExpander.hh>

#include <PEL.hh>
#include <PEL_Context.hh>
#include <PEL_Exec.hh>
#include <PEL_ExtractionExp.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>

#include <fstream>
#include <iostream>


PEL_ModuleExpander const*
PEL_ModuleExpander::PROTOTYPE = new PEL_ModuleExpander() ;

//----------------------------------------------------------------------
PEL_ModuleExpander:: PEL_ModuleExpander( void )
//----------------------------------------------------------------------
   : PEL_Application( "PEL_ModuleExpander" )
   , BASE( "" )
   , OUTPUT( "" )
   , INPUT( "" )
{
}

//----------------------------------------------------------------------
PEL_ModuleExpander* 
PEL_ModuleExpander:: create_replica(
              PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExpander:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   std::string submod = ( exp->has_entry( "submodule_to_expand" ) ?
                          exp->string_data( "submodule_to_expand" ) :
                          "" ) ;
   
   PEL_ModuleExpander* result =
      new PEL_ModuleExpander( a_owner,
                              exp->string_data( "skeleton_file" ),
                              exp->string_data( "input_file" ),
                              exp->string_data( "expanded_file" ),
                              submod ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExpander* 
PEL_ModuleExpander:: create_replica_from_args(
                         PEL_Object* a_owner, stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExpander:: create_replica_from_args" ) ;
   
   if( args.size() != 3 ) notify_error_in_arguments() ;
   
   PEL_ModuleExpander* result =
      new PEL_ModuleExpander( a_owner, args(0), args(1), args(2), "" ) ;
   args.remove_at(0) ;
   args.remove_at(0) ;
   args.remove_at(0) ;

   PEL_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ModuleExpander:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExpander:: run" ) ;

   // Read input file:
   PEL_Module* mod = PEL_Module::create( 0, "Root", INPUT,
                                         PEL_Exec::execution_context() ) ;
   PEL_ModuleIterator* it = mod->create_module_iterator( mod ) ;
   it->start() ;
   PEL_Module* input_mod = it->item()->create_clone( this ) ;
   if( !SUBMOD.empty() )
   {
      if( !input_mod->has_module( SUBMOD ) )
      {
         PEL_Error::object()->raise_plain(
            "Submodule " + SUBMOD + " is not a submodule name in " + INPUT ) ;
      }
      input_mod = input_mod->module( SUBMOD ) ;
   }
   
   PEL_Module const* expanded_module =
                           create_expanded_module( mod, input_mod, BASE ) ;

   std::ofstream out( OUTPUT.c_str() ) ;
   out.close() ;
   expanded_module->write( OUTPUT, "text" ) ;

   mod->destroy() ; mod = 0 ; it = 0 ; expanded_module = 0 ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExpander:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << usage_title( "pelsdd" )  ;
   PEL::out() << "<base_dir> <input_file> <expanded_file>"
              << std::endl << std::endl ;
   PEL::out() << "     Create the full data deck associated to a simplified"
              << std::endl
              << "     one and a data base."
              << std::endl ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExpander:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << operands_title() ;
   PEL::out() << "     <skeleton_file>" << std::endl
              << "          skeleton data deck"
              << std::endl << std::endl ;
   PEL::out() << "     <input_file>" << std::endl
              << "          simplified data deck"
              << std::endl << std::endl ;
   PEL::out() << "     <expanded_file>" << std::endl
              << "          full data deck created"
              << std::endl << std::endl ;
}

//----------------------------------------------------------------------
PEL_ModuleExpander:: PEL_ModuleExpander( PEL_Object* a_owner,
                                         std::string const& skeleton_file,
                                         std::string const& input_file,
                                         std::string const& expanded_file,
                                         std::string const& submod_name  )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , BASE( skeleton_file )
   , OUTPUT( expanded_file )
   , INPUT( input_file )
   , SUBMOD( submod_name )
{
   PEL_LABEL( "PEL_ModuleExpander:: PEL_ModuleExpander" ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExpander:: ~PEL_ModuleExpander( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
PEL_Module const*
PEL_ModuleExpander:: create_expanded_module(
                                 PEL_Object* a_owner,
                                 PEL_Module const* input_mod,
                                 std::string const& skeleton_file_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExpander:: create_expanded_module" ) ;
   PEL_CHECK_PRE( input_mod != 0 ) ;
   PEL_CHECK_PRE( ! skeleton_file_name.empty() ) ;

   // Initialize data-base:
   PEL_ExtractionExp::initialize( input_mod ) ;

   // Read skeleton file:
   PEL_Module* result = create_skeleton_module( a_owner, skeleton_file_name ) ;

   PEL_ExtractionExp::reset() ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_ModuleExpander:: create_skeleton_module( PEL_Object* a_owner,
                                             std::string const& file_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExpander:: create_skeleton_module" ) ;
   PEL_CHECK( !file_name.empty() ) ;
   PEL_CHECK( PEL_ExtractionExp::is_initialized() ) ;
   
   PEL_Module* mod = PEL_Module::create( a_owner, "Root", file_name,
                                         PEL_Exec::execution_context() ) ;
   PEL_ModuleIterator* it = mod->create_module_iterator( 0 ) ;
   it->start() ;
   PEL_ASSERT( it->is_valid() ) ;
   PEL_Module* result = it->item() ;
   it->go_next() ;
   PEL_ASSERT( !it->is_valid() ) ;
   it->destroy() ; it = 0 ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( a_owner ) ) ;
   return( result ) ;
}

