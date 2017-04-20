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

#include <PEL_Application.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>
#include <PEL_assertions.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

struct PEL_Application_ERROR
{
   static void n0( stringVector const& args ) ;
   static void n1( std::string const& class_name ) ;
} ;

//----------------------------------------------------------------------
PEL_Application*
PEL_Application:: make( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   PEL_Application const* proto =
      static_cast<PEL_Application const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PEL_Application* result = proto->create_replica( a_owner, exp ) ;
   result->initialize_objects_storage( exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Application*
PEL_Application:: make( PEL_Object* a_owner,
                        stringVector& args )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: make" ) ;
   PEL_CHECK_PRE( args.size() != 0 ) ;

   if( args( 0 ) != "-A" ) PEL_Application_ERROR::n0( args ) ;
   args.remove_at( 0 ) ;
   string name = args( 0 ) ;
   args.remove_at( 0 ) ;
   PEL_Application const* proto =
      static_cast<PEL_Application const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PEL_Application* result = proto->create_replica_from_args( a_owner, args ) ;
   if( args.size() != 0 ) result->print_usage_then_exit() ;
   
   result->initialize_objects_storage( 0 ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Application:: PEL_Application( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , SAVER( 0 )
   , persistent_objects( 0 )
   , persistent_objects_it( 0 )
{
   persistent_objects = PEL_ListIdentity::create( this ) ;
   persistent_objects_it = PEL_ListIterator::create( this,
                                                     persistent_objects ) ;

   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_Application:: PEL_Application( std::string const& name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , SAVER( 0 )
   , persistent_objects( 0 )
   , persistent_objects_it( 0 )
{
   PEL_LABEL( "PEL_Application:: PEL_Application" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_under_ownership_of( plugins_map() ) ) ;
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_Application:: ~PEL_Application( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Application:: register_storable_objects( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: register_storable_objects" ) ;

   add_storable_objects( persistent_objects ) ;
}


//----------------------------------------------------------------------
void
PEL_Application:: add_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: add_storable_objects" ) ;
}


//----------------------------------------------------------------------
void
PEL_Application:: write_storable_objects( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: write_storable_objects" ) ;

   if( SAVER != 0 )
   {
      SAVER->start_cycle() ;

//?????????????????? si une certaine verbosite
      PEL::out() << "\n---Saving (cycle:"
                 << SAVER->cycle_number()
                 << ")---\n" ;

      SAVER->start_new_object( "PEL_Application" ) ;

      for( persistent_objects_it->start();
	   persistent_objects_it->is_valid() ;
	   persistent_objects_it->go_next() )
      {
	 persistent_objects_it->item()->save_state( SAVER ) ;
      }


      SAVER->finalize_object() ;

      SAVER->terminate_cycle() ;
   }
}

//----------------------------------------------------------------------
void
PEL_Application:: restore_registered_objects( PEL_ObjectReader* ret ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: restore_registered_objects" ) ;

   ret->start_object_retrieval( "PEL_Application" ) ;

   for( persistent_objects_it->start();
        persistent_objects_it->is_valid() ;
        persistent_objects_it->go_next() )
   {
      persistent_objects_it->item()->restore_state( ret ) ;
   }

   ret->end_object_retrieval() ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Application:: object_writer_module_name( void ) const
//----------------------------------------------------------------------
{
   static std::string const result = "PEL_ObjectWriter" ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Application:: create_replica_PRE( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Application:: create_replica_POST( PEL_Application const* result,
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Application:: create_replica_from_args_POST( PEL_Application const* result,
                                                 PEL_Object* a_owner,
                                                 stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Application:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_Application*
PEL_Application:: create_replica_from_args( PEL_Object* a_owner,
                                            stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Application:: create_replica_from_args" ) ;
   
   PEL_Application* result = 0 ;

   PEL_Application_ERROR::n1( type_name() ) ;

   PEL_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Application:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
void
PEL_Application:: notify_error_in_arguments( void ) const
//----------------------------------------------------------------------
{
   print_usage_then_exit( -1 ) ;
}


//----------------------------------------------------------------------
void
PEL_Application:: print_usage( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Application:: print_options( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Application:: print_operands( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Application:: print_exit_status( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string
PEL_Application:: usage_title( std::string const& name ) const
//----------------------------------------------------------------------
{
   std::string result = "USAGE\n" ;
   result += "     <exe> -A " + name ;
   result += " [-h] " ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Application:: options_title( void ) const
//----------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "OPTIONS" << endl ;
   mesg << "     -h   Display help and exit" << endl << endl ;
   return mesg.str() ;
}

//----------------------------------------------------------------------
std::string
PEL_Application:: operands_title( void ) const
//----------------------------------------------------------------------
{
   std::string result = "OPERANDS\n" ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Application:: exit_status_title( void ) const
//----------------------------------------------------------------------
{
   std::string result = "EXIT STATUS\n" ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Application:: print_usage_then_exit( int exit_status ) const
//----------------------------------------------------------------------
{
   print_usage() ;
   PEL::out() << endl <<"OPTIONS" << endl ;
   PEL::out() << "     -h   Display help and exit" << endl << endl ;
   print_options() ;
   print_operands() ;
   print_exit_status() ;
   PEL::out() << endl ;
   PEL_System::exit( exit_status ) ;
}

//----------------------------------------------------------------------
void
PEL_Application:: initialize_objects_storage( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   if( ( exp != 0 ) && exp->has_module( object_writer_module_name() ) )
   {
      PEL_ModuleExplorer const* e =
         exp->create_subexplorer( 0, object_writer_module_name() ) ;
      SAVER = PEL_ObjectWriter::create( this, e, exp ) ;   
      e->destroy() ;
   }
   register_storable_objects() ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_Application:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
            PEL_ObjectRegister::create( PEL_Root::object(),
                                        "PEL_Application descendant" ) ;
   return( result ) ;
}

//internal---------------------------------------------------------------
void
PEL_Application_ERROR:: n0( stringVector const& args )
//internal---------------------------------------------------------------
{
   std::ostringstream os ;
   os << "invalid command-line arguments : " << endl << "  ";
   for( size_t i=0 ; i<args.size() ; ++i ) os << " " << args(i) ;
   PEL_Error::object()->raise_plain( os.str() ) ;
}

//internal---------------------------------------------------------------
void
PEL_Application_ERROR:: n1( std::string const& class_name )
//internal---------------------------------------------------------------
{
   std::ostringstream os ;
   os << "Class " << class_name << " derived from PEL_Application :" << endl ;
   os << "   instantiation from the command_line is impossible" << endl ;
   os << "   (\"create_replica_from_args\" has not been overridden)." ;
   PEL_Error::object()->raise_plain( os.str() ) ;
}
