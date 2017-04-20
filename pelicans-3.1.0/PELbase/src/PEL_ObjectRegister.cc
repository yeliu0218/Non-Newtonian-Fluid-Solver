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

#include <PEL_ObjectRegister.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Map.hh>
#include <PEL_MapIterator.hh>
#include <PEL_String.hh>

#include <stringVector.hh>

#include <sstream>

struct PEL_ObjectRegister_ERROR
{
   static void n0( std::string const& register_name,
                   std::string const& item_name ) ;
   static void n1( std::string const& register_name,
                   std::string const& item_name,
                   PEL_Map const* items ) ;
   static void n3( std::string const& register_name,
                   std::string const& item_name ) ;
   static void n4( std::string const& register_name ) ;
} ;

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_ObjectRegister:: create( PEL_Object* a_owner,
                             std::string const& a_register_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: create" ) ;
   PEL_CHECK_PRE( !a_register_name.empty() ) ;

   PEL_ObjectRegister* result =
                    new PEL_ObjectRegister( a_owner, a_register_name ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->register_name() == a_register_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister:: PEL_ObjectRegister( PEL_Object* a_owner,
                                         std::string const& a_register_name )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , REGISTER( PEL_Map::create( this ) )
   , NAME( a_register_name )
{
}

//----------------------------------------------------------------------
PEL_ObjectRegister:: ~PEL_ObjectRegister( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const&
PEL_ObjectRegister:: register_name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//----------------------------------------------------------------------
bool
PEL_ObjectRegister:: has( std::string const& a_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: has" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_String const* nn = PEL_String::create( 0, a_name ) ;
   bool result = REGISTER->has_key( nn ) ;
   nn->destroy() ; nn = 0 ;

   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Object*
PEL_ObjectRegister:: item( std::string const& a_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: item" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_String const* nn = PEL_String::create( 0, a_name ) ;
   PEL_Object* result = REGISTER->item_at( nn ) ;
   nn->destroy() ; nn = 0 ;
   if( result == 0 )
   {
      PEL_ObjectRegister_ERROR::n1( register_name(), a_name, REGISTER ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Iterator*
PEL_ObjectRegister:: create_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: create_iterator" ) ;

   PEL_Iterator* result = REGISTER->create_iterator( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ObjectRegister:: register_item( std::string const& a_name,
                                    PEL_Object* an_item )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: register_item" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( an_item != 0 ) ;

   PEL_String* nn = PEL_String::create( an_item, a_name ) ;
   if( REGISTER->has_key( nn ) )
   {
      PEL_ObjectRegister_ERROR::n0( register_name(), a_name ) ;
   }
   REGISTER->set_item_at( nn, an_item ) ;

   PEL_CHECK_POST( has( a_name ) ) ;
   PEL_CHECK_POST( item( a_name ) == an_item ) ;
}

//----------------------------------------------------------------------
void
PEL_ObjectRegister:: unregister_item( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: unregister_item" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   if( !has( a_name ) )
   {
      PEL_ObjectRegister_ERROR::n3( register_name(), a_name ) ;
   }
   
   PEL_String const* nn = PEL_String::create( 0, a_name ) ;
   REGISTER->remove_at( nn ) ;
   nn->destroy() ; nn = 0 ;

   PEL_CHECK_POST( !has( a_name ) ) ;
}

//----------------------------------------------------------------------
void
PEL_ObjectRegister:: unregister_item( PEL_Object* an_item )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectRegister:: unregister_item" ) ;
   PEL_CHECK_PRE( an_item != 0 ) ;

   PEL_Object const* nn = 0 ;
   PEL_MapIterator* it = REGISTER->create_iterator( 0 ) ;
   for( it->start() ; nn == 0 && it->is_valid() ; it->go_next() )
   {
      if( it->item() == an_item ) nn = it->key() ;
   }
   it->destroy() ; it = 0 ;
   
   if( nn == 0 )
   {
      PEL_ObjectRegister_ERROR::n4( register_name() ) ;
   }
   REGISTER->remove_at( nn ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectRegister_ERROR:: n0( std::string const& register_name,
                               std::string const& item_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to register : "   << register_name << std::endl ;
   msg << "           of  name : \"" << item_name << "\"" << std::endl ;
   msg << "which has already been registered" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectRegister_ERROR:: n1( std::string const& register_name,
                               std::string const& item_name ,
                               PEL_Map const* items )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Request for non registered : "   << register_name << std::endl ;
   msg << "                  of  name : \"" << item_name << "\""
       << std::endl << std::endl ;
   if( items->count() == 0 )
   {
      msg << "No " + register_name + " registered" << std::endl ;
   }
   else
   {  
      msg << items->count() << " registered " + register_name + "(s) of name : " 
          << std::endl ;
      stringVector ref_table( items->count() ) ;
      PEL_MapIterator* it = items->create_iterator( 0 ) ;
      for( size_t i=0 ; it->is_valid() ; it->go_next(), ++i )
      {
         PEL_String const* ss = dynamic_cast<PEL_String const*>( it->key() ) ;
         ref_table(i) = ss->to_string() ;
      }
      it->destroy() ; it = 0 ;
      ref_table.sort() ;
      for( size_t i=0 ; i<ref_table.size() ; ++i )
      {
         msg << "   - \"" << ref_table(i) << "\"" << std::endl ;
      }
   }
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectRegister_ERROR:: n3( std::string const& register_name,
                               std::string const& item_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to unregister : "   << register_name << std::endl ;
   msg << "             of  name : \"" << item_name << "\"" << std::endl ;
   msg << "which is not registered" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectRegister_ERROR:: n4( std::string const& register_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to unregister : "   << register_name << std::endl ;
   msg << "a no registered object" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
