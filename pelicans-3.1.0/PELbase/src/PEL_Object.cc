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

#include <PEL_Object.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_Module.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>

#include <typeinfo>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <map>

using std::map ;
using std::string ;

static map<PEL_Object const*,size_t>* _PEL_OBJECT_LIST = 0 ;
static size_t _RANK = 0 ;
static PEL_Object const* catched_object = 0 ;
static size_t catched_object_rank = (size_t) ~0 ;

static map<PEL_Object const*,size_t>* _PEL_OBJECT_LIST2 = 0 ;
static bool trace_allocation = false ;

//----------------------------------------------------------------------
size_t PEL_Object:: ALLOCATED = 0 ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
PEL_Object:: PEL_Object( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : MY_OWNER( 0 ),
     POSSESSIONS( 0 )   
{
   
   // Construct the object
   ALLOCATED++ ;
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      static bool prem = true ;
      if( prem )
      {
         _PEL_OBJECT_LIST2 = new map<PEL_Object const*,size_t> ;
         _PEL_OBJECT_LIST = new map<PEL_Object const*,size_t> ;
         _RANK = 0 ;
         prem = false ;
      }
      if( a_owner==0 )
      {
         _PEL_OBJECT_LIST->operator[]( this ) = _RANK ;
      }
      if( catched_object==this )
      {
         PEL_Error::object()->trace( "Catched object creation" ) ;
      }
      if( _RANK==catched_object_rank )
      {
         PEL_Error::object()->trace( "Catched object creation" ) ;
         catched_object=this ;
      }
      if( trace_allocation )
      {
         _PEL_OBJECT_LIST2->operator[]( this ) = _RANK ;
      }
      _RANK++ ;
   }
   if( a_owner != 0 )
   {
      a_owner->insert_possession( this ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Object:: ~PEL_Object( void )
//----------------------------------------------------------------------
{
   if( MY_OWNER!=0 &&  owner()->POSSESSIONS!=0 )
   {
      MY_OWNER->POSSESSIONS->remove( this ) ;
      MY_OWNER=0 ;
   }
   
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      if( owner() == 0 )
      {
         _PEL_OBJECT_LIST->erase( this ) ;
      }
      if( catched_object==this )
      {
         PEL_Error::object()->trace( "Catched object deletion" ) ;
      }
      if( trace_allocation )
      {
         _PEL_OBJECT_LIST2->erase( this ) ;
      }
   }
   if( POSSESSIONS != 0 )
   {
      PEL_ListIterator *it = POSSESSIONS->create_iterator( 0 ) ;
      for( ; it->is_valid() ; it->go_next() )
      {
         PEL_Object* child = it->item() ;
         child->MY_OWNER = 0 ;
         delete child ;
      }
      it->destroy() ;
      POSSESSIONS->destroy() ; POSSESSIONS=0 ;
   }
   ALLOCATED-- ;
}

//----------------------------------------------------------------------
PEL_Object*
PEL_Object:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: create_clone" ) ;
   PEL_Error::object()->raise_not_implemented( this, "create_clone" ) ;

   PEL_Object* result = 0 ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: update( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: update" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Error::object()->raise_not_implemented( this, "update" ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: destroy( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: destroy" ) ;
   PEL_CHECK_PRE( owner()==0 ) ;
   
   delete this ;
}

//----------------------------------------------------------------------
void
PEL_Object:: destroy_possession( PEL_Object const* a_possession )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: destroy_possession" ) ;
   PEL_CHECK_PRE( a_possession->owner()==this ) ;

   delete a_possession ;
}

//----------------------------------------------------------------------
size_t
PEL_Object:: address( void ) const
//----------------------------------------------------------------------
{
   return( (size_t)(long) this  ) ; //????????????????????????????????
}

bool valid_first( char c )
{
   return isalpha(c)!=0 ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Object:: type_name( void ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: type_name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static string result ;

   string nn = typeid(*this).name() ;

   // on supprime quelques chiffres qui sont accolés au début
   // du nom de la classe par le compilateur gcc
   string::iterator it = std::find_if( nn.begin(), nn.end(), valid_first ) ;
   if( it == nn.end() )
   {
      result = nn ;
   }
   else
   {
      result.assign( it, nn.end() ) ;
   }
   // on supprime egalement le "class" ajoute par le compilateur VC++
   size_t idx ;
   if( ( idx = result.find( "class " ) ) < result.length() )
   {
      result = result.substr( idx+6 ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Object const*
PEL_Object:: owner( void ) const
//----------------------------------------------------------------------
{
   return( MY_OWNER ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: is_under_ownership_of( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: is_under_ownership_of" ) ;
   PEL_CHECK_PRE( other!=this ) ;
   PEL_CHECK_INV( invariant() ) ;
   bool result = false ;
   if( MY_OWNER==0 )
   {
      result = ( other==0 ) ;
   }
   else if( MY_OWNER==other )
   {
      result = true ;
   }
   else
   {
      result = MY_OWNER->is_under_ownership_of( other ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: same_type( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: same_type" ) ;
   PEL_CHECK_PRE( other != 0 ) ;

   bool result = true ;
// ( typeid(*this) == typeid(*other) ) ;

   PEL_CHECK_POST( 
       FORMAL( EQUIVALENT( result, typeid(*this) == typeid(*other ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: set_owner( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: set_owner" ) ;
   PEL_CHECK_PRE( a_owner != 0 ) ;
   PEL_CHECK_PRE( owner() == 0 ) ;

   a_owner->insert_possession( this ) ;

   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      _PEL_OBJECT_LIST->erase( this ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: change_owner( PEL_Object* a_owner,
                           PEL_Object* a_possession )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: change_owner" ) ;
   PEL_CHECK_PRE( a_possession->is_under_ownership_of(this) ) ;

   PEL_Object* father = a_possession->MY_OWNER ;
   
   father->POSSESSIONS->remove( a_possession ) ;
   a_possession->MY_OWNER = 0 ;
   if( a_owner!=0 )
   {
      a_possession->set_owner( a_owner ) ;
   }

   PEL_CHECK_POST( a_possession->owner()==a_owner ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: comparable( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   return( same_type( other ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: is_equal( PEL_Object const* other ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = has_same_address( other ) ;

   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
PEL_Object:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "three_way_comparison" ) ;
   int result = -1 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Object:: hash_code( void ) const 
//----------------------------------------------------------------------
{
   size_t result = address() ;
   return( result  ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: has_same_address( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   return( other == this ) ;
}

//----------------------------------------------------------------------
int
PEL_Object:: address_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   int result = 0 ;
   if( this > other )
   {
      result = 1 ;
   }
   else if( this < other )
   {
      result = -1 ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Error::object()->raise_not_implemented( this, "save_state" ) ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Error::object()->raise_not_implemented( this, "restore_state" ) ;
   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: update_for_restore_state( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: update_for_restore_state" ) ;
   PEL_CHECK_INV( invariant() ) ;
   update() ;
}

//----------------------------------------------------------------------
void
PEL_Object:: display_info( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: display_info" ) ;
   std::string const s( indent_width, ' ' ) ;
   os << s << "type    : " << type_name() << std::endl ;
   if( owner()==0 )
   {
      os << s << "owner   : nil" << std::endl ;
   }
   else
   {
      os << s << "owner   : " << std::endl ;
      os << s << "   type    : " << owner()->type_name() << std::endl ;
      os << s << "   address : "
         << std::setiosflags( std::ios::hex )
         << std::setiosflags( std::ios::showbase )
         << owner()
         << std::resetiosflags( std::ios::showbase )
         << std::resetiosflags( std::ios::hex )
         << std::endl ;
   }
   if( POSSESSIONS!=0 && POSSESSIONS->count()!=0 )
   {
      os << s << "owned objects : " << std::endl ;
      for( size_t i=0 ; i<POSSESSIONS->count() ; i++ )
      {
         POSSESSIONS->at(i)->display_info(os,indent_width+3) ;
      }
   }
   os << s << "address : "
      << std::setiosflags( std::ios::hex )
      << std::setiosflags( std::ios::showbase )
      << this
      << std::resetiosflags( std::ios::showbase )
      << std::resetiosflags( std::ios::hex )
      << std::endl ;
}

//----------------------------------------------------------------------
void
PEL_Object:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
int
PEL_Object:: GetNumberOf_PEL_objects( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: GetNumberOf_PEL_objects" ) ;
   return( (int)ALLOCATED ) ;
}

//----------------------------------------------------------------------
std::ostream&  
PEL_Object:: TraceRemainingObjects( std::ostream& out  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Object:: TraceRemainingObjects" ) ;
   out << std::endl ;
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      for( map<PEL_Object const*,size_t>::const_iterator it=
              _PEL_OBJECT_LIST->begin() ;
           it!=_PEL_OBJECT_LIST->end() ;
           ++it ) 
      {
         out << "[" << (*it).second << "] :" << std::endl ;
         (*it).first->display_info( out, 3 ) ;
         (*it).first->print( out, 3 ) ;
         out << std::endl ;
      }
   }
   else
   {
      out << "To enable traceRemainingObjects facilities do run \n"
          << " application with -Cobjects command line option.\n" ;
   }
   return out ;
}

//----------------------------------------------------------------------
void
PEL_Object:: catch_object( PEL_Object const* obj )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( PEL_Assertion::is_handling_check( PEL_Assertion::Objects )) ;
   catched_object = obj ;
}

//----------------------------------------------------------------------
void
PEL_Object:: catch_object_by_rank( size_t rank )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( PEL_Assertion::is_handling_check( PEL_Assertion::Objects )) ;
   catched_object_rank = rank ;
}

//----------------------------------------------------------------------
void
PEL_Object:: start_trace_allocating( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( !trace_allocating() ) ;
   trace_allocation = true ;
   PEL_CHECK_POST( trace_allocating() ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: trace_allocating( void )
//----------------------------------------------------------------------
{
   return( trace_allocation ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: stop_trace_allocating( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( trace_allocating() ) ;
   trace_allocation = false ;
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      _PEL_OBJECT_LIST2->clear() ;
   }
   PEL_CHECK_POST( !trace_allocating() ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: trace_not_destroyed_object( std::ostream& out )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( trace_allocating() ) ;
   if( PEL_Assertion::is_handling_check( PEL_Assertion::Objects ) )
   {
      for( map<PEL_Object const*,size_t>::const_iterator it=
                       _PEL_OBJECT_LIST2->begin() ;
           it!=_PEL_OBJECT_LIST2->end() ;
           ++it )
      {
         out << "[" << it->second << "] :" << std::endl ;
         it->first->display_info( out, 3 ) ;
         it->first->print( out, 3 ) ;
         out << std::endl ;
      }
   }
}

//----------------------------------------------------------------------
bool
PEL_Object:: invariant( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: create_clone_POST( PEL_Object const* result,
                                PEL_Object const* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: save_state_PRE( PEL_ObjectWriter const* writer ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( writer != 0 ) ;
   PEL_ASSERT( writer->has_an_opened_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: save_state_POST( PEL_ObjectWriter const* writer ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( writer->has_an_opened_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: restore_state_PRE( PEL_ObjectReader const* reader ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( reader != 0 ) ;
   PEL_ASSERT( reader->positioned_in_a_valid_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: restore_state_POST( PEL_ObjectReader const* reader ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( reader->positioned_in_a_valid_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: is_equal_PRE( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( comparable( other ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: is_equal_POST( bool result, PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result, hash_code()==other->hash_code() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: three_way_comparison_PRE( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( comparable( other ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Object:: three_way_comparison_POST( int result, 
                                        PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( EQUIVALENT( result==0, is_equal(other) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
PEL_Object:: insert_possession( PEL_Object* obj )
//----------------------------------------------------------------------
{
   PEL_CHECK( obj!=0 ) ;

   if( POSSESSIONS == 0 )
   {
       POSSESSIONS = PEL_ListIdentity:: create( 0 ) ;  
   }
   POSSESSIONS->prepend( obj ) ;
   obj->MY_OWNER = this ;

   PEL_CHECK_POST( obj->has_same_address( POSSESSIONS->item( obj ) ) ) ;
}

