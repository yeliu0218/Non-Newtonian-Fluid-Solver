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

#include <PEL_Expression.hh>

#include <PEL_assertions.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_List.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>

#include <stringVector.hh>

#include <iostream>
#include <sstream>
#include <set>

//----------------------------------------------------------------------
PEL_Expression*
PEL_Expression:: create( PEL_Object* a_owner,
                         std::string const& a_name,
                         PEL_Sequence const* argument_list,
                         std::string const& a_comment ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( argument_list!=0 ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<argument_list->index_limit() ; ++i ),
         dynamic_cast<PEL_Data*>( argument_list->at(i) ) != 0 ) ) ;
   PEL_CHECK_PRE( valid_arguments_of( a_name, argument_list ) ) ;

   PEL_Expression const* proto =
      static_cast<PEL_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
   
   PEL_Expression* result = proto->create_replica( a_owner, argument_list ) ;
   result->COMMENT = a_comment ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->external_brackets_are_set() ) ;
   PEL_CHECK_POST( result->comment() == a_comment ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_Expression:: registered_expressions( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: registered_expressions" ) ;

   return( plugins_name() ) ;
}

//----------------------------------------------------------------------
PEL_Expression:: PEL_Expression( PEL_Object* a_owner,
                                 std::string const& a_name,
                                 PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Data( a_owner )
   , NAME( a_name )
   , HAS_BRACKETS( true )
   , COMMENT( "" )
   , ARGUMENTS( argument_list )
   , IT( argument_list->create_iterator( this ) )
{
   PEL_LABEL( "PEL_Expression:: PEL_Expression" ) ;
   PEL_CHECK_PRE( argument_list != 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_Expression:: PEL_Expression( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Data( plugins_map() )
   , NAME( a_name )
   , HAS_BRACKETS( true )
   , COMMENT( "" )
   , ARGUMENTS( 0 )
   , IT( 0 )
{
   PEL_LABEL( "PEL_Expression:: PEL_Expression" ) ;
   
   plugins_map()->register_item( a_name, this ) ;
   plugins_name().append( a_name ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_Expression*
PEL_Expression:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: create_clone" ) ;
   PEL_Expression* result = 0 ;
   if( is_a_prototype() )
   {
      result = const_cast<PEL_Expression*>(this) ;
   }
   else
   {
      PEL_Sequence const* args = ARGUMENTS ;
      PEL_Sequence* l = args->create_clone( 0 ) ;
      l->clear() ;
      for( size_t i=0 ; i<args->index_limit() ; ++i )
      {
         PEL_Object* obj = args->at(i) ;
         l->append( obj->create_clone( l ) ) ;
      }
      result = PEL_Expression::create( a_owner, name(), l, comment() ) ;
      l->set_owner( result ) ;
      if( HAS_BRACKETS )
         result->set_external_brackets() ;
      else
         result->unset_external_brackets() ;
   }
   PEL_CHECK_POST( IMPLIES( is_a_prototype(), result==this ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_a_prototype(),
                            create_clone_POST( result, a_owner ) ) ) ;
   PEL_CHECK_POST( result->name()==name() ) ;
   PEL_CHECK_POST( result->comment()==comment() ) ;
   PEL_CHECK_POST( result->external_brackets_are_set()==external_brackets_are_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Expression:: ~PEL_Expression( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: ~PEL_Expression" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Expression:: name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( NAME ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Expression:: usage_of( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: usage_of" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   
   PEL_Expression const* proto =
      static_cast<PEL_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   return( proto->usage() ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: valid_arguments_of( std::string const& a_name,
                                     PEL_Sequence const* argument_list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: valid_arguments_of" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( argument_list!=0 ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<argument_list->index_limit() ; ++i ),
         dynamic_cast<PEL_Data*>( argument_list->at(i) ) != 0 ) ) ;
   
   PEL_Expression const* proto =
      static_cast<PEL_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   return( proto->valid_arguments( argument_list ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: declare( PEL_List * lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
      data->declare( lst ) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: context_has_required_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: context_has_required_variables" ) ;
   PEL_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = true ;
   for( IT->start() ; result && IT->is_valid() ; IT->go_next() )
   {
      PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
      result = data->context_has_required_variables(ct) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }

   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_arguments() ; ++i ),
              !result || arg(i)->context_has_required_variables(ct) ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Expression:: nb_arguments( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: nb_arguments" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   return( ARGUMENTS->count() ) ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_Expression:: arg( size_t idx ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: arguments" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK( idx<nb_arguments() ) ;

   PEL_Data const* result = extract_arg( ARGUMENTS, idx) ;

   PEL_CHECK_POST( result != 0 ) ;   
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_Expression:: extract_arg( PEL_Sequence const* some_arguments,
                              size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: extract_arg" ) ;
   PEL_CHECK( idx<some_arguments->count() ) ;
   PEL_CHECK_PRE( dynamic_cast<PEL_Data*>( some_arguments->at(idx) ) != 0 ) ;

   PEL_Data const* result = static_cast<PEL_Data*>( some_arguments->at(idx) ) ;

   PEL_CHECK_POST( result != 0 ) ;   
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: raise_error( std::string const& message ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: raise_error" ) ;
   PEL_CHECK_PRE( !message.empty() ) ;

   std::ostringstream msg ;
   msg << "*** Expression " << name() << " error:" << std::endl ;
   msg << "    usage: " << usage() << std::endl ;
   msg << "    " << message << std::endl ;
   if( !COMMENT.empty() ) msg << "    " << COMMENT << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: value_can_be_evaluated " ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
   bool result = true ;
   for( IT->start() ; IT->is_valid() && result ; IT->go_next() )
   {
      PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
      result = data->value_can_be_evaluated( ct ) ;
      if( !IT->is_valid() )
      {
         result = false ; // Cycle
         break ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_Expression:: undefined_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: undefined_variables" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
   
   static stringVector result(0) ;
   stringVector undef_var(0) ;
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
      stringVector const& s = data->undefined_variables( ct ) ;
      for( size_t i=0 ; i<s.size() ; ++i )
      {
         undef_var.extend( s(i) ) ;
      }
      if( !IT->is_valid() )
      {
         break ; // Cycle
      }
   }
   result = undef_var ;
   result.sort() ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Expression:: create_derivative( PEL_Object* a_owner,
                                    PEL_Variable const* var,
                                    PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Error::object()->raise_plain(
      "No create_derivative method implemented for expression "+name() ) ;
   PEL_Expression* result = 0 ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: external_brackets_are_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_BRACKETS ) ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: set_external_brackets( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: set_external_brackets" ) ;
   HAS_BRACKETS = true ;
   PEL_CHECK_POST( external_brackets_are_set() ) ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: unset_external_brackets( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: unset_external_brackets" ) ;
   HAS_BRACKETS = false ;
   PEL_CHECK_POST( !external_brackets_are_set() ) ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << name() << "(" ;
   if( ARGUMENTS!=0 )
   {
      bool prem = true ;
      for( IT->start() ; IT->is_valid() ; IT->go_next() )
      {
         if( !prem )
         {
            os << ", " ;
         }
         prem = false ;
         PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
         data->print( os, 0 ) ;
         if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
      }
   }
   else
   {
      os << "prototype" ;
   }
   os << ")" ;
}

//----------------------------------------------------------------------
void
PEL_Expression:: print_prototypes( std::ostream& os, size_t indent_width )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: print_prototypes" ) ;

   std::string space( indent_width, ' ' ) ;
   PEL_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   std::set< std::string > names ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_Expression const* expr = 
                     static_cast< PEL_Expression* >( it->item() ) ;
      names.insert( expr->name() ) ;
   }
   std::set< std::string >::const_iterator itn = names.begin() ;
   for( ; itn != names.end() ; ++itn ) 
   {
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         PEL_Expression const* expr = 
                     static_cast< PEL_Expression* >( it->item() ) ;
         if( expr->name() == (*itn) )
         {
            os << "     (\"" << expr->usage() << "\" \"\" \"\")" << std::endl ;
            break ;
         }
      }
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( ARGUMENTS==0 ) ;  
}

//----------------------------------------------------------------------
std::string const&
PEL_Expression:: comment( void ) const
//----------------------------------------------------------------------
{
   return( COMMENT ) ;  
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Expression:: create_non_const_simplification( PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: create_non_const_simplification" ) ;
   PEL_CHECK_PRE( create_non_const_simplification_PRE( a_owner ) ) ;

   PEL_List* new_list = PEL_List::create( 0 ) ;
   for( IT->start(); IT->is_valid() ; IT->go_next() )
   {
      PEL_Data const* data = static_cast<PEL_Data const*>( IT->item() ) ;
      new_list->append( data->create_simplification( new_list ) ) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }
   PEL_Expression* tmp = create_replica( 0, new_list ) ;
   new_list->set_owner( tmp ) ;
   
   PEL_Data* result =0 ;
   
   if( tmp->is_constant() )
   {
      result = tmp->create_simplification( a_owner ) ;
   }
   else
   {
      result = tmp->create_operator_simplification( a_owner ) ;
   }
   if( result!=tmp )
   {
      tmp->destroy() ;
   }
   else
   {
      tmp->set_owner( a_owner ) ;
   }
   
   PEL_CHECK_POST( create_non_const_simplification_POST( a_owner, result ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Expression:: create_operator_simplification( PEL_Object* a_owner )  
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: create_operator_simplification" ) ;
   PEL_Data* result = this ;
   PEL_CHECK_POST( create_operator_simplification_POST( a_owner, result ) ) ;
   return result ;
}


//----------------------------------------------------------------------
bool
PEL_Expression:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression:: is_raw_data" ) ;

   PEL_CHECK_POST( is_raw_data_POST( false ) ) ;
   
   return false ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: create_replica_PRE( PEL_Object const* a_owner,
                                     PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( argument_list != 0 ) ;
   PEL_ASSERT( argument_list->count() == argument_list->index_limit() ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<argument_list->count() ; ++i ),
              dynamic_cast<PEL_Data const*>( argument_list->at(i) ) != 0 ) ) ;
   PEL_ASSERT( valid_arguments( argument_list ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: create_replica_POST( PEL_Expression const* result,
                                      PEL_Object const* a_owner,
                                      PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->name() == name() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: valid_arguments_PRE( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( some_arguments != 0 ) ;
   PEL_ASSERT( some_arguments->count() == some_arguments->index_limit() ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<some_arguments->count() ; ++i ),
              dynamic_cast<PEL_Data const*>( some_arguments->at(i) ) != 0 ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: create_operator_simplification_POST(
               PEL_Object const* a_owner, PEL_Data const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result==this || result->owner()==a_owner ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Expression:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Data::invariant() ) ;
   return true ;  
}

//----------------------------------------------------------------------
stringVector&
PEL_Expression:: plugins_name( void )
//----------------------------------------------------------------------
{
   static stringVector result(0) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_Expression:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
             PEL_ObjectRegister::create( PEL_Root::object(),
                                         "PEL_Expression descendant" ) ;
   return( result ) ;
}
