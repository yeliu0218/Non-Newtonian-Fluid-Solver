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

#include <PEL_DataWithContextExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_String.hh>
#include <PEL_Variable.hh>
#include <stringVector.hh>

#include <iostream>

PEL_DataWithContextExp const*
PEL_DataWithContextExp::PROTO =
                     new PEL_DataWithContextExp( "data_with_context" ) ;

//----------------------------------------------------------------------
PEL_DataWithContextExp:: PEL_DataWithContextExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_TransferExp( a_name )
   , DATA( 0 )
{
   PEL_LABEL( "PEL_DataWithContextExp:: PEL_DataWithContextExp" ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContextExp:: PEL_DataWithContextExp(
                                    PEL_Object* a_owner,
                                    std::string const& a_name,
                                    PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_TransferExp( a_owner, a_name, argument_list )
   , DATA( 0 )
{
   PEL_LABEL( "PEL_DataWithContextExp:: PEL_DataWithContextExp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Data const* d = arg(0) ;
   PEL_ContextSimple* ct = PEL_ContextSimple::create( this ) ;
   for( size_t i=1 ; i<argument_list->count() ; )
   {
      std::string const& n = arg(i++)->to_string() ;
      PEL_Variable const* v = PEL_Variable::object( n ) ;
      ct->extend( v, arg(i++)->create_clone( ct ) ) ;
   }
   DATA = PEL_DataWithContext::create( this, d, ct ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContextExp:: ~PEL_DataWithContextExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: ~PEL_DataWithContextExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContextExp*
PEL_DataWithContextExp:: create(
                           PEL_Object* a_owner,
                           PEL_Data const* data, PEL_Context const* ct )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: create" ) ;
   PEL_CHECK_PRE( data != 0 ) ;

   PEL_List* argument_list = PEL_List::create( 0 ) ;
   argument_list->append( data->create_clone( argument_list ) ) ;
   if( ct != 0 )
   {
      for( size_t i=0 ; i<ct->nb_variables() ; ++i )
      {
         PEL_Variable const* v = ct->variable(i) ;
         argument_list->append( PEL_String::create( argument_list, v->name() ) ) ;
         argument_list->append( ct->value( v )->create_clone( argument_list ) ) ;
      }
   }
   
   PEL_DataWithContextExp* result =
                       PROTO->create_replica( a_owner, argument_list ) ;

   argument_list->set_owner( result ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContextExp*
PEL_DataWithContextExp:: create_replica(
                             PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_DataWithContextExp* result =
      new PEL_DataWithContextExp( a_owner, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_DataWithContextExp:: declare( PEL_List* lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;
   DATA->declare( lst ) ;
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContextExp:: context_has_required_variables(
                                           PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp::context_has_required_variables " ) ;
   PEL_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   return( DATA->context_has_required_variables( ct ) ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_DataWithContextExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: data_type" ) ;
   return( DATA->data_type() ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContextExp:: is_constant( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContextExp:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: value_can_be_evaluated" ) ;
   return( DATA->value_can_be_evaluated( ct ) ) ;
}

//----------------------------------------------------------------------
stringVector const& 
PEL_DataWithContextExp:: undefined_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: undefined_variables" ) ;
   return( DATA->undefined_variables( ct ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContextExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   size_t const n = some_arguments->count() ;
   bool result = ( n%2 == 1 ) ;
   for( size_t i=1 ; result && i<n ; i += 2 )
   {
      result = ( extract_arg( some_arguments, i )->data_type() == String ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_DataWithContextExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: usage" ) ;
   static std::string result =
      name() + "(<expression>[,SS,<value>])" ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_DataWithContextExp:: data( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContextExp:: data" ) ;
   PEL_CHECK( data_PRE( ct ) ) ;

   PEL_Data const* result = DATA ;

   PEL_CHECK_POST( data_POST( result, ct ) ) ;
   return( result ) ;
}
