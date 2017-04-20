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

#include <PEL_ContextSimple.hh>

#include <PEL_assertions.hh>
#include <PEL_Iterator.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>
#include <PEL_Variable.hh>

//----------------------------------------------------------------------
PEL_ContextSimple*
PEL_ContextSimple:: create( PEL_Object* a_owner ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: create" ) ;
   PEL_ContextSimple* result = new PEL_ContextSimple( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ContextSimple:: PEL_ContextSimple( PEL_Object* a_owner ) 
//----------------------------------------------------------------------
   : PEL_Context( a_owner )
   , VALUES( PEL_Vector::create( this, 0 ) )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ContextSimple*
PEL_ContextSimple:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: create_clone" ) ;

   PEL_ContextSimple* result = new PEL_ContextSimple( a_owner ) ;
   
   PEL_Iterator* it = VALUES->create_iterator( 0 ) ;
   for( size_t i=0 ; i<VALUES->index_limit() ; i++  )
   {
      PEL_Data const* dat = static_cast<PEL_Data const*>( VALUES->at(i) ) ;
      if( dat!=0 )
      {
         PEL_Variable const*var=PEL_Variable::object(i) ;
         result->extend( var, dat->create_clone( result ) ) ;
      }
   }
   it->destroy() ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ContextSimple:: ~PEL_ContextSimple( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   notify_observers_of_my_destruction() ;
}


//----------------------------------------------------------------------
size_t
PEL_ContextSimple:: nb_variables( void ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: nb_variables" ) ;
   return( VALUES->count() ) ;
}

//----------------------------------------------------------------------
PEL_Variable const* 
PEL_ContextSimple:: variable( size_t i ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: variable" ) ;
   PEL_CHECK_PRE( variable_PRE(i) ) ;
   
   PEL_Variable const* result = 0 ;
   for( size_t j=0 ; j<VALUES->index_limit() ; j++ )
   {
      if( VALUES->at(j)!=0 )
      {
         if( i==0 )
         {
            result = PEL_Variable::object(j) ;
            break ;
         }
         i-- ;
      }
   }
   
   PEL_CHECK_POST( variable_POST( result ) ) ;
   return( result ) ;
      
}

//----------------------------------------------------------------------
bool
PEL_ContextSimple:: has_variable( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: has_variable" ) ;
   PEL_CHECK_PRE( has_variable_PRE( var ) ) ;

   size_t idx = var->id_number() ;
   bool result = ( VALUES->index_limit()>idx && VALUES->at( idx )!=0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data* 
PEL_ContextSimple:: value( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: value" ) ;
   PEL_CHECK_PRE( value_PRE( var ) ) ;
   
   PEL_Data* result = static_cast<PEL_Data* >(
      VALUES->at( var->id_number() ) ) ;

   PEL_CHECK_POST( value_POST( result, var ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_ContextSimple:: has_circular_definition( PEL_Variable const* var,
                                             PEL_Data const* a_value ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: has_circular_definition" ) ;
   PEL_CHECK_PRE( var != 0 ) ;
   PEL_CHECK_PRE( a_value != 0 ) ;
   
   PEL_List* dummy_context = PEL_List::create( 0 ) ;
   a_value->declare( dummy_context ) ;
   bool result = dummy_context->has( var ) ;
   dummy_context->destroy() ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ContextSimple:: extend( PEL_Variable const* var,
                            PEL_Data const* a_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: extend" ) ;
   PEL_CHECK_PRE( var != 0 ) ;
   PEL_CHECK_PRE( a_value != 0 ) ;
   PEL_CHECK_PRE( !has_circular_definition( var, a_value ) ) ;
   PEL_CHECK_PRE( var->data_type()==a_value->data_type() ) ;
   PEL_CHECK_PRE( a_value->is_under_ownership_of( this ) ) ;
   
   size_t idx = var->id_number() ;
   if( VALUES->index_limit()<=idx )
   {
      VALUES->resize( idx+10 ) ;
   }
   VALUES->set_at( idx, const_cast<PEL_Data*>(a_value) ) ;
   update_observers() ;
   PEL_CHECK_POST( has_variable( var ) ) ;
   PEL_CHECK_POST( value(var) == a_value ) ;
}

//----------------------------------------------------------------------
void
PEL_ContextSimple:: extend( PEL_Context const* other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: extend( PEL_Context const* )" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   
   for( size_t i=0 ; i<other->nb_variables() ; i++ )
   {
      PEL_Variable const* var = other->variable(i) ;
      extend( other->variable(i), other->value(var)->create_clone( this ) ) ;
   }
   update_observers() ;
}

//----------------------------------------------------------------------
void
PEL_ContextSimple:: set_value_of( PEL_Variable const* var,
                                  PEL_Data const* a_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextSimple:: set_value_of" ) ;
   PEL_CHECK_PRE( var !=0 ) ;
   PEL_CHECK_PRE( has_variable( var ) ) ;
   PEL_CHECK_PRE( a_value !=0 ) ;
   PEL_CHECK_PRE( var->data_type()==a_value->data_type() ) ;
   PEL_CHECK_PRE( !has_circular_definition(var,a_value ) ) ;
   PEL_CHECK_PRE( a_value->is_under_ownership_of( this ) ) ;
   
   VALUES->set_at( var->id_number(), const_cast<PEL_Data*>(a_value) ) ;
   update_observers() ;

   PEL_CHECK_POST( value( var )==a_value ) ;
   
}

//----------------------------------------------------------------------
void
PEL_ContextSimple:: update( void )
//----------------------------------------------------------------------
{
   update_observers() ;
}

//----------------------------------------------------------------------
void
PEL_ContextSimple:: update_for_destruction_of( PEL_Context const* subject )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
PEL_ContextSimple:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( VALUES!=0 ) ;
   return( true ) ;
}




