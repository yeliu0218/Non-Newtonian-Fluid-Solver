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

#include <PEL_Variable.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Double.hh>
#include <PEL_Context.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Root.hh>
#include <PEL_List.hh>
#include <stringVector.hh>

#include <iostream>

//----------------------------------------------------------------------
PEL_Variable:: PEL_Variable( PEL_Object* a_owner,
                             std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Data( a_owner )
   , NAME( a_name )
   , KIND( PEL_Variable::data_type( a_name ) )
   , EVALUATING( false )
{
   variable_list()->append( this ) ;
   ID = variable_list()->count()-1 ;
  
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Variable*
PEL_Variable:: create_clone( PEL_Object * a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: create_clone" ) ;
   PEL_Variable* result = new PEL_Variable( a_owner, this ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Variable:: PEL_Variable( PEL_Object* a_owner,
                             PEL_Variable const* other ) 
//----------------------------------------------------------------------
      : PEL_Data( a_owner )
      , NAME( other->NAME )
      , KIND( other->KIND )
      , ID( other->ID )
      , EVALUATING( false )
{
   PEL_LABEL( "PEL_Variable:: PEL_Variable" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Variable:: ~PEL_Variable( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Variable:: nb_objects( void )
//----------------------------------------------------------------------
{
   return variable_list()->index_limit() ;
}

//----------------------------------------------------------------------
PEL_Variable const *
PEL_Variable:: object( std::string const& a_name ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: object" ) ;

   // Check name :
   PEL_Variable:: data_type( a_name ) ;
   
   PEL_Variable const* result = 0 ;
   PEL_Iterator* it = variable_list()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_Variable const* var = static_cast<PEL_Variable const*>(it->item()) ;
      if( var->name()==a_name )
      {
         result=var ;
      }
   }
   it->destroy() ; it = 0 ;
   if( result==0 )
   {
      result = new PEL_Variable( variable_list(), a_name ) ;
   }
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Variable const*
PEL_Variable:: object( size_t id )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: object" ) ;
   PEL_CHECK_PRE( id < nb_objects() ) ;
   
   PEL_Variable const* result=
             static_cast<PEL_Variable const*>( variable_list()->at( id ) ) ;

   PEL_CHECK_POST( result->id_number() == id ) ;
   return result ;
}

//----------------------------------------------------------------------
size_t
PEL_Variable:: id_number( void ) const
//----------------------------------------------------------------------
{
   return ID ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Variable:: name( void ) const
//----------------------------------------------------------------------
{
   return NAME ;
}

//----------------------------------------------------------------------
void
PEL_Variable:: declare( PEL_List* lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;
   
   lst->extend( const_cast<PEL_Variable*>(this) ) ;
   
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_Variable:: context_has_required_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: context_has_required_variables" ) ;
   PEL_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;

   bool result = ct->has_variable( this ) ;

   PEL_CHECK_POST( EQUIVALENT( result, ct->has_variable(this) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_Variable:: data_type( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: data_type(a_name)" ) ;
   
   PEL_Data::Type result = PEL_Data::Undefined ;

   if( a_name.length() >= 2 )
   {
      std::string typ = a_name.substr( 0, 2 ) ;
      if( typ=="IS" )
      {
         result = Int ;
      }
      else if( typ=="IV" )
      {
         result = IntVector ;
      }
      else if( typ=="IA" )
      {
         result = IntArray2D ;
      }
      else if( typ=="BA" )
      {
         result = BoolArray2D ;
      }
      else if( typ=="SA" )
      {
         result = StringArray2D ;
      }
      else if( typ=="DS" )
      {
         result = Double ;
      }
      else if( typ=="DV" )
      {
         result = DoubleVector ;
      }
      else if( typ=="DA" )
      {
         result = DoubleArray2D ;
      }
      else if( typ=="BS" )
      {
         result = Bool ;
      }
      else if( typ=="BV" )
      {
         result = BoolVector ;
      }
      else if( typ=="SS" )
      {
         result = String ;
      }
      else if( typ=="SV" )
      {
         result = StringVector ;
      }
   }
   if( result == PEL_Data::Undefined )
   {
      std::string msg = "\""+a_name+"\" is not a valid variable name\n" ;
      msg += "A valid name is \"XY_name\"\n" ;
      msg += "   where \"X\" is the scalar type of the variable :\n" ;
      msg += "       - \"I\" : integer\n" ;
      msg += "       - \"D\" : double\n" ;
      msg += "       - \"B\" : boolean\n" ;
      msg += "       - \"S\" : string\n" ;
      msg += "   and \"Y\" defined its dimension :\n" ;
      msg += "       - \"S\" : simple (only one element)\n" ;
      msg += "       - \"V\" : vector\n" ;
      msg += "       - \"A\" : array2D\n" ;
      msg += "Examples : \"DV_coordinates\", \"SS_name\", \"IA_connectivity\"\n" ;
      PEL_Error::object()->raise_plain( msg ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_Variable:: data_type( void ) const
//----------------------------------------------------------------------
{
   return KIND ;
}

//----------------------------------------------------------------------
bool
PEL_Variable:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: value_can_be_evaluated" ) ;

   bool result = ! EVALUATING ;
   if( result )
   {
      EVALUATING = true ;
      result =  ( ct!=0 &&
                  ct->has_variable(this) &&
                  ct->value(this)!=0 &&
                  ct->value(this)->value_can_be_evaluated(ct) ) ;
      EVALUATING = false ;
   }

   PEL_CHECK_POST(
      IMPLIES(
         result,
         ct!=0 &&
            ct->has_variable(this) &&
            ct->value(this) != 0 &&
            ct->value(this)->value_can_be_evaluated(ct) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_Variable::  undefined_variables( PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: undefined_variables" ) ;

   static stringVector result(0) ;
   result.re_initialize(0) ;
   if( EVALUATING ||
       ct == 0 || !ct->has_variable(this) || ct->value(this) == 0 )
   {
      result.append( name() ) ;
   }
   else
   {
      EVALUATING = true ;
      result = ct->value(this)->undefined_variables( ct ) ;
      EVALUATING = false ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Variable:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return data( ct )->to_bool( ct ) ;
}

//----------------------------------------------------------------------
double
PEL_Variable:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   return data( ct )->to_double( ct ) ;
}

//----------------------------------------------------------------------
int
PEL_Variable:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;   
   return data( ct )->to_int( ct ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Variable:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE( ct ) ) ;
   return data( ct )->to_string( ct ) ;
}

//----------------------------------------------------------------------
doubleVector const& 
PEL_Variable:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return data( ct )->to_double_vector( ct ) ;
}

//----------------------------------------------------------------------
intVector const&
PEL_Variable:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return data( ct )->to_int_vector( ct ) ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_Variable:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return data( ct )->to_string_vector(ct ) ;
}

//----------------------------------------------------------------------
boolVector const&
PEL_Variable:: to_bool_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_bool_vector" ) ;
   PEL_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return data( ct )->to_bool_vector( ct ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_Variable:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return data( ct )->to_double_array2D( ct ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
PEL_Variable:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return data( ct )->to_bool_array2D( ct ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
PEL_Variable:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return data( ct )->to_string_array2D( ct ) ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_Variable:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return data( ct )->to_int_array2D( ct ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PEL_Variable:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_double_array3D" ) ;
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return data( ct )->to_double_array3D( ct ) ;
}

//----------------------------------------------------------------------
intArray3D const&
PEL_Variable:: to_int_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: to_int_array3D" ) ;
   PEL_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return data( ct )->to_int_array3D( ct ) ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_Variable:: data( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( ct!=0 && ct->has_variable(this) ) ;

   PEL_Data const* result = ct->value( this ) ;

   PEL_CHECK( result!=0 ) ;   
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_Variable:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: is_raw_data" ) ;

   PEL_CHECK_POST( is_raw_data_POST( false ) ) ;
   
   return false ;
}

//----------------------------------------------------------------------
void
PEL_Variable:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << "$" << NAME ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Variable:: create_derivative( PEL_Object* a_owner,
                                  PEL_Variable const* var,
                                  PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Variable:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Data* result =
      ( var->name()==name() ?
        PEL_Double::create( a_owner, 1.0 ) :
        ct->value(this)->create_derivative( a_owner, var, ct ) ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_List*
PEL_Variable:: variable_list( void ) 
//----------------------------------------------------------------------
{
   static PEL_List* result = PEL_List::create( PEL_Root::object() ) ;
   return result ;
}
