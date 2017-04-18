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

#include <PEL_Data.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_BoolArray2D.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_IntVector.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_IntArray3D.hh>
#include <PEL_List.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleArray3D.hh>
#include <PEL_String.hh>
#include <PEL_StringArray2D.hh>
#include <PEL_StringVector.hh>
#include <PEL_Variable.hh>
#include <boolVector.hh>
#include <boolArray2D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <stringArray2D.hh>

#include <sstream>
#include <iostream>

//----------------------------------------------------------------------
PEL_Data:: PEL_Data( PEL_Object* a_owner ) 
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   PEL_LABEL( "PEL_Data:: PEL_Data" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data:: ~PEL_Data( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: ~PEL_Data" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PEL_Data:: declare( PEL_List* lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: is_constant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: is_constant" ) ;
   PEL_List* lst = PEL_List::create( 0 ) ;
   declare( lst ) ;
   bool result = lst->count()==0 ;
   lst->destroy() ;
   
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: is_raw_data" ) ;

   PEL_CHECK_POST( is_raw_data_POST( true ) ) ;
   
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: is_raw_data_POST( bool result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result, is_constant() ) ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: context_has_required_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: context_has_required_variables" ) ;
   PEL_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Data:: type_name( Type kind ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: type_name" ) ;
   std::string result ;
   switch( kind )
   {
      case Double :
         result = "Double" ;
         break ;
      case Int :
         result = "Int" ;
         break ;
      case Bool :
         result = "Bool" ;
         break ;
      case String :
         result = "String" ;
         break ;
      case DoubleVector :
         result = "DoubleVector" ;
         break ;
      case IntVector :
         result = "IntVector" ;
         break ;
      case BoolVector :
         result = "BoolVector" ;
         break ;
      case DoubleArray2D :
         result = "DoubleArray2D" ;
         break ;
      case BoolArray2D :
         result = "BoolArray2D" ;
         break ;
      case StringArray2D :
         result = "StringArray2D" ;
         break ;
      case DoubleArray3D :
         result = "DoubleArray3D" ;
         break ;
      case IntArray2D :
         result = "IntArray2D" ;
         break ;
      case IntArray3D :
         result = "IntArray3D" ;
         break ;
      case StringVector :
         result = "StringVector" ;
         break ;
      default :
         PEL_Error::object()->raise_plain( "Bad kind" ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   return true ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_Data:: undefined_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: undefined_variables" ) ;
   static stringVector result(0) ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   exitWithError( "to_bool" ) ;
   return 0 ;
}

//----------------------------------------------------------------------
double
PEL_Data:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   exitWithError( "to_double" ) ;
   return 0 ;
}

//----------------------------------------------------------------------
int
PEL_Data:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;
   exitWithError( "to_int" ) ;
   return 0 ;
}

//----------------------------------------------------------------------
std::string
const&
PEL_Data:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE( ct ) ) ;
   exitWithError( "to_string" ) ;
   static std::string dummy ;
   return dummy ;
}

//----------------------------------------------------------------------
doubleVector
const& 
PEL_Data:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   exitWithError( "to_double_vector" ) ;
   static doubleVector ret(0) ;
   return ret ;
}

//----------------------------------------------------------------------
intVector
const&
PEL_Data:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   exitWithError( "to_int_vector" ) ;
   static intVector ret(0) ;
   return ret ;
}

//----------------------------------------------------------------------
stringVector const& 
PEL_Data:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   exitWithError( "to_string_vector" ) ;
   static stringVector ret( 0 ) ;
   return ret ;     
}

//----------------------------------------------------------------------
boolVector const&
PEL_Data:: to_bool_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_bool_vector" ) ;
   PEL_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   exitWithError( "to_bool_vector" ) ;
   static boolVector ret( 0 ) ;
   return ret ;  
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_Data:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   exitWithError( "to_double_array2D" ) ;
   static doubleArray2D ret(0,0) ;  
   return ret ;
}

//----------------------------------------------------------------------
boolArray2D const&
PEL_Data:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   exitWithError( "to_bool_array2D" ) ;
   static boolArray2D ret(0,0) ;  
   return ret ;
}

//----------------------------------------------------------------------
stringArray2D const&
PEL_Data:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   exitWithError( "to_string_array2D" ) ;
   static stringArray2D ret(0,0) ;  
   return ret ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_Data:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   exitWithError( "to_int_array2D" ) ;
   static intArray2D ret(0,0) ;  
   return( ret ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PEL_Data:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_double_array3D" ) ;
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   exitWithError( "to_double_array3D" ) ;
   static doubleArray3D ret(0,0,0) ;  
   return ret ;
}

//----------------------------------------------------------------------
intArray3D const&
PEL_Data:: to_int_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: to_int_array3D" ) ;
   PEL_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   exitWithError( "to_int_array3D" ) ;
   static intArray3D ret(0,0,0) ;  
   return( ret ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Data:: value_as_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: value_as_string" ) ;
   PEL_CHECK_PRE( value_can_be_evaluated( ct ) ) ;
   
   std::ostringstream strout ;
   strout.precision(16) ;
   
   if( data_type()==String ) 
   {
      strout << "\"" << to_string( ct ) << "\"" ;
   }
   else if( data_type()==Bool ) 
   {
      strout <<  ( to_bool(ct) ? "true" : "false" ) ;
   }
   else if( data_type()==Int ) 
   {
      strout << to_int(ct) ;
   }
   else if( data_type()==Double )
   {
      PEL::print_double(strout, to_double(ct)) ;
   }
   else if( data_type()==StringVector )
   {
      stringVector const& val = to_string_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << "\"" << val(i) << "\" "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==BoolVector )
   {
      boolVector const& val = to_bool_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << ( val(i) ? "true" : "false" ) <<  " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==IntVector )
   {
      intVector const& val = to_int_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << val(i) << " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==DoubleVector )
   {
      doubleVector const& val = to_double_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         PEL::print_double(strout, val(i) ) ;
         strout << " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==BoolArray2D )
   {
      boolArray2D const& val = to_bool_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << ( val(i,j) ? "true" : "false" ) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==IntArray2D )
   {
      intArray2D const& val = to_int_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << val(i,j) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==StringArray2D )
   {
      stringArray2D const& val = to_string_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << "\"" << val(i,j) << "\"" ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==DoubleArray2D )
   {
      doubleArray2D const& val = to_double_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            PEL::print_double(strout, val(i,j) ) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else
   {
      PEL_Error::object()->raise_internal(
         "*** PEL_Data:: value_as_string error\n"
         "    unexpected data type \""+type_name( data_type() )+"\"" ) ;
   }      
   return( strout.str() ) ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Data:: create_derivative( PEL_Object* a_owner,
                              PEL_Variable const* var,
                              PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Error::object()->raise_plain(
      "No create_derivative method implemented for expression of class "
      + PEL_Object::type_name() ) ;
   PEL_Data* result = 0 ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: create_derivative_PRE( PEL_Object* a_owner,
                                  PEL_Variable const* var,
                                  PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( var!=0 ) ;
   PEL_ASSERT( data_type()==Double || data_type()==DoubleVector ) ;
   PEL_ASSERT( var->data_type()==Double ) ;
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: create_derivative_POST( PEL_Object* a_owner,
                                   PEL_Variable const* var,
                                   PEL_Data const* result ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->owner()==a_owner ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Data:: create_simplification( PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: create_simplification" ) ;
   PEL_CHECK_PRE( data_type()!=Undefined ) ;

   PEL_Data* result = 0 ;
   
   if( is_constant() )
   {
      switch( data_type() ) 
      {
         case Int :
            result = PEL_Int::create( a_owner, to_int() ) ;
            break ;
         case IntVector :
            result = PEL_IntVector::create( a_owner, to_int_vector() ) ;
            break ;
         case IntArray2D :
            result = PEL_IntArray2D::create( a_owner, to_int_array2D() ) ;
            break ;
         case BoolArray2D :
            result = PEL_BoolArray2D::create( a_owner, to_bool_array2D() ) ;
            break ;
         case StringArray2D :
            result = PEL_StringArray2D::create( a_owner, to_string_array2D() ) ;
            break ;
         case IntArray3D :
            result = PEL_IntArray3D::create( a_owner, to_int_array3D() ) ;
            break ;
         case Double :
            result = PEL_Double::create( a_owner, to_double() ) ;
            break ;
         case DoubleVector :
            result = PEL_DoubleVector::create( a_owner, to_double_vector() ) ;
            break ;
         case DoubleArray2D :
            result = PEL_DoubleArray2D::create( a_owner,to_double_array2D() ) ;
            break ;
         case DoubleArray3D :
            result = PEL_DoubleArray3D::create( a_owner,to_double_array3D() ) ;
            break ;
         case String :
            result = PEL_String::create( a_owner, to_string() ) ;
            break ;
         case StringVector :
            result = PEL_StringVector::create( a_owner, to_string_vector() ) ;
            break ;
         case Bool :
            result = PEL_Bool::create( a_owner, to_bool() ) ;
            break ;
        default :
            PEL_Error::object()->raise_internal( "Bad type" ) ;
      }
   }
   else
   {
      result = create_non_const_simplification( a_owner ) ;
   }
   
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_Data:: create_non_const_simplification( PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data:: create_non_const_simplification" ) ;
   PEL_CHECK_PRE( create_non_const_simplification_PRE( a_owner ) ) ;

   PEL_Data* result = create_clone( a_owner ) ;

   PEL_CHECK_POST( create_non_const_simplification_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: declare_PRE( PEL_List const* lst ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( lst!=0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: declare_POST( PEL_List const* lst ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<lst->count(); i++),
                       dynamic_cast<PEL_Variable*>(lst->at(i))!=0 ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: context_has_required_variables_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( ct != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_double_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==Double ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_int_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==Int ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_bool_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==Bool ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_string_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==String ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_double_vector_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==DoubleVector ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_int_vector_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==IntVector ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_bool_vector_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==BoolVector ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_string_vector_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==StringVector ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_double_array2D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==DoubleArray2D ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_bool_array2D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==BoolArray2D ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_string_array2D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==StringArray2D ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_int_array2D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==IntArray2D ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_double_array3D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==DoubleArray3D ) ;
   return true ;  
}

//----------------------------------------------------------------------
bool
PEL_Data:: to_int_array3D_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   PEL_ASSERT( data_type()==IntArray3D ) ;
   return true ;  
}

//----------------------------------------------------------------------
void
PEL_Data:: exitWithError( std::string const& mess ) const
//----------------------------------------------------------------------
{
   PEL::out() << "Current lexical object : " << std::endl ;
   display_info( PEL::out(), 5 ) ;
   PEL::out() << "can't be converted through " << mess << " method ! " ;
   PEL_Error::object()->raise_plain(
      "Error in PEL_Data derived class in method "+mess ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: create_non_const_simplification_PRE( PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( data_type()!=Undefined ) ;
   PEL_ASSERT( !is_constant() ) ;
   
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: create_non_const_simplification_POST( PEL_Object* a_owner,
                                                 PEL_Data const* result ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Data:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return true ;  
}

