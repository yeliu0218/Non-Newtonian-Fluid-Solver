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

#include <PEL_StringExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Vector.hh>
#include <PEL_Sequence.hh>
#include <PEL_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>
#include <sstream>

PEL_StringExp const*
PEL_StringExp:: PROTOTYPE_EMPTY = new PEL_StringExp( "empty" ) ;

PEL_StringExp const*
PEL_StringExp:: PROTOTYPE_TO_STRING = new PEL_StringExp( "to_string" ) ;


//----------------------------------------------------------------------
PEL_StringExp:: PEL_StringExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
{
   PEL_LABEL( "PEL_StringExp:: PEL_StringExp" ) ;
}

//----------------------------------------------------------------------
PEL_StringExp:: PEL_StringExp( PEL_Object* a_owner,
                               std::string const& a_name,
                               PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
{
   PEL_LABEL( "PEL_StringExp:: PEL_StringExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_StringExp:: ~PEL_StringExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: ~PEL_StringExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_StringExp*
PEL_StringExp:: create_replica( PEL_Object* a_owner,
                                PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_StringExp* result = 
                  new PEL_StringExp( a_owner, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_StringExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: data_type" ) ;
   PEL_Data::Type result = ( name()=="empty" ? Bool : String ) ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_StringExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   if( name()=="to_string" )
   {
      result = ( some_arguments->count() == 1 ) ;
      if( result )
      {
         PEL_Data::Type t = extract_arg( some_arguments, 0 )->data_type() ;
         result = result && ( t==Int || t==Double ) ;
      }
   }
   else if( name()=="empty" )
   {
      result = ( some_arguments->count() == 1 ) ;
      if( result )
      {
         result = result &&
         ( extract_arg( some_arguments, 0 )->data_type() == String ) ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
std::string const&
PEL_StringExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: usage" ) ;

   static std::string result ;
   if( name()=="to_string" )
   {
      result = "to_string(DS|IS)" ;
   }
   else if( name()=="empty" )
   {
      result = "empty(SS)" ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_StringExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   bool result = false ;
   
   if( name()=="empty" )
   {
      std::string const& str = arg(0)->to_string( ct ) ;
      result = str.empty() ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
std::string const&
PEL_StringExp:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_StringExp:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE(ct) ) ;
   if( name()=="to_string" )
   {
      std::ostringstream os ;
      PEL_Data::Type t = arg(0)->data_type() ;
      if( t==Double )
      {
         os << std::scientific << arg(0)->to_double( ct ) ;
      }
      else
      {
         os << arg(0)->to_int( ct ) ;
      }
      STR_RESULT = os.str() ;
   }
   
   return STR_RESULT ;
}
