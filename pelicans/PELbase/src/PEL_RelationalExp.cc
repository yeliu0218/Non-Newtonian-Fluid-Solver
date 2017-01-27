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

#include <PEL_RelationalExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <iostream>

PEL_RelationalExp const* 
PEL_RelationalExp::PROTOTYPE_LE = new PEL_RelationalExp( "<=", LE ) ;

PEL_RelationalExp const* 
PEL_RelationalExp::PROTOTYPE_GE = new PEL_RelationalExp( ">=", GE ) ;

PEL_RelationalExp const* 
PEL_RelationalExp::PROTOTYPE_LT = new PEL_RelationalExp( "<", LT ) ;

PEL_RelationalExp const*
PEL_RelationalExp::PROTOTYPE_GT = new PEL_RelationalExp( ">", GT ) ;

PEL_RelationalExp const* 
PEL_RelationalExp::PROTOTYPE_EQ = new PEL_RelationalExp( "=", EQ ) ;

PEL_RelationalExp const* 
PEL_RelationalExp::PROTOTYPE_NEQ = new PEL_RelationalExp( "!=", NEQ ) ;

//----------------------------------------------------------------------
PEL_RelationalExp*
PEL_RelationalExp:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RelationalExp:: create_replica" ) ;
   return new PEL_RelationalExp( a_owner, name(), argument_list, OP ) ;
}

//----------------------------------------------------------------------
PEL_RelationalExp:: PEL_RelationalExp( std::string const& a_name,
                                       ComparisonOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP(a_op)
   , ARG0(0)
   , ARG1(0)
{
   PEL_LABEL( "PEL_RelationalExp:: PEL_RelationalExp" ) ;
}

//----------------------------------------------------------------------
PEL_RelationalExp:: PEL_RelationalExp( PEL_Object* a_owner,
                                      std::string const& a_name,
                                      PEL_Sequence const* argument_list,
                                      ComparisonOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP(a_op)
   , ARG0( arg(0) )
   , ARG1( arg(1) )
{
   PEL_LABEL( "PEL_RelationalExp:: PEL_RelationalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_RelationalExp:: ~PEL_RelationalExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RelationalExp:: ~PEL_RelationalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_RelationalExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( PEL_Data::Bool ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_RelationalExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   if( OP == EQ || OP == NEQ )
   {
      result = "IS|DS|SS|BS " + name() + " <same type>" ;
   }
   else
   {
      result = "IS|DS " + name() + " <same type>" ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_RelationalExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RelationalExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      PEL_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      PEL_Data::Type k1 =  extract_arg( some_arguments, 1 )->data_type() ;
      result = ( k0==k1 ) ;
      if( result )
      {
         if( OP == EQ || OP == NEQ )
         {
            result = k0==Double || k0==Int || k0==String || k0== Bool ;
         }
         else
         {
            result = k0==Double || k0==Int ;
         }
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_RelationalExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RelationalExp:: to_bool" ) ;
   
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   bool result = true ;
   PEL_Data::Type k0 = ARG0->data_type() ;
   
   if( OP==EQ && k0==String )
   {
      result = ( ARG0->to_string( ct ) == ARG1->to_string( ct ) ) ;
   }
   else if( OP==EQ && k0==Bool )
   {
      result = ( ARG0->to_bool( ct ) == ARG1->to_bool( ct ) ) ;
   }
   else if( OP==NEQ && k0==String )
   {
      result = ( ARG0->to_string( ct ) != ARG1->to_string( ct ) ) ;
   }
   else if( OP==NEQ && k0==Bool )
   {
      result = ( ARG0->to_bool( ct ) != ARG1->to_bool( ct ) ) ;
   }
   else
   {
      double v1 ;
      double v2 ;
      if( ARG0->data_type()==Double )
      {
         v1 = ARG0->to_double( ct ) ;
      }
      else
      {
         v1 = ARG0->to_int( ct ) ;
      }
      if( ARG1->data_type()==Double )
      {
         v2 = ARG1->to_double( ct ) ;
      }
      else
      {
         v2 = ARG1->to_int( ct ) ;
      }
      switch(OP)
      {
         case LE : result = ( v1 <= v2 ) ;
            break ;
         case GE : result = ( v1 >= v2 ) ;
            break ;
         case LT : result = ( v1 < v2 ) ;
            break ;
         case GT : result = ( v1 > v2 ) ;
            break ;
         case EQ : result = ( v1 == v2 ) ;
            break ;
         case NEQ : result = ( v1 != v2 ) ;
            break ;
         default :
            PEL_Error::object()->raise_plain( "Internal error" ) ;
      }
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_RelationalExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   if( external_brackets_are_set() ) os << "(" ;
   ARG0->print( os, 0 ) ;
   os << name() ;
   ARG1->print( os, 0 ) ;
   if( external_brackets_are_set() ) os << ")" ;
}
