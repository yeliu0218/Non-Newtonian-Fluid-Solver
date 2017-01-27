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

#include <PEL_BooleanExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <iostream>

PEL_BooleanExp const*  
PEL_BooleanExp::PROTOTYPE_or  = new PEL_BooleanExp( OR, "||" ) ;

PEL_BooleanExp const*
PEL_BooleanExp::PROTOTYPE_and = new PEL_BooleanExp( AND, "&&" ) ;

PEL_BooleanExp const*
PEL_BooleanExp::PROTOTYPE_not = new PEL_BooleanExp( NOT, "!" ) ;

//----------------------------------------------------------------------
PEL_BooleanExp:: PEL_BooleanExp( BoolExp exp_id, 
                                 std::string const& a_name  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( exp_id )
   , ARG0( 0 )
   , ARG1( 0 )
{
   PEL_LABEL( "PEL_BooleanExp:: PEL_BooleanExp" ) ;
}

//----------------------------------------------------------------------
PEL_BooleanExp*
PEL_BooleanExp:: create_replica( PEL_Object* a_owner,
                                 PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BooleanExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_BooleanExp* result = new PEL_BooleanExp( a_owner, 
                                                OP,
                                                name(),
                                                argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_BooleanExp:: PEL_BooleanExp( PEL_Object* a_owner,
                                 BoolExp exp_id,
                                 std::string const& a_name,
                                 PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list ),
     OP( exp_id ),
     ARG0( arg(0) ),
     ARG1( nb_arguments()>1 ? arg(1) : 0 )
{
   PEL_LABEL( "PEL_BooleanExp:: PEL_BooleanExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_BooleanExp:: ~PEL_BooleanExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BooleanExp:: ~PEL_BooleanExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_BooleanExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case OR :
         result = "BS || BS" ;
	 break ;
      case AND :
         result = "BS && BS" ;
	 break ;
      case NOT :
         result = "! BS" ;
	 break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_BooleanExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BooleanExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( OP )
   {
      case OR  :
      case AND :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
	    result = result && ( t0==Bool && t1==Bool ) ;
	 }
	 break ;
      case NOT :
         result = ( some_arguments->count() == 1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = result && ( t0 == Bool ) ;
         }
	 break ;
      default : result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_BooleanExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return PEL_Data::Bool ;
}

//----------------------------------------------------------------------
bool
PEL_BooleanExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BooleanExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result ;
   switch( OP )
   {
      case AND :
         result = ( ARG0->to_bool(ct) && ARG1->to_bool(ct) ) ;
         break ;
      case OR :
         result = ( ARG0->to_bool(ct) || ARG1->to_bool(ct) ) ;
         break ;
      case NOT :
         result = !( ARG0->to_bool( ct ) ) ;
         break ;
      default :
	PEL_Error::object()->raise_internal( "Bad operator" ) ;
	result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
PEL_BooleanExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   if( external_brackets_are_set() ) os << "(" ;
   switch( OP )
   {
      case AND :
         ARG0->print( os, 0 ) ;
         os << " && " ;
         ARG1->print( os, 0 ) ;
         break ;
      case OR :
         ARG0->print( os, 0 ) ;
         os << " || " ;
         ARG1->print( os, 0 ) ;
         break ;
      case NOT :
         os << "! " ;
         ARG0->print( os, 0 ) ;
         break ;
      default :
	PEL_Error::object()->raise_internal( "Bad operator" ) ;
   }
   if( external_brackets_are_set() ) os << ")" ;
}
