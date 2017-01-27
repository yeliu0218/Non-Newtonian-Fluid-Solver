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

#include <PEL_ComparisonExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <iostream>

PEL_ComparisonExp const*
PEL_ComparisonExp::PROTOTYPE_min = new PEL_ComparisonExp( "min", min_op ) ;

PEL_ComparisonExp const*
PEL_ComparisonExp::PROTOTYPE_max = new PEL_ComparisonExp( "max", max_op ) ;

//----------------------------------------------------------------------
PEL_ComparisonExp:: PEL_ComparisonExp(
   std::string const& a_name,
   Function a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , FIRST( 0 ) 
   , SECOND( 0 )
{
   PEL_LABEL( "PEL_ComparisonExp:: PEL_ComparisonExp" ) ;
}

//----------------------------------------------------------------------
PEL_ComparisonExp:: PEL_ComparisonExp( PEL_Object* a_owner,
				       std::string const& a_name,
				       PEL_Sequence const* argument_list,
				       Function a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP(a_op)
   , FIRST( arg(0) ) 
   , SECOND( arg(1) )
{
   PEL_LABEL( "PEL_ComparisonExp:: PEL_ComparisonExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ComparisonExp:: ~PEL_ComparisonExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ComparisonExp:: ~PEL_ComparisonExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ComparisonExp*
PEL_ComparisonExp:: create_replica( PEL_Object* a_owner,
				    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ComparisonExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ComparisonExp* result = new PEL_ComparisonExp( a_owner, 
						      name(), 
						      argument_list, 
						      OP ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_ComparisonExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = name() + "(DS,DS) " ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ComparisonExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ComparisonExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      result &=
         ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
         ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ComparisonExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return Double ;
}

//----------------------------------------------------------------------
double
PEL_ComparisonExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ComparisonExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   
   double result = PEL::bad_double() ;
   double v1 = FIRST->to_double( ct ) ;
   double v2 = SECOND->to_double( ct ) ;
   if( OP == min_op ) result = PEL::min( v1, v2 ) ;
   else if ( OP == max_op ) result = PEL::max( v1, v2 ) ;
   return( result ) ;
}
