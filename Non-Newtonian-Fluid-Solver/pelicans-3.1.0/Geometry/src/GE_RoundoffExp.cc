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

#include <GE_RoundoffExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>
#include <fstream>

GE_RoundoffExp const* 
GE_RoundoffExp::PROTOTYPE = new GE_RoundoffExp( "default_roundoff" ) ;

//----------------------------------------------------------------------
GE_RoundoffExp:: GE_RoundoffExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
{
   PEL_LABEL( "GE_RoundoffExp:: GE_RoundoffExp" ) ;
}

//----------------------------------------------------------------------
PEL_Expression*
GE_RoundoffExp:: create_replica( PEL_Object* a_owner,
                                 PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   GE_RoundoffExp* result = new GE_RoundoffExp( a_owner,
                                                name(), 
                                                argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_RoundoffExp:: GE_RoundoffExp( PEL_Object* a_owner,
                                 std::string const& a_name,
                                 PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
{
   PEL_LABEL( "GE_RoundoffExp:: GE_RoundoffExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_RoundoffExp:: ~GE_RoundoffExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: ~GE_RoundoffExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
GE_RoundoffExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: data_type" ) ;

   PEL_Data::Type result = PEL_Data::Double ;

   return( result ) ;
}

//----------------------------------------------------------------------
double
GE_RoundoffExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;

   double   xx = arg( 0 )->to_double( ct ) ;
   int    nsig = arg( 1 )->to_int( ct ) ;
   double xmin = arg( 2 )->to_double( ct ) ;

   double result = PEL::bad_double() ;

   double a = PEL::abs( xx ) ;
   if( a > xmin )
   {
      int n = static_cast<int>( PEL::log10( a ) ) ;
      double q = PEL::pow( 10.0, (double)(nsig-n) ) ;
      result = (int)( a*q + 0.5 )/q ;
   }
   else
   {
      result = 0.0 ;
   }
   if( xx < 0.0 ) result = -result ;

   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
GE_RoundoffExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: usage" ) ;

   static std::string result =  "default_roundoff(DV,IS,DS)" ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_RoundoffExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RoundoffExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = ( some_arguments->count() == 3 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Int )  &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) ;
   return( result ) ;
}

