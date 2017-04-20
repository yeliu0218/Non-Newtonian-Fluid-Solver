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

#include <RS_SolidBodyRotationEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <PEL_Context.hh>
#include <PEL_Data.hh>

#include <doubleVector.hh>

#include <iostream>

RS_SolidBodyRotationEXP const* 
RS_SolidBodyRotationEXP::PROTOTYPE_U = 
   new RS_SolidBodyRotationEXP( "SolidBodyRotation_value", U ) ;

RS_SolidBodyRotationEXP const* 
RS_SolidBodyRotationEXP::PROTOTYPE_B = 
   new RS_SolidBodyRotationEXP( "SolidBodyRotation_velocity", B ) ;

//----------------------------------------------------------------------
RS_SolidBodyRotationEXP:: RS_SolidBodyRotationEXP( 
   std::string const& a_name,
   Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
   , BETA( 0 )
{
}

//----------------------------------------------------------------------
RS_SolidBodyRotationEXP*
RS_SolidBodyRotationEXP:: create_replica( PEL_Object* a_owner,
                               PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_SolidBodyRotationEXP* result = 
      new RS_SolidBodyRotationEXP( a_owner, 
				   name(), 
				   argument_list, 
				   EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_SolidBodyRotationEXP:: RS_SolidBodyRotationEXP( 
   PEL_Object* a_owner,
   std::string const& a_name,
   PEL_Sequence const* argument_list,
   Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
   , BETA( 2 )
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: RS_SolidBodyRotationEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_SolidBodyRotationEXP:: ~RS_SolidBodyRotationEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: ~RS_SolidBodyRotationEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_SolidBodyRotationEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = Double ; 
         break ;
      case B : result = DoubleVector ; 
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
std::string const& 
RS_SolidBodyRotationEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = 
	 "SolidBodyRotation_value($DV_X,$DS_T)" ; 
         break ;
      case B : result = 
	 "SolidBodyRotation_velocity($DV_X,$DS_T)" ; 
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
RS_SolidBodyRotationEXP:: valid_arguments( 
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: valid_arguments" ) ;

   bool result = ( some_arguments->count() == 2 ) &&
    ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector )
    && ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) ;
   return result ;
}

//----------------------------------------------------------------------
double
RS_SolidBodyRotationEXP:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double result = PEL::max_double() ;

   doubleVector const& res = arg( 0 )->to_double_vector( ct ) ;
   double x = res( 0 ) ;
   double y = res( 1 ) ; 
   double t = arg( 1 )->to_double( ct ) ;

   // Rotation
   double xr = x*PEL::cos(2.*t)-y*PEL::sin(2.*t) ;
   double yr = x*PEL::sin(2.*t)+y*PEL::cos(2.*t) ;
   double rr = PEL::sqrt( PEL::sqr(xr+0.45)+PEL::sqr(yr) ) ;

   // Value
   double uu = PEL::max_double() ;
   if( rr<=0.35 )
   {
      uu = 1.-rr/0.35 ;
   }
   else if( xr>=0.1 && xr<=0.6 && yr>=-0.25 && yr<=0.25 )
   {
      uu = 1. ;
   }
   else
   {
      uu= 0. ;
   }
   result = uu ;

   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_SolidBodyRotationEXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SolidBodyRotationEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   doubleVector& result = BETA ;

   doubleVector const& res = arg( 0 )->to_double_vector( ct ) ;
   double x = res( 0 ) ;
   double y = res( 1 ) ; 

   result( 0 ) =  2.*y ;
   result( 1 ) = -2.*x ;

   return result ;
}
