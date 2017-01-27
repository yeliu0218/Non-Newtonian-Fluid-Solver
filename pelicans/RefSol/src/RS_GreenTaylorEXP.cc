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

#include <RS_GreenTaylorEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_GreenTaylorEXP const* 
RS_GreenTaylorEXP::PROTOTYPE_U = 
             new RS_GreenTaylorEXP( "GreenTaylor_velocity", U ) ;

RS_GreenTaylorEXP const* 
RS_GreenTaylorEXP::PROTOTYPE_P = 
             new RS_GreenTaylorEXP( "GreenTaylor_pressure", P ) ;

RS_GreenTaylorEXP const* 
RS_GreenTaylorEXP::PROTOTYPE_UD = 
             new RS_GreenTaylorEXP( "GreenTaylor_grad_velocity", UD ) ;

//----------------------------------------------------------------------
RS_GreenTaylorEXP:: RS_GreenTaylorEXP( std::string const& a_name,
                                       Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
}

//----------------------------------------------------------------------
RS_GreenTaylorEXP*
RS_GreenTaylorEXP:: create_replica( 
                             PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_GreenTaylorEXP* result = new RS_GreenTaylorEXP( a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_GreenTaylorEXP:: RS_GreenTaylorEXP( PEL_Object* a_owner,
				       std::string const& a_name,
				       PEL_Sequence const* argument_list,
				       Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
   PEL_LABEL( "RS_GreenTaylorEXP:: RS_GreenTaylorEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_GreenTaylorEXP:: ~RS_GreenTaylorEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: ~RS_GreenTaylorEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_GreenTaylorEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : 
         result = DoubleVector ; 
         break ;
      case P : 
         result = DoubleVector ;
         break ;
      case UD : 
         result = DoubleArray2D ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_GreenTaylorEXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( EXPR == P ? DV_result_1 : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
   double x = res( 0 ) * 2. * PEL::pi() ;
   double y = res( 1 ) * 2. * PEL::pi() ; 
   double t = arg(1)->to_double(ct) ;
   double mu = arg(2)->to_double(ct) ;
   double expt = PEL::exp( -8.*mu*PEL::pi()*PEL::pi()*t ) ;

   if( EXPR == P )
   {
      result( 0 ) = -0.25*( PEL::cos(2.*x)+PEL::cos(2.*y) )*expt*expt ;
   }
   else if( EXPR == U )
   { 
      result( 0 ) = -PEL::cos(x)*PEL::sin(y)*expt ;
      result( 1 ) =  PEL::sin(x)*PEL::cos(y)*expt ;
   }

   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_GreenTaylorEXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ; 

   doubleArray2D& result = doubleArray2D_result ;

   if( EXPR == UD )
   {
      doubleVector const& res = arg(0)->to_double_vector(ct) ;
      double x = res( 0 ) * 2. * PEL::pi() ;
      double y = res( 1 ) * 2. * PEL::pi() ; 
      double t = arg(1)->to_double(ct) ;
      double mu = arg(2)->to_double(ct) ;
      double expt = PEL::exp( -8.*mu*PEL::pi()*PEL::pi()*t ) ;

      result( 0, 0 ) =  2.*PEL::pi()* PEL::sin(x)*PEL::sin(y)*expt ;
      result( 0, 1 ) = -2.*PEL::pi()* PEL::cos(x)*PEL::cos(y)*expt ;
      result( 1, 0 ) =  2.*PEL::pi()* PEL::cos(x)*PEL::cos(y)*expt ;
      result( 1, 1 ) = -2.*PEL::pi()* PEL::sin(x)*PEL::sin(y)*expt ;
   }

   return result ;
}

//----------------------------------------------------------------------
std::string const& 
RS_GreenTaylorEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : 
         result = "GreenTaylor_velocity($DV_X,$DS_T,$DS_NU)" ; 
         break ;
      case P : 
         result = "GreenTaylor_pressure($DV_X,$DS_T,$DS_NU)" ;
         break ;
      case UD : 
         result = "GreenTaylor_grad_velocity($DV_X,$DS_T,$DS_NU)" ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
RS_GreenTaylorEXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_GreenTaylorEXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = ( some_arguments->count() == 3 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )  &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
