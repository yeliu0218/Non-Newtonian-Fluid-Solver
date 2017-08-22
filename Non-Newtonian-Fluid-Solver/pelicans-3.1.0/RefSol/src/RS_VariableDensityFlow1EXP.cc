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

#include <RS_VariableDensityFlow1EXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_U = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_velocity", U ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_UD = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_grad_velocity", UD ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_P = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_pressure", P ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_RHO = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_rho", RHO ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_RHODT = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_rhodt", RHODT ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_RHOD = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_grad_rho", RHOD ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_RHOU = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_rho_velocity", 
                                   RHOU ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_RHOUD = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_grad_rho_velocity", 
                                   RHOUD ) ;

RS_VariableDensityFlow1EXP const* 
RS_VariableDensityFlow1EXP:: PROTOTYPE_F = 
   new RS_VariableDensityFlow1EXP( "VariableDensityFlow1_rhsf", F ) ;

//----------------------------------------------------------------------
RS_VariableDensityFlow1EXP:: RS_VariableDensityFlow1EXP( 
                                                   std::string const& a_name,
                			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
}

//----------------------------------------------------------------------
RS_VariableDensityFlow1EXP*
RS_VariableDensityFlow1EXP:: create_replica( 
                                PEL_Object* a_owner,
				PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_VariableDensityFlow1EXP* result = new RS_VariableDensityFlow1EXP( 
                                                      a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_VariableDensityFlow1EXP:: RS_VariableDensityFlow1EXP(
                                   PEL_Object* a_owner,
   			           std::string const& a_name,
			           PEL_Sequence const* argument_list,
			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: RS_VariableDensityFlow1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_VariableDensityFlow1EXP:: ~RS_VariableDensityFlow1EXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: ~RS_VariableDensityFlow1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_VariableDensityFlow1EXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = DoubleVector ; 
         break ;
      case F : result = DoubleVector ; 
         break ;
      case P : result = DoubleVector ;
	 break ;
      case RHO : result = DoubleVector ;
	 break ;
      case RHOD : result = DoubleArray2D ;
	 break ;
      case RHODT : result = DoubleVector ;
	 break ;
      case RHOU : result = DoubleVector ;
	 break ;
      case UD : result = DoubleArray2D ;
         break ;
      case RHOUD : result = DoubleArray2D ;
         break ;
   }
   return result ;      
}

//----------------------------------------------------------------------
doubleVector const&
RS_VariableDensityFlow1EXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( ( EXPR == P || EXPR == RHO || EXPR == RHODT ) ? 
                               DV_result_1 
                            : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
   double t     = arg(1)->to_double(ct) ;
   double alpha = arg(2)->to_double(ct) ;
   double beta  = arg(3)->to_double(ct) ;

   if( EXPR == U )
   {
      double y = res( 1 ) ; 
      result( 0 ) = alpha*y*(1.-y)*(2.+ PEL::cos(beta*t)) ;
      result( 1 ) = 0.0 ;
   }
   else if( EXPR == P )
   { 
      double x = res( 0 ) ; 
      double mu = arg(4)->to_double(ct) ;
      double xmax = arg(8)->to_double(ct) ;
      result( 0 ) = -2.*mu*alpha*(2.+ PEL::cos(beta*t))*( x-0.5*xmax ) ;
   }
   else if( EXPR == RHO )
   { 
      double x = res( 0 ) ;
      double y = res( 1 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      result( 0 ) = cc *
                   ((xint-x0)*(xint-x0)*(-2.*xint+3.*x1-x0)/PEL::pow(h,3)+1.) ;
   }
   else if( EXPR == RHODT )
   {
      double x = res( 0 ) ;
      double y = res( 1 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      double ux = alpha*y*(1.-y)*(2.+ PEL::cos(beta*t)) ;
      result( 0 ) = 6.*cc* ux * (xint-x0)*(x1-xint)/PEL::pow( h, 3 ) ;
     }
   else if( EXPR == RHOU )
   {
      double x = res( 0 ) ;
      double y = res( 1 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      double ux = alpha*y*(1.-y)*(2.+ PEL::cos(beta*t)) ;
      result( 0 )= ux * cc *
                  ((xint-x0)*(xint-x0)*(-2.*xint+3.*x1-x0)/PEL::pow(h,3)+1.);
      result( 1 )= 0. ;
   }
   else if( EXPR == F )
   { 
      double x = res( 0 ) ;
      double y = res( 1 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      result( 0 ) = cc * alpha*y*(1.-y)*(-beta*PEL::sin(beta*t))
	       *( (xint-x0)*(xint-x0)*(-2.*xint+3.*x1-x0)/PEL::pow(h,3)+1.) ;
      result( 1 ) = 0. ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_VariableDensityFlow1EXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   doubleArray2D& result = doubleArray2D_result ;
   doubleVector const& res = arg(0)->to_double_vector(ct) ;
   double y = res( 1 ) ; 
   double t = arg(1)->to_double(ct) ;
   double alpha = arg(2)->to_double(ct) ;
   double beta  = arg(3)->to_double(ct) ;
   
   if( EXPR == UD )
   {
      result( 0, 0 ) = 0 ;
      result( 0, 1 ) = alpha*(1.-2.*y )*(2+PEL::cos(beta*t))  ;
      result( 1, 0 ) = 0.0 ;
      result( 1, 1 ) = 0.0 ;
   }
   else if( EXPR == RHOUD )
   {  
      double x = res( 0 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      double ux = alpha*y*(1.-y)*(2.+ PEL::cos(beta*t)) ;     
      result( 0, 0 ) =  ux*6.*cc*(xint-x0)*(x1-xint)/PEL::pow( h, 3 ) ;
      result( 0, 1 ) =  ux*6.*cc*alpha*(2.*y-1.)*(2.*t+PEL::sin(beta*t)/beta)
                         *(xint-x0)*(x1-xint)/PEL::pow( h, 3 )
                       + cc*((xint-x0)*(xint-x0)
                         *(-2.*xint+3.*x1-x0)/PEL::pow(h,3)+1.)
                         *alpha*(1.-2.*y)*(2.+ PEL::cos(beta*t)) ;
      result( 1, 0 ) = 0.0 ;
      result( 1, 1 ) = 0.0 ;
   }
   else if( EXPR == RHOD )
   {
      double x = res( 0 ) ;
      double x0 = arg(5)->to_double(ct) ;
      double x1 = arg(6)->to_double(ct) ;
      double h = x1 - x0 ;
      double cc = arg(7)->to_double(ct) ;
      double xint = x-alpha*y*(1.-y)*(2.*t+ PEL::sin(beta*t)/beta);
      result( 0, 0 ) = 6.* cc *(xint-x0)*(x1-xint)/PEL::pow( h, 3 ) ; 
      result( 0, 1 ) = 6.* cc * alpha*(2.*y-1.)*(2.*t+PEL::sin(beta*t)/beta)*
	            (xint-x0)*(x1-xint)/PEL::pow( h, 3 )  ;
   }

   return result ;
}


//----------------------------------------------------------------------
std::string const& 
RS_VariableDensityFlow1EXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "VariableDensityFlow1_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"; 
         break ;
      case F : result =  "VariableDensityFlow1_rhsf($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)";
	 break ;
      case P : result = "VariableDensityFlow1_pressure($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"; 
	 break ;
      case RHO : result = "VariableDensityFlow1_rho($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L) ";
	 break ;
      case RHOD : result = "VariableDensityFlow1_grad_rho($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"; 
	 break ;
      case RHODT : result = "VariableDensityFlow1_rhsg($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)";
	 break ;
      case RHOU : result = "VariableDensityFlow1_rho_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"; 
	 break ;
      case UD : result = "VariableDensityFlow1_grad_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"  ;
         break ;
      case RHOUD : result = "VariableDensityFlow1_grad_rho_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"  ; 
         break ;
   }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_VariableDensityFlow1EXP:: valid_arguments( 
                                   PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_VariableDensityFlow1EXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =  
 ( some_arguments->count() == 9 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 8 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
