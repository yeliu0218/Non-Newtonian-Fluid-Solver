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

#include <RS_AvulaEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>
using std::cout ;
using std::endl ; 

RS_AvulaEXP const* 
RS_AvulaEXP::PROTOTYPE_U_AXI = 
             new RS_AvulaEXP( "Avula_velocity_axi", U_AXI ) ;
RS_AvulaEXP const* 
RS_AvulaEXP::PROTOTYPE_P_AXI = 
             new RS_AvulaEXP( "Avula_pressure_axi", P_AXI ) ;
RS_AvulaEXP const* 
RS_AvulaEXP::PROTOTYPE_U_3D = 
             new RS_AvulaEXP( "Avula_velocity_3D", U_3D ) ;
RS_AvulaEXP const* 
RS_AvulaEXP::PROTOTYPE_P_3D = 
             new RS_AvulaEXP( "Avula_pressure_3D", P_3D ) ;

//--------------------------------------------------------------------------
RS_AvulaEXP:: RS_AvulaEXP( std::string const& a_name,
                           Func an_expr) 
//--------------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR(an_expr)
   , DV_result_1( 1 )
   , DV_result_2( 2 )
   , DV_result_3( 3 )
{
   PEL_LABEL( "RS_AvulaEXP:: RS_AvulaEXP" ) ;
}

//--------------------------------------------------------------------------
RS_AvulaEXP*
RS_AvulaEXP:: create_replica( PEL_Object* a_owner,
                              PEL_Sequence const* argument_list ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_AvulaEXP* result = new RS_AvulaEXP( a_owner, 
                                          name(), 
                                          argument_list, 
                                          EXPR ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
RS_AvulaEXP:: RS_AvulaEXP( PEL_Object* a_owner,
                           std::string const& a_name,
                           PEL_Sequence const* argument_list,
                           Func an_expr ) 
//--------------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 )
   , DV_result_3( 3 )
{
   PEL_LABEL( "RS_AvulaEXP:: RS_AvulaEXP" ) ;

//   compute_more_J0_roots( 2000 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//--------------------------------------------------------------------------
RS_AvulaEXP:: ~RS_AvulaEXP( void ) 
//--------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: ~RS_AvulaEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;   
}

//--------------------------------------------------------------------------
PEL_Data::Type
RS_AvulaEXP:: data_type( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: data_type" ) ;
   PEL_Data::Type my_type = PEL_Data::Undefined ;
   switch( EXPR )
   {
      case U_AXI : my_type = DoubleVector ; 
         break ;
      case P_AXI : my_type = DoubleVector ; 
         break ;
      case U_3D : my_type = DoubleVector ;
         break ;
      case P_3D : my_type = DoubleVector ; 
         break ;
      default :
	 PEL_Error::object()->raise_plain( "Choose:SolAvula_velocity_axi, SolAvula_pressure_axi, SolAvula_velocity_3D or SolAvula_pressure_3D " ) ;
   }
   return my_type ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_AvulaEXP:: to_double_vector( PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   doubleVector const& xx  = arg( 0 )->to_double_vector( ct ) ; 
   double t                = arg( 1 )->to_double( ct ) ;
   double radius           = arg( 2 )->to_double( ct ) ; 
   double length           = arg( 3 )->to_double( ct ) ;
   double p1               = arg( 4 )->to_double( ct ) ;
   double p0               = arg( 5 )->to_double( ct ) ; 
   double period           = arg( 6 )->to_double( ct ) ; 
   double t_cut            = arg( 7 )->to_double( ct ) ;
   double rho              = arg( 8 )->to_double( ct ) ; 
   double mu               = arg( 9 )->to_double( ct ) ;  
   double tol              = arg( 10 )->to_double( ct ) ;

   double tau = 2.0 * PEL::pi() * t / period ;
   double tau_cut = 2.0 * PEL::pi() * t_cut / period ;
   double u_ref = 2.0 * PEL::pi() * radius / period ;
   double reynolds = rho * u_ref * radius / mu ;

   double yyy = rho * u_ref * u_ref * length / radius ;
   double a = p1 / yyy ;
   double c = p0 / yyy ;

   if( EXPR==P_AXI || EXPR==P_3D )
   {
      doubleVector& result = DV_result_1 ;
      double p_in = 0.0 ;
      double z = ( EXPR==P_AXI ? xx( 1 ) : xx( 2 ) ) ;
      z = xx( 1 ) ;
      if( tau <= tau_cut ) 
      {
         p_in = p1 * PEL::sin( tau ) + p0 ;
      }
      else 
      {
         p_in = p1 * PEL::sin( tau_cut ) + p0 ;
      }
      result( 0 ) = p_in * ( length - z ) / length  ;
      return result ;
   }
   else if( EXPR == U_AXI )
   { 
      doubleVector& result = DV_result_2 ;
      double eta = xx( 0 ) / radius ;
      result ( 0 ) = 0 ;
      result ( 1 ) = u_ref * uz_adim( eta, tau, 
                                      tau_cut, a, c, reynolds, tol ) ;
      return result ;
   }
   else if( EXPR == U_3D )
   { 
      doubleVector& result = DV_result_3 ;
      double eta = PEL::sqrt( xx( 0 )*xx( 0 ) + xx( 1 )*xx( 1 ) ) / radius ;
      result ( 0 ) = 0 ;
      result ( 1 ) = 0 ;
      result ( 2 ) = u_ref * uz_adim( eta, tau, 
                                      tau_cut, a, c, reynolds, tol ) ;
      return result ;
   }
}

//----------------------------------------------------------------------
std::string const& 
RS_AvulaEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: usage" ) ;

   static std::string result ;
   switch(EXPR)
   {
      case U_AXI : 
	 result = "Avula_velocity_axi($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)" ; 
         break ; 
      case U_3D : 
	 result = "Avula_velocity_3D($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)" ; 
         break ;
      case P_AXI : 
	 result = "Avula_pressure_axi($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)" ; 
         break ;
      case P_3D : 
	 result = "Avula_pressure_3D($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)" ; 
         break ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_AvulaEXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: valid_arguments" ) ;
   bool result = ( some_arguments->count() == 11 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) && 
 ( extract_arg( some_arguments, 8 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 9 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 10 )->data_type() == PEL_Data::Double ) ;

   return result ; 
}

//----------------------------------------------------------------------------
double
RS_AvulaEXP:: uz_adim( double eta, double tau, double tau_cut, 
                       double a, double c, double reynolds, double tol ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: uz_adim" ) ;

   double f_tau_cut, alpha_n, gamma_n, h_n, ee, du ;
   double result = 0.0 ;

   if( tau > tau_cut ) f_tau_cut = a * PEL::sin( tau_cut ) + c ;

   size_t n = 0 ;
   do
   {
      alpha_n = root_of_J0( n ) ;
      gamma_n = alpha_n * alpha_n / reynolds ;
      if( tau <= tau_cut )
      {
         h_n = hconvol( tau, a, c, gamma_n ) ;
      }
      else
      {
         h_n = hconvol( tau_cut, a, c, gamma_n ) ;
         ee = PEL::exp( -gamma_n * ( tau - tau_cut ) ) ;
         h_n = ( 1.0 - ee ) / gamma_n * f_tau_cut + ee * h_n ;
      }
      du = PEL::j0( alpha_n * eta ) / alpha_n / PEL::j1( alpha_n ) * h_n ;
      result += du ;
      n++ ;
   }
   while( PEL::abs(du) > tol ) ;

   result *= 2.0 ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
RS_AvulaEXP:: hconvol( double tau, double a, double c, double gamma_n ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: hconvol" ) ;

   double ee = PEL::exp( - gamma_n * tau ) ;
   double result = a/(gamma_n*gamma_n+1.0) * 
                   ( gamma_n*PEL::sin(tau) - PEL::cos(tau) + ee ) 
                 + c/gamma_n* (1.0 - ee ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
RS_AvulaEXP:: root_of_J0( size_t n ) const 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "RS_AvulaEXP:: root_of_J0" ) ;

   static doubleVector J0_ROOTS( 0 ) ;

   if( n >= J0_ROOTS.size() )
   {
      PEL_ASSERT( n == J0_ROOTS.size() ) ;
      double x = ( J0_ROOTS.size()==0 ?  
                      1.0 : 
                   J0_ROOTS( J0_ROOTS.size()-1 ) + 3.0 ) ;
      double tolerance = 1.e-10 ;
      double f, df, dx ;
      for( size_t i=0 ; i<100 ; ++i )
      {
         do
         {
            f  =  PEL::j0( x ) ;
            df = -PEL::j1( x ) ;
            PEL_CHECK( PEL::abs( df ) > tolerance ) ;
            dx = - f/df ;
            x += dx ;
         }
         while( PEL::abs( dx/x) > tolerance ) ;
         J0_ROOTS.append( x ) ;
         x += 3.0 ;
      }
   }
   double result = J0_ROOTS( n ) ;
   return( result ) ;
}
