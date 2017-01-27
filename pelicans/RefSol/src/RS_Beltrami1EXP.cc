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

#include <RS_Beltrami1EXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_Beltrami1EXP const* 
RS_Beltrami1EXP::PROTOTYPE_U = 
             new RS_Beltrami1EXP( "Beltrami1_velocity", U ) ;

RS_Beltrami1EXP const* 
RS_Beltrami1EXP::PROTOTYPE_P = 
             new RS_Beltrami1EXP( "Beltrami1_pressure", P ) ;

RS_Beltrami1EXP const*
RS_Beltrami1EXP::PROTOTYPE_UD =
             new RS_Beltrami1EXP( "Beltrami1_grad_velocity", UD ) ;


//----------------------------------------------------------------------
RS_Beltrami1EXP:: RS_Beltrami1EXP( std::string const& a_name, Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
   , DV_result_1( 1 ) 
   , DV_result_2( 3 )
   , doubleArray2D_result( 3, 3 )
{
}

//----------------------------------------------------------------------
RS_Beltrami1EXP*
RS_Beltrami1EXP:: create_replica( PEL_Object* a_owner,
                                  PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Beltrami1EXP* result = new RS_Beltrami1EXP( a_owner, 
                                                  name(), 
                                                  argument_list, 
                                                  EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Beltrami1EXP:: RS_Beltrami1EXP( PEL_Object* a_owner,
				   std::string const& a_name,
				   PEL_Sequence const* argument_list,
				   Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
   , DV_result_1( 1 ) 
   , DV_result_2( 3 )   
   , doubleArray2D_result( 3, 3 )
{
   PEL_LABEL( "RS_Beltrami1EXP:: RS_Beltrami1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Beltrami1EXP:: ~RS_Beltrami1EXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: ~RS_Beltrami1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Beltrami1EXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = DoubleVector ; 
         break ;
      case P : result = DoubleVector ;
	 break ;  
      case UD : result = DoubleArray2D ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_Beltrami1EXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( EXPR == P ? DV_result_1 : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector( ct ) ;
   double t = arg(1)->to_double( ct ) ;
   doubleVector const& wave = arg(2)->to_double_vector( ct ) ;
   double a = arg(3)->to_double( ct ) ;
   double nu = arg(4)->to_double( ct ) ;

   double kx = res( 0 ) * wave ( 0 ) ;
   double ly = res( 1 ) * wave ( 1 ) ;
   double mz = res( 2 ) * wave ( 2 ) ;
   double mk = wave( 2 ) * wave( 0 ) ;
   double ml = wave( 2 ) * wave( 1 ) ;

   double ll = wave( 0 ) * wave( 0 ) + wave( 1 ) * wave( 1 ) 
             + wave( 2 ) * wave( 2 ) ;
   double lambda = PEL::sqrt( ll ) ;
   double laml = lambda * wave( 1 ) ;
   double lamk = lambda * wave( 0 ) ;

   double tt = -nu * ll * t ;
   double expt = PEL::exp( tt ) ;

   double b = wave( 0 ) * wave( 0 ) + wave( 1 ) * wave( 1 ) ;
   double coef = a / b ;

   double ckx = PEL::cos( kx ) ;
   double skx = PEL::sin( kx ) ;
   double cly = PEL::cos( ly ) ;
   double sly = PEL::sin( ly ) ;
   double cmz = PEL::cos( mz ) ;
   double smz = PEL::sin( mz ) ;

   double uu = -coef * ( laml * ckx * sly * smz + mk * skx * cly * cmz ) 
                     * expt ;
   double vv = coef * ( lamk * skx * cly * smz - ml * ckx * sly * cmz ) 
                    * expt ;
   double ww = a * ckx * cly * smz * expt ;

   if( EXPR == U )
   {
      result( 0 ) = uu ;
      result( 1 ) = vv ;
      result( 2 ) = ww ;
   }
   else if( EXPR == P )
   {
      result( 0 ) = -0.5 * ( uu * uu + vv * vv + ww * ww )  ; 
   }

   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_Beltrami1EXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ; 

   doubleArray2D& result = doubleArray2D_result ;

   PEL_ASSERT( EXPR == UD ) ;

   doubleVector const& res = arg(0)->to_double_vector( ct ) ;
   double t = arg(1)->to_double( ct ) ;
   doubleVector const& wave = arg(2)->to_double_vector( ct ) ;
   double a = arg(3)->to_double( ct ) ;
   double nu = arg(4)->to_double( ct ) ;

   double k = wave( 0 ) ;
   double l = wave( 1 ) ;
   double m = wave( 2 ) ;
   double kx = res( 0 ) * k ;
   double ly = res( 1 ) * l ;
   double mz = res( 2 ) * m ;
   double mk = m * k ;
   double ml = m * l ;

   double ll = wave( 0 ) * wave( 0 ) + wave( 1 ) * wave( 1 ) 
      + wave( 2 ) * wave( 2 ) ;
   double lambda = PEL::sqrt( ll ) ;
   double laml = lambda * wave( 1 ) ;
   double lamk = lambda * wave( 0 ) ;

   double tt = -nu * ll * t ;
   double expt = PEL::exp( tt ) ;
 
   double b = k * k + l * l ;
   double coef = a / b ;

   double ckx = PEL::cos( kx ) ;
   double skx = PEL::sin( kx ) ;
   double cly = PEL::cos( ly ) ;
   double sly = PEL::sin( ly ) ;
   double cmz = PEL::cos( mz ) ;
   double smz = PEL::sin( mz ) ;


   result( 0, 0 ) = -coef * ( -laml * k * skx * sly * smz 
			      + mk * k * ckx * cly * cmz ) * expt ;
   result( 0, 1 ) = -coef * ( laml * l * ckx * cly * smz 
			      - mk * l * skx * sly * cmz ) * expt ;
   result( 0, 2 ) = -coef * ( laml * m * ckx * sly * cmz 
			      - mk * m * skx * cly * smz ) * expt ;
   result( 1, 0 ) = coef * ( lamk * k * ckx * cly * smz 
			     + ml * k * skx * sly * cmz ) * expt ;
   result( 1, 1 ) = coef * ( -lamk * l * skx * sly * smz 
			     - ml * l * ckx * cly * cmz ) * expt ;
   result( 1, 2 ) = coef * ( lamk * m * skx * cly * cmz
			     + ml * m * ckx * sly * smz ) * expt ;
   result( 2, 0 ) = -a * k * skx * cly * smz * expt ;
   result( 2, 1 ) = -a * l * ckx * sly * smz * expt ;
   result( 2, 2 ) = a * m * ckx * cly * cmz * expt ;
 
   return result ;
}

//----------------------------------------------------------------------
std::string const& 
RS_Beltrami1EXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {        
      case U : result = 
         "Beltrami1_velocity($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)" ;
         break ;
      case P : result = 
         "Beltrami1_pressure($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)" ;
	 break ;
      case UD : result = 
         "Beltrami1_grad_velocity($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)" ; 
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
RS_Beltrami1EXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Beltrami1EXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = ( some_arguments->count() == 5 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )  &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )  &&
 ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
