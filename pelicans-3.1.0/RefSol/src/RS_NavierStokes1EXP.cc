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

#include <RS_NavierStokes1EXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_NavierStokes1EXP const* 
RS_NavierStokes1EXP::PROTOTYPE_U = 
             new RS_NavierStokes1EXP( "NavierStokes1_velocity", U ) ;

RS_NavierStokes1EXP const* 
RS_NavierStokes1EXP::PROTOTYPE_P = 
             new RS_NavierStokes1EXP( "NavierStokes1_pressure", P ) ;

RS_NavierStokes1EXP const* 
RS_NavierStokes1EXP::PROTOTYPE_F = 
             new RS_NavierStokes1EXP( "NavierStokes1_force", F ) ;

RS_NavierStokes1EXP const* 
RS_NavierStokes1EXP::PROTOTYPE_UD = 
             new RS_NavierStokes1EXP( "NavierStokes1_grad_velocity", UD ) ;

//----------------------------------------------------------------------
RS_NavierStokes1EXP:: RS_NavierStokes1EXP( std::string const& a_name,
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
RS_NavierStokes1EXP*
RS_NavierStokes1EXP:: create_replica( 
                             PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_NavierStokes1EXP* result = new RS_NavierStokes1EXP( a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_NavierStokes1EXP:: RS_NavierStokes1EXP( PEL_Object* a_owner,
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
   PEL_LABEL( "RS_NavierStokes1EXP:: RS_NavierStokes1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_NavierStokes1EXP:: ~RS_NavierStokes1EXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: ~RS_NavierStokes1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_NavierStokes1EXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = DoubleVector ; 
         break ;
      case P : result = DoubleVector ;
	 break ;
      case F : result = DoubleVector ;
	 break ;
      case UD : result = DoubleArray2D ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_NavierStokes1EXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( EXPR == P ? DV_result_1 : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
   double dens = arg(1)->to_double(ct) ;
   double visc = arg(2)->to_double(ct) ;
   double x = res( 0 ) ;
   double x2 = x * x ;
   double x3 = x * x2 ;
   double x4 = x * x3 ;
   double y = res( 1 ) ;
   double y2 = y * y ;
   double y3 = y * y2 ;
   double y4 = y * y3 ;

   double ux =  1000.0 * ( x2 - 2. * x3 + x4 )
                       * ( 2. * y - 6. * y2 + 4. * y3 ) ;
   double uy = - 1000.0 * ( y2 - 2. * y3 + y4 )
                        * ( 2. * x - 6. * x2 + 4. * x3 ) ;
   if( EXPR == P )
   {
      result( 0 ) =  100.0 * ( x2 + y2 - 2.0/3.0 ) ;
   }
   else if( EXPR == U )
   { 
      result( 0 ) = ux ;
      result( 1 ) = uy ;
   }
   else if( EXPR == F )
   {
      double dux_dx = 1000.0 * ( 2. * x - 6. * x2 + 4. * x3 )
                             * ( 2. * y - 6. * y2 + 4. * y3 ) ;
      double dux_dy = 1000.0 * ( x2 - 2. * x3 + x4 )
                             * ( 2. - 12. * y + 12. * y2 ) ;
      double duy_dx = - 1000.0 * ( y2 - 2. * y3 + y4 )
                               * ( 2. - 12. * x + 12. * x2 ) ;
      double duy_dy = - 1000.0 * ( 2. * y - 6. * y2 + 4. * y3 )
                               * ( 2. * x - 6. * x2 + 4. * x3 ) ;

      double d2ux_dx2 = 1000.0 * ( 2. - 12. * x + 12. * x2 )
                               * ( 2. * y - 6. * y2 + 4. * y3 ) ;
      double d2ux_dy2 = 1000.0 * ( x2 - 2. * x3 + x4 )
                               * ( -12. + 24. * y ) ;
      double d2uy_dy2 = - 1000.0 * ( 2. - 12. * y + 12. * y2 )
                                 * ( 2. * x - 6. * x2 + 4. * x3 ) ;
      double d2uy_dx2 = - 1000.0 * ( y2 - 2. * y3 + y4 )
                                 * ( -12. + 24. * x ) ;
      double dp_dx = 200.0 * x ;
      double dp_dy = 200.0 * y ;
      result( 0 ) = dens*( ux*dux_dx + uy*dux_dy )
                  - visc*( d2ux_dx2 + d2ux_dy2 ) 
                  + dp_dx ;
      result( 1 ) = dens*( ux*duy_dx + uy*duy_dy )
                  - visc*( d2uy_dx2 + d2uy_dy2 ) 
                  + dp_dy ;   
   }

   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_NavierStokes1EXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ; 

   doubleArray2D& result = doubleArray2D_result ;

   if( EXPR == UD )
   {
      doubleVector const& res = arg(0)->to_double_vector(ct) ;
      double x = res( 0 ) ;
      double x2 = x * x ;
      double x3 = x * x2 ;
      double x4 = x * x3 ;
      double y = res( 1 ) ;
      double y2 = y * y ;
      double y3 = y * y2 ;
      double y4 = y * y3 ;

      double dux_dx = 1000.0 * ( 2. * x - 6. * x2 + 4. * x3 )
                             * ( 2. * y - 6. * y2 + 4. * y3 ) ;
      double dux_dy = 1000.0 * ( x2 - 2. * x3 + x4 )
                             * ( 2. - 12. * y + 12. * y2 ) ;
      double duy_dx = - 1000.0 * ( y2 - 2. * y3 + y4 )
                               * ( 2. - 12. * x + 12. * x2 ) ;
      double duy_dy = - 1000.0 * ( 2. * y - 6. * y2 + 4. * y3 )
                               * ( 2. * x - 6. * x2 + 4. * x3 ) ;

      result( 0, 0 ) = dux_dx ;
      result( 0, 1 ) = dux_dy ;
      result( 1, 0 ) = duy_dx ;
      result( 1, 1 ) = duy_dy ;
   }

   return result ;
}

//----------------------------------------------------------------------
std::string const& 
RS_NavierStokes1EXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "NavierStokes1_velocity($DV_X,$DS_ALPHA,$DS_MU)" ; 
         break ;
      case P : result = "NavierStokes1_pressure($DV_X,$DS_ALPHA,$DS_MU)" ;
	 break ;
      case F : result = "NavierStokes1_force($DV_X,$DS_ALPHA,$DS_MU)" ;
	 break ;
      case UD : result = "NavierStokes1_grad_velocity($DV_X,$DS_ALPHA,$DS_MU)" ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
RS_NavierStokes1EXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_NavierStokes1EXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = ( some_arguments->count() == 3 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )       &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
