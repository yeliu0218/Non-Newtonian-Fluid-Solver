#include <AP_FSIsolutionEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>

#include <iostream>

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_UU = 
             new AP_FSIsolutionEXP( "FSIsolution_uu", UU ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_PP = 
             new AP_FSIsolutionEXP( "FSIsolution_pp", PP ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_MU = 
             new AP_FSIsolutionEXP( "FSIsolution_mu", MU ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_DD = 
             new AP_FSIsolutionEXP( "FSIsolution_dd", DD ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_SF = 
             new AP_FSIsolutionEXP( "FSIsolution_sf", SF ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_SS = 
             new AP_FSIsolutionEXP( "FSIsolution_ss", SS ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_UD = 
             new AP_FSIsolutionEXP( "FSIsolution_grad_u", UD ) ;

AP_FSIsolutionEXP const* 
AP_FSIsolutionEXP::PROTOTYPE_SD = 
             new AP_FSIsolutionEXP( "FSIsolution_grad_s", SD ) ;

//----------------------------------------------------------------------
AP_FSIsolutionEXP:: AP_FSIsolutionEXP( std::string const& a_name,
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
AP_FSIsolutionEXP*
AP_FSIsolutionEXP:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   AP_FSIsolutionEXP* result = new AP_FSIsolutionEXP( a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_FSIsolutionEXP:: AP_FSIsolutionEXP( PEL_Object* a_owner,
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
   PEL_LABEL( "AP_FSIsolutionEXP:: AP_FSIsolutionEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
AP_FSIsolutionEXP:: ~AP_FSIsolutionEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: ~AP_FSIsolutionEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
AP_FSIsolutionEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case UU : result = PEL_Data::DoubleVector ; 
         break ;
      case PP : result = PEL_Data::DoubleVector ; 
         break ;
      case MU : result = PEL_Data::DoubleVector ; 
         break ;
      case DD : result = PEL_Data::DoubleVector ; 
         break ;
      case SF : result = PEL_Data::DoubleVector ; 
         break ;
      case SS : result = PEL_Data::DoubleVector ; 
         break ;
      case UD : result = PEL_Data::DoubleArray2D ; 
         break ;
      case SD : result = PEL_Data::DoubleArray2D ;
         break ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
AP_FSIsolutionEXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( ( EXPR == PP || EXPR == MU ) ? DV_result_1 
                                                         : DV_result_2 ) ;

   double alpha = PEL::pi()/4. ;

   doubleVector const& xy = arg(0)->to_double_vector( ct ) ;
   double x = xy( 0 ) ;
   double y = xy( 1 ) ;
   double ri = arg(1)->to_double( ct ) ;
   double re = arg(2)->to_double( ct ) ;

   double r = PEL::sqrt( x*x + y*y ) ;
   double rr = r * r ;

   double sinx = y/r ;
   double cosy = x/r ;
   double sinA = PEL::sin( alpha ) ;
   double cosA = PEL::cos( alpha );
   
   if( EXPR == UU )
   {
      double arg1 = r - ri ;
      double arg2 = re - r ;
      double arg3 = re - ri ;

      double vit_ex = PEL::tanh( arg1 )+PEL::tanh( arg2 )-PEL::tanh( arg3 ) ;

      result( 0 ) = -sinx * vit_ex * sinA ;
      result( 1 ) =  cosy * vit_ex * sinA ;
   }
   if( EXPR == PP )
   {
      result( 0 ) = -4.0 * sinA * sinA * ( 1.0 - 4.0 * sinA * sinA ) ;
   }
   else if( EXPR == MU )
   {
      double tan2 = PEL::tanh( 2.0 ) * PEL::tanh( 2.0 ) ;
    
      result( 0 ) = 2.0 * cosA * ( 8.0 * sinA * sinA + 1.0 )/tan2 ; 
   }
   else if( EXPR == DD )
   {
      result( 0 ) = -2.0 * sinA * re * re * ( sinA * x + cosA * y ) / rr ;
      result( 1 ) = -2.0 * sinA * re * re * ( sinA * y - cosA * x ) / rr ;
   }
   else if( EXPR == SF )
   {
      double arg1 = r - ri ;
      double arg2 = re - r ;
      double arg3 = re - ri ;

      double vit_ex = PEL::tanh( arg1 )+PEL::tanh( arg2 )-PEL::tanh( arg3 ) ;
      double tan2 = PEL::tanh( 2.0 ) * PEL::tanh( 2.0 ) ;
      double mu =  2.0 * cosA * ( 8.0 * sinA * sinA + 1.0 )/tan2 ;
      double theta_0 = PEL::tanh( arg2 )*PEL::tanh( arg2 ) 
                     - PEL::tanh( arg1 )*PEL::tanh( arg1 ) ;
      double theta_1 = PEL::tanh( arg2 )
                     * ( 1.0 - PEL::tanh( arg2 ) * PEL::tanh( arg2 ) ) ;
      double theta_2 = PEL::tanh( arg1 )
                     * ( 1.0 - PEL::tanh( arg1 ) * PEL::tanh( arg1 ) ) ;

      double f_r = -sinA * vit_ex/r ;
      double f_theta = -mu * sinA * ( theta_0/r - 2.0 *( theta_1 + theta_2 ) 
                                      - vit_ex/rr ) ;

      result( 0 ) = cosy * f_r - sinx * f_theta ;
      result( 1 ) = sinx * f_r + cosy * f_theta ;
   }
   else if( EXPR == SS )
   { 
      double r3 = rr * r ;
      double r4 = rr * rr ;
      double r5 = rr * r3 ;
      
      double sinA2 = sinA * sinA ;
      double cosA2 = cosA * cosA ;
      double sinA3 = sinA2 * sinA ;
      double sinA4 = sinA2 * sinA2 ;

      double s_r = 8.0 * sinA2 
                 *( 1.0 + 36.0 * sinA2/rr - 64.0 * cosA2/rr + 160.0 * sinA4/r4
                  - 384.0 * sinA2 * cosA2/r4 + 16.0 * sinA2 * cosA2/rr )/r3 ;
      double s_theta = -128.0 * sinA3 * cosA * ( 1.0 + sinA2 )/r5 ;

      result( 0 ) = cosy * s_r - sinx * s_theta ;
      result( 1 ) = sinx * s_r + cosy * s_theta ;
   }

   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
AP_FSIsolutionEXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ; 

   doubleArray2D& result = doubleArray2D_result ;

   double alpha = PEL::pi()/4. ;

   doubleVector const& xy = arg(0)->to_double_vector( ct ) ;
   double x = xy( 0 ) ;
   double y = xy( 1 ) ;
   double ri = arg(1)->to_double( ct ) ;
   double re = arg(2)->to_double( ct ) ;

   double r = PEL::sqrt( x*x + y*y ) ;
   double rr = r * r ;

   double sinA = PEL::sin( alpha ) ;
   double cosA = PEL::cos( alpha );

   if( EXPR == UD )
   {
      double arg1 = r - ri ;
      double arg2 = re - r ;
      double arg3 = re - ri ;

      double vit_ex_1 = PEL::tanh( arg1 )+PEL::tanh( arg2 )-PEL::tanh( arg3 ) ;
      double vit_ex_2 = PEL::tanh( arg2 )*PEL::tanh( arg2 )
                      - PEL::tanh( arg1 )*PEL::tanh( arg1 ) ;

      result( 0, 0 ) =  x * y * sinA * ( vit_ex_1/r - vit_ex_2 )/rr ;
      result( 0, 1 ) = -sinA * ( x * x * vit_ex_1/r + y * y * vit_ex_2 )/rr ;
      result( 1, 0 ) =  sinA * ( y * y * vit_ex_1/r + x * x * vit_ex_2 )/rr ;
      result( 1, 1 ) = -x * y * sinA * ( vit_ex_1/r - vit_ex_2 )/rr ;
   }
   else if( EXPR == SD )
   {
      result( 0, 0 ) = -2.*sinA*re*re*( sinA*(y*y-x*x) - 2.*x*y*cosA )/rr/rr ;
      result( 0, 1 ) = -2.*sinA*re*re*( cosA*(x*x-y*y) - 2.*x*y*sinA )/rr/rr ;
      result( 1, 0 ) = -2.*sinA*re*re*( cosA*(x*x-y*y) - 2.*x*y*sinA )/rr/rr ;
      result( 1, 1 ) = -2.*sinA*re*re*( sinA*(x*x-y*y) + 2.*x*y*cosA )/rr/rr ;
   }

   return result ;
}

//----------------------------------------------------------------------
std::string const& 
AP_FSIsolutionEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case UU : result = "FSIsolution_uu(<doubleVector>,<double>,<double>)" ; 
         break ;
      case PP : result = "FSIsolution_pp(<doubleVector>,<double>,<double>)" ; 
         break ;
      case MU : result = "FSIsolution_mu(<doubleVector>,<double>,<double>)" ; 
         break ;
      case DD : result = "FSIsolution_dd(<doubleVector>,<double>,<double>)" ; 
         break ;
      case SF : result = "FSIsolution_sf(<doubleVector>,<double>,<double>)" ; 
         break ;
      case SS : result = "FSIsolution_ss(<doubleVector>,<double>,<double>)" ; 
         break ;
      case UD : result = "FSIsolution_grad_u(<doubleVector>,<double>,<double>)" ; 
         break ;
      case SD : result = "FSIsolution_grad_s(<doubleVector>,<double>,<double>)" ; 
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
AP_FSIsolutionEXP:: valid_arguments( 
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSIsolutionEXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   result =  ( some_arguments->count() == 3 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) ; 

   return result ;
}
