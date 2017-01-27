#include <RS_PotentialField.hh>

#include <RS_Bingham.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <math.h>
#include <iostream>

RS_PotentialField const* 
RS_PotentialField::PROTOTYPE_ONE 
    = new RS_PotentialField( "ConstantPhi_ON_OFF", ONE ) ;

//RS_PotentialField const* 
//RS_PotentialField::PROTOTYPE_TWO 
//    = new RS_PotentialField( "Step_Function_inner", TWO ) ;
    
//----------------------------------------------------------------------
RS_PotentialField:: RS_PotentialField( 
                               std::string const& a_name, Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_PotentialField" ) ;
}

//----------------------------------------------------------------------
RS_PotentialField*
RS_PotentialField:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_PotentialField:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_PotentialField* result =
               new RS_PotentialField( a_owner, name(), argument_list, EXPR ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_PotentialField:: RS_PotentialField(
         PEL_Object* a_owner, std::string const& a_name,
		 PEL_Sequence const* argument_list, Func an_expr  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_PotentialField:: RS_PotentialField" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_PotentialField:: ~RS_PotentialField( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_PotentialField:: ~RS_PotentialField" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_PotentialField:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}
//----------------------------------------------------------------------
double
RS_PotentialField:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_PotentialField:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double result = 0.0;
   switch( EXPR )
   {
      case ONE:
	  {
         double const time = arg(0)->to_double(ct) ;
         double const Tinterval_ON  = arg(1)->to_double(ct) ;
         double const Tinterval_OFF = arg(2)->to_double(ct) ;
         double const valuePhi = arg(3)->to_double(ct) ;
         
         double const Tinterval_total_length = Tinterval_ON + Tinterval_OFF ; 
         
         int N; // number of periods
         
         N = int( time / Tinterval_total_length );
         
         double tmp = N * Tinterval_total_length + Tinterval_ON ;
         
         if ( time < tmp )
            result = valuePhi ;
         
	     break;
	  }
      case TWO:
	  {
	      break;
      }
   }
	
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
RS_PotentialField:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_PotentialField:: usage" ) ;

   static std::string result ;

   switch( EXPR )
   {
      case ONE :
         result = "ConstantPhi_ON_OFF($DS_T, $DS_Tinterval_ON, $DS_Tinterval_OFF, $DS_const_valuePhi)" ;
         break ;
      case TWO :
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_PotentialField:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_PotentialField:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( EXPR )
   {
      case ONE :
	   result = ( some_arguments->count() == 4 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) ;
      break ;
      case TWO :
          
      break ;
   }

   return( result ) ;
}
