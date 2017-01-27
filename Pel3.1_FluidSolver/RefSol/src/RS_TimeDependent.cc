/*
 *  Copyright : 
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#include <RS_TimeDependent.hh>

#include <RS_Bingham.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <math.h>
#include <iostream>

RS_TimeDependent const* 
RS_TimeDependent::PROTOTYPE_ONE 
    = new RS_TimeDependent( "TimeDependentProfile", ONE ) ;

RS_TimeDependent const* 
RS_TimeDependent::PROTOTYPE_TWO 
    = new RS_TimeDependent( "Step_Function_inner", TWO ) ;
    
RS_TimeDependent const* 
	RS_TimeDependent::PROTOTYPE_THREE 
	= new RS_TimeDependent( "Step_Function_outer", THREE ) ;

RS_TimeDependent const*
RS_TimeDependent::PROTOTYPE_FOUR
    = new RS_TimeDependent( "StartAtTimePoint", FOUR ) ;

//----------------------------------------------------------------------
RS_TimeDependent:: RS_TimeDependent( 
                               std::string const& a_name, Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_TimeDependent" ) ;
}

//----------------------------------------------------------------------
RS_TimeDependent*
RS_TimeDependent:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TimeDependent:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_TimeDependent* result =
               new RS_TimeDependent( a_owner, name(), argument_list, EXPR ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_TimeDependent:: RS_TimeDependent(
                 PEL_Object* a_owner, std::string const& a_name,
		 PEL_Sequence const* argument_list, Func an_expr  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_TimeDependent:: RS_TimeDependent" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_TimeDependent:: ~RS_TimeDependent( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TimeDependent:: ~RS_TimeDependent" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_TimeDependent:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}
//----------------------------------------------------------------------
double
RS_TimeDependent:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TimeDependent:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double result = 0.0;
   switch( EXPR )
   {
	case ONE:
	{
	   double const U = arg(0)->to_double(ct) ;
	   double const T = arg(1)->to_double(ct) ;
	   double const A = arg(2)->to_double(ct) ;
	   double const B = arg(3)->to_double(ct) ;

	   result = U + A *PEL::sin(B*T);

	   break;
	}
        case TWO:
	{

           double const T1 = arg(0)->to_double(ct) ;
           double const DT1 = arg(1)->to_double(ct) ;
           double const U1 = arg(2)->to_double(ct) ;
           if ((int(T1)/int(DT1))%90==0) 
	   {
              result = U1;
	   }
           else 
	   {
              result = 0.0;
	   }	

	   break;
	}
	case THREE :
	{
           double const T2 = arg(0)->to_double(ct) ;
           double const DT2 = arg(1)->to_double(ct) ;
           double const U2 = arg(2)->to_double(ct) ;
           if ((int(T2)/int(DT2))%90==0) 
	   {
              result = 0.0;
	   }
           else 
	   {
              result  = U2;
	   }

  	   break;
	}
	case FOUR:
	{
	   double const T = arg(0)->to_double(ct) ;
	   double const pt = arg(1)->to_double(ct) ;
	   if (T>pt)
	      result = 1.0;

	   break;
	}


   }
	

   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
RS_TimeDependent:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TimeDependent:: usage" ) ;

   static std::string result ;

   switch( EXPR )
   {
      case ONE :
         result = "TimeDependentProfile($DS_T, $DS_A, $DS_B, $DS_C)" ;
         break ;
      case TWO :
         result = "Step_Function_inner($DS_T1, $DS_DT1, $DS_U1)" ;
         break ;
      case THREE:
	 result = "Step_Function_outer($DS_T2, $DS_DT2, $DS_U2)" ;
	 break ;
      case FOUR :
	 result = "StartAtTimePoint($DS_T, $DS_StartPoint)" ;
         break ;	   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_TimeDependent:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TimeDependent:: valid_arguments" ) ;
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
           result = ( some_arguments->count() == 3 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double );
          
           break ;
      case THREE :
	   result = ( some_arguments->count() == 3 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) ;
	   break ;
      case FOUR :
           result = ( some_arguments->count() == 2 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) ;
           break ;
   }

   return( result ) ;
}
