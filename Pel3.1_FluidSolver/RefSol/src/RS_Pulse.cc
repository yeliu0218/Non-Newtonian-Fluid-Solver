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

#include <RS_Pulse.hh>

#include <RS_Bingham.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <math.h>
#include <iostream>

RS_Pulse const* 
RS_Pulse::PROTOTYPE_ONE 
    = new RS_Pulse( "Periodic_Inphase", ONE ) ;

RS_Pulse const* 
RS_Pulse::PROTOTYPE_TWO 
    = new RS_Pulse( "Periodic_Outphase", TWO ) ;
    
RS_Pulse const* 
	RS_Pulse::PROTOTYPE_THREE 
	= new RS_Pulse( "Step_Inphase", THREE ) ;

RS_Pulse const* 
	RS_Pulse::PROTOTYPE_FOUR 
	= new RS_Pulse( "Step_Outphase", FOUR ) ;

//----------------------------------------------------------------------
RS_Pulse:: RS_Pulse( 
                               std::string const& a_name, Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_Pulse" ) ;
}

//----------------------------------------------------------------------
RS_Pulse*
RS_Pulse:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Pulse:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Pulse* result =
               new RS_Pulse( a_owner, name(), argument_list, EXPR ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Pulse:: RS_Pulse(
                 PEL_Object* a_owner, std::string const& a_name,
		 PEL_Sequence const* argument_list, Func an_expr  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_Pulse:: RS_Pulse" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Pulse:: ~RS_Pulse( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Pulse:: ~RS_Pulse" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Pulse:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}
//----------------------------------------------------------------------
double
RS_Pulse:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Pulse:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double result = 0.0;
   switch( EXPR )
   {
	case ONE:
{          double const r = arg(0)->to_double(ct) ;
	   double const Um1 = arg(1)->to_double(ct) ;
	   double const A1 = arg(2)->to_double(ct) ;
	   double const B1 = arg(3)->to_double(ct) ;
	   double const Um2 = arg(4)->to_double(ct) ;
           double const A2 = arg(5)->to_double(ct) ;
	   double const B2 = arg(6)->to_double(ct) ;
	   double const T = arg(7)->to_double(ct) ;
	   double const Rin = arg(8)->to_double(ct) ;
           double const Rout = arg(9)->to_double(ct) ;
           
           if (r >=0.0 && r<=Rin)
             {
               result = Um1 + A1 *PEL::sin(B1*T);
               }
           if (r >Rin && r<=1)
             {
               result = Um2 + A2 *PEL::sin(B2*T);
              }
           

	   break;
}
	case TWO:
{          double const r = arg(0)->to_double(ct) ;
	   double const Um1 = arg(1)->to_double(ct) ;
	   double const A1 = arg(2)->to_double(ct) ;
	   double const B1 = arg(3)->to_double(ct) ;
	   double const Um2 = arg(4)->to_double(ct) ;
           double const A2 = arg(5)->to_double(ct) ;
	   double const B2 = arg(6)->to_double(ct) ;
	   double const T = arg(7)->to_double(ct) ;
	   double const Rin = arg(8)->to_double(ct) ;
           //double const Rout = arg(9)->to_double(ct) ;
          
           if (r >=0.0 && r<=Rin)
             {
               result = Um1 + A1 *PEL::sin(B1*T);
               }
           if (r >Rin && r<=1)
             {
               result = Um2 + A2 *PEL::cos(B2*T);
              }
           

	   break;
}
	case THREE:
{          double const r = arg(0)->to_double(ct) ;
	   double const Um1 = arg(1)->to_double(ct) ;
	   double const Um2 = arg(2)->to_double(ct) ;           
	   double const T = arg(3)->to_double(ct) ;
           double const TI1 = arg(4)->to_double(ct) ;
           double const TI2 = arg(5)->to_double(ct) ;
	   double const Rin = arg(6)->to_double(ct) ;
           //double const Rout = arg(7)->to_double(ct) ;
           
           if (r >=0.0 && r<=Rin)
             {
                     if(int(floor(T/TI1))%2==0)
                       {
                        result = 0.0;
                       }
                        else 
                       {
                       result  = Um1;
                       }
               }
           if (r >Rin && r<=1)
              {
                     if(int(floor(T/TI2))%2==0)
                       {
                        result = 0.0;
                       }
                        else 
                       {
                       result  = Um2;
                       }
               }
           

	   break;
}

	case FOUR:
{          double const r = arg(0)->to_double(ct) ;
	   double const Um1 = arg(1)->to_double(ct) ;
	   double const Um2 = arg(2)->to_double(ct) ;           
	   double const T = arg(3)->to_double(ct) ;
           double const TI1 = arg(4)->to_double(ct) ;
           double const TI2 = arg(5)->to_double(ct) ;
	   double const Rin = arg(6)->to_double(ct) ;
           //double const Rout = arg(7)->to_double(ct) ;
          
           if (r >=0.0 && r<=Rin)
             {
                     if(int(floor(T/TI1))%2==0)
                       {
                        result = 0.0;
                       }
                        else 
                       {
                       result  = Um1;
                       }
               }
           if (r >Rin && r<=1)
              {
                     if(int(floor(T/TI2))%2==0)
                       {
                        result = Um2;
                       }
                        else 
                       {
                       result  = 0.0;
                       }
               }
           

	   break;
}

   }
	

   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
RS_Pulse:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Pulse:: usage" ) ;

   static std::string result ;

   switch( EXPR )
   {
      case ONE :
         result = "Periodic_Inphase($DS_r,$DS_Um1,$DS_A1,$DS_B1,$DS_Um2,$DS_A2,$DS_B2,$DS_T,$DS_Rin,$DS_Rout)" ;
         break ;
      case TWO :
        result = "Periodic_Outphase($DS_r,$DS_Um1,$DS_A1,$DS_B1,$DS_Um2,$DS_A2,$DS_B2,$DS_T,$DS_Rin,$DS_Rout)" ;
        break ;
	case THREE :
	  result = "Step_Inphase($DS_r,$DS_Um1,$DS_Um2,$DS_T,$DS_TI1,$DS_TI2,$DS_Rin,$DS_Rout)" ;
	  break ;
	   case FOUR :
	  result = "Step_Outphase($DS_r,$DS_Um1,$DS_Um2,$DS_T,$DS_TI1,$DS_TI2,$DS_Rin,$DS_Rout)" ;
	  break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_Pulse:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Pulse:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( EXPR )
   {
      case ONE :
	   result = ( some_arguments->count() == 10 ) &&
                   ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
                   ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double )&&
                   ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 8 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 9 )->data_type() == PEL_Data::Double ) ;
	   break ;
      case TWO :
         result = ( some_arguments->count() == 10 ) &&
                   ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
                   ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double )&&
                   ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 8 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 9 )->data_type() == PEL_Data::Double ) ;
          
         break ;
	case THREE :
	   result = ( some_arguments->count() == 8 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) ;
	   break ;
        case FOUR :
	      result = ( some_arguments->count() == 8 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double )&&
	   ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) ;
	   break ;
 
   }

   return( result ) ;
}
