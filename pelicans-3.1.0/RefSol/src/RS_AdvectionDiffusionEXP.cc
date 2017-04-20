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

#include <RS_AdvectionDiffusionEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh> 

#include <iostream>

RS_AdvectionDiffusionEXP const* 
RS_AdvectionDiffusionEXP::PROTOTYPE_ONE 
    = new RS_AdvectionDiffusionEXP( "AdvectionDiffusion1_value", ONE ) ;

RS_AdvectionDiffusionEXP const* 
RS_AdvectionDiffusionEXP::PROTOTYPE_TWO 
    = new RS_AdvectionDiffusionEXP( "AdvectionDiffusion2_value", TWO ) ;

//----------------------------------------------------------------------
RS_AdvectionDiffusionEXP:: RS_AdvectionDiffusionEXP( 
                               std::string const& a_name, Func an_expr ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: RS_AdvectionDiffusionEXP" ) ;
}

//----------------------------------------------------------------------
RS_AdvectionDiffusionEXP*
RS_AdvectionDiffusionEXP:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_AdvectionDiffusionEXP* result =
               new RS_AdvectionDiffusionEXP( a_owner,
                                             name(), 
                                             argument_list,
                                             EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_AdvectionDiffusionEXP:: RS_AdvectionDiffusionEXP(
                 PEL_Object* a_owner, std::string const& a_name,
		 PEL_Sequence const* argument_list, Func an_expr  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: RS_AdvectionDiffusionEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_AdvectionDiffusionEXP:: ~RS_AdvectionDiffusionEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: ~RS_AdvectionDiffusionEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_AdvectionDiffusionEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}

//----------------------------------------------------------------------
double
RS_AdvectionDiffusionEXP:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double result = PEL::max_double() ;
   if( EXPR == ONE )
   {
      doubleVector const& po = arg( 0 )->to_double_vector( ct ) ;
      double const av = arg( 1 )->to_double( ct ) ;
      PEL_ASSERT( arg( 2 )->to_double( ct )>0. ) ;
      double const df = arg( 2 )->to_double( ct ) ;
      double const ui = arg( 3 )->to_double( ct ) ;
      double const uo = arg( 4 )->to_double( ct ) ;
      double const Pe = av/df ;
      result = ui+(uo-ui)*(PEL::exp( Pe*po(0) )-1.)/(PEL::exp( Pe )-1.) ;
   }
   else if( EXPR == TWO )
   {
      doubleVector const& po = arg( 0 )->to_double_vector( ct ) ;
      double const t = arg( 1 )->to_double( ct ) ;
      doubleVector const& av = arg( 2 )->to_double_vector( ct ) ;
      doubleVector const& ip = arg( 3 )->to_double_vector( ct ) ;
      PEL_ASSERT( arg( 4 )->to_double_vector( ct )( 0 )>0. ) ;
      PEL_ASSERT( arg( 4 )->to_double_vector( ct )( 1 )>0. ) ;
      doubleVector const& df = arg( 4 )->to_double_vector( ct ) ;

      double const tt = 4.*t+1. ;
      result=PEL::exp( -PEL::pow(po(0)-av(0)*t-ip(0),2.)/df(0)/tt
		       -PEL::pow(po(1)-av(1)*t-ip(1),2.)/df(1)/tt )/tt ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
RS_AdvectionDiffusionEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case ONE :
         result = "AdvectionDiffusion1_value($DV_X,$DS_a,$DS_k,$DS_uin,$DS_out)" ;
         break ;
      case TWO :
	 result = "AdvectionDiffusion2_value($DV_X,$DS_T,$DV_a,$DV_X0,$DV_K)" ;
	 break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_AdvectionDiffusionEXP:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_AdvectionDiffusionEXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( EXPR )
   {
      case ONE :
         result =  ( some_arguments->count() == 5 ) &&
	    ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector )&&
	    ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
	    ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
	    ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
	    ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) ;
         break ;
      case TWO :
	 result = ( some_arguments->count() == 5 ) &&
	    ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector )&&
	    ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
	    ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::DoubleVector )&&
	    ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::DoubleVector )&&
	    ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::DoubleVector ) ;
	 break ;
   }
   
   return( result ) ;
}
