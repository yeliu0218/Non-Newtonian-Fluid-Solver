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

#include <RS_SmoothedTubeEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_SmoothedTubeEXP const* 
RS_SmoothedTubeEXP:: PROTOTYPE_U = 
   new RS_SmoothedTubeEXP( "SmoothedTube", U ) ;

RS_SmoothedTubeEXP const* 
RS_SmoothedTubeEXP:: PROTOTYPE_DU = 
   new RS_SmoothedTubeEXP( "SmoothedTubeForce", DU ) ;

//----------------------------------------------------------------------
RS_SmoothedTubeEXP:: RS_SmoothedTubeEXP( std::string const& a_name,
                                   Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , FIRST(0)
   , SECOND(0)
   , THIRD(0)
   , FOURTH(0)
   , EXPR( an_exp )
   , DV_result( 1 ) 
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: RS_SmoothedTubeEXP" ) ;
}

//----------------------------------------------------------------------
RS_SmoothedTubeEXP*
RS_SmoothedTubeEXP:: create_replica( PEL_Object* a_owner,
                                     PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_SmoothedTubeEXP* result = new RS_SmoothedTubeEXP( a_owner, 
                                                        name(), 
                                                        argument_list, 
                                                        EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_SmoothedTubeEXP:: RS_SmoothedTubeEXP( PEL_Object* a_owner,
                                         std::string const& a_name,
                                         PEL_Sequence const* argument_list,
                                         Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , FIRST( arg( 0 ) )
   , SECOND( arg( 1 ) )
   , THIRD( arg( 2 ) )
   , FOURTH( arg( 3 ) )
   , EXPR( an_exp )
   , DV_result( 1 ) 
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: RS_SmoothedTubeEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_SmoothedTubeEXP:: ~RS_SmoothedTubeEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: ~RS_SmoothedTubeEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_SmoothedTubeEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: data_type" ) ;
   
   return DoubleVector ;      
}

//----------------------------------------------------------------------
std::string const& 
RS_SmoothedTubeEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch(EXPR)
   {
      case U : 
         result = "SmoothedTube( $DV_X, $DV_Center, $DS_Delta, $DS_epsilon )" ; 
         break ;
      case DU : 
         result = "SmoothedTubeForce( $DV_X, $DV_Center, $DS_Delta, $DS_epsilon )" ; 
	 break ;
   }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_SmoothedTubeEXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = ( some_arguments->count() == 4 ) ;
   result = ( extract_arg(some_arguments,0)->data_type() == DoubleVector ) &&
            ( extract_arg(some_arguments,1)->data_type() == DoubleVector ) &&
            ( extract_arg(some_arguments,2)->data_type() == Double )       &&
            ( extract_arg(some_arguments,3)->data_type() == Double ) ;
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
RS_SmoothedTubeEXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedTubeEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = DV_result ;
   result( 0 ) = 0.0 ;

   doubleVector const& xx = FIRST->to_double_vector( ct ) ;
   doubleVector const& cc = SECOND->to_double_vector( ct ) ;
   double dd = THIRD->to_double( ct ) ;
   double ee = FOURTH->to_double( ct ) ;

   PEL_ASSERT( xx.size() == 3 ) ; //????????????
   PEL_ASSERT( cc.size() == 3 ) ; //????????????
   double phi = ( xx(1) - cc(1) ) * ( xx(1) - cc(1) ) +
                ( xx(2) - cc(2) ) * ( xx(2) - cc(2) ) ;
   phi = PEL::sqrt( phi ) ;
   double aa = PEL::pi()*(dd-phi)/ee ;

   if( EXPR == U )
   {
      if( PEL::abs(dd-phi) < ee )
      {
         result( 0 ) = 0.5*( 1.0 + (dd-phi)/ee + PEL::sin(aa)/PEL::pi() ) ;
      }
      else if( dd-phi >= ee )
      {
         result( 0 ) = 1.0 ;
      }
   }
   else if( EXPR == DU )
   { 
      if( PEL::abs(dd-phi) < ee )
      {
         double Hp  = ( 1.0 + PEL::cos(aa) )/2.0/ee ;
         double Hpp = - PEL::sin(aa)*PEL::pi()/2.0/ee/ee ;
         result( 0 ) = -( Hpp - Hp/phi ) ;
      }
   }
   return result ;
}
