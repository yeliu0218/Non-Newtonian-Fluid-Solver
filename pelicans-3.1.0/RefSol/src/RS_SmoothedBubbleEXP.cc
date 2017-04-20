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

#include <RS_SmoothedBubbleEXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>

#include <iostream>

RS_SmoothedBubbleEXP const* 
RS_SmoothedBubbleEXP:: PROTOTYPE_U = 
   new RS_SmoothedBubbleEXP( "SmoothedBubble", U ) ;

RS_SmoothedBubbleEXP const* 
RS_SmoothedBubbleEXP:: PROTOTYPE_GRADU = 
   new RS_SmoothedBubbleEXP( "SmoothedBubbleGrad", GRADU ) ;

RS_SmoothedBubbleEXP const* 
RS_SmoothedBubbleEXP:: PROTOTYPE_RHS = 
   new RS_SmoothedBubbleEXP( "SmoothedBubbleForce", RHS ) ;


//----------------------------------------------------------------------
RS_SmoothedBubbleEXP:: RS_SmoothedBubbleEXP( std::string const& a_name,
                                             Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result( 1 ) 
   , doubleArray2D_result( 1, 2 )
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: RS_SmoothedBubbleEXP" ) ;
}

//----------------------------------------------------------------------
RS_SmoothedBubbleEXP*
RS_SmoothedBubbleEXP:: create_replica( PEL_Object* a_owner,
                                    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_SmoothedBubbleEXP* result = new RS_SmoothedBubbleEXP( a_owner, 
                                                            name(), 
                                                            argument_list,
                                                            EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_SmoothedBubbleEXP:: RS_SmoothedBubbleEXP(  PEL_Object* a_owner,
                                              std::string const& a_name,
                                              PEL_Sequence const* argument_list,
                                              Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_exp )
   , DV_result( 1 ) 
   , doubleArray2D_result( 1, 2 )
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: RS_SmoothedBubbleEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_SmoothedBubbleEXP:: ~RS_SmoothedBubbleEXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: ~RS_SmoothedBubbleEXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_SmoothedBubbleEXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: data_type" ) ;
   
   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : 
         result = DoubleVector ; 
         break ;
      case GRADU : 
         result = DoubleArray2D ;
         break ;
      case RHS : 
         result = DoubleVector ;
         break ;
   }
   return result ;      
}

//----------------------------------------------------------------------
std::string const& 
RS_SmoothedBubbleEXP:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( EXPR )
   {
      case U : 
         result = "SmoothedBubble( $DV_X, $DV_Center, $DS_Delta, $DS_epsilon )" ; 
         break ;
      case GRADU : 
         result = "SmoothedBubbleGrad( $DV_X, $DV_Center, $DS_Delta, $DS_epsilon )" ; 
         break ;
      case RHS : 
         result = "SmoothedBubbleForce( $DV_X, $DV_Center, $DS_Delta, $DS_epsilon )" ; 
         break ;
   }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_SmoothedBubbleEXP:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: valid_arguments" ) ;
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
RS_SmoothedBubbleEXP:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = DV_result ;
   result( 0 ) = 0.0 ;

   doubleVector const& xx = arg(0)->to_double_vector( ct ) ;
   doubleVector const& cc = arg(1)->to_double_vector( ct ) ;
   double dd = arg(2)->to_double( ct ) ;
   double ee = arg(3)->to_double( ct ) ;

   size_t nb_sp_dims = xx.size() ;
   PEL_ASSERT( cc.size() == nb_sp_dims ) ; //????????????
   double phi = 0.0 ;
   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      phi += (xx(d)-cc(d))*(xx(d)-cc(d)) ;
   }
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
   else if( EXPR == RHS )
   { 
      if( PEL::abs(dd-phi) < ee )
      {
         double Hp  = ( 1.0 + PEL::cos(aa) )/2.0/ee ;
         double Hpp = - PEL::sin(aa)*PEL::pi()/2.0/ee/ee ;
         result( 0 ) = -( Hpp - ( (double)nb_sp_dims - 1.0 )*Hp/phi ) ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_SmoothedBubbleEXP:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_SmoothedBubbleEXP:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ; 

   doubleArray2D& result = doubleArray2D_result ;
   
   doubleVector const& xx = arg(0)->to_double_vector( ct ) ;
   doubleVector const& cc = arg(1)->to_double_vector( ct ) ;
   double dd = arg(2)->to_double( ct ) ;
   double ee = arg(3)->to_double( ct ) ;

   size_t nb_sp_dims = xx.size() ;
   result.re_initialize( 1, nb_sp_dims ) ;
   PEL_ASSERT( cc.size() == nb_sp_dims ) ; //????????????
   double phi = 0.0 ;
   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      phi += (xx(d)-cc(d))*(xx(d)-cc(d)) ;
   }
   phi = PEL::sqrt( phi ) ;
   double aa = PEL::pi()*(dd-phi)/ee ;

   if( EXPR == GRADU )
   {
      if( PEL::abs(dd-phi) < ee )
      {
         double Hp  = ( 1.0 + PEL::cos(aa) )/2.0/ee ;
         for( size_t d=0 ; d<nb_sp_dims ; ++d )
         {
            result( 0, d ) = - Hp * ( xx(d) - cc(d) )/phi ;            
         }
      }
   }

   return result ;
}


