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

#include <GE_SegmentSegment2_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Point.hh>

GE_SegmentSegment2_INT const*
GE_SegmentSegment2_INT::PROTOTYPE = new GE_SegmentSegment2_INT() ;

//----------------------------------------------------------------------
GE_SegmentSegment2_INT*
GE_SegmentSegment2_INT:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: create_clone" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   GE_SegmentSegment2_INT* result =
                    new GE_SegmentSegment2_INT( a_owner,
                                                ALPHA_BETA_EPSILON,
                                                DET_EPSILON,
                                                COORDS_EPSILON ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment2_INT*
GE_SegmentSegment2_INT:: create( PEL_Object* a_owner,
                                 double a_alpha_beta_epsilon,
                                 double a_determinant_epsilon,
                                 double a_coordinates_epsilon )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: create" ) ;
   PEL_CHECK_PRE( a_alpha_beta_epsilon>0. ) ;
   PEL_CHECK_PRE( a_determinant_epsilon>0. ) ;
   PEL_CHECK_PRE( a_coordinates_epsilon>0. ) ;
   
   GE_SegmentSegment2_INT* result =
                    new GE_SegmentSegment2_INT( a_owner,
                                                a_alpha_beta_epsilon,
                                                a_determinant_epsilon,
                                                a_coordinates_epsilon ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment2_INT:: GE_SegmentSegment2_INT( void )
//----------------------------------------------------------------------
   : GE_SegmentSegment_INT( "GE_SegmentSegment2_INT" )
   , ALPHA_BETA_EPSILON( -PEL::max_double() )
   , DET_EPSILON( -PEL::max_double() )
   , COORDS_EPSILON( -PEL::max_double() )
   , INTER_TYPE( no_intersection )
   , ALPHA( -PEL::max_double() )
   , BETA( -PEL::max_double() )
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: GE_SegmentSegment2_INT" ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
GE_SegmentSegment2_INT:: GE_SegmentSegment2_INT(
                              PEL_Object* a_owner,
                              double a_alpha_beta_epsilon,
                              double a_determinant_epsilon,
                              double a_coordinates_epsilon )
//----------------------------------------------------------------------
   : GE_SegmentSegment_INT( a_owner )
   , ALPHA_BETA_EPSILON( a_alpha_beta_epsilon )
   , DET_EPSILON( a_determinant_epsilon )
   , COORDS_EPSILON( a_coordinates_epsilon )
   , INTER_TYPE( no_intersection )
   , ALPHA( -PEL::max_double() )
   , BETA( -PEL::max_double() )
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: GE_SegmentSegment2_INT" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment2_INT:: ~GE_SegmentSegment2_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: ~GE_SegmentSegment2_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentSegment2_INT:: has_intersection( GE_Point const* P1,
                                           GE_Point const* P2,
                                           GE_Point const* Q1,
                                           GE_Point const* Q2 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: has_intersection" ) ;
   PEL_CHECK_PRE( has_intersection_PRE( P1, P2, Q1, Q2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double aX = P1->coordinate( 0 ) ;
   double aZ = P1->coordinate( 1 ) ;
   double bX = P2->coordinate( 0 ) ;
   double bZ = P2->coordinate( 1 ) ;
   double cX = Q1->coordinate( 0 ) ;
   double cZ = Q1->coordinate( 1 ) ;
   double dX = Q2->coordinate( 0 ) ;
   double dZ = Q2->coordinate( 1 ) ;

   double num_AB = (dX-cX)*(aZ-cZ)-(dZ-cZ)*(aX-cX) ;
   double num_CD = (bX-aX)*(aZ-cZ)-(bZ-aZ)*(aX-cX) ;
   double denomi = (bX-aX)*(dZ-cZ)-(bZ-aZ)*(dX-cX) ;

   bool result = false ;
   
   if( PEL::abs( denomi )>DET_EPSILON )
   {
      // Denomi non nul ... the sides are not parallel
      double aa = num_AB/denomi ;
      double bb = num_CD/denomi ;
      if( aa>-ALPHA_BETA_EPSILON && aa-1.<ALPHA_BETA_EPSILON && 
          bb>-ALPHA_BETA_EPSILON && bb-1.<ALPHA_BETA_EPSILON )
      {
         result = true ;
      }
      else
      {
         result = false ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
GE_SegmentSegment2_INT:: compute_intersection( GE_Point const* P1,
                                               GE_Point const* P2,
                                               GE_Point const* Q1,
                                               GE_Point const* Q2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: compute_intersection" ) ;
   PEL_CHECK_PRE( compute_intersection_PRE( P1, P2, Q1, Q2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   INTER_TYPE = no_intersection ;
   ALPHA = -PEL::max_double() ;
   BETA = -PEL::max_double() ;

   double aX = P1->coordinate( 0 ) ;
   double aZ = P1->coordinate( 1 ) ;
   double bX = P2->coordinate( 0 ) ;
   double bZ = P2->coordinate( 1 ) ;
   double cX = Q1->coordinate( 0 ) ;
   double cZ = Q1->coordinate( 1 ) ;
   double dX = Q2->coordinate( 0 ) ;
   double dZ = Q2->coordinate( 1 ) ;

   double num_AB = (dX-cX)*(aZ-cZ)-(dZ-cZ)*(aX-cX) ;
   double num_CD = (bX-aX)*(aZ-cZ)-(bZ-aZ)*(aX-cX) ;
   double denomi = (bX-aX)*(dZ-cZ)-(bZ-aZ)*(dX-cX) ;
   
   if( PEL::abs( denomi )>DET_EPSILON )
   {
      // Denomi non nul ... the sides are not parallel
      ALPHA = num_AB/denomi ;
      BETA = num_CD/denomi ;
      if( ALPHA>-ALPHA_BETA_EPSILON && ALPHA-1.<ALPHA_BETA_EPSILON && 
          BETA>-ALPHA_BETA_EPSILON && BETA-1.<ALPHA_BETA_EPSILON )
      {
         INTER_TYPE = one_intersection ;
         ALPHA = PEL::max( 0., PEL::min( 1., ALPHA ) ) ;
         BETA = PEL::max( 0., PEL::min( 1., BETA ) ) ;
      }
      else
      {
         INTER_TYPE = no_intersection ;
         ALPHA = -PEL::max_double() ;
         BETA = -PEL::max_double() ;
      }
   }
   else // Sides are parallel
   {
      if( PEL::abs( num_AB )>DET_EPSILON &&
          PEL::abs( num_CD )>DET_EPSILON )
      {
         // Sides are strictly parallel
	 INTER_TYPE = parallel ;
      }
      else // Sides are colinear
      {
	 if( PEL::abs( aX - bX )>COORDS_EPSILON )
            // [AB] is not vertical projection on the horizontal axis
	 {
	    double abXmin = PEL::min( aX, bX ), abXmax = PEL::max( aX, bX ) ;
	    double cdXmin = PEL::min( cX, dX ), cdXmax = PEL::max( cX, dX ) ;
	    if( cdXmax-abXmin<-COORDS_EPSILON ||
                abXmax-cdXmin<-COORDS_EPSILON )
            {
               // Sides do not overlap
	       INTER_TYPE = colinear_disjoint ;
            }
	    else if( (abXmin<=cdXmin&&cdXmax<=abXmax) || 
                     (cdXmin<=abXmin&&abXmax<=cdXmax) )
            {
               // One side lies enterely in the other one
	       INTER_TYPE = colinear_one_in_the_other ;
            }
	    else
	    {
	       // each side share one enpoint of the other one
               INTER_TYPE = colinear ;
	    }
	 }
	 else // [AB] is not horizontal projection on the vertical axis
	 {
	    double abZmin = PEL::min( aZ, bZ ), abZmax = PEL::max( aZ, bZ ) ;
	    double cdZmin = PEL::min( cZ, dZ ), cdZmax = PEL::max( cZ, dZ ) ;
	    if( cdZmax-abZmin<-COORDS_EPSILON ||
                abZmax-cdZmin<-COORDS_EPSILON )
            {
               // Sides do not overlap
	       INTER_TYPE = colinear_disjoint ;
            }
	    else if( (abZmin<=cdZmin&&cdZmax<=abZmax) || 
                     (cdZmin<=abZmin&&abZmax<=cdZmax) )
            {
               // One side lies enterely in the other one
	       INTER_TYPE = colinear_one_in_the_other ;
            }
	    else
	    {
	       // each side share one enpoint of the other one
               INTER_TYPE = colinear ;
            }
         }
      }
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( compute_intersection_POST( P1, P2, Q1, Q2 ) ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment_INT::IntersectionType
GE_SegmentSegment2_INT:: intersection_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: intersection_type" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( INTER_TYPE ) ;
}

//----------------------------------------------------------------------
double
GE_SegmentSegment2_INT:: alpha( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: alpha" ) ;
   PEL_CHECK_PRE( alpha_PRE() ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   double result = ALPHA ;
   PEL_CHECK_POST( alpha_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
GE_SegmentSegment2_INT:: beta( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: beta" ) ;
   PEL_CHECK_PRE( beta_PRE() ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   double result = BETA ;
   PEL_CHECK_POST( beta_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment_INT*
GE_SegmentSegment2_INT:: create_replica(
               PEL_Object* a_owner,
               PEL_ModuleExplorer const* a_mod_exp,
               GE_PointPoint_INT const* a_pt_pt_intersector,
               GE_PointSegment_INT const* a_pt_seg_intersector ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment2_INT:: create_replica" ) ;
   PEL_CHECK( a_mod_exp!=0 ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const a_alpha_beta_epsilon =
                           a_mod_exp->double_data( "alpha_beta_epsilon" ) ;
   a_mod_exp->test_data( "alpha_beta_epsilon", "alpha_beta_epsilon>0." ) ;

   double const a_determinant_epsilon =
                           a_mod_exp->double_data( "determinant_epsilon" ) ;
   a_mod_exp->test_data( "determinant_epsilon", "determinant_epsilon>0." ) ;

   double const a_coordinates_epsilon =
                           a_mod_exp->double_data( "coordinates_epsilon" ) ;
   a_mod_exp->test_data( "coordinates_epsilon", "coordinates_epsilon>0." ) ;
   
   GE_SegmentSegment_INT* result =
      new GE_SegmentSegment2_INT( a_owner,
                                  a_alpha_beta_epsilon,
                                  a_determinant_epsilon,
                                  a_coordinates_epsilon ) ;

   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentSegment2_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( GE_SegmentSegment_INT::invariant() ) ;
   return( true ) ;
}
