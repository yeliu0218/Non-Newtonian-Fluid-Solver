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

#include <GE_SegmentSegment1_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Point.hh>
#include <GE_PointSegment_INT.hh>

GE_SegmentSegment1_INT const*
GE_SegmentSegment1_INT::PROTOTYPE = new GE_SegmentSegment1_INT() ;

//----------------------------------------------------------------------
GE_SegmentSegment1_INT*
GE_SegmentSegment1_INT:: create( PEL_Object* a_owner,
                                 GE_PointSegment_INT const* a_pt_seg_intersector,
                                 double a_epsilon )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: create" ) ;
   PEL_CHECK_PRE( a_pt_seg_intersector!=0 ) ;
   PEL_CHECK_PRE( a_epsilon>0. ) ;
   
   GE_SegmentSegment1_INT* result =
                      new GE_SegmentSegment1_INT( a_owner,
                                                  a_pt_seg_intersector,
                                                  a_epsilon ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment1_INT*
GE_SegmentSegment1_INT:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: create_clone" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   
   GE_SegmentSegment1_INT* result =
                      new GE_SegmentSegment1_INT( a_owner,
                                                  PT_SEG_INT,
                                                  EPSILON ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment1_INT:: GE_SegmentSegment1_INT( void )
//----------------------------------------------------------------------
   : GE_SegmentSegment_INT( "GE_SegmentSegment1_INT" )
   , EPSILON( -PEL::max_double() )
   , PT_SEG_INT( 0 )
   , INTER_TYPE( no_intersection )
   , ALPHA( -PEL::max_double() )
   , BETA( -PEL::max_double() )
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: GE_SegmentSegment1_INT" ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment1_INT:: GE_SegmentSegment1_INT(
                              PEL_Object* a_owner,
                              GE_PointSegment_INT const* a_pt_seg_intersector,
                              double a_epsilon )
//----------------------------------------------------------------------
   : GE_SegmentSegment_INT( a_owner )
   , EPSILON( a_epsilon )
   , PT_SEG_INT( a_pt_seg_intersector->create_clone( this ) )
   , INTER_TYPE( no_intersection )
   , ALPHA( -PEL::max_double() )
   , BETA( -PEL::max_double() )
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: GE_SegmentSegment1_INT" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment1_INT:: ~GE_SegmentSegment1_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: ~GE_SegmentSegment1_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentSegment1_INT:: has_intersection( GE_Point const* P1,
                                           GE_Point const* P2,
                                           GE_Point const* Q1,
                                           GE_Point const* Q2 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: has_intersection" ) ;
   PEL_CHECK_PRE( has_intersection_PRE( P1, P2, Q1, Q2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool const Q1inP2P1 = PT_SEG_INT->point_in_segment( Q1, P1, P2 ) ;
   bool const Q2inP2P1 = PT_SEG_INT->point_in_segment( Q2, P1, P2 ) ;
   bool const P1inQ2Q1 = PT_SEG_INT->point_in_segment( P1, Q1, Q2 ) ;
   bool const P2inQ2Q1 = PT_SEG_INT->point_in_segment( P2, Q1, Q2 ) ;
   bool result = !( ( Q1inP2P1 && Q2inP2P1 ) || ( P1inQ2Q1 && P2inQ2Q1 ) ) ;
   if( result )
   {
      double const Q1x = Q1->coordinate( 0 ) ;
      double const Q1z = Q1->coordinate( 1 ) ;
      double const Q2x = Q2->coordinate( 0 ) ;
      double const Q2z = Q2->coordinate( 1 ) ;

      double const P1x = P1->coordinate( 0 ) ;
      double const P1z = P1->coordinate( 1 ) ;
      double const P2x = P2->coordinate( 0 ) ;
      double const P2z = P2->coordinate( 1 ) ;

      double vecP2P1_Q1P1 = (P2x-P1x)*(Q1z-P1z)-(P2z-P1z)*(Q1x-P1x) ;
      double vecP2P1_Q2P1 = (P2x-P1x)*(Q2z-P1z)-(P2z-P1z)*(Q2x-P1x) ;
      double vecQ2Q1_P1Q1 = (Q2x-Q1x)*(P1z-Q1z)-(Q2z-Q1z)*(P1x-Q1x) ;
      double vecQ2Q1_P2Q1 = (Q2x-Q1x)*(P2z-Q1z)-(Q2z-Q1z)*(P2x-Q1x) ;
      double vecQ2Q1_P2P1 = (Q2x-Q1x)*(P2z-P1z)-(Q2z-Q1z)*(P2x-P1x) ;
      bool const alpha1 = Q1inP2P1 || Q2inP2P1 ||
                                       vecP2P1_Q1P1*vecP2P1_Q2P1<0. ;
      bool const alpha2 =  P1inQ2Q1 ||  P2inQ2Q1 ||
                                       vecQ2Q1_P1Q1*vecQ2Q1_P2Q1<0. ;      
      result = ( alpha1 && alpha2 && PEL::abs( vecQ2Q1_P2P1 )>EPSILON ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
GE_SegmentSegment1_INT:: compute_intersection( GE_Point const* P1,
                                               GE_Point const* P2,
                                               GE_Point const* Q1,
                                               GE_Point const* Q2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: compute_intersection" ) ;
   PEL_CHECK_PRE( compute_intersection_PRE( P1, P2, Q1, Q2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   INTER_TYPE = no_intersection ;
   ALPHA = -PEL::max_double() ;
   BETA = -PEL::max_double() ;
   
   bool const Q1inP2P1 = PT_SEG_INT->point_in_segment( Q1, P1, P2 ) ;
   bool const Q2inP2P1 = PT_SEG_INT->point_in_segment( Q2, P1, P2 ) ;
   bool const P1inQ2Q1 = PT_SEG_INT->point_in_segment( P1, Q1, Q2 ) ;
   bool const P2inQ2Q1 = PT_SEG_INT->point_in_segment( P2, Q1, Q2 ) ;
   
   if( ( Q1inP2P1 && Q2inP2P1 ) || ( P1inQ2Q1 && P2inQ2Q1 ) )
   {
      INTER_TYPE = colinear_one_in_the_other ;
   }
   else
   {
      double const Q1x = Q1->coordinate( 0 ) ;
      double const Q1z = Q1->coordinate( 1 ) ;
      double const Q2x = Q2->coordinate( 0 ) ;
      double const Q2z = Q2->coordinate( 1 ) ;

      double const P1x = P1->coordinate( 0 ) ;
      double const P1z = P1->coordinate( 1 ) ;
      double const P2x = P2->coordinate( 0 ) ;
      double const P2z = P2->coordinate( 1 ) ;

      double vecP2P1_Q1P1 = (P2x-P1x)*(Q1z-P1z)-(P2z-P1z)*(Q1x-P1x) ;
      double vecP2P1_Q2P1 = (P2x-P1x)*(Q2z-P1z)-(P2z-P1z)*(Q2x-P1x) ;
      double vecQ2Q1_P1Q1 = (Q2x-Q1x)*(P1z-Q1z)-(Q2z-Q1z)*(P1x-Q1x) ;
      double vecQ2Q1_P2Q1 = (Q2x-Q1x)*(P2z-Q1z)-(Q2z-Q1z)*(P2x-Q1x) ;
      double vecQ2Q1_P2P1 = (Q2x-Q1x)*(P2z-P1z)-(Q2z-Q1z)*(P2x-P1x) ;
      bool const alpha1 = Q1inP2P1 || Q2inP2P1 ||
                                       vecP2P1_Q1P1*vecP2P1_Q2P1<0. ;
      bool const alpha2 =  P1inQ2Q1 ||  P2inQ2Q1 ||
                                       vecQ2Q1_P1Q1*vecQ2Q1_P2Q1<0. ;
      
      if( alpha1 && alpha2 && PEL::abs( vecQ2Q1_P2P1 )>EPSILON ) // intersection
      {
         ALPHA = PEL::max( 0., PEL::min( 1., -vecQ2Q1_P1Q1/vecQ2Q1_P2P1 ) ) ;
         BETA = PEL::max( 0., PEL::min( 1.,  vecP2P1_Q1P1/vecQ2Q1_P2P1 ) ) ;
         INTER_TYPE = one_intersection  ;
      }
      else // no intersection
      {
         if( PEL::abs( vecP2P1_Q1P1 )>EPSILON &&
             PEL::abs( vecQ2Q1_P1Q1 )>EPSILON )
         {
            // no colinear
            if( PEL::abs( vecQ2Q1_P2P1 )<EPSILON )
            {
               INTER_TYPE = parallel ; // Sides are strictly parallel
            }
            else
            {
               INTER_TYPE = no_intersection ;
            }
         }
         else // Sides are colinear
         {
            if( !Q1inP2P1 && !Q2inP2P1 && !P1inQ2Q1 && !P2inQ2Q1 )
            {
               // Sides do not overlap
               INTER_TYPE = colinear_disjoint ; // Sides do not overlap
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
GE_SegmentSegment1_INT:: intersection_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: intersection_type" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( INTER_TYPE ) ;
}

//----------------------------------------------------------------------
double
GE_SegmentSegment1_INT:: alpha( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: alpha" ) ;
   PEL_CHECK_PRE( alpha_PRE() ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   double result = ALPHA ;
   PEL_CHECK_POST( alpha_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
GE_SegmentSegment1_INT:: beta( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: beta" ) ;
   PEL_CHECK_PRE( beta_PRE() ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   double result = BETA ;
   PEL_CHECK_POST( beta_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentSegment_INT*
GE_SegmentSegment1_INT:: create_replica(
              PEL_Object* a_owner,
              PEL_ModuleExplorer const* a_mod_exp,
              GE_PointPoint_INT const* a_pt_pt_intersector,
              GE_PointSegment_INT const* a_pt_seg_intersector ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentSegment1_INT:: create_replica" ) ;
   PEL_CHECK( a_mod_exp!=0 ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const a_epsilon = a_mod_exp->double_data( "epsilon" ) ;
   a_mod_exp->test_data( "epsilon", "epsilon>0." ) ;

   if( a_pt_seg_intersector==0 )
   {
      PEL_Error::object()->raise_module_error(
         a_mod_exp,
         "A point-segment intersector is expected\n"
         "for the building of \"GE_SegmentSegment1_INT\" object" ) ;
   }
   
   GE_SegmentSegment_INT* result =
      new GE_SegmentSegment1_INT( a_owner,
                                  a_pt_seg_intersector,
                                  a_epsilon ) ;
   
   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentSegment1_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( GE_SegmentSegment_INT::invariant() ) ;
   return( true ) ;
}
