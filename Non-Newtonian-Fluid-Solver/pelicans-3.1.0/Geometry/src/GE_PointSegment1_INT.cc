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

#include <GE_PointSegment1_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Point.hh>
#include <GE_PointPoint_INT.hh>

GE_PointSegment1_INT const*
GE_PointSegment1_INT::PROTOTYPE = new GE_PointSegment1_INT() ;

//----------------------------------------------------------------------
GE_PointSegment1_INT*
GE_PointSegment1_INT:: create( PEL_Object* a_owner,
                               GE_PointPoint_INT const* a_pt_pt_intersector,
                               double a_epsilon )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointSegment1_INT:: create" ) ;
   PEL_CHECK_PRE( a_pt_pt_intersector!=0 ) ;
   PEL_CHECK_PRE( a_epsilon>0. ) ;
   
   GE_PointSegment1_INT* result =
      new GE_PointSegment1_INT( a_owner,
                                a_pt_pt_intersector,
                                a_epsilon ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointSegment1_INT*
GE_PointSegment1_INT:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointSegment1_INT:: create_clone" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   
   GE_PointSegment1_INT* result =
      new GE_PointSegment1_INT( a_owner,
                                PT_PT_INTERSECTOR,
                                EPSILON ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointSegment1_INT:: GE_PointSegment1_INT( void )
//----------------------------------------------------------------------
   : GE_PointSegment_INT( "GE_PointSegment1_INT" )
   , EPSILON( -PEL::max_double() )
   , PT_PT_INTERSECTOR( 0 )
{
   PEL_LABEL( "GE_PointSegment1_INT:: GE_PointSegment1_INT" ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointSegment1_INT:: GE_PointSegment1_INT(
                               PEL_Object* a_owner,
                               GE_PointPoint_INT const* a_pt_pt_intersector,
                               double a_epsilon )
//----------------------------------------------------------------------
   : GE_PointSegment_INT( a_owner )
   , EPSILON( a_epsilon )
   , PT_PT_INTERSECTOR( a_pt_pt_intersector->create_clone( this ) )
{
   PEL_LABEL( "GE_PointSegment1_INT:: GE_PointSegment1_INT" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointSegment1_INT:: ~GE_PointSegment1_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointSegment1_INT:: ~GE_PointSegment1_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_PointSegment1_INT:: point_in_segment( GE_Point const* P,
                                         GE_Point const* Q1,
                                         GE_Point const* Q2 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointSegment1_INT:: point_in_segment" ) ;
   PEL_CHECK_PRE( point_in_segment_PRE( P, Q1, Q2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = PT_PT_INTERSECTOR->points_are_close( P, Q1 ) ||
                 PT_PT_INTERSECTOR->points_are_close( P, Q2 ) ;
   if( !result )
   {
      double const dQ1Q2 = Q1->distance( Q2 ) ;
      double const dPQ1  = P->distance( Q1 ) ;
      double const dPQ2  = P->distance( Q2 ) ;
      result = ( PEL::abs( dQ1Q2-(dPQ1+dPQ2) )<EPSILON ) ;
   }
   return( result )  ;
}

//----------------------------------------------------------------------
GE_PointSegment_INT*
GE_PointSegment1_INT:: create_replica(
                 PEL_Object* a_owner,
                 PEL_ModuleExplorer const* a_mod_exp,
                 GE_PointPoint_INT const* a_pt_pt_intersector ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointSegment1_INT:: create_replica" ) ;
   PEL_CHECK( a_mod_exp!=0 ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const a_epsilon = a_mod_exp->double_data( "epsilon" ) ;
   a_mod_exp->test_data( "epsilon", "epsilon>0." ) ;

   if( a_pt_pt_intersector==0 )
   {
      PEL_Error::object()->raise_module_error(
         a_mod_exp,
         "A point-point intersector is expected\n"
         "for the building of \"GE_PointSegment1_INT\" object" ) ;
   }
   
   GE_PointSegment_INT* result =
      new GE_PointSegment1_INT( a_owner,
                                a_pt_pt_intersector,
                                a_epsilon ) ;
   
   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_PointSegment1_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( GE_PointSegment_INT::invariant() ) ;
   return( true )  ;
}
