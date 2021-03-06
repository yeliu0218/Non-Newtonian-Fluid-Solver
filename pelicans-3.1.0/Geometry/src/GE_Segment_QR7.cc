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

#include <GE_Segment_QR7.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_Product_QR.hh>
#include <GE_ReferenceCube.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_RefinedCube_QR.hh>
#include <GE_RefinedSegment_QR.hh>
#include <GE_RefinedSquare_QR.hh>

GE_Segment_QR7 const* GE_Segment_QR7::REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Segment_QR7:: GE_Segment_QR7( void )
//----------------------------------------------------------------------
  : GE_QuadratureRule( "GE_Segment_QR7", GE_ReferenceSegment::object(), 7 )
{
   double const a = PEL::sqrt((3.-2.*PEL::sqrt(1.2))/7.) ;
   double const b = PEL::sqrt((3.+2.*PEL::sqrt(1.2))/7.) ;
   double const wa = ( 0.5 + 1./(6.*PEL::sqrt(1.2)))/2 ;
   double const wb = ( 0.5 - 1./(6.*PEL::sqrt(1.2)))/2 ;

   append_point( GE_Point::create( this, (1.+ a)/2. ), wa ) ;
   append_point( GE_Point::create( this, (1.- a)/2. ),  wa ) ;
   append_point( GE_Point::create( this, (1.+ b)/2. ), wb ) ;
   append_point( GE_Point::create( this, (1.- b)/2. ),  wb ) ;
 
   set_sum_of_weights( 1. ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Segment_QR7:: ~GE_Segment_QR7( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Segment_QR7 const*
GE_Segment_QR7:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Segment_QR7 const* rule1D = 0 ;
   if( rule1D == 0 )
   {
      rule1D = new GE_Segment_QR7() ;

      GE_RefinedSegment_QR const* rrule1D =
                GE_RefinedSegment_QR::create( "GE_Segment2R_QR7", rule1D ) ;

      GE_Product_QR const* rule2D = GE_Product_QR::create( 
                                               "GE_Segment_QR7_to_Square", 
                                               rule1D, 
                                               GE_ReferenceSquare::object() ) ;

      GE_RefinedSquare_QR const* rrule2D = GE_RefinedSquare_QR::create( 
                                           "GE_Segment_QR7_to_Square4R", 
                                           rule2D ) ;

      GE_Product_QR const* rule3D = GE_Product_QR::create( 
                                               "GE_Segment_QR7_to_Cube", 
                                               rule1D, 
                                               GE_ReferenceCube::object() ) ;
      
      GE_RefinedCube_QR const* rrule3D = GE_RefinedCube_QR::create( 
                                           "GE_Segment_QR7_to_Cube8R", 
                                           rule3D ) ;

      rrule1D = 0 ;
      rule2D = 0 ; rrule2D = 0 ;
      rule3D = 0 ; rrule3D = 0 ;
   }
   return( rule1D ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Segment_QR7:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment_QR7:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( nb_points()==4 ) ;
   return( true ) ;
}
