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

#include <GE_Triangle_QR9.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_RefinedTriangle_QR.hh>

GE_Triangle_QR9 const* 
GE_Triangle_QR9:: REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Triangle_QR9:: GE_Triangle_QR9( void )
//----------------------------------------------------------------------
   : GE_QuadratureRule( "GE_Triangle_QR9", GE_ReferenceTriangle::object(), 9 )
{
   double const r1 = (1.-1./3.*PEL::sqrt(5.+4.*PEL::sqrt(5./14.)))/2. ;
   double const r2 = (1.-1./3.*PEL::sqrt(5.-4.*PEL::sqrt(5./14.)))/2. ;
   double const r3 = 1./2. ;
   double const r4 = (1.+1./3.*PEL::sqrt(5.-4.*PEL::sqrt(5./14.)))/2. ;
   double const r5 = (1.+1./3.*PEL::sqrt(5.+4.*PEL::sqrt(5./14.)))/2. ;
   double const A1 = (161./450.-13./(180.*PEL::sqrt(5./14.)))/2. ;
   double const A2 = (161./450.+13./(180.*PEL::sqrt(5./14.)))/2. ;
   double const A3 = 128./(2.*225.) ;
   double const A4 = (161./450.+13./(180.*PEL::sqrt(5./14.)))/2. ;
   double const A5 = (161./450.-13./(180.*PEL::sqrt(5./14.)))/2. ;
   double const s1 = 0.0398098571 ;
   double const s2 = 0.1980134179 ;
   double const s3 = 0.4379748102 ;
   double const s4 = 0.6954642734 ;
   double const s5 = 0.9014649142 ;
   double const B1 = 0.1007941926*(1.-s1) ;
   double const B2 = 0.2084506672*(1.-s2) ;
   double const B3 = 0.2604633916*(1.-s3) ;
   double const B4 = 0.2426935942*(1.-s4) ;
   double const B5 = 0.1598203766*(1.-s5) ;

   append_point( GE_Point::create( this, s1, r1*(1.-s1) ), A1*B1 ) ;
   append_point( GE_Point::create( this, s2, r1*(1.-s2) ), A1*B2 ) ;
   append_point( GE_Point::create( this, s3, r1*(1.-s3) ), A1*B3 ) ;
   append_point( GE_Point::create( this, s4, r1*(1.-s4) ), A1*B4 ) ;
   append_point( GE_Point::create( this, s5, r1*(1.-s5) ), A1*B5 ) ;
   append_point( GE_Point::create( this, s1, r2*(1.-s1) ), A2*B1 ) ;
   append_point( GE_Point::create( this, s2, r2*(1.-s2) ), A2*B2 ) ;
   append_point( GE_Point::create( this, s3, r2*(1.-s3) ), A2*B3 ) ;
   append_point( GE_Point::create( this, s4, r2*(1.-s4) ), A2*B4 ) ;
   append_point( GE_Point::create( this, s5, r2*(1.-s5) ), A2*B5 ) ;
   append_point( GE_Point::create( this, s1, r3*(1.-s1) ), A3*B1 ) ;
   append_point( GE_Point::create( this, s2, r3*(1.-s2) ), A3*B2 ) ;
   append_point( GE_Point::create( this, s3, r3*(1.-s3) ), A3*B3 ) ;
   append_point( GE_Point::create( this, s4, r3*(1.-s4) ), A3*B4 ) ;
   append_point( GE_Point::create( this, s5, r3*(1.-s5) ), A3*B5 ) ;
   append_point( GE_Point::create( this, s1, r4*(1.-s1) ), A4*B1 ) ;
   append_point( GE_Point::create( this, s2, r4*(1.-s2) ), A4*B2 ) ;
   append_point( GE_Point::create( this, s3, r4*(1.-s3) ), A4*B3 ) ;
   append_point( GE_Point::create( this, s4, r4*(1.-s4) ), A4*B4 ) ;
   append_point( GE_Point::create( this, s5, r4*(1.-s5) ), A4*B5 ) ;
   append_point( GE_Point::create( this, s1, r5*(1.-s1) ), A5*B1 ) ;
   append_point( GE_Point::create( this, s2, r5*(1.-s2) ), A5*B2 ) ;
   append_point( GE_Point::create( this, s3, r5*(1.-s3) ), A5*B3 ) ;
   append_point( GE_Point::create( this, s4, r5*(1.-s4) ), A5*B4 ) ;
   append_point( GE_Point::create( this, s5, r5*(1.-s5) ), A5*B5 ) ;

   set_sum_of_weights( 0.5 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Triangle_QR9:: ~GE_Triangle_QR9( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Triangle_QR9 const*
GE_Triangle_QR9:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Triangle_QR9 const* rule = 0 ;
   if( rule == 0 )
   {
      rule = new GE_Triangle_QR9() ;

      GE_RefinedTriangle_QR const* rrule = 
                 GE_RefinedTriangle_QR::create( "GE_Triangle4R_QR9", rule ) ;
      
      rrule = 0 ;
   }
   return( rule ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle_QR9:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle_QR9:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( IMPLIES(REGISTRATOR!=0,GE_QuadratureRule::nb_points()==25) ) ;

   return( true ) ;
}
