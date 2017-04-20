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

#include <GE_Square_QR5.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_RefinedSquare_QR.hh>

GE_Square_QR5 const* 
GE_Square_QR5:: REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Square_QR5:: GE_Square_QR5( void )
//----------------------------------------------------------------------
   : GE_QuadratureRule( "GE_Square_QR5", GE_ReferenceSquare::object(), 5 )
{
   double const q = 0.0 ;
   double const wq = 2.0/7.0 ;
   double const a1 = 0.0 ;
   double const a2 = PEL::sqrt(14.0/15.0) ;
   double const wa = 5.0/63.0 ;
   double const b1 = PEL::sqrt(3.0/5.0) ;
   double const b2 = PEL::sqrt(1.0/3.0) ;
   double const wb = 5.0/36.0 ;

   append_point( GE_Point::create( this, 0.5*(  q+1.0), 0.5*(  q+1.0) ), wq ) ;
   append_point( GE_Point::create( this, 0.5*( a1+1.0), 0.5*( a2+1.0) ), wa ) ;
   append_point( GE_Point::create( this, 0.5*( a1+1.0), 0.5*(-a2+1.0) ), wa ) ;
   append_point( GE_Point::create( this, 0.5*( b1+1.0), 0.5*( b2+1.0) ), wb ) ;
   append_point( GE_Point::create( this, 0.5*( b1+1.0), 0.5*(-b2+1.0) ), wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b1+1.0), 0.5*( b2+1.0) ), wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b1+1.0), 0.5*(-b2+1.0) ), wb ) ;

   set_sum_of_weights( 1.0 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Square_QR5:: ~GE_Square_QR5( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Square_QR5 const*
GE_Square_QR5:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Square_QR5 const* rule = 0 ;
   if( rule == 0 )
   {
      rule = new GE_Square_QR5() ;

      GE_RefinedSquare_QR const* rrule = 
                GE_RefinedSquare_QR::create( "GE_Square4R_QR5", rule ) ;
      
      rrule = 0 ;
   }
   return( rule ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Square_QR5:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Square_QR5:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( IMPLIES(REGISTRATOR!=0,GE_QuadratureRule::nb_points()==7) ) ;

   return( true ) ;
}
