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

#include <GE_Cube_QR5.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceCube.hh>
#include <GE_RefinedCube_QR.hh>

GE_Cube_QR5 const* 
GE_Cube_QR5:: REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Cube_QR5:: GE_Cube_QR5( void )
//-----------------------------------------------------------------------------
   : GE_QuadratureRule( "GE_Cube_QR5", GE_ReferenceCube::object(), 5 )
{
   double const a = PEL::sqrt(19.0/30.0) ;
   double const b = PEL::sqrt(19.0/33.0) ;
   double const wa = 320.0/361.0/8.0 ;
   double const wb = 121.0/361.0/8.0 ;

   append_point( GE_Point::create( this, 0.5*( a+1.),0.5,0.5), wa ) ;
   append_point( GE_Point::create( this, 0.5*(-a+1.),0.5,0.5), wa ) ;
   append_point( GE_Point::create( this, 0.5,0.5*( a+1.),0.5), wa ) ;
   append_point( GE_Point::create( this, 0.5,0.5*(-a+1.),0.5), wa ) ;
   append_point( GE_Point::create( this, 0.5,0.5,0.5*( a+1.)), wa ) ;
   append_point( GE_Point::create( this, 0.5,0.5,0.5*(-a+1.)), wa ) ;
   append_point( GE_Point::create( this, 0.5*(b+1.),0.5*(b+1.),0.5*( b+1.)   ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(b+1.),0.5*(b+1.),0.5*(-b+1.)   ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(b+1.),0.5*(-b+1.),0.5*(b+1.)   ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(b+1.),0.5*(-b+1.),0.5*(-b+1.)  ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b+1.),0.5*(b+1.),0.5*(b+1.)   ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b+1.),0.5*(b+1.),0.5*(-b+1.)  ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b+1.),0.5*(-b+1.),0.5*(b+1.)  ),
                 wb ) ;
   append_point( GE_Point::create( this, 0.5*(-b+1.),0.5*(-b+1.),0.5*(-b+1.) ),
                 wb ) ;

   set_sum_of_weights( 1.0 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Cube_QR5:: ~GE_Cube_QR5( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Cube_QR5 const*
GE_Cube_QR5:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Cube_QR5 const* rule = 0 ;
   if( rule == 0 )
   {
      rule = new GE_Cube_QR5() ;

      GE_RefinedCube_QR const* rrule = 
                GE_RefinedCube_QR::create( "GE_Cube8R_QR5", rule ) ;
      
      rrule = 0 ;
   }
   return( rule ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cube_QR5:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cube_QR5:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( IMPLIES(REGISTRATOR!=0,GE_QuadratureRule::nb_points()==14) ) ;

   return( true ) ;
}
