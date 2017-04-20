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

#include <GE_ReferenceTriangle.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_Vector.hh>

GE_ReferenceTriangle*
GE_ReferenceTriangle:: UNIQUE_INSTANCE = 0 ;

//---------------------------------------------------------------------
GE_ReferenceTriangle const*
GE_ReferenceTriangle:: object( void )
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangle:: object" ) ;

   if( UNIQUE_INSTANCE == 0 )
   {
      UNIQUE_INSTANCE = new GE_ReferenceTriangle() ;
   }
   GE_ReferenceTriangle const* result = UNIQUE_INSTANCE ;
   
   PEL_CHECK_PRE( result!=0 ) ;
   PEL_CHECK_PRE( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_ReferenceTriangle:: GE_ReferenceTriangle( void )
//----------------------------------------------------------------------
   : GE_ReferencePolyhedron( "GE_ReferenceTriangle", 3, 3, 2, 0.5 )
{
   set_vertex( 0, GE_Point::create( this, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 0.0, 1.0 ) ) ;

   doubleVector n(2) ;

   size_t i_face = 0 ;
   append_face_vertex( i_face, 0 ) ;
   append_face_vertex( i_face, 1 ) ;
   n(0) =  0. ; n(1) = -1. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   i_face = 1 ;
   append_face_vertex( i_face, 1 ) ;
   append_face_vertex( i_face, 2 ) ;
   n(0) =  1./PEL::sqrt(2.) ; n(1) = 1./PEL::sqrt(2.) ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   i_face = 2 ;
   append_face_vertex( i_face, 0 ) ;
   append_face_vertex( i_face, 2 ) ;
   n(0) = -1. ; n(1) = 0. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_ReferenceTriangle:: ~GE_ReferenceTriangle( void )
//----------------------------------------------------------------------
{
   UNIQUE_INSTANCE = 0 ;
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------
bool
GE_ReferenceTriangle:: contains( GE_Point const* pt_ref,
                                 double tol ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangle:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt_ref, tol ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const x = pt_ref->coordinate( 0 ) ;
   double const y = pt_ref->coordinate( 1 ) ;

   bool const result = ( x+tol>0. && y+tol>0. && x+y-tol<1. ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
bool
GE_ReferenceTriangle:: face_contains( size_t i_face,
                                      GE_Point const* pt_ref,
                                      double tol ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangle:: face_contains" ) ;
   PEL_CHECK_PRE( face_contains_PRE( i_face, pt_ref, tol ) ) ;

   double const x = pt_ref->coordinate( 0 ) ;
   double const y = pt_ref->coordinate( 1 ) ;

   bool result = false ;
   if( i_face == 0 )
   {
      result = x+tol>0. && x-tol<1. &&        // 0<=x<=1
               y+tol>0. && y-tol<0. ;         // y=0
   }
   else if( i_face == 1 )
   {
      result = x+y+tol>1. && x+y-tol<1. &&    // x+y=1
               x+tol>0.                 &&    // x>=0
               y+tol>0.                 ;     // y>=0
   }
   else if( i_face == 2 )
   {
      result = x+tol>0. && x-tol<0. &&        // x=0
               y+tol>0. && y-tol<1. ;         // 0<=y<=1
   }
   return( result ) ;
}

//---------------------------------------------------------------------
void
GE_ReferenceTriangle:: project( GE_Point* pt_ref ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangle:: project" ) ;
   PEL_CHECK_PRE( project_PRE( pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   double xi  = PEL::max( 0., pt_ref->coordinate( 0 ) ) ;
   double eta = PEL::max( 0., pt_ref->coordinate( 1 ) ) ;
   if( xi+eta > 1. )
   {
      eta = PEL::max( 0., PEL::min( 1., 0.5*( 1. + eta - xi ) ) ) ;
      xi  = 1. - eta ;
   }
   pt_ref->set_coordinate( 0,  xi ) ;
   pt_ref->set_coordinate( 1, eta ) ;

   PEL_CHECK_POST( project_POST( pt_ref ) ) ;
}

//----------------------------------------------------------------------
void
GE_ReferenceTriangle:: build_neighbor( GE_Point const* pt_ref,
                                       size_t ic,
                                       GE_Point* neighbor,
                                       double& delta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangle:: build_neighbor" ) ;
   PEL_CHECK_PRE( build_neighbor_PRE( pt_ref, ic, neighbor ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   delta = 1.E-4 ;

   for( size_t i=0 ; i<dimension() ; ++i )
   {
      neighbor->set_coordinate( i, pt_ref->coordinate(i) ) ;
   }
   if( neighbor->coordinate(0)+neighbor->coordinate(1)+delta>1. )
   {
      if( neighbor->coordinate(ic)<delta )
      {
         size_t jc = ic+1 ;
         if( jc==2 )
         {
            jc = 0 ;
         }
         double const x_new = neighbor->coordinate(jc)-delta ;
         neighbor->set_coordinate( jc, x_new ) ;
      }
      else
      {
         delta = -delta ;
      }
   }
   neighbor->set_coordinate( ic, neighbor->coordinate(ic)+delta ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( build_neighbor_POST( pt_ref, ic, neighbor, delta ) ) ;
}
