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

#include <GE_ReferenceSquare.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_Vector.hh>

GE_ReferenceSquare*
GE_ReferenceSquare:: UNIQUE_INSTANCE = 0 ;

//---------------------------------------------------------------------
GE_ReferenceSquare const*
GE_ReferenceSquare:: object( void )
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquare:: object" ) ;

   if( UNIQUE_INSTANCE == 0 )
   {
      UNIQUE_INSTANCE = new GE_ReferenceSquare() ;
   }
   GE_ReferenceSquare const* result = UNIQUE_INSTANCE ;
   
   PEL_CHECK_PRE( result!=0 ) ;
   PEL_CHECK_PRE( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_ReferenceSquare:: GE_ReferenceSquare( void )
//----------------------------------------------------------------------
   : GE_ReferencePolyhedron( "GE_ReferenceSquare", 4, 4, 2, 1.0 )
{
   set_vertex( 0, GE_Point::create( this, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 1.0, 1.0 ) ) ;
   set_vertex( 3, GE_Point::create( this, 0.0, 1.0 ) ) ;

   doubleVector n(2) ;

   size_t i_face = 0 ;
   append_face_vertex( i_face, 0 ) ;
   append_face_vertex( i_face, 1 ) ;
   n(0) =  0. ; n(1) = -1. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   i_face = 1 ;
   append_face_vertex( i_face, 1 ) ;
   append_face_vertex( i_face, 2 ) ;
   n(0) =  1. ; n(1) =  0. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   i_face = 2 ;
   append_face_vertex( i_face, 2 ) ;
   append_face_vertex( i_face, 3 ) ;
   n(0) =  0. ; n(1) =  1. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   i_face = 3 ;
   append_face_vertex( i_face, 0 ) ;
   append_face_vertex( i_face, 3 ) ;
   n(0) = -1. ; n(1) =  0. ;
   set_face_normal( i_face, GE_Vector::create( this, n ) ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_ReferenceSquare:: ~GE_ReferenceSquare( void )
//----------------------------------------------------------------------
{
   UNIQUE_INSTANCE = 0 ;
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------
bool
GE_ReferenceSquare:: contains( GE_Point const* pt_ref, double tol ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquare:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt_ref, tol ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = true ;
   for( size_t i=0 ; result && i<dimension() ; ++i )
   {
      double x = pt_ref->coordinate(i) ;
      result = ( x+tol>0. && x-tol<1. ) ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------
bool
GE_ReferenceSquare:: face_contains( size_t i_face,
                                    GE_Point const* pt_ref,
                                    double tol ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquare:: face_contains" ) ;
   PEL_CHECK_PRE( face_contains_PRE( i_face, pt_ref, tol ) ) ;

   double const x = pt_ref->coordinate( 0 ) ;
   double const y = pt_ref->coordinate( 1 ) ;

   bool result = false ;
   if( i_face == 0 )
   {
      result = x+tol>0. && x-tol<1. &&     // 0<=x<=1
               y+tol>0. && y-tol<0. ;      // y=0
   }
   else if( i_face == 1 )
   {
      result = x+tol>1. && x-tol<1. &&     // x=1
               y+tol>0. && y-tol<1. ;      // 0<=y<=1
   }
   else if( i_face == 2 )
   {
      result = x+tol>0. && x-tol<1. &&     // 0<=x<=1
               y+tol>1. && y-tol<1. ;      // y=1
   }
   else if( i_face == 3 )
   {
      result = x+tol>0. && x-tol<0. &&     // x=0
               y+tol>0. && y-tol<1. ;      // 0<=y<=1
   }
   return( result ) ;
}

//---------------------------------------------------------------------
void
GE_ReferenceSquare:: project( GE_Point* pt_ref ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquare:: project" ) ;
   PEL_CHECK_PRE( project_PRE( pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<dimension() ; ++i )
   {
      pt_ref->set_coordinate( i,
                              PEL::max(
                                 0.,
                                 PEL::min( 1.,
                                           pt_ref->coordinate( i ) ) ) ) ;
   }

   PEL_CHECK_POST( project_POST( pt_ref ) ) ;
}

//----------------------------------------------------------------------
void
GE_ReferenceSquare:: build_neighbor( GE_Point const* pt_ref,
                                     size_t ic,
                                     GE_Point* neighbor,
                                     double& delta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquare:: build_neighbor" ) ;
   PEL_CHECK_PRE( build_neighbor_PRE( pt_ref, ic, neighbor ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<dimension() ; ++i )
   {
      neighbor->set_coordinate( i, pt_ref->coordinate(i) ) ;
   }

   delta = 1.E-4 ;
   if( pt_ref->coordinate(ic)+delta>1. )
   {
      delta = -delta ;
   }
   neighbor->set_coordinate( ic, pt_ref->coordinate(ic)+delta ) ;   

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( build_neighbor_POST( pt_ref, ic, neighbor, delta ) ) ;
}
