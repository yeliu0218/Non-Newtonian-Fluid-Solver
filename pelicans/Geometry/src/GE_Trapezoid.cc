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

#include <GE_Trapezoid.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_SetOfPoints.hh>

#include <doubleArray3D.hh>

#include <iostream>

GE_Trapezoid const* GE_Trapezoid::MODEL = new GE_Trapezoid() ;

//----------------------------------------------------------------------------
GE_Trapezoid:: GE_Trapezoid( void )
//----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Trapezoid" )
   , MEASURE( 0 )
{
   PEL_LABEL( "GE_Trapezoid:: GE_Trapezoid" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
GE_Trapezoid*
GE_Trapezoid:: create_replica( PEL_Object* a_owner,
                               PEL_Vector* vertices ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Trapezoid* result = new GE_Trapezoid( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Trapezoid:: GE_Trapezoid( PEL_Object* a_owner, PEL_Vector* vertices )
//----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Trapezoid:: GE_Trapezoid" ) ;
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
GE_Trapezoid:: ~GE_Trapezoid( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: ~GE_Trapezoid" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//----------------------------------------------------------------------------
std::string const&
GE_Trapezoid:: name( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static std::string const result = "GE_Trapezoid" ;

   PEL_CHECK_POST( result == "GE_Trapezoid" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
GE_Trapezoid:: is_consistent( std::ostream& os, bool verbose ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool ok = ( nb_space_dimensions()==2 ) ;

   if( ok )
   {
      ok = PEL::toler( vertex(0)->coordinate(0)-vertex(3)->coordinate(0) ) &&
           PEL::toler( vertex(1)->coordinate(0)-vertex(2)->coordinate(0) ) ;
   }

   // IL MANQUE Qqchose sur les faces qui ne doivent pas se croiser ...
   return( ok ) ;
}

//-----------------------------------------------------------------------------
double
GE_Trapezoid:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   double result = MEASURE ;
   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Trapezoid:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = 0 ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
GE_Trapezoid:: contains( GE_Point const* pt ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   static GE_Point* pt_ref =
                         GE_Point::create( PEL_Root::object(), (size_t) 2 ) ;
   PEL_CHECK( pt_ref!=0 &&
              pt_ref->nb_coordinates()==dimension() ) ;
   apply_inverse_mapping( pt, pt_ref ) ;
   bool is_in = reference_polyhedron()->contains( pt_ref ) ;
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Trapezoid:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceSquare::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: apply_mapping( GE_Point const* pt_ref,
                              GE_Point* pt  ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double const x_ref = pt_ref->coordinate( 0 ) ;
   double const y_ref = pt_ref->coordinate( 1 ) ;

   double const x_v0 = vertex(0)->coordinate(0) ;
   double const y_v0 = vertex(0)->coordinate(1) ;

   double const x_real = x_ref*( vertex(1)->coordinate(0)-x_v0 )+x_v0 ;
   double const y_low= (1.-x_ref)*y_v0+x_ref*vertex(1)->coordinate(1) ;
   double const y_high =
       (1.-x_ref)*vertex(3)->coordinate(1)+x_ref*vertex(2)->coordinate(1) ;
   double const y_real = y_ref*(y_high-y_low)+y_low ;

   pt->set_coordinate( (size_t)0, x_real ) ;
   pt->set_coordinate( (size_t)1, y_real ) ;

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: apply_inverse_mapping( GE_Point const* pt,
                                      GE_Point* pt_ref ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double const x_real = pt->coordinate(0) ;
   double const y_real = pt->coordinate(1) ;

   double const x_v0 = vertex(0)->coordinate(0) ;
   double const y_v0 = vertex(0)->coordinate(1) ;

   double const x_ref = ( x_real-x_v0 )/( vertex(1)->coordinate(0)-x_v0 ) ;
   double const y_low= (1.-x_ref)*y_v0+x_ref*vertex(1)->coordinate(1) ;
   double const y_high =
       (1.-x_ref)*vertex(3)->coordinate(1)+x_ref*vertex(2)->coordinate(1) ;
   double const y_ref = (y_real-y_low)/(y_high-y_low) ;

   pt_ref->set_coordinate( (size_t)0, x_ref ) ;
   pt_ref->set_coordinate( (size_t)1, y_ref ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: build_mapping_derivative( GE_Point const* pt_ref,
                                         GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   double hx = vertex( 1 )->coordinate( 0 ) - vertex( 0 )->coordinate( 0 ) ;

   jac->set_item( 0, 0, hx ) ;
   jac->set_item( 0, 1, 0.0 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double dh1y = vertex( 1 )->coordinate( 1 ) - vertex( 0 )->coordinate( 1 ) ;
   double dh2y = vertex( 2 )->coordinate( 1 ) - vertex( 3 )->coordinate( 1 ) ;

   jac->set_item( 1, 0, y*( dh2y - dh1y ) + dh1y ) ;
   jac->set_item( 1, 1,  vertex( 3 )->coordinate( 1 ) + x * dh2y
                       - vertex( 0 )->coordinate( 1 ) - x * dh1y ) ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                            GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   double hx = vertex( 1 )->coordinate( 0 ) - vertex( 0 )->coordinate( 0 ) ;

   tjac->set_item( 0, 0, hx ) ;
   tjac->set_item( 1, 0, 0.0 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double dh1y = vertex( 1 )->coordinate( 1 ) - vertex( 0 )->coordinate( 1 ) ;
   double dh2y = vertex( 2 )->coordinate( 1 ) - vertex( 3 )->coordinate( 1 ) ;

   tjac->set_item( 0, 1, y*( dh2y - dh1y ) + dh1y ) ;
   tjac->set_item( 1, 1,  vertex( 3 )->coordinate( 1 ) + x * dh2y
                        - vertex( 0 )->coordinate( 1 ) - x * dh1y ) ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: build_mapping_hessian( GE_Point const* pt_ref,
                                      doubleArray3D* hessian,
                                      bool& nonzero_hessian ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: build_mapping_hessian" ) ;
   PEL_CHECK_PRE( build_mapping_hessian_PRE( pt_ref,
                                             hessian, nonzero_hessian ) ) ;

   nonzero_hessian = true ;
   hessian->set( 0.0 ) ;

   double dh1y = vertex( 1 )->coordinate( 1 ) - vertex( 0 )->coordinate( 1 ) ;
   double dh2y = vertex( 2 )->coordinate( 1 ) - vertex( 3 )->coordinate( 1 ) ;

   double xx = dh2y - dh1y ;
   (*hessian)(1,0,1) = xx ;
   (*hessian)(1,1,0) = xx ;
}

//----------------------------------------------------------------------------
void
GE_Trapezoid:: update_internal( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
}

//----------------------------------------------------------------------------
double
GE_Trapezoid:: computed_measure( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Trapezoid:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;

   double norm01 = 0., norm03 = 0., scal0103 = 0. ;
   for( size_t d=0; d<2 ; ++d )
   {
      double contr01 = vertex(1)->coordinate(d) - vertex(0)->coordinate(d) ;
      double contr03 = vertex(3)->coordinate(d) - vertex(0)->coordinate(d) ;
      norm01   += contr01*contr01 ;
      norm03   += contr03*contr03 ;
      scal0103 += contr01*contr03 ;
   }

   double area =
            0.5*PEL::sqrt( PEL::abs( norm01*norm03 - scal0103*scal0103 ) ) ;

   double norm21 = 0., norm23 = 0., scal2123 = 0. ;
   for( size_t d=0; d<2 ; ++d )
   {
      double contr21 = vertex(1)->coordinate(d) - vertex(2)->coordinate(d) ;
      double contr23 = vertex(3)->coordinate(d) - vertex(2)->coordinate(d) ;
      norm21   += contr21*contr21 ;
      norm23   += contr23*contr23 ;
      scal2123 += contr21*contr23 ;
   }
   area += 0.5*PEL::sqrt( PEL::abs( norm21*norm23 - scal2123*scal2123 ) ) ;

   return( area ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Trapezoid:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Trapezoid:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceSquare::object() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
GE_Trapezoid:: invariant( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT( PEL::equal( MEASURE, computed_measure() ) ) ;
   }
   return( true ) ;
}
