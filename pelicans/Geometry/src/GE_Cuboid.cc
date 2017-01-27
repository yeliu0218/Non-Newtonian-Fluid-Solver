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

#include <GE_Cuboid.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceCube.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>

#include <iostream>

GE_Cuboid const* GE_Cuboid::MODEL = new GE_Cuboid() ;

//-----------------------------------------------------------------------------
GE_Cuboid:: GE_Cuboid( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Cuboid" )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Cuboid:: GE_Cuboid" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Cuboid*
GE_Cuboid:: create_replica( PEL_Object* a_owner,
                            PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Cuboid* result = new GE_Cuboid( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Cuboid:: GE_Cuboid( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Cuboid:: GE_Cuboid" ) ;
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Cuboid:: ~GE_Cuboid( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: ~GE_Cuboid" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Cuboid:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static std::string const result = "GE_Cuboid" ;

   PEL_CHECK_POST( result == "GE_Cuboid" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Cuboid:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double const result = MEASURE ;

   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Cuboid:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = center() ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cuboid:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   static GE_Point* pt_ref =
                 GE_Point:: create( PEL_Root::object(), (size_t) 3 ) ;
   PEL_CHECK( pt_ref!=0 && pt_ref->nb_coordinates()==dimension() ) ;
   apply_inverse_mapping( pt, pt_ref ) ;
   bool const is_in = reference_polyhedron()->contains( pt_ref ) ;
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Cuboid:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceCube::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Cuboid:: apply_mapping( GE_Point const* pt_ref,
                           GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   double umxmymz = 1.0 - x - y - z ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, umxmymz * V0->coordinate( d ) +
                                   x * V1->coordinate( d ) +
                                   y * V3->coordinate( d ) +
                                   z * V4->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Cuboid:: apply_inverse_mapping( GE_Point const* pt,
                                   GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;

   double dot10P=0., dot30P=0., dot40P=0., norm01=0., norm03=0., norm04=0. ;
   for( size_t d=0; d<nb_space_dimensions(); d++ )
   {
      double v0d  = V0->coordinate( d)  ;
      double v0v1 = V1->coordinate( d ) - v0d ;
      double v0v3 = V3->coordinate( d ) - v0d ;
      double v0v4 = V4->coordinate( d ) - v0d ;
      norm01 += v0v1 * v0v1 ;
      norm03 += v0v3 * v0v3 ;
      norm04 += v0v4 * v0v4 ;

      double x = pt->coordinate( d ) ;
      dot10P += v0v1 * ( x - v0d ) ;
      dot30P += v0v3 * ( x - v0d ) ;
      dot40P += v0v4 * ( x - v0d ) ;
   }

   pt_ref->set_coordinate( (size_t)0, dot10P/norm01 ) ;
   pt_ref->set_coordinate( (size_t)1, dot30P/norm03 ) ;
   pt_ref->set_coordinate( (size_t)2, dot40P/norm04 ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Cuboid:: build_mapping_derivative( GE_Point const* pt_ref,
                                      GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double const v0d = V0->coordinate( d)  ;
      jac->set_item( d, 0, V1->coordinate( d ) - v0d ) ;
      jac->set_item( d, 1, V3->coordinate( d ) - v0d ) ;
      jac->set_item( d, 2, V4->coordinate( d ) - v0d ) ;
   }
}

//----------------------------------------------------------------------------
void
GE_Cuboid:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                         GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double const v0d = V0->coordinate( d)  ;
      tjac->set_item( 0, d, V1->coordinate( d ) - v0d ) ;
      tjac->set_item( 1, d, V3->coordinate( d ) - v0d ) ;
      tjac->set_item( 2, d, V4->coordinate( d ) - v0d ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Cuboid:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
}

//-----------------------------------------------------------------------------
double
GE_Cuboid:: computed_measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;
   double const vol =
      vertex( 0 )->distance( vertex( 1 ) )*
      vertex( 0 )->distance( vertex( 3 ) )*
      vertex( 0 )->distance( vertex( 4 ) ) ;
   return( vol ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cuboid:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Cuboid:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool ok = ( nb_space_dimensions()>=3 ) ;

   // "Bottom" quadrilateral is a rectangle orthogonal to a "vertical" edge
   GE_Vector* a = GE_Vector::create( 0, vertex( 1 ), vertex( 0 ) ) ;
   GE_Vector* b = GE_Vector::create( 0, vertex( 3 ), vertex( 0 ) ) ;
   ok = PEL::toler( a->dot_product( b ) ) ;
   GE_Vector* c = GE_Vector::create( 0, vertex( 3 ), vertex( 2 ) ) ;
   ok = ok && PEL::toler( b->dot_product( c ) ) ;
   ok = ok && PEL::equal( a->norm(), c->norm() ) ;
   c->re_initialize( vertex( 4 ), vertex( 0 ) ) ;
   ok = ok && PEL::toler( a->dot_product( c ) ) ;
   ok = ok && PEL::toler( b->dot_product( c ) ) ;

   // "Top" quadrilateral is a rectangle orthogonal to a "vertical" edge
   a->re_initialize( vertex( 5 ), vertex( 4 ) ) ;
   b->re_initialize( vertex( 7 ), vertex( 4 ) ) ;
   ok = ok && PEL::toler( a->dot_product( b ) ) ;
   ok = ok && PEL::toler( a->dot_product( c ) ) ;
   ok = ok && PEL::toler( b->dot_product( c ) ) ;
   c->re_initialize( vertex( 6 ), vertex( 7 ) ) ;
   ok = ok && PEL::toler( b->dot_product( c ) ) ;
   ok = ok && PEL::equal( a->norm(), c->norm() ) ;

   a->destroy() ;
   b->destroy() ;
   c->destroy() ;

   return( ok ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cuboid:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == center() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cuboid:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceCube::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Cuboid:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT( PEL::equal( MEASURE, computed_measure() ) ) ;
   }
   return( true ) ;
}
