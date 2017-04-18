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

#include <GE_Tetrahedron.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTetrahedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>

#include <iostream>

GE_Tetrahedron const*
GE_Tetrahedron::MODEL = new GE_Tetrahedron() ;

//-----------------------------------------------------------------------------
GE_Tetrahedron:: GE_Tetrahedron( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Tetrahedron" )
   , ORIENT( 0. )
   , MEASURE( 0. )
   , FV_CENTER( 0 )
{
   PEL_LABEL( "GE_Tetrahedron:: GE_Tetrahedron" ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Tetrahedron*
GE_Tetrahedron:: create_replica( PEL_Object* a_owner,
                                 PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Tetrahedron* result = new GE_Tetrahedron( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Tetrahedron:: GE_Tetrahedron( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , ORIENT( 0. )
   , MEASURE( 0. )
   , FV_CENTER( 0 )
{
   PEL_LABEL( "GE_Tetrahedron:: GE_Tetrahedron" ) ;
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Tetrahedron:: ~GE_Tetrahedron( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: ~GE_Tetrahedron" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Tetrahedron:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static const std::string result = "GE_Tetrahedron" ;

   PEL_CHECK_POST( result == "GE_Tetrahedron" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Tetrahedron:: measure( GE_Point const* v0,
                          GE_Point const* v1,
                          GE_Point const* v2,
                          GE_Point const* v3 )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: measure" ) ;
   PEL_CHECK_PRE( v0 != 0 && v0->nb_coordinates() == 3 ) ;
   PEL_CHECK_PRE( v1 != 0 && v1->nb_coordinates() == 3 ) ;
   PEL_CHECK_PRE( v2 != 0 && v2->nb_coordinates() == 3 ) ;
   PEL_CHECK_PRE( v3 != 0 && v3->nb_coordinates() == 3 ) ;

   double const result = PEL::abs( signed_measure( v0, v1, v2, v3 ) ) ;

   PEL_CHECK_POST( result>=0. ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Tetrahedron:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( !is_prototype() && !is_updating() ) ;

   double result = MEASURE ;

   PEL_CHECK_POST( measure_POST( result ) ) ;
   PEL_CHECK_POST(
      PEL::equal( result,
                  PEL::abs( signed_measure( vertex(0), vertex(1),
                                            vertex(2), vertex(3) ) ) ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Tetrahedron:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   if( FV_CENTER == 0 )
   {
      FV_CENTER = GE_Point::create( const_cast<GE_Tetrahedron*>( this ),
                                    (size_t)3 ) ; // necessarily 3D
      compute_finite_volume_center( vertex(0), vertex(1),
                                    vertex(2), vertex(3),
                                    FV_CENTER ) ;
   }

   GE_Point const* result = FV_CENTER ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Tetrahedron:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( !is_prototype() && !is_updating() ) ;

   static GE_Point* pt_ref =
                 GE_Point:: create( PEL_Root::object(), (size_t) 3 ) ;
   apply_inverse_mapping( pt, pt_ref ) ;
   bool result = reference_polyhedron()->contains( pt_ref ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Tetrahedron:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result =
                                            GE_ReferenceTetrahedron::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Tetrahedron:: apply_mapping( GE_Point const* pt_ref,
                                GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   double umxmymz = 1.0 - x - y - z ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, umxmymz * V0->coordinate( d ) +
                                   x * V1->coordinate( d ) +
                                   y * V2->coordinate( d ) +
                                   z * V3->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Tetrahedron:: apply_inverse_mapping( GE_Point const* pt,
                                        GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( !is_prototype() && !is_updating() ) ;

   static GE_Vector* p = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;
   static GE_Vector* a = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;
   static GE_Vector* b = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;
   static GE_Vector* c = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   p->re_initialize( pt, V0 ) ;
   a->re_initialize( V2, V0 ) ;
   b->re_initialize( V3, V0 ) ;
   c->set_as_cross_product( a, b ) ;
   pt_ref->set_coordinate( 0, ORIENT*c->dot_product( p )/MEASURE/6. ) ;

   a->re_initialize( V1, V0 ) ;
   c->set_as_cross_product( a, b ) ;
   pt_ref->set_coordinate( 1, -1.*ORIENT*c->dot_product( p )/MEASURE/6. ) ;
   b->re_initialize( V2, V0 ) ;
   c->set_as_cross_product( a, b ) ;
   pt_ref->set_coordinate( 2,  ORIENT*c->dot_product( p )/MEASURE/6. ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Tetrahedron:: build_mapping_derivative( GE_Point const* pt_ref,
                                           GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      jac->set_item( d, 0, V1->coordinate( d ) - v0d ) ;
      jac->set_item( d, 1, V2->coordinate( d ) - v0d ) ;
      jac->set_item( d, 2, V3->coordinate( d ) - v0d ) ;
   }
}

//----------------------------------------------------------------------------
void
GE_Tetrahedron:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                              GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      tjac->set_item( 0, d, V1->coordinate( d ) - v0d ) ;
      tjac->set_item( 1, d, V2->coordinate( d ) - v0d ) ;
      tjac->set_item( 2, d, V3->coordinate( d ) - v0d ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Tetrahedron:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   double const vol = signed_measure( vertex(0), vertex(1),
                                      vertex(2), vertex(3) ) ;

   MEASURE = PEL::abs( vol ) ;
   ORIENT = vol/PEL::max( MEASURE, 1.E-14 ) ;

   // le calcul de FV_CENTER n'est fait que si il avait deja été fait avant
   if( FV_CENTER != 0 )
      compute_finite_volume_center( vertex(0), vertex(1),
                                    vertex(2), vertex(3),
                                    FV_CENTER ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Tetrahedron:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool ok = ( nb_space_dimensions()>=3 ) ;

   return( ok ) ;
}

//-----------------------------------------------------------------------------
double
GE_Tetrahedron:: signed_measure( GE_Point const* v0,
                                 GE_Point const* v1,
                                 GE_Point const* v2,
                                 GE_Point const* v3 )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: signed_measure" ) ;
   PEL_CHECK( v0 != 0 && v0->nb_coordinates() == 3 ) ;
   PEL_CHECK( v1 != 0 && v1->nb_coordinates() == 3 ) ;
   PEL_CHECK( v2 != 0 && v2->nb_coordinates() == 3 ) ;
   PEL_CHECK( v3 != 0 && v3->nb_coordinates() == 3 ) ;

   static GE_Vector* N1 = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;
   static GE_Vector* N2 = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;
   static GE_Vector* N3 = GE_Vector::create( PEL_Root::object(), (size_t)3 ) ;

   N1->re_initialize( v1, v0 ) ;
   N2->re_initialize( v2, v0 ) ;
   N3->set_as_cross_product( N1, N2 ) ;
   N1->re_initialize( v3, v0 ) ;

   return( N3->dot_product( N1 )/6. ) ;
}

//-----------------------------------------------------------------------------
void
GE_Tetrahedron:: compute_finite_volume_center( GE_Point const* v0,
                                               GE_Point const* v1,
                                               GE_Point const* v2,
                                               GE_Point const* v3,
                                               GE_Point* result )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Tetrahedron:: compute_finite_volume_center" ) ;
   PEL_CHECK( v0 != 0 && v0->nb_coordinates() == 3 ) ;
   PEL_CHECK( v1 != 0 && v1->nb_coordinates() == 3 ) ;
   PEL_CHECK( v2 != 0 && v2->nb_coordinates() == 3 ) ;
   PEL_CHECK( v3 != 0 && v3->nb_coordinates() == 3 ) ;
   PEL_CHECK( result != 0 && result->nb_coordinates() == 3 ) ;

   // Reference:
   // http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
   //
   // Note: the following implementation is different from that of GE_Triangle.
   // it may be optimized by suppressing multiply done calculations.

   GE_Point const* a = v0 ;
   GE_Point const* b = v1 ;
   GE_Point const* c = v2 ;
   GE_Point const* d = v3 ;

   // Use coordinates relative to point `a' of the tetrahedron.
   double xba = b->coordinate( 0 ) - a->coordinate( 0 ) ;
   double yba = b->coordinate( 1 ) - a->coordinate( 1 ) ;
   double zba = b->coordinate( 2 ) - a->coordinate( 2 ) ;
   double xca = c->coordinate( 0 ) - a->coordinate( 0 ) ;
   double yca = c->coordinate( 1 ) - a->coordinate( 1 ) ;
   double zca = c->coordinate( 2 ) - a->coordinate( 2 ) ;
   double xda = d->coordinate( 0 ) - a->coordinate( 0 ) ;
   double yda = d->coordinate( 1 ) - a->coordinate( 1 ) ;
   double zda = d->coordinate( 2 ) - a->coordinate( 2 ) ;

   // Squares of lengths of the edges incident to `a'.
   double balength = xba * xba + yba * yba + zba * zba ;
   double calength = xca * xca + yca * yca + zca * zca ;
   double dalength = xda * xda + yda * yda + zda * zda ;

   // Cross products of these edges.
   double xcrosscd = yca * zda - yda * zca ;
   double ycrosscd = zca * xda - zda * xca ;
   double zcrosscd = xca * yda - xda * yca ;
   double xcrossdb = yda * zba - yba * zda ;
   double ycrossdb = zda * xba - zba * xda ;
   double zcrossdb = xda * yba - xba * yda ;
   double xcrossbc = yba * zca - yca * zba ;
   double ycrossbc = zba * xca - zca * xba ;
   double zcrossbc = xba * yca - xca * yba ;

   // Calculate the denominator of the formulae.
   // Take your chances with floating-point roundoff.
   double denominator =
          0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd) ;

   // Calculate offset (from `a') of circumcenter.
   double xcirca =
          (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
           denominator ;
   double ycirca =
          (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
           denominator ;
   double zcirca =
          (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
           denominator ;
   result->set_coordinate( 0, a->coordinate( 0 ) + xcirca ) ;
   result->set_coordinate( 1, a->coordinate( 1 ) + ycirca ) ;
   result->set_coordinate( 2, a->coordinate( 2 ) + zcirca ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Tetrahedron:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result != 0 &&
               result->owner() == this &&
               result->nb_coordinates() == nb_space_dimensions() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Tetrahedron:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceTetrahedron::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Tetrahedron:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT(
         PEL::equal(
            MEASURE,
            PEL::abs( signed_measure( vertex(0), vertex(1),
                                      vertex( 2), vertex( 3 ) ) ) ) ) ;
      PEL_ASSERT( PEL::abs( PEL::abs(ORIENT)-1. ) < 1.e-8 ) ;
   }
   return( true ) ;
}
