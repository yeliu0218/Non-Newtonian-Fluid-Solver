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

#include <GE_Triangle.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>

#include <iostream>

GE_Triangle const* GE_Triangle::MODEL = new GE_Triangle() ;

GE_Vector* GE_Triangle::V0V1 = 0 ;
GE_Vector* GE_Triangle::V0V2 = 0 ;
GE_Vector* GE_Triangle::N    = 0 ;

//-----------------------------------------------------------------------------
GE_Triangle:: GE_Triangle( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Triangle" )
   , UNIT_NORMAL( 0 )
   , MEASURE( PEL::bad_double() )
   , ORIENT( false )
   , FV_CENTER( 0 )
{
   PEL_LABEL( "GE_Triangle:: GE_Triangle" ) ;

   V0V1 = GE_Vector::create( this, (size_t) 2 ) ;
   V0V2 = GE_Vector::create( this, (size_t) 2 ) ;
   N = GE_Vector::create( this, (size_t) 3 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Triangle*
GE_Triangle:: create_replica( PEL_Object* a_owner,
                              PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Triangle* result = new GE_Triangle( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Triangle:: GE_Triangle( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , UNIT_NORMAL( 0 )
   , MEASURE( PEL::bad_double() )
   , ORIENT( false )
   , FV_CENTER( 0 )
{
   PEL_LABEL( "GE_Triangle:: GE_Triangle" ) ;

   if( nb_space_dimensions() == 3 )
   {
      UNIT_NORMAL = GE_Vector::create( this, 3 ) ;
   }

   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Triangle:: ~GE_Triangle( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: ~GE_Triangle" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
      V0V1 = 0 ;
      V0V2 = 0 ;
      N = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Triangle:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: name" ) ;

   static std::string const result = "GE_Triangle" ;

   PEL_CHECK_POST( result == "GE_Triangle" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Triangle:: measure( GE_Point const* v0,
                       GE_Point const* v1,
                       GE_Point const* v2 )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: measure" ) ;
   PEL_CHECK_PRE( v0 != 0 ) ;
   PEL_CHECK_PRE( v1 != 0 && v1->nb_coordinates() == v0->nb_coordinates() ) ;
   PEL_CHECK_PRE( v2 != 0 && v2->nb_coordinates() == v0->nb_coordinates() ) ;

   double const result = PEL::abs( signed_measure( v0, v1, v2 ) ) ;

   PEL_CHECK_POST( result>=0. ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Triangle:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   PEL_CHECK( PEL::equal( MEASURE,
                          signed_measure( vertex(0), vertex(1), vertex(2) ) ) ) ;
   
   double result = PEL::abs( MEASURE ) ;
   
   PEL_CHECK_POST( measure_POST( result ) ) ;
   PEL_CHECK_POST(
      PEL::equal( result, measure( vertex(0), vertex(1), vertex(2) ) ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Triangle:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   if( FV_CENTER == 0 )
   {
      FV_CENTER = GE_Point::create( const_cast<GE_Triangle*>( this ),
                                    (size_t)2 ) ; // necessarily 2D
      V0V1->re_initialize( vertex(1), vertex(0) ) ;
      V0V2->re_initialize( vertex(2), vertex(0) ) ;
      compute_finite_volume_center( vertex(0), V0V1, V0V2, FV_CENTER ) ;
   }

   GE_Point const* result = FV_CENTER ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   //??? this implementation is questionable
   //???   ==> see preferably that of GE_Quadrilateral
   bool is_in = ( nb_space_dimensions()==2 ) ;
   if( !is_in )
   {
      PEL_CHECK( nb_space_dimensions()==3 ) ;
      double normP_Pproject = 0. ;
      for( size_t d=0 ; d<3 ; ++d )
      {
         normP_Pproject +=
                 UNIT_NORMAL->component(d)
                           *( pt->coordinate(d) - vertex(0)->coordinate(d) ) ;
      }
      double const epsilon_real =
              reference_polyhedron()->epsilon()*PEL::sqrt( measure() ) ;
      is_in = PEL::toler( normP_Pproject, epsilon_real ) ;
   }
   if( is_in )
   {
      static GE_Point* pt_ref =
                         GE_Point::create( PEL_Root::object(), (size_t) 2 ) ;
      PEL_CHECK( pt_ref!=0 &&
                 pt_ref->nb_coordinates()==dimension() ) ;
      apply_inverse_mapping( pt, pt_ref ) ;
      is_in = reference_polyhedron()->contains( pt_ref ) ;
   }
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_Triangle:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Vector const* result = UNIT_NORMAL ;

   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Triangle:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceTriangle::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Triangle:: apply_mapping( GE_Point const* pt_ref,
                             GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double umxmy = 1.0 - x - y ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, umxmy * V0->coordinate( d ) +
                                 x * V1->coordinate( d ) +
                                 y * V2->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
void
GE_Triangle:: apply_inverse_mapping( GE_Point const* v0,
                                     GE_Point const* v1,
                                     GE_Point const* v2,
                                     double v0v1v2_measure,
                                     GE_Point const* pt,
                                     GE_Point* pt_ref )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==2 || pt->nb_coordinates()==3 ) ;
   PEL_CHECK_PRE( v0 != 0 && v0->nb_coordinates()==pt->nb_coordinates() ) ;
   PEL_CHECK_PRE( v1 != 0 && v1->nb_coordinates()==pt->nb_coordinates() ) ;
   PEL_CHECK_PRE( v2 != 0 && v2->nb_coordinates()==pt->nb_coordinates() ) ;
   PEL_CHECK_PRE( v0v1v2_measure>0. ) ;
   PEL_CHECK_PRE( pt_ref !=0 && pt_ref->nb_coordinates()==2 ) ;

   double norm01=0., norm02=0., dot102=0., dot10P=0., dot20P=0. ;
   size_t const nb_sp_dims = pt->nb_coordinates() ;
   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double x = pt->coordinate(d) ;
      double v0d  = v0->coordinate( d ) ;
      double v0v1 = v1->coordinate( d ) - v0d ;
      double v0v2 = v2->coordinate( d ) - v0d ;
      norm01 += v0v1 * v0v1 ;
      norm02 += v0v2 * v0v2 ;
      dot10P += v0v1 * ( x - v0d ) ;
      dot20P += v0v2 * ( x - v0d ) ;
      dot102 += v0v1 * v0v2 ;
   }

   double const coef = 4.*v0v1v2_measure*v0v1v2_measure ;
   pt_ref->set_coordinate( (size_t)0,
                          (norm02*dot10P-dot102*dot20P)/coef ) ;
   pt_ref->set_coordinate( (size_t)1,
                          (norm01*dot20P-dot102*dot10P)/coef) ;
}

//-----------------------------------------------------------------------------
void
GE_Triangle:: apply_inverse_mapping( GE_Point const* pt,
                                     GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   apply_inverse_mapping( vertex( 0 ), vertex( 1 ), vertex( 2 ),
                          measure(),
                          pt, pt_ref ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//??? the inversion of a linear system whose matrix is jac is immediate
//??? --> may be used to save time
//----------------------------------------------------------------------------
void
GE_Triangle:: build_mapping_derivative( GE_Point const* pt_ref,
                                        GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      jac->set_item( d, 0, V1->coordinate( d ) - v0d ) ;
      jac->set_item( d, 1, V2->coordinate( d ) - v0d ) ;
   }
}

//----------------------------------------------------------------------------
void
GE_Triangle:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                           GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      tjac->set_item( 0, d, V1->coordinate( d ) - v0d ) ;
      tjac->set_item( 1, d, V2->coordinate( d ) - v0d ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Triangle:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   bool const first = ( MEASURE == PEL::bad_double() ) ;

   MEASURE = signed_measure( vertex(0), vertex(1), vertex(2) ) ;

   V0V1->re_initialize( vertex(1), vertex(0) ) ;
   V0V2->re_initialize( vertex(2), vertex(0) ) ;

   if( nb_space_dimensions()==3 )
   {
      UNIT_NORMAL->set_as_cross_product( V0V1, V0V2 ) ;
      double const n = 1./UNIT_NORMAL->norm() ;
      UNIT_NORMAL->scale( n ) ;
   }

   // Initial orientation:
   if( first ) ORIENT = ( MEASURE<0. ? false : true ) ;

   // le calcul de FV_CENTER n'est fait que si il avait deja été fait avant
   if( FV_CENTER != 0 )
      compute_finite_volume_center( vertex(0), V0V1, V0V2, FV_CENTER ) ;
}

//-----------------------------------------------------------------------------
double
GE_Triangle:: signed_measure( GE_Point const* v0,
                              GE_Point const* v1,
                              GE_Point const* v2 )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: signed_measure" ) ;
   PEL_CHECK( v0 != 0 ) ;
   PEL_CHECK( v1 != 0 && v1->nb_coordinates() == v0->nb_coordinates() ) ;
   PEL_CHECK( v2 != 0 && v2->nb_coordinates() == v0->nb_coordinates() ) ;

   size_t nb_dims = v0->nb_coordinates() ;
   V0V1->re_initialize( v1, v0 ) ;
   V0V2->re_initialize( v2, v0 ) ;
   N->set_as_cross_product( V0V1, V0V2 ) ;

   double const result = ( nb_dims==2 ? N->component(2)/2. : N->norm()/2. ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Triangle:: compute_finite_volume_center( GE_Point const* a,
                                            GE_Vector const* ab,
                                            GE_Vector const* ac,
                                            GE_Point* result )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: compute_finite_volume_center" ) ;
   PEL_CHECK( a != 0 && a->nb_coordinates() == 2 ) ;
   PEL_CHECK( ab != 0 && ab->nb_components() == 2 ) ;
   PEL_CHECK( ac != 0 && ac->nb_components() == 2 ) ;
   PEL_CHECK( result != 0 && result->nb_coordinates() == 2 ) ;

   // Let a, b, c be vectors in R^2.
   // Let ax, ay, be the components of a, and likewise for b, c.
   // Let |a| denote the Euclidean norm of a
   // The circumcenter m of triangle a,b,c is given by:
   //
   //           | by-ay  |b-a|^2 |
   //           | cy-ay  |c-a|^2 |
   // mx = ax - ------------------
   //             | bx-ax  by-ay |
   //           2 | cx-ax  cy-ay |
   //
   //           | bx-ax  |b-a|^2 |
   //           | cx-ax  |c-a|^2 |
   // my = ay + ------------------.
   //             | bx-ax  by-ay |
   //           2 | cx-ax  cy-ay |
   //
   // Reference: http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html

   double denom = ab->component( 0 ) * ac->component( 1 ) -
                  ac->component( 0 ) * ab->component( 1 ) ;

   double ab2 = ab->component( 0 ) * ab->component( 0 ) +
                ab->component( 1 ) * ab->component( 1 ) ;

   double ac2 = ac->component( 0 ) * ac->component( 0 ) +
                ac->component( 1 ) * ac->component( 1 ) ;

   double xx = ab->component( 1 ) * ac2 - ac->component( 1 ) * ab2 ;
   double yy = ab->component( 0 ) * ac2 - ac->component( 0 ) * ab2 ;

   result->set_coordinate( 0, a->coordinate( 0 ) - xx/2.0/denom ) ;
   result->set_coordinate( 1, a->coordinate( 1 ) + yy/2.0/denom ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool const o = ( MEASURE<0. ? false : true ) ;
   bool ok = ( ORIENT == o && nb_space_dimensions()>=2 ) ;

   return( ok ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle:: finite_volume_center_POST( GE_Point const* result ) const
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
GE_Triangle:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceTriangle::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      if( nb_space_dimensions()==3 )
      {
         PEL_ASSERT( UNIT_NORMAL!=0 ) ;
         PEL_ASSERT( UNIT_NORMAL->nb_components()==3 ) ;
         double dot01 = 0., dot02 = 0. ;
         for( size_t d=0; d<3; d++ )
         {
            dot01 += UNIT_NORMAL->component(d)*
               ( vertex(1)->coordinate(d)-vertex(0)->coordinate(d) ) ;
            dot02 += UNIT_NORMAL->component(d)*
               ( vertex(2)->coordinate(d)-vertex(0)->coordinate(d) ) ;
         }
         PEL_ASSERT( PEL::toler( dot01 ) ) ;
         PEL_ASSERT( PEL::toler( dot02 ) ) ;
      }
      else
      {
         PEL_ASSERT( UNIT_NORMAL == 0 ) ;
      }
   }
   return( true ) ;
}
