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

#include <GE_Rectangle.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <iostream>

GE_Rectangle const* GE_Rectangle::MODEL = new GE_Rectangle() ;

//----------------------------------------------------------------------------
GE_Rectangle:: GE_Rectangle( void )
//----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Rectangle" )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0 )
{
   PEL_LABEL( "GE_Rectangle:: GE_Rectangle" ) ;
}

//----------------------------------------------------------------------------
GE_Rectangle*
GE_Rectangle:: create_replica( PEL_Object* a_owner,
                               PEL_Vector* vertices ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Rectangle* result = new GE_Rectangle( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Rectangle:: GE_Rectangle( PEL_Object* a_owner, PEL_Vector* vertices )
//----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Rectangle:: GE_Rectangle" ) ;

   if( nb_space_dimensions()==3 )
   {
      UNIT_NORMAL = GE_Vector::create( this, 3 ) ;
   }
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
GE_Rectangle:: ~GE_Rectangle( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: ~GE_Rectangle" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//----------------------------------------------------------------------------
std::string const&
GE_Rectangle:: name( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static std::string const result = "GE_Rectangle" ;

   PEL_CHECK_POST( result == "GE_Rectangle" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Rectangle:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double result = MEASURE ;

   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Rectangle:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = center() ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
GE_Rectangle:: contains( GE_Point const* pt ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: contains" ) ;
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
              reference_polyhedron()->epsilon()*PEL::sqrt( MEASURE ) ;
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
GE_Rectangle:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   GE_Vector const* result = UNIT_NORMAL ;
   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Rectangle:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceSquare::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
GE_Rectangle:: apply_mapping( GE_Point const* pt_ref,
                              GE_Point* pt  ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double umxmy = 1.0 - x - y ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, umxmy * V0->coordinate( d ) +
                                 x * V1->coordinate( d ) +
		   	                     y * V3->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_Rectangle:: apply_inverse_mapping( GE_Point const* pt,
                                      GE_Point* pt_ref ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   double dot10P=0., dot30P=0., norm01=0., norm03=0. ;
   size_t dim = nb_space_dimensions() ;
   for( size_t d=0; d<dim; ++d )
   {
      double v0d  = V0->coordinate( d ) ;
      double v0v1 = V1->coordinate( d ) - v0d ;
      double v0v3 = V3->coordinate( d ) - v0d ;
      norm01 += v0v1 * v0v1 ;
      norm03 += v0v3 * v0v3 ;

      double x = pt->coordinate( d ) ;
      dot10P += v0v1 * ( x - v0d ) ;
      dot30P += v0v3 * ( x - v0d ) ;
   }

   pt_ref->set_coordinate( (size_t)0, dot10P/norm01 ) ;
   pt_ref->set_coordinate( (size_t)1, dot30P/norm03 ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//??? the inversion of a linear system whose matrix is jac is immediate
//??? --> may be used to save time
//----------------------------------------------------------------------------
void
GE_Rectangle:: build_mapping_derivative( GE_Point const* pt_ref,
                                         GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      jac->set_item( d, 0, V1->coordinate( d ) - v0d ) ;
      jac->set_item( d, 1, V3->coordinate( d ) - v0d ) ;
   }
}

//??? the inversion of a linear system whose matrix is tjac is immediate
//??? --> may be used to save time
//----------------------------------------------------------------------------
void
GE_Rectangle:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                            GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      tjac->set_item( 0, d, V1->coordinate( d ) - v0d ) ;
      tjac->set_item( 1, d, V3->coordinate( d ) - v0d ) ;
   }
}

//----------------------------------------------------------------------------
void
GE_Rectangle:: update_internal( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
   if( nb_space_dimensions()==3 )
   {
      GE_Vector* V0V1 = GE_Vector::create( 0, vertex(1), vertex(0) ) ;
      GE_Vector* V0V3 = GE_Vector::create( 0, vertex(3), vertex(0) ) ;
      UNIT_NORMAL->set_as_cross_product( V0V1, V0V3 ) ;
      double const inv_norm = 1./UNIT_NORMAL->norm() ;
      UNIT_NORMAL->scale( inv_norm ) ;
      V0V1->destroy() ; V0V1 = 0 ;
      V0V3->destroy() ; V0V3 = 0 ;
   }
}

//----------------------------------------------------------------------------
double
GE_Rectangle:: computed_measure( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;

   double const area =
      vertex(0)->distance( vertex(1) )*vertex(0)->distance( vertex(3) ) ;
   return( area ) ;
}

//----------------------------------------------------------------------------
bool
GE_Rectangle:: is_consistent( std::ostream& os, bool verbose ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Rectangle:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool success = ( nb_space_dimensions()>=2 ) ;
   if( verbose )
      display_check( os, "nb space dimensions greater than 2", success ) ;
   bool ok = success ;

   GE_Vector* V0V1 = GE_Vector::create( 0, vertex(1), vertex(0) ) ;
   double norm_V0V1 = V0V1->norm() ;
   success = ( norm_V0V1 > 10.e-50 ) ;
   if( verbose )
      display_check( os, "length of V0V1 greater than 1.E-50", success ) ;
   ok = ok && success ;

   GE_Vector* V0V3 = GE_Vector::create( 0, vertex(3), vertex(0) ) ;
   double norm_V0V3 = V0V3->norm() ;
   success = ( norm_V0V3 > 10.e-50 ) ;
   if( verbose )
      display_check( os, "length of V0V3 greater than 1.E-50", success ) ;
   ok = ok && success ;

   GE_Vector* V2V3 = GE_Vector::create( 0, vertex(3), vertex(2) ) ;
   double norm_V2V3 = V2V3->norm() ;
   success = ( norm_V2V3 > 10.e-50 ) ;
   if( verbose )
      display_check( os, "length of V2V3 greater than 1.E-50", success ) ;
   ok = ok && success ;

   double xx = V0V1->dot_product( V0V3 ) / norm_V0V3 / norm_V0V1 ;
   success = ( PEL::abs( xx ) < 1.0e-8 ) ;
   if( verbose ) display_check( os, "V0V1 orthogonal to V0V3", success ) ;
   ok = ok && success ;

   xx = V0V3->dot_product( V2V3 ) / norm_V0V3 / norm_V2V3 ;
   success = ( PEL::abs( xx ) < 1.0e-8 ) ;
   if( verbose ) display_check( os, "V2V3 orthogonal to V0V3", success ) ;
   ok = ok && success ;

   xx = norm_V0V1 / norm_V2V3 ;
   success = ( ( xx > (1.0-1.0e-8) ) && ( xx < (1.0+1.0e-8) ) ) ;
   if( verbose ) display_check( os, "V0V1 and V2V3 of same length", success ) ;
   ok = ok && success ;

   if( nb_space_dimensions()==3 )
   {
      static GE_Vector* cros = GE_Vector::create( PEL_Root::object(),
                                                  (size_t)3 ) ;
      cros->set_as_cross_product( V0V1, V0V3 ) ;
      xx = cros->dot_product( V2V3 ) / cros->norm() / norm_V2V3 ;
      success = ( PEL::abs( xx ) < 1.0e-8 ) ;
      if( verbose )
         display_check( os, "Flatness", success ) ;
      ok = ok && success ;
   }

   V0V1->destroy() ;
   V0V3->destroy() ;
   V2V3->destroy() ;

   return( ok ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Rectangle:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == center() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Rectangle:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceSquare::object() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
GE_Rectangle:: invariant( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT( PEL::equal( MEASURE, computed_measure() ) ) ;
      if( nb_space_dimensions()==3 )
      {
         GE_Vector* V0V1 =
            GE_Vector::create( 0, vertex( 1 ), vertex( 0 ) ) ;
         GE_Vector* V0V3 =
            GE_Vector::create( 0, vertex( 3 ), vertex( 0 ) ) ;
         PEL_ASSERT( UNIT_NORMAL!=0 &&
                     UNIT_NORMAL->nb_components()==3 ) ;
         PEL_ASSERT( PEL::toler( UNIT_NORMAL->dot_product( V0V1 ) ) ) ;
         PEL_ASSERT( PEL::toler( UNIT_NORMAL->dot_product( V0V3 ) ) ) ;
         V0V3->destroy() ;
         V0V1->destroy() ;
      }
      else
      {
         PEL_ASSERT( UNIT_NORMAL==0 ) ;
      }
   }
   return( true ) ;
}
