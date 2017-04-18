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

#include <GE_Hexahedron.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_ReferenceCube.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <doubleArray2D.hh>
#include <doubleArray3D.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;

GE_Hexahedron const* GE_Hexahedron::MODEL = new GE_Hexahedron() ;

//-----------------------------------------------------------------------------
GE_Hexahedron:: GE_Hexahedron( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Hexahedron" )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Hexahedron:: GE_Hexahedron" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Hexahedron*
GE_Hexahedron:: create_replica( PEL_Object* a_owner,
                                PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Hexahedron* result = new GE_Hexahedron( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Hexahedron:: GE_Hexahedron( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Hexahedron:: GE_Hexahedron" ) ;
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Hexahedron:: ~GE_Hexahedron( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: ~GE_Hexahedron" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Hexahedron:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static std::string const result = "GE_Hexahedron" ;

   PEL_CHECK_POST( result == "GE_Hexahedron" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Hexahedron:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   PEL_CHECK( PEL::equal( MEASURE, computed_measure() ) ) ;

   double const result = MEASURE ;

   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Hexahedron:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = 0 ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Hexahedron:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   static GE_Point* pt_ref =
                 GE_Point:: create( PEL_Root::object(), (size_t) 3 ) ;
   PEL_CHECK( pt_ref!=0 && pt_ref->nb_coordinates()==3 ) ;

   //??? The Newton algorithm can not converge if the point is outside the
   //??? hexahedron - make first a ball inclusion test.
   bool is_in = pt->distance( center() )<inter_vertices_maximum_distance() ;
   if ( is_in )
   {
      apply_inverse_mapping( pt, pt_ref ) ;
      is_in = reference_polyhedron()->contains( pt_ref ) ;
   }
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Hexahedron:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceCube::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Hexahedron:: apply_mapping( GE_Point const* pt_ref,
                               GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;
   GE_Point const* V5 = vertex( 5 ) ;
   GE_Point const* V6 = vertex( 6 ) ;
   GE_Point const* V7 = vertex( 7 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   double Ngeom_0 = (1.0-x)*(1.0-y)*(1.0-z) ;  // VERTEX 0
   double Ngeom_1 =      x *(1.0-y)*(1.0-z) ;  // VERTEX 1
   double Ngeom_2 =      x *     y *(1.0-z) ;  // VERTEX 2
   double Ngeom_3 = (1.0-x)*     y *(1.0-z) ;  // VERTEX 3
   double Ngeom_4 = (1.0-x)*(1.0-y)*     z  ;  // VERTEX 4
   double Ngeom_5 =      x *(1.0-y)*     z  ;  // VERTEX 5
   double Ngeom_6 =      x *     y *     z  ;  // VERTEX 6
   double Ngeom_7 = (1.0-x)*     y *     z  ;  // VERTEX 7

   for( size_t d=0 ; d<3; ++d )
   {
      pt->set_coordinate( d, Ngeom_0 * V0->coordinate( d ) +
                             Ngeom_1 * V1->coordinate( d ) +
                             Ngeom_2 * V2->coordinate( d ) +
                             Ngeom_3 * V3->coordinate( d ) +
                             Ngeom_4 * V4->coordinate( d ) +
                             Ngeom_5 * V5->coordinate( d ) +
                             Ngeom_6 * V6->coordinate( d ) +
                             Ngeom_7 * V7->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Hexahedron:: apply_inverse_mapping( GE_Point const* pt,
                                       GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double error_bound_required = 0.1*reference_polyhedron()->epsilon() ;
   size_t const newton_max_iteration = 500 ;

   static size_t nb_sp_dims = 3 ;
   static GE_Matrix* mat = GE_Matrix::create( PEL_Root::object(),
                                              nb_sp_dims, nb_sp_dims ) ;
   static doubleArray2D vec1( 1, nb_sp_dims ) ;
   static doubleArray2D vec2( 1, nb_sp_dims ) ;
   static GE_Point* pt_new = GE_Point::create( PEL_Root::object(), nb_sp_dims ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;

   // Initial guess: same as in GE_Cuboid
   double dot10P=0., dot30P=0., dot40P=0., norm01=0., norm03=0., norm04=0. ;
   for( size_t d=0 ; d<nb_sp_dims ; d++ )
   {
      double v0d  = V0->coordinate( d ) ;
      double v0v1 = V1->coordinate( d ) - v0d ;
      double v0v3 = V3->coordinate( d ) - v0d ;
      double v0v4 = V4->coordinate( d ) - v0d ;
      norm01 += v0v1 * v0v1 ;
      norm03 += v0v3 * v0v3 ;
      norm04 += v0v4 * v0v4 ;

      double x = pt->coordinate(d) ;
      dot10P += v0v1 * ( x - v0d ) ;
      dot30P += v0v3 * ( x - v0d ) ;
      dot40P += v0v4 * ( x - v0d ) ;
   }

   pt_ref->set_coordinate( (size_t)0, dot10P/norm01 ) ;
   pt_ref->set_coordinate( (size_t)1, dot30P/norm03 ) ;
   pt_ref->set_coordinate( (size_t)2, dot40P/norm04 ) ;

   // Newton iterations start here
   size_t iter = 0 ;
   bool cv = false ;
   while( !cv )
   {
      //??? it may be faster if the mapping and the mapping derivative
      //??? are computed in the same function
      apply_mapping( pt_ref, pt_new ) ;
      build_mapping_derivative( pt_ref, mat ) ;

      for( size_t d=0 ; d<nb_sp_dims ; ++d )
      {
         vec1(0,d) = pt->coordinate( d ) - pt_new->coordinate( d ) ;
      }

      // Solve for a new increment
      mat->compute_determinant() ;
      mat->invert( vec1, vec2 ) ;

      // Update
      double err = 0. ;
      for( size_t iDim=0 ; iDim<nb_sp_dims ; ++iDim )
      {
         err = PEL::max( err, PEL::abs( vec2(0,iDim) ) ) ;
         double const x = pt_ref->coordinate(iDim)+vec2(0,iDim) ;
         pt_ref->set_coordinate( iDim, x ) ;
      }

      // Check convergence
      cv = err<error_bound_required ;
      if( iter>=newton_max_iteration )
      {
         std::ostringstream message ;
         message << "Unable to apply inverse mapping." << std::endl ;
         print( message, 0 ) ;
         message << std::endl ;
         message << "pre_image: " ;
         pt->print( message, 0 ) ;
         message << std::endl ;
         PEL_Error::object()->raise_internal( message.str() ) ;
      }
      ++iter ;
   }

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Hexahedron:: build_mapping_derivative( GE_Point const* pt_ref,
                                          GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   doubleArray2D const& dN = dNgeom( pt_ref ) ;

   size_t const ref_dim = dimension() ;
   size_t const nb_sp_dims = nb_space_dimensions() ;

   // the access to the d-th coordinate of the i-th vertex is slow

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;
   GE_Point const* V5 = vertex( 5 ) ;
   GE_Point const* V6 = vertex( 6 ) ;
   GE_Point const* V7 = vertex( 7 ) ;

   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      double v1d = V1->coordinate( d ) ;
      double v2d = V2->coordinate( d ) ;
      double v3d = V3->coordinate( d ) ;
      double v4d = V4->coordinate( d ) ;
      double v5d = V5->coordinate( d ) ;
      double v6d = V6->coordinate( d ) ;
      double v7d = V7->coordinate( d ) ;

      for( size_t j=0 ; j<ref_dim ; ++j )
      {
         jac->set_item( d, j, dN( 0, j ) * v0d +
                              dN( 1, j ) * v1d +
                              dN( 2, j ) * v2d +
                              dN( 3, j ) * v3d +
                              dN( 4, j ) * v4d +
                              dN( 5, j ) * v5d +
                              dN( 6, j ) * v6d +
                              dN( 7, j ) * v7d ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
GE_Hexahedron:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   doubleArray2D const& dN = dNgeom( pt_ref ) ;

   size_t const ref_dim = dimension() ;
   size_t const nb_sp_dims = nb_space_dimensions() ;

   // the access to the d-th coordinate of the i-th vertex is slow

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;
   GE_Point const* V4 = vertex( 4 ) ;
   GE_Point const* V5 = vertex( 5 ) ;
   GE_Point const* V6 = vertex( 6 ) ;
   GE_Point const* V7 = vertex( 7 ) ;

   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      double v1d = V1->coordinate( d ) ;
      double v2d = V2->coordinate( d ) ;
      double v3d = V3->coordinate( d ) ;
      double v4d = V4->coordinate( d ) ;
      double v5d = V5->coordinate( d ) ;
      double v6d = V6->coordinate( d ) ;
      double v7d = V7->coordinate( d ) ;

      for( size_t j=0 ; j<ref_dim ; ++j )
      {
         tjac->set_item( j, d, dN( 0, j ) * v0d +
                               dN( 1, j ) * v1d +
                               dN( 2, j ) * v2d +
                               dN( 3, j ) * v3d +
                               dN( 4, j ) * v4d +
                               dN( 5, j ) * v5d +
                               dN( 6, j ) * v6d +
                               dN( 7, j ) * v7d ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
GE_Hexahedron:: build_mapping_hessian( GE_Point const* pt_ref,
                                       doubleArray3D* hessian,
                                       bool& nonzero_hessian ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: build_mapping_hessian" ) ;
   PEL_CHECK_PRE( build_mapping_hessian_PRE( pt_ref,
                                             hessian, nonzero_hessian ) ) ;

   //??? it would be faster to implement it directly, as for build_mapping
   //??? and build_(tr_)mapping_derivative
   geobfs_2_mapping_hessian( d2Ngeom( pt_ref ), hessian ) ;
   nonzero_hessian = true ;
}

//-----------------------------------------------------------------------------
void
GE_Hexahedron:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
}

//-----------------------------------------------------------------------------
double
GE_Hexahedron:: computed_measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;
   
   static GE_Matrix* mat = GE_Matrix::create( PEL_Root::object(), 3, 3 ) ;
   static GE_QuadratureRule const* qr =
                               GE_QuadratureRule::object( "GE_Cube_QR3" ) ;
   double vol = 0. ;
   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      build_tr_mapping_derivative( qr->point( i ), mat ) ;
      mat->compute_determinant() ;
      vol += qr->weight(i)*PEL::abs( mat->determinant() ) ;
   }
   return( vol ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Hexahedron:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Hexahedron:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   bool ok = ( nb_space_dimensions()>=3 ) ;

   //??? IL MANQUE DES PROPRIETES SUR LES SOMMETS :
   //???    - les faces doivent etre valides ;
   //???    - elles ne doivent pas se croiser ;
   //???    - ...
   return( ok ) ;
}

//---------------------------------------------------------------------
doubleArray2D const&
GE_Hexahedron:: dNgeom( GE_Point const* pt_ref ) const
//---------------------------------------------------------------------
{
   static doubleArray2D result( 8, 3 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   // VERTEX 0 :
   result( 0, 0 )  =  -(1.0-y)*(1.0-z) ;   // derivative with respect to x
   result( 0, 1 )  =  -(1.0-x)*(1.0-z) ;   // derivative with respect to y
   result( 0, 2 )  =  -(1.0-x)*(1.0-y) ;   // derivative with respect to z

   // VERTEX 1 :
   result( 1, 0 )  =   (1.0-y)*(1.0-z) ;   // derivative with respect to x
   result( 1, 1 )  =  -x*(1.0-z) ;         // derivative with respect to y
   result( 1, 2 )  =  -x*(1.0-y) ;         // derivative with respect to z

   // VERTEX 2 :
   result( 2, 0 )  =   y*(1.0-z) ;         // derivative with respect to x
   result( 2, 1 )  =   x*(1.0-z) ;         // derivative with respect to y
   result( 2, 2 )  =  -x*y ;               // derivative with respect to z

   // VERTEX 3 :
   result( 3, 0 )  =  -y*(1.0-z) ;         // derivative with respect to x
   result( 3, 1 ) =   (1.0-x)*(1.0-z) ;    // derivative with respect to y
   result( 3, 2 ) =  -(1.0-x)*y ;          // derivative with respect to z

   // VERTEX 4 :
   result( 4, 0 ) =  -(1.0-y)*z ;          // derivative with respect to x
   result( 4, 1 ) =  -(1.0-x)*z ;          // derivative with respect to y
   result( 4, 2 ) =   (1.0-x)*(1.0-y) ;    // derivative with respect to z

   // VERTEX 5 :
   result( 5, 0 ) =   (1.0-y)*z ;          // derivative with respect to x
   result( 5, 1 ) =  -x*z ;                // derivative with respect to y
   result( 5, 2 ) =   x*(1.0-y) ;          // derivative with respect to z

   // VERTEX 6 :
   result( 6, 0 ) =   y*z ;                // derivative with respect to x
   result( 6, 1 ) =   x*z ;                // derivative with respect to y
   result( 6, 2 ) =   x*y ;                // derivative with respect to z

   // VERTEX 7 :
   result( 7, 0 ) =  -y*z ;                // derivative with respect to x
   result( 7, 1 ) =   (1.0-x)*z ;          // derivative with respect to y
   result( 7, 2 ) =   (1.0-x)*y ;          // derivative with respect to z

   return( result ) ;
}

//---------------------------------------------------------------------
doubleArray3D const&
GE_Hexahedron:: d2Ngeom( GE_Point const* pt_ref ) const
//---------------------------------------------------------------------
{
   static doubleArray3D result( 8, 3, 3 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   // VERTEX 0 :
   result( 0, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 0, 0, 1 ) =  (1.0-z) ;   // derivative with respect to x and y
   result( 0, 0, 2 ) =  (1.0-y) ;   // derivative with respect to x and z

   result( 0, 1, 0 ) =  (1.0-z) ;   // derivative with respect to y and x
   result( 0, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 0, 1, 2 ) =  (1.0-x) ;   // derivative with respect to y and z

   result( 0, 2, 0 ) =  (1.0-y) ;   // derivative with respect to z and x
   result( 0, 2, 1 ) =  (1.0-x) ;   // derivative with respect to z and y
   result( 0, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 1 :
   result( 1, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 1, 0, 1 ) = -(1.0-z) ;   // derivative with respect to x and y
   result( 1, 0, 2 ) = -(1.0-y) ;   // derivative with respect to x and z

   result( 1, 1, 0 ) = -(1.0-z) ;   // derivative with respect to y and x
   result( 1, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 1, 1, 2 ) =  x       ;   // derivative with respect to y and z

   result( 1, 2, 0 ) = -(1.0-y) ;   // derivative with respect to z and x
   result( 1, 2, 1 ) =  x       ;   // derivative with respect to z and y
   result( 1, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 2 :
   result( 2, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 2, 0, 1 ) =  (1.0-z) ;   // derivative with respect to x and y
   result( 2, 0, 2 ) = -y       ;   // derivative with respect to x and z

   result( 2, 1, 0 ) =  (1.0-z) ;   // derivative with respect to y and x
   result( 2, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 2, 1, 2 ) = -x       ;   // derivative with respect to y and z

   result( 2, 2, 0 ) = -y       ;   // derivative with respect to z and x
   result( 2, 2, 1 ) = -x       ;   // derivative with respect to z and y
   result( 2, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 3 :
   result( 3, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 3, 0, 1 ) = -(1.0-z) ;   // derivative with respect to x and y
   result( 3, 0, 2 ) =  y       ;   // derivative with respect to x and z

   result( 3, 1, 0 ) = -(1.0-z) ;   // derivative with respect to y and x
   result( 3, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 3, 1, 2 ) = -(1.0-x) ;   // derivative with respect to y and z

   result( 3, 2, 0 ) =  y       ;   // derivative with respect to z and x
   result( 3, 2, 1 ) = -(1.0-x) ;   // derivative with respect to z and y
   result( 3, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 4 :
   result( 4, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 4, 0, 1 ) =  z       ;   // derivative with respect to x and y
   result( 4, 0, 2 ) = -(1.0-y) ;   // derivative with respect to x and z

   result( 4, 1, 0 ) =  z       ;   // derivative with respect to y and x
   result( 4, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 4, 1, 2 ) = -(1.0-x) ;   // derivative with respect to y and z

   result( 4, 2, 0 ) = -(1.0-y) ;   // derivative with respect to z and x
   result( 4, 2, 1 ) = -(1.0-x) ;   // derivative with respect to z and y
   result( 4, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 5 :
   result( 5, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 5, 0, 1 ) = -z       ;   // derivative with respect to x and y
   result( 5, 0, 2 ) =  (1.0-y) ;   // derivative with respect to x and z

   result( 5, 1, 0 ) = -z       ;   // derivative with respect to y and x
   result( 5, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 5, 1, 2 ) = -x       ;   // derivative with respect to y and z

   result( 5, 2, 0 ) =  (1.0-y) ;   // derivative with respect to z and x
   result( 5, 2, 1 ) = -x       ;   // derivative with respect to z and y
   result( 5, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 6 :
   result( 6, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 6, 0, 1 ) =  z       ;   // derivative with respect to x and y
   result( 6, 0, 2 ) =  y       ;   // derivative with respect to x and z

   result( 6, 1, 0 ) =  z       ;   // derivative with respect to y and x
   result( 6, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 6, 1, 2 ) =  x       ;   // derivative with respect to y and z

   result( 6, 2, 0 ) =  y       ;   // derivative with respect to z and x
   result( 6, 2, 1 ) =  x       ;   // derivative with respect to z and y
   result( 6, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   // VERTEX 7 :
   result( 7, 0, 0 ) =  0.0     ;   // derivative with respect to x and x
   result( 7, 0, 1 ) = -z       ;   // derivative with respect to x and y
   result( 7, 0, 2 ) = -y       ;   // derivative with respect to x and z

   result( 7, 1, 0 ) = -z       ;   // derivative with respect to y and x
   result( 7, 1, 1 ) =  0.0     ;   // derivative with respect to y and y
   result( 7, 1, 2 ) =  (1.0-x) ;   // derivative with respect to y and z

   result( 7, 2, 0 ) = -y       ;   // derivative with respect to z and x
   result( 7, 2, 1 ) =  (1.0-x) ;   // derivative with respect to z and y
   result( 7, 2, 2 ) =  0.0     ;   // derivative with respect to z and z

   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Hexahedron:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Hexahedron:: reference_polyhedron_POST(
                          GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceCube::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Hexahedron:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   return( true ) ;
}

