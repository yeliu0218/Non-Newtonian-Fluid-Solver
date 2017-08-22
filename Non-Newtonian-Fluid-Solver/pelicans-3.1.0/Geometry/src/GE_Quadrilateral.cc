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

#include <GE_Quadrilateral.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
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

GE_Quadrilateral const* GE_Quadrilateral::MODEL = new GE_Quadrilateral() ;

//-----------------------------------------------------------------------------
GE_Quadrilateral:: GE_Quadrilateral( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Quadrilateral" )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0 )
{
   PEL_LABEL( "GE_Quadrilateral:: GE_Quadrilateral" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Quadrilateral*
GE_Quadrilateral:: create_replica( PEL_Object* a_owner,
                                   PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Quadrilateral* result = new GE_Quadrilateral( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Quadrilateral:: GE_Quadrilateral( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0 )
{
   PEL_LABEL( "GE_Quadrilateral:: GE_Quadrilateral" ) ;

   if( nb_space_dimensions()==3 )
   {
      UNIT_NORMAL = GE_Vector::create( this, 3 ) ;
   }
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Quadrilateral:: ~GE_Quadrilateral( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: ~GE_Quadrilateral" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Quadrilateral:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static std::string const result = "GE_Quadrilateral" ;

   PEL_CHECK_POST( result == "GE_Quadrilateral" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Quadrilateral:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   PEL_CHECK( PEL::equal( MEASURE, computed_measure() ) ) ;

   double result = MEASURE ;

   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Quadrilateral:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = 0 ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_ReferencePolyhedron const* ref_poly = reference_polyhedron() ;

   // The Newton algorithm may not converge if the point is outside the
   // hexahedron - make first a ball inclusion test.
   double diam = inter_vertices_maximum_distance() ;
   bool result = pt->distance( center() )<diam ;
   if( result )
   {
      static GE_Point* pt_ref =
                       GE_Point::create( PEL_Root::object(), (size_t) 2 ) ;
      apply_inverse_mapping( pt, pt_ref ) ;
      result = ref_poly->contains( pt_ref ) ;
      if( result && nb_space_dimensions() == 3 )
      {
         static GE_Point* pt_new =
                          GE_Point::create( PEL_Root::object(), (size_t) 3 ) ;
         apply_mapping( pt_ref, pt_new ) ;
         result = ( pt->distance( pt_new )/diam ) < ref_poly->epsilon() ;
      }
   }

   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_Quadrilateral:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Vector const* result = UNIT_NORMAL ;

   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Quadrilateral:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;

   static GE_ReferencePolyhedron const* result = GE_ReferenceSquare::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Quadrilateral:: apply_mapping( GE_Point const* pt_ref,
                                  GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double Ngeom_0 = (1.0-x)*(1.0-y) ;  //--- VERTEX 0
   double Ngeom_1 = x*(1.0-y)       ;  //--- VERTEX 1
   double Ngeom_2 = x*y             ;  //--- VERTEX 2
   double Ngeom_3 = (1.0-x)*y       ;  //--- VERTEX 3

   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, Ngeom_0 * V0->coordinate( d ) +
                             Ngeom_1 * V1->coordinate( d ) +
                             Ngeom_2 * V2->coordinate( d ) +
                             Ngeom_3 * V3->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Quadrilateral:: apply_inverse_mapping( GE_Point const* pt,
                                          GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   size_t nbsp = nb_space_dimensions() ;

   if ( nbsp == 2 )
   {
      apply_inverse_mapping_2D( pt, pt_ref ) ;
   }
   else
   {
      PEL_CHECK( nb_space_dimensions()==3 ) ;
      apply_inverse_mapping_3D( pt, pt_ref ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Quadrilateral:: build_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;

   doubleArray2D const& dN = dNgeom( pt_ref ) ;

   size_t const ref_dim = dimension() ;
   size_t const nb_sp_dims = nb_space_dimensions() ;

   // the access to the d-th coordinate of the i-th vertex is slow

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      double v1d = V1->coordinate( d ) ;
      double v2d = V2->coordinate( d ) ;
      double v3d = V3->coordinate( d ) ;

      for( size_t j=0 ; j<ref_dim ; ++j )
      {
         jac->set_item( d, j, dN( 0, j ) * v0d +
                              dN( 1, j ) * v1d +
                              dN( 2, j ) * v2d +
                              dN( 3, j ) * v3d ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
GE_Quadrilateral:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;

   doubleArray2D const& dN = dNgeom( pt_ref ) ;

   size_t const ref_dim = dimension() ;
   size_t const nb_sp_dims = nb_space_dimensions() ;

   // the access to the d-th coordinate of the i-th vertex is slow

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V2 = vertex( 2 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double v0d = V0->coordinate( d ) ;
      double v1d = V1->coordinate( d ) ;
      double v2d = V2->coordinate( d ) ;
      double v3d = V3->coordinate( d ) ;

      for( size_t j=0 ; j<ref_dim ; ++j )
      {
         tjac->set_item( j, d, dN( 0, j ) * v0d +
                               dN( 1, j ) * v1d +
                               dN( 2, j ) * v2d +
                               dN( 3, j ) * v3d ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
GE_Quadrilateral:: build_mapping_hessian( GE_Point const* pt_ref,
                                          doubleArray3D* hessian,
                                          bool& nonzero_hessian ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: build_mapping_hessian" ) ;
   PEL_CHECK_PRE( build_mapping_hessian_PRE( pt_ref,
                                             hessian, nonzero_hessian ) ) ;

   //??? it would be faster to implement it directly, as for build_mapping
   //??? and build_(tr_)mapping_derivative
   geobfs_2_mapping_hessian( d2Ngeom( pt_ref ), hessian ) ;
   nonzero_hessian = true ;
}

//-----------------------------------------------------------------------------
void
GE_Quadrilateral:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
   if( nb_space_dimensions()==3 )
   {
      static GE_Point* c = GE_Point::create( PEL_Root::object(), (size_t)3 ) ;
      static GE_Vector* v0 = GE_Vector::create( PEL_Root::object(), 3 ) ;
      static GE_Vector* v1 = GE_Vector::create( PEL_Root::object(), 3 ) ;
      static GE_Vector* v =  GE_Vector::create( PEL_Root::object(), 3 ) ;
      for( size_t i=0 ; i<3 ; ++i )
      {
         double x = vertex(0)->coordinate(i) ;
         for( size_t j=1 ; j<3 ; ++j ) x += vertex(j)->coordinate(i) ;
         c->set_coordinate( i, 0.25*x ) ;
      }
      v0->re_initialize( c, vertex(0) ) ;
      v1->re_initialize( c, vertex(1) ) ;
      UNIT_NORMAL->set_as_cross_product( v0, v1 ) ;
      v0->re_initialize( c, vertex(2) ) ;
      v->set_as_cross_product( v1, v0 ) ;
      UNIT_NORMAL->sum( 1., v ) ;
      v1->re_initialize( c, vertex(3) ) ;
      v->set_as_cross_product( v0, v1 ) ;
      UNIT_NORMAL->sum( 1., v ) ;
      v0->re_initialize( c, vertex(0) ) ;
      v->set_as_cross_product( v1, v0 ) ;
      UNIT_NORMAL->sum( 1., v ) ;
      double const inv_norm = 1./UNIT_NORMAL->norm() ;
      UNIT_NORMAL->scale( inv_norm ) ;
   }
}

//-----------------------------------------------------------------------------
double
GE_Quadrilateral:: computed_measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;

   double norm01 = 0., norm03 = 0., scal0103 = 0. ;
   for( size_t d=0; d<nb_space_dimensions(); ++d )
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
   for( size_t d=0; d<nb_space_dimensions(); ++d )
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
GE_Quadrilateral:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: is_consistent" ) ;
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

//??????
//   double xx = 0. ;
//
//   if( nb_space_dimensions()==3 )
//   {
//      GE_Vector* cros = GE_Vector::create( 0, V0V1->nb_components() ) ;
//      cros->set_as_cross_product( V0V1, V0V3 ) ;
//      xx = cros->dot_product( V2V3 ) / cros->norm() / norm_V2V3 ;
//
//      success = ( xx < 1.E-8 ) ;
//      if( verbose )
//         display_check( os, "Flatness", success ) ;
//      ok = ok && success ;
//      cros->destroy() ;
//   }
   //??????????????????????????

   V0V1->destroy() ;
   V0V3->destroy() ;
   V2V3->destroy() ;

   if ( nb_space_dimensions() == 2 )
      success = is_convex2D() ;
   else if ( nb_space_dimensions() == 3 )
      success = is_convex3D() ;
   if( verbose )
      display_check( os, "Convexity", success ) ;
   ok = ok && success ;

   return( ok ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: is_convex2D( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: is_convex2D" ) ;
   PEL_CHECK( !is_prototype() ) ;
   PEL_CHECK( nb_space_dimensions()==2 ) ;

   double const x0 = vertex( 0 )->coordinate( 0 ) ;
   double const y0 = vertex( 0 )->coordinate( 1 ) ;
   double const x1 = vertex( 1 )->coordinate( 0 ) ;
   double const y1 = vertex( 1 )->coordinate( 1 ) ;
   double const x2 = vertex( 2 )->coordinate( 0 ) ;
   double const y2 = vertex( 2 )->coordinate( 1 ) ;
   double const x3 = vertex( 3 )->coordinate( 0 ) ;
   double const y3 = vertex( 3 )->coordinate( 1 ) ;
   // Mapping(x,y) = Sxy + Ux + Qy + K
   double const xs = x0+x2-x1-x3 ;
   double const ys = y0+y2-y1-y3 ;
   double const xu = x1-x0 ;
   double const yu = y1-y0 ;
   double const xq = x3-x0 ;
   double const yq = y3-y0 ;

   // Jac = [(Sy+U) (Sx+Q)]
   // det Jac = (Sy+U)^(Sx+Q)
   double const u_vec_q = xu*yq - yu*xq ;
   double const s_vec_q = xs*yq - ys*xq ;
   double const u_vec_s = xu*ys - yu*xs ;
   double const eps = 1.0e-12 ;
   bool result ;
   if( PEL::abs(u_vec_q)>eps )
   {
      // det J=0 <=> ax + by = 1
      double const a = -u_vec_s/u_vec_q ;
      double const b = -s_vec_q/u_vec_q ;
      result = a<1. && b<1. && (a+b)<1. ;
   }
   else
   {
      double const u_scal_s = xu*xq + yu*yq ;
      result = u_scal_s<0. ;
   }


   PEL_CHECK_INV( invariant() ) ;
   return result ;
}

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: is_convex3D( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: is_convex3D" ) ;
   PEL_CHECK( !is_prototype() ) ;
   PEL_CHECK( nb_space_dimensions()==3 ) ;

   GE_Vector* vecAB = GE_Vector::create( 0, vertex( 1 ), vertex( 0 ) ) ;
   GE_Vector* vecAD = GE_Vector::create( 0, vertex( 3 ), vertex( 0 ) ) ;
   GE_Vector* vecCD = GE_Vector::create( 0, vertex( 3 ), vertex( 2 ) ) ;
   GE_Vector* vecCB = GE_Vector::create( 0, vertex( 1 ), vertex( 2 ) ) ;
   GE_Vector* ABxAD = GE_Vector::create( 0, 3 ) ;
   GE_Vector* CDxCB = GE_Vector::create( 0, 3 ) ;
   GE_Vector* ADxCD = GE_Vector::create( 0, 3 ) ;
   GE_Vector* CBxAB = GE_Vector::create( 0, 3 ) ;

   ABxAD->set_as_cross_product( vecAB, vecAD ) ;
   ABxAD->scale( 1./ABxAD->norm() ) ;
   CDxCB->set_as_cross_product( vecCD, vecCB ) ;
   CDxCB->scale( 1./CDxCB->norm() ) ;
   ADxCD->set_as_cross_product( vecAD, vecCD ) ;
   ADxCD->scale( 1./ADxCD->norm() ) ;
   CBxAB->set_as_cross_product( vecCB, vecAB ) ;
   CBxAB->scale( 1./CBxAB->norm() ) ;

   bool result ;

   result = (    ( ABxAD->dot_product( UNIT_NORMAL ) * CBxAB->dot_product( UNIT_NORMAL ) >= -1.E-08 )
         && ( CBxAB->dot_product( UNIT_NORMAL ) * CDxCB->dot_product( UNIT_NORMAL ) >= -1.E-08 )
         && ( CDxCB->dot_product( UNIT_NORMAL ) * ADxCD->dot_product( UNIT_NORMAL ) >= -1.E-08 ) ) ;

   vecAB->destroy() ;
   vecAD->destroy() ;
   vecCD->destroy() ;
   vecCB->destroy() ;
   ABxAD->destroy() ;
   CDxCB->destroy() ;
   ADxCD->destroy() ;
   CBxAB->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
   return result ;
}

//---------------------------------------------------------------------
doubleArray2D const&
GE_Quadrilateral:: dNgeom( GE_Point const* pt_ref ) const
//---------------------------------------------------------------------
{
   static doubleArray2D result( 4, 2 ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   //-- VERTEX 0
   result( 0, 0 )  = -(1.0-y)  ;    // derivative with respect to x
   result( 0, 1 )  = -(1.0-x)  ;    // derivative with respect to y

   //-- VERTEX 1
   result( 1, 0 )  =  (1.0-y)  ;    // derivative with respect to x
   result( 1, 1 )  = -x        ;    // derivative with respect to y

   //-- VERTEX 2
   result( 2, 0 )  =  y        ;    // derivative with respect to x
   result( 2, 1 )  =  x        ;    // derivative with respect to y

   //-- VERTEX 3
   result( 3, 0 )  = -y        ;    // derivative with respect to x
   result( 3, 1 )  =  (1.0-x)  ;    // derivative with respect to y

   return( result ) ;
}

//---------------------------------------------------------------------
doubleArray3D const&
GE_Quadrilateral:: d2Ngeom( GE_Point const* pt_ref ) const
//---------------------------------------------------------------------
{
   static doubleArray3D result( 4, 2, 2 ) ;

   //-- VERTEX 0
   result( 0, 0, 0 ) =  0.0 ;    // derivative with respect to x and x
   result( 0, 0, 1 ) =  1.0 ;    // derivative with respect to x and x

   result( 0, 1, 0 ) =  1.0 ;    // derivative with respect to y and x
   result( 0, 1, 1 ) =  0.0 ;    // derivative with respect to y and y

   //-- VERTEX 1
   result( 1, 0, 0 ) =  0.0 ;    // derivative with respect to x and x
   result( 1, 0, 1 ) = -1.0 ;    // derivative with respect to x and y

   result( 1, 1, 0 ) = -1.0 ;    // derivative with respect to y and x
   result( 1, 1, 1 ) =  0.0 ;    // derivative with respect to y and y

   //-- VERTEX 2
   result( 2, 0, 0 ) =  0.0 ;    // derivative with respect to x and x
   result( 2, 0, 1 ) =  1.0 ;    // derivative with respect to x and y

   result( 2, 1, 0 ) =  1.0 ;    // derivative with respect to y and x
   result( 2, 1, 1 ) =  0.0 ;    // derivative with respect to y and y

   //-- VERTEX 3
   result( 3, 0, 0 ) =  0.0 ;    // derivative with respect to x and x
   result( 3, 0, 1 ) = -1.0 ;    // derivative with respect to x and y

   result( 3, 1, 0 ) = -1.0 ;    // derivative with respect to y and x
   result( 3, 1, 1 ) =  0.0 ;    // derivative with respect to y and y

   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Quadrilateral:: apply_inverse_mapping_2D( GE_Point const* pt,
                                             GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: apply_inverse_mapping_2D" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double const error_bound_required = 0.1*reference_polyhedron()->epsilon() ;

   size_t const newton_max_iteration = 50 ;

   double const x = pt->coordinate(0) ;
   double const y = pt->coordinate(1) ;

   double const x0 = vertex( 0 )->coordinate( 0 ) ;
   double const y0 = vertex( 0 )->coordinate( 1 ) ;
   double const x1 = vertex( 1 )->coordinate( 0 ) ;
   double const y1 = vertex( 1 )->coordinate( 1 ) ;
   double const x2 = vertex( 2 )->coordinate( 0 ) ;
   double const y2 = vertex( 2 )->coordinate( 1 ) ;
   double const x3 = vertex( 3 )->coordinate( 0 ) ;
   double const y3 = vertex( 3 )->coordinate( 1 ) ;

   double xi = 0.5, dxi = 0., eta = 0.5 , deta = 0. ;
   size_t iter = 0 ;
   bool fail = false ;
   do
   {
      double const ax_xi  = (1.-eta)*(x1-x0)+eta*(x2-x3) ;
      double const ay_xi  = (1.-eta)*(y1-y0)+eta*(y2-y3) ;
      double const ax_eta = (1. -xi)*(x3-x0)+ xi*(x2-x1) ;
      double const ay_eta = (1. -xi)*(y3-y0) +xi*(y2-y1) ;
      double const det = ( ax_xi*ay_eta-ax_eta*ay_xi ) ;
      double const fx =
         x - (1-xi-eta+xi*eta)*x0-(xi-xi*eta)*x1-xi*eta*x2-(eta-xi*eta)*x3 ;
      double const fy =
         y - (1-xi-eta+xi*eta)*y0-(xi-xi*eta)*y1-xi*eta*y2-(eta-xi*eta)*y3 ;
      dxi = ( fx*ay_eta - fy*ax_eta )/det ;
      deta = -( fx*ay_xi - fy*ax_xi )/det ;
      xi += dxi ; eta += deta ;
      iter++ ;
      fail = iter>newton_max_iteration ;
   } while( !fail &&
       PEL::abs(dxi)>error_bound_required &&
       PEL::abs(deta)>error_bound_required ) ;

   bool found = false ;
   if( !fail )
   {
      pt_ref->set_coordinate( 0, xi ) ;
      pt_ref->set_coordinate( 1, eta ) ;
      found = reference_polyhedron()->contains( pt_ref ) ;
   }
}

//---------------------------------------------------------------------
void
GE_Quadrilateral:: apply_inverse_mapping_3D( GE_Point const* pt,
                                               GE_Point* pt_ref ) const
//---------------------------------------------------------------------
{
   PEL_LABEL( "GE_Quadrilateral:: apply_inverse_mapping_3D" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   PEL_CHECK( dimension()<=nb_space_dimensions() ) ;
   PEL_CHECK( pt!=0 &&
              pt->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_CHECK( pt_ref!=0 &&
              pt_ref->nb_coordinates()==dimension() ) ;

   double error_bound_required = 0.1*reference_polyhedron()->epsilon() ;
   size_t newton_max_iteration = 500 ;

   static size_t nb_sp_dims = 3 ;
   static GE_Matrix* mat = GE_Matrix::create( PEL_Root::object(), 2, 2 ) ;
   static GE_Matrix* jac = GE_Matrix::create( PEL_Root::object(), nb_sp_dims, 2 ) ;
   static GE_Point* pt_new = GE_Point::create( PEL_Root::object(), (size_t)3 ) ;
   static doubleArray2D vec1( 1, 2 ) ;
   static doubleArray2D vec2( 1, 2 ) ;
   static doubleVector u1( nb_sp_dims ) ;
   static doubleVector u2( nb_sp_dims ) ;

   GE_Point const* V0 = vertex( 0 ) ;
   GE_Point const* V1 = vertex( 1 ) ;
   GE_Point const* V3 = vertex( 3 ) ;

   // Initial guess: same as in GE_Rectangle
   double dot10P=0., dot30P=0., norm01=0., norm03=0. ;
   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      double v0d  = V0->coordinate( d ) ;
      double v0v1 = V1->coordinate( d ) - v0d ;
      double v0v3 = V3->coordinate( d ) - v0d ;
      norm01 += v0v1 * v0v1 ;
      norm03 += v0v3 * v0v3 ;

      double x = pt->coordinate( d ) ;
      dot10P += v0v1 * ( x - v0d ) ;
      dot30P += v0v3 * ( x - v0d ) ;

      u1( d ) = v0v1 ;
      u2( d ) = v0v3 ;
   }
   pt_ref->set_coordinate( (size_t)0, dot10P/norm01 ) ;
   pt_ref->set_coordinate( (size_t)1, dot30P/norm03 ) ;
   for( size_t d=0 ; d<nb_sp_dims ; ++d )
   {
      u1( d ) /= norm01 ;
      u2( d ) /= norm03 ;
   }

   // Newton iterations start here
   size_t iter = 0 ;
   bool cv = false ;
   while( !cv )
   {
      //??? it may be faster if the mapping and the mapping derivative
      //??? are computed in the same function
      apply_mapping( pt_ref, pt_new ) ;
      build_mapping_derivative( pt_ref, jac ) ;
      for( size_t j=0 ; j<2 ; ++j )
      {
         double xx_0=0.0, xx_1=0.0 ;
         for( size_t d=0 ; d<nb_sp_dims ; ++d )
         {
            xx_0 += jac->item( d, j ) * u1( d ) ;
            xx_1 += jac->item( d, j ) * u2( d ) ;
         }
         mat->set_item( 0, j, xx_0 ) ;
         mat->set_item( 1, j, xx_1 ) ;
      }
      double bb_0=0.0, bb_1=0.0 ;
      for( size_t d=0 ; d<nb_sp_dims ; ++d )
      {
         bb_0 -= ( pt_new->coordinate(d) - pt->coordinate(d) ) * u1( d ) ;
         bb_1 -= ( pt_new->coordinate(d) - pt->coordinate(d) ) * u2( d ) ;
      }
      vec1( 0, 0 ) = bb_0 ;
      vec1( 0, 1 ) = bb_1 ;

      // Solve for a new increment
      mat->compute_determinant() ;
      mat->invert( vec1, vec2 ) ;

      // Update
      double err = 0. ;
      for( size_t iDim=0 ; iDim<2 ; ++iDim )
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

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: reference_polyhedron_POST(
                                   GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceSquare::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Quadrilateral:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      if( nb_space_dimensions()==3 )
      {
         PEL_ASSERT( UNIT_NORMAL!=0 &&
                     UNIT_NORMAL->nb_components()==3 ) ;
         PEL_ASSERT( PEL::equal( UNIT_NORMAL->norm(), 1. ) ) ;
      }
      else
      {
         PEL_ASSERT( UNIT_NORMAL==0 ) ;
      }
   }
   return( true ) ;
}
