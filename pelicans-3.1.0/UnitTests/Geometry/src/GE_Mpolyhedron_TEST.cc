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

#include <GE_Mpolyhedron_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_SimplePolygon2D.hh>
#include <GE_Vector.hh>
#include <GE_SetOfPoints.hh>

#include <iostream>

using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
GE_Mpolyhedron_TEST*
GE_Mpolyhedron_TEST:: REGISTRATOR = new GE_Mpolyhedron_TEST() ;
//-------------------------------------------------------------------------


//----------------------------------------------------------------------------
GE_Mpolyhedron_TEST:: GE_Mpolyhedron_TEST( void )
//----------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Mpolyhedron", "GE_Mpolyhedron_TEST" )
{
}

//----------------------------------------------------------------------------
GE_Mpolyhedron_TEST:: ~GE_Mpolyhedron_TEST( void )
//----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   GE_Mpolyhedron const* poly = create_polyhedron( this, exp ) ;
   adapt_to_polyhedron( poly ) ;

   test_generalities( exp, poly ) ;

   test_mapping( poly ) ;

   double d_eps = exp->double_data( "dbl_epsilon" ) ;
   double d_min = exp->double_data( "dbl_minimum" ) ;
   double dx = exp->double_data( "dx" ) ;
   test_mapping_derivative( poly, dx, d_eps, d_min ) ;

   test_contains( poly ) ;

   if( exp->has_entry( "test_if_fv_center_is_circumcenter" ) )
   {
      bool bb = exp->bool_data( "test_if_fv_center_is_circumcenter" ) ;
      if( bb ) test_if_fv_center_is_circumcenter( poly ) ;
   }
}

//-------------------------------------------------------------------------
GE_Mpolyhedron const*
GE_Mpolyhedron_TEST:: create_polyhedron( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "create_polyhedron_then_reset" ) ;

   size_t nb_sp_dims = exp->int_data( "nb_space_dimensions" ) ;
   PEL_ASSERT( nb_sp_dims>0 && nb_sp_dims<4 ) ;

   doubleVector const& vec = exp->doubleVector_data( "vertices_list" ) ;
   doubleVector coord( nb_sp_dims ) ;
   GE_SetOfPoints* vm = GE_SetOfPoints:: create( a_owner, nb_sp_dims ) ;
   for( size_t i=0 ; i<vec.size() ; i+=nb_sp_dims )
   {
      for( size_t j=0; j<nb_sp_dims; ++j )
      {
         coord( j ) = vec(i+j) ;
      }
      GE_Point* pt = GE_Point::create( 0, coord ) ;
      vm->append( pt ) ;
      pt->destroy() ; pt = 0 ;
   }
   size_t_vector vertices( vm->nb_points() ) ;
   for( size_t i=0 ; i<vm->nb_points() ; ++i )
   {
      vertices(i) = i ;
   }
   GE_Mpolyhedron const* poly =
                      GE_Mpolyhedron:: create( exp->string_data("name"),
                                               vm, vertices ) ;
   return( poly ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: adapt_to_polyhedron( GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------------
{
   size_t nb_sp_dims = poly->nb_space_dimensions() ;
   ORIGIN = GE_Point::origin( nb_sp_dims ) ;

   if( PT_REF == 0 ) PT_REF = GE_Point::create( this, poly->dimension() ) ;
   else PT_REF->copy( PT_REF->origin( poly->dimension() ) ) ;

   if( PT==0 ) PT = GE_Point::create( this, nb_sp_dims ) ;
   else PT->copy( ORIGIN ) ;

   if( PT2==0 ) PT2 = GE_Point::create( this, nb_sp_dims ) ;
   else PT2->copy( ORIGIN ) ;

   if( VECTOR==0 ) VECTOR = GE_Vector::create(  this, nb_sp_dims ) ;
   else VECTOR->re_initialize( ORIGIN, ORIGIN ) ;

   if( JACOB != 0 ) destroy_possession( JACOB ) ;
   JACOB = GE_Matrix::create( this, nb_sp_dims, poly->dimension() ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: test_generalities( PEL_ModuleExplorer const* exp,
                                         GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------------
{
   string mesg = exp->name() + ":" ;

   bool ok =  exp->string_data("name")==poly->name() ;
   notify_one_test_result( mesg+"name", ok ) ;

   ok = PEL::toler( poly->measure()-estimate_measure( poly ), 1.E-10 ) ;
   notify_one_test_result( "measure", ok ) ;
   if( !ok )
      out() <<"measure : "<<poly->measure()
            <<" estimated : " <<estimate_measure( poly ) <<std::endl ;

   ok = PEL::toler( poly->inter_vertices_maximum_distance()-
			  compute_inter_vert_max_dist( poly ), 1.E-10 ) ;
   notify_one_test_result( "inter_vertices_maximum_distance", ok ) ;

   ok = PEL::toler( poly->equivalent_ball_diameter()-
			  compute_equi_ball_diameter( poly ), 1.E-10 ) ;
   notify_one_test_result( "equivalent_ball_diameter", ok ) ;

   for( size_t i=0 ; i<poly->nb_space_dimensions() ; ++i )
   {
      ok = ok && PEL::toler( poly->inter_vertices_maximum_distance(i)-
                             compute_inter_vert_max_dist( poly, i ), 1.E-10 ) ;
   }

   notify_one_test_result( "inter_vertices_maximum_distance(dim)", ok ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: test_mapping( GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron_TEST:: test_mapping" ) ;

   GE_ReferencePolyhedron const* ref_poly = poly->reference_polyhedron() ;
   size_t nb_verts = poly->nb_vertices() ;

   bool ok = true ;
   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      poly->apply_mapping( ref_poly->vertex( i ), PT ) ;
      ok = ( PT->distance( poly->vertex( i ) ) < 1.e-10 ) ;
   }

   notify_one_test_result( "mapping(vertices)", ok ) ;

   ok = true ;
   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      poly->apply_inverse_mapping( poly->vertex( i ), PT_REF ) ;
      ok &= ( PT_REF->distance( ref_poly->vertex( i ) ) < 1.e-8 ) ;
   }

   notify_one_test_result( "inverse_mapping(vertices)", ok ) ;

   poly->apply_inverse_mapping( poly->center(), PT_REF ) ;
   ok = ( PT_REF->distance( ref_poly->center() ) < 1.e-8 ) ;
   notify_one_test_result( "inverse_mapping(center)", ok ) ;

   ok = true ;
   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_5" ) ;
   GE_QuadratureRule const* qr = qrp->quadrature_rule( ref_poly ) ;
   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      poly->apply_mapping( qr->point( i ), PT ) ;
      poly->apply_inverse_mapping( PT, PT_REF ) ;
      ok &= ( PT_REF->distance( qr->point( i ) ) < 1.e-8 ) ;
   }
   notify_one_test_result( "inverse_mapping/mapping", ok ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: test_mapping_derivative( GE_Mpolyhedron const* poly,
                                               double dx,
                                               double d_eps, double d_min )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron_TEST:: test_mapping_derivative" ) ;

   size_t nb_sp_dims = poly->nb_space_dimensions() ;
   size_t ref_dim = poly->dimension() ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_5" ) ;
   GE_Matrix* jac  = GE_Matrix::create( 0, nb_sp_dims, ref_dim ) ;
   GE_Matrix* jac2 = GE_Matrix::create( 0, nb_sp_dims, ref_dim ) ;
   GE_Matrix* tjac = GE_Matrix::create( 0, ref_dim, nb_sp_dims ) ;
   doubleArray3D hessian( nb_sp_dims, ref_dim, ref_dim ) ;
   bool nonzero_hessian ;

   GE_ReferencePolyhedron const* ref_poly = poly->reference_polyhedron() ;
   GE_QuadratureRule const* qr = qrp->quadrature_rule( ref_poly ) ;

   bool eq = false ;
   bool ok_1 = true ;
   bool ok_2 = true ;
   bool ok_3 = true ;
   bool ok_4 = true ;
   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      poly->apply_mapping( qr->point( i ), PT ) ;
      poly->build_mapping_derivative( qr->point( i ), jac ) ;
      poly->build_tr_mapping_derivative( qr->point( i ), tjac ) ;
      poly->build_mapping_hessian( qr->point( i ), &hessian, nonzero_hessian ) ;

      for( size_t jj=0 ; jj<ref_dim ; ++jj )
      {
         PT_REF->set( qr->point( i ) ) ;
         PT_REF->set_coordinate( jj, PT_REF->coordinate( jj ) + dx ) ;
         poly->apply_mapping( PT_REF, PT2 ) ;
         poly->build_mapping_derivative( PT_REF, jac2 ) ;
         for( size_t d=0 ; d<nb_sp_dims ; ++d )
         {
            // *** mapping % (finite differences)
            double dF = ( PT2->coordinate( d ) - PT->coordinate( d ) ) / dx ;
            eq = PEL::double_equality( jac->item( d, jj ), dF, d_eps, d_min ) ;
            if( !eq )
            {
               PEL::out() << "dF_{" << d << "}/dx_{" << jj << "}: " ;
               PEL::out() << "  exact=" << jac->item( d, jj ) ;
               PEL::out() << "  diff=" << dF << endl ;
            }
            ok_1 &= eq ;

            // *** transpose mapping % mapping
            eq = PEL::double_equality( jac->item( d, jj ),
                                       tjac->item( jj, d ), 1.e-8, 1.e-12 ) ;
            ok_2 &= eq ;

            // ***  hessian % (finite differences)
            for( size_t j=0 ; j<ref_dim ; ++j )
            {
               double d2F = ( jac2->item( d, j ) - jac->item( d, j ) ) / dx ;
               double exact = ( nonzero_hessian ? hessian( d, j, jj ) : 0.0 ) ;
               eq = PEL::double_equality( exact, d2F, d_eps, d_min ) ;
               if( !eq )
               {
                  PEL::out() << "d2F_{" << d << "}/dx_{" << j
                             << "}dx_{" << jj << "}: " ;
                  PEL::out() << "  exact=" << exact ;
                  PEL::out() << "  diff=" << d2F << endl ;
               }
               ok_3 &= eq ;
            }
         }
      }

      for( size_t d=0 ; d<nb_sp_dims ; ++d )
         for( size_t j=0 ; j<ref_dim ; ++j )
            for( size_t k=0 ; k<ref_dim ; ++k )
               ok_4 &= PEL::double_equality( hessian( d, j, k ),
                                             hessian( d, k, j ),
                                             1.e-14, 1.e-50 ) ;


   }

   notify_one_test_result( "mapping_derivative", ok_1 ) ;
   notify_one_test_result( "tr_mapping_derivative/mapping_derivative", ok_2 ) ;
   notify_one_test_result( "mapping_hessian", ok_3 ) ;
   notify_one_test_result( "symmetry of mapping_hessian", ok_4 ) ;

   jac->destroy() ;
   jac2->destroy() ;
   tjac->destroy() ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: test_contains( GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron_TEST:: test_contains" ) ;

   size_t nb_sp_dims = poly->nb_space_dimensions() ;
   double h = poly->inter_vertices_maximum_distance() ;

   bool ok = true ;
   for( size_t iv=0 ; iv<poly->nb_vertices() ; ++iv )
   {
      ok &= poly->contains( poly->vertex( iv ) ) ;
   }
   notify_one_test_result( "contains(vertices)", ok ) ;

   ok = poly->contains( poly->center() ) ;
   notify_one_test_result( "contains(center)", ok ) ;

   for( size_t d=0 ; d<nb_sp_dims; ++d )
   {
      PT->set_coordinate( d, poly->center()->coordinate(d) +
                             100.*((double)(d+1))*h) ;
   }
   ok = !poly->contains( PT ) ;
   notify_one_test_result( "!contains(far point)", ok ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron_TEST:: test_if_fv_center_is_circumcenter(
                                                  GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------------
{
   bool ok = true ;
   GE_Point const* fvc = poly->finite_volume_center() ;
   PEL_ASSERT( fvc != 0 ) ;
   double dd = fvc->distance( poly->vertex( 0 ) ) ;

   for( size_t iv=1 ; iv<poly->nb_vertices() ; ++iv )
   {
      double dd2 = fvc->distance( poly->vertex( iv ) ) ;
      bool eq = PEL::double_equality( dd, dd2, 1.e-10, 1.e-100 ) ;
      ok = ok && eq ;
      if( !eq ) out() << dd << "  " << dd2 << endl ;
   }
   notify_one_test_result( "fv_center_is_circumcenter", ok ) ;
}

//----------------------------------------------------------------------------
double
GE_Mpolyhedron_TEST:: estimate_measure( GE_Mpolyhedron const* poly ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron_TEST:: estimate_measure" ) ;

   double measure = 0. ;

   // Case same dimension for the polyhedron and its reference use of a change
   // of variable.
   if( poly->dimension() == poly->nb_space_dimensions() )
   {
      GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_5" ) ;

      GE_QuadratureRule const* qr =
         qrp->quadrature_rule( poly->reference_polyhedron() ) ;

      for( size_t i=0; i<qr->nb_points(); i++ )
      {
         poly->build_mapping_derivative( qr->point( i ), JACOB ) ;
         double contrib = 1.0 ;
         JACOB->compute_determinant() ;
         measure += contrib*PEL::abs( JACOB->determinant() )*qr->weight( i ) ;
      }
   }
   else if( poly->name()=="GE_Segment" )
      measure = poly->vertex( 0 )->distance( poly->vertex( 1 ) ) ;
   else if( poly->name()=="GE_Triangle" )
   {
      GE_Vector* V1 = GE_Vector::create( 0, poly->vertex( 1 ), poly->vertex( 0 ) ) ;
      GE_Vector* V2 = GE_Vector::create( 0, poly->vertex( 2 ), poly->vertex( 0 ) ) ;
      measure = PEL::sqrt( PEL::pow( V1->norm(), (size_t)2 )*PEL::pow( V2->norm(), (size_t)2 )
                      -PEL::pow( V1->dot_product( V2 ), (size_t)2 ) )/2. ;
      V1->destroy() ;
      V2->destroy() ;
   }
   else if( poly->name()=="GE_Rectangle" || poly->name()=="GE_Quadrilateral"|| poly->name()=="GE_Trapezoid"  )
   {
      GE_Vector* V1 = GE_Vector::create( 0, poly->vertex( 1 ), poly->vertex( 0 ) ) ;
      GE_Vector* V2 = GE_Vector::create( 0, poly->vertex( 3 ), poly->vertex( 0 ) ) ;
      measure = PEL::sqrt( PEL::pow( V1->norm(), (size_t)2 )*PEL::pow( V2->norm(), (size_t)2 )
                      -PEL::pow( V1->dot_product( V2 ), (size_t)2 ) )/2. ;
      V1->re_initialize( poly->vertex( 1 ), poly->vertex( 2 ) ) ;
      V2->re_initialize( poly->vertex( 3 ), poly->vertex( 2 ) ) ;
      measure += PEL::sqrt( PEL::pow( V1->norm(), (size_t)2 )*PEL::pow( V2->norm(), (size_t)2 )
                      -PEL::pow( V1->dot_product( V2 ), (size_t)2 ) )/2. ;
      V1->destroy() ;
      V2->destroy() ;
   }
   else
   {
      PEL_Error::object()->raise_internal(
                                    "Unknown polyhedron "+poly->name() ) ;
   }
   return( measure ) ;
}

//----------------------------------------------------------------------------
double
GE_Mpolyhedron_TEST:: compute_inter_vert_max_dist(
                                            GE_Mpolyhedron const* poly ) const
//----------------------------------------------------------------------------
{
   size_t nb_verts = poly->nb_vertices() ;
   double max_dist = 0.0 ;
   for( size_t iV1=0 ; iV1<nb_verts ; iV1++ )
   {
      for( size_t iV2=iV1+1 ; iV2<nb_verts ; iV2++ )
      {
         double d = poly->vertex(iV1)->distance( poly->vertex(iV2) ) ;
         if( d>max_dist )
            max_dist = d ;
      }
   }

   return( max_dist ) ;
}

//----------------------------------------------------------------------------
double
GE_Mpolyhedron_TEST:: compute_equi_ball_diameter(
                                            GE_Mpolyhedron const* poly ) const
//----------------------------------------------------------------------------
{
   size_t dd = poly->nb_space_dimensions() ;
   double diam = PEL::pow( PEL::pow( 2, dd-1 )*dd*poly->measure()/PEL::pi()/
                    PEL::max( (size_t)1, (size_t)(dd-1) ), 1.0/dd ) ;

   return( diam ) ;
}

//----------------------------------------------------------------------------
double
GE_Mpolyhedron_TEST:: compute_inter_vert_max_dist( GE_Mpolyhedron const* poly,
                                                   size_t dim ) const
//----------------------------------------------------------------------------
{
   size_t nb_verts = poly->nb_vertices() ;
   double max_dist = 0.0 ;
   for( size_t iV1=0 ; iV1<nb_verts ; iV1++ )
   {
      double x1 = poly->vertex(iV1)->coordinate(dim) ;
      for( size_t iV2=iV1+1 ; iV2<nb_verts ; iV2++ )
      {
         double x2 = poly->vertex(iV2)->coordinate(dim) ;
         double d = PEL::abs( x1-x2 ) ;
         if( d>max_dist ) max_dist = d ;
      }
   }
   return( max_dist ) ;
}
