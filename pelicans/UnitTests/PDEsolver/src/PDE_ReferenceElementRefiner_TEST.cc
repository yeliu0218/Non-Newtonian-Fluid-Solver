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

#include <PDE_ReferenceElementRefiner_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <boolVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <GE_ReferencePolyhedronRefiner.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ReferenceElementRefiner.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::endl ;
using std::ios_base ;
using std::ostringstream ;
using std::setw ;
using std::string ;

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner_TEST*
PDE_ReferenceElementRefiner_TEST:: REGISTRATOR
                                     = new PDE_ReferenceElementRefiner_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner_TEST:: PDE_ReferenceElementRefiner_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_ReferenceElementRefiner",
                     "PDE_ReferenceElementRefiner_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner_TEST:: ~PDE_ReferenceElementRefiner_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: process_one_test(
                                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: process_one_test" ) ;

   TEST_NAME = exp->name() ;
   DIM = exp->int_data( "nb_space_dimensions" ) ;
   if( exp->has_entry( "verbose_level" ) )
      VERB_LEVEL = exp->int_data( "verbose_level" ) ;
   else VERB_LEVEL = 0 ;
   MY_EPS = exp->double_data( "dbl_epsilon" ) ;
   MY_MIN = exp->double_data( "dbl_minimum" ) ;

   std::string nn = exp->string_data( "reference_element_refiner" ) ;
   ELRF = PDE_ReferenceElementRefiner::object( nn ) ;
   PDE_ReferenceElementRefiner::set_one_level_difference_rule(
                   PDE_ReferenceElementRefiner::Parents ) ;

   GE_ReferencePolyhedronRefiner const* crf = ELRF->cell_refiner() ;
   PDE_ReferenceElement const* elm = ELRF->reference_element() ;

   GE_SetOfPoints* verts = GE_SetOfPoints:: create( 0, DIM ) ;
   CPOLY = create_coarse_polyhedron( verts, exp ) ;

   bool ok_ref_poly =  ( CPOLY->reference_polyhedron() ==
                         elm->reference_polyhedron() ) ;
   ok_ref_poly = ok_ref_poly && ( crf->subcell_reference_polyhedron() ==
                                  elm->reference_polyhedron() ) ;

   notify_one_test_result( TEST_NAME+" reference polyhedron ", ok_ref_poly ) ;

   build_new_geometry( verts ) ;

   bool check_coefs = true ;
   if( exp->has_entry( "check_refinement_coefs" ) )
       check_coefs = exp->bool_data( "check_refinement_coefs" ) ;

   if( check_coefs ) check_refi_coefficients() ;

   nn =  exp->string_data( "quadrature_rule_provider" ) ;
   GE_QRprovider const* qrp = GE_QRprovider::object( nn ) ;
   GE_QuadratureRule const* qr = qrp->quadrature_rule(
                                      elm->reference_polyhedron() ) ;

   check_refi_equation( qr ) ;

   verts->destroy() ;
}

//-------------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_ReferenceElementRefiner_TEST:: create_coarse_polyhedron(
                                           GE_SetOfPoints* verts,
                                           PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: create_coarse_polyhedron" ) ;

   doubleVector const& vec = exp->doubleVector_data( "vertices_list" ) ;
   doubleVector coord( DIM ) ;
   for( size_t i=0; i<vec.size(); i=i+DIM )
   {
      for( size_t j=0; j<(size_t) DIM; j++ )
      {
         coord( j ) = vec(i+j) ;
      }
      GE_Point* pt = GE_Point::create( 0, coord ) ;
      verts->append( pt ) ;
      pt->destroy() ; pt = 0 ;
   }
   size_t_vector vertices( verts->nb_points() ) ;
   for( size_t i=0 ; i<verts->nb_points() ; ++i )
   {
      vertices(i) = i ;
   }
   GE_Mpolyhedron const* poly =
                  GE_Mpolyhedron:: create( exp->string_data("polyhedron_name"),
                                           verts, vertices ) ;
   return( poly ) ;
}

//-----------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: build_new_geometry( GE_SetOfPoints* verts )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: build_new_geometry" ) ;

   GE_ReferencePolyhedronRefiner const* crf = ELRF->cell_refiner() ;

   // numero local (dans crf) de vertex --> numero global
   size_t_vector glob_vert( crf->nb_vertices() ) ;

   // *** new vertices

   GE_Point* pt = GE_Point::create( 0, DIM ) ;

   bool newvert = false ;
   if( VERB_LEVEL > 1 ) out() << "   vertices involved :" << endl ;

   for( size_t i=0 ; i<crf->nb_vertices() ; ++i )
   {
      CPOLY->apply_mapping( crf->vertex( i ), pt ) ;
      if( VERB_LEVEL > 1 ) newvert = !verts->has( pt ) ;
      verts->extend( pt ) ;
      glob_vert( i ) = verts->index( pt ) ;
      if( VERB_LEVEL > 1 )
      {
         out() << setw( 9 ) << glob_vert( i ) ;
         pt->print( out(), 2 ) ;
         if( newvert ) out() << "   NEW" ;
         out() << endl ;
      }
   }

   pt->destroy() ;

   // *** new polyhedra

   if( VERB_LEVEL > 1 ) out() << "   subcells" << endl ;

   RPOLYS.clear() ;
   for( size_t ic=0 ; ic<crf->nb_subcells() ; ++ic )
   {
      if( VERB_LEVEL > 1 ) out() << "      vertices: " ;
      size_t nbvs = crf->subcell_reference_polyhedron()->nb_vertices() ;
      size_t_vector idx( nbvs ) ;
      for( size_t iv=0 ; iv<nbvs ; ++iv )
      {
         idx( iv ) = glob_vert( crf->subcell_vertex( ic, iv ) ) ;
         if( VERB_LEVEL > 1 ) out() << setw( 4 ) << idx( iv ) << " " ;
      }
      GE_Mpolyhedron* rpoly =
         GE_Mpolyhedron::create( CPOLY->name(), verts, idx ) ;
      RPOLYS.push_back( rpoly ) ;
      if( VERB_LEVEL > 1 ) out() << endl ;
   }
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: check_refi_coefficients( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: check_refi_coefficients" ) ;

   GE_Point* pt  = GE_Point::create( 0, DIM ) ;
   GE_Point* pt2 = GE_Point::create( 0, DIM ) ;
   GE_Point* pt_ref_c = GE_Point::create( 0, DIM ) ;

   GE_ReferencePolyhedronRefiner const* crf = ELRF->cell_refiner() ;
   PDE_ReferenceElement const* elm = ELRF->reference_element() ;

   boolVector ok_parent( elm->nb_nodes() ) ;
   ok_parent.set( true ) ;

   for( size_t ic=0 ; ic<crf->nb_subcells() ; ++ic )
   {
      GE_Mpolyhedron const* rpoly = RPOLYS[ic] ;
      for( size_t pn=0 ; pn<elm->nb_nodes() ; ++pn )
      {
         CPOLY->apply_mapping( elm->node_location( pn ), pt2 ) ;
         size_t leading_n = ELRF->leading_child_node( pn, ic ) ;
         for( size_t locn=0 ; locn<elm->nb_nodes() ; ++locn )
         {
            GE_Point const* pt_ref_r = elm->node_location( locn ) ;
            rpoly->apply_mapping( pt_ref_r, pt ) ;
            bool same_location = ( pt2->distance( pt ) < 1.e-8 ) ;

            bool ok_lead = ( leading_n == locn &&  same_location ) ||
                           ( leading_n != locn && !same_location ) ;
            if( !ok_lead )
            {
               out() << "subcell " << ic << " : invalid leading info" << endl ;
               CPOLY->print( out(), 0 ) ;
               out() << "parent node : " << pn << "   " ;
               pt2->print( out(), 0 ) ; out() << endl ;
               out() << "child node  : " << locn << "   " ;
               pt->print( out(), 0 ) ; out() << endl ;
               out() << "leading child node : " << leading_n << endl ;
            }
            ok_parent( pn ) = ok_parent( pn ) && ok_lead ;

            CPOLY->apply_inverse_mapping( pt, pt_ref_c ) ;
            double bf = elm->N_local( pn, pt_ref_c ) ;
            if( PEL::abs( bf ) > 1.e-6 )
            {
               bool found = false ;
               for( size_t ii=0 ; ii<ELRF->nb_childs( pn, ic ) ; ++ii )
               {
                  if( ELRF->child_node( pn, ic, ii ) == locn )
                  {
                     found = true ;
                     break ;
                  }
               }
               ok_parent( pn ) = ok_parent( pn ) && found ;
               if( !found )
               {
                  out() << "parent node " << pn
                       << " : forgotten child node " <<  locn
                       << " in subcell " << ic << endl ;
               }
            }
         }
      }
      for( size_t rn=0 ; rn<elm->nb_nodes() ; ++rn )
      {
         GE_Point const* pt_ref_r = elm->node_location( rn ) ;
         rpoly->apply_mapping( pt_ref_r, pt ) ;
         CPOLY->apply_inverse_mapping( pt, pt_ref_c ) ;
         for( size_t i = 0 ; i<ELRF->nb_parents( rn, ic ) ; ++i )
         {
            size_t pn = ELRF->parent_node( rn, ic, i ) ;
            double bf = elm->N_local( pn, pt_ref_c ) ;
            double r_coef = ELRF->refinement_coef( rn, pn, ic ) ;
            bool eq = PEL::double_equality( bf, r_coef, MY_EPS, MY_MIN ) ;
            ok_parent( pn ) = ok_parent( pn ) && eq ;
            if( !eq ) display_error_coef( ic, rn, pn, bf, r_coef ) ;
         }
      }
   }

   for( size_t pn=0 ; pn<elm->nb_nodes() ; ++pn )
   {
      ostringstream mesg ;
      mesg <<  TEST_NAME << " refinement coefs (" << pn << ")" ;
      notify_one_test_result( mesg.str(), ok_parent( pn ) ) ;
   }

   pt->destroy() ;
   pt2->destroy() ;
   pt_ref_c->destroy() ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: check_refi_equation(
                                              GE_QuadratureRule const* qr )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: check_refi_equation" ) ;

   GE_Point* pt = GE_Point::create( 0, DIM ) ;
   GE_Point* pt_ref_c = GE_Point::create( 0, DIM ) ;

   GE_ReferencePolyhedronRefiner const* crf = ELRF->cell_refiner() ;
   PDE_ReferenceElement const* elm = ELRF->reference_element() ;

   boolVector ok_parent( elm->nb_nodes() ) ;
   ok_parent.set( true ) ;
   for( size_t ic=0 ; ic<crf->nb_subcells() ; ++ic )
   {
      GE_Mpolyhedron const* rpoly = RPOLYS[ic] ;
      for( size_t ip=0 ; ip<qr->nb_points() ; ++ip )
      {
         GE_Point const* pt_ref_r = qr->point( ip ) ;
         rpoly->apply_mapping( pt_ref_r, pt ) ;
         CPOLY->apply_inverse_mapping( pt, pt_ref_c ) ;
         for( size_t pn=0 ; pn<elm->nb_nodes() ; ++pn )
         {
            double val_c = elm->N_local( pn, pt_ref_c ) ;
            double val_r = 0.0 ;
            for( size_t i=0 ; i<ELRF->nb_childs( pn, ic ) ; ++i )
            {
               size_t cn = ELRF->child_node( pn, ic, i ) ;

               double coef = ELRF->refinement_coef( cn, pn, ic ) ;
               val_r += coef * elm->N_local( cn, pt_ref_r ) ;
            }
            bool eq = PEL::double_equality( val_c, val_r, MY_EPS, MY_MIN ) ;
            ok_parent( pn ) = ok_parent( pn ) && eq ;
            if( !eq ) display_error( pn, val_c, val_r ) ;
         }
      }
   }

   for( size_t pn=0 ; pn<elm->nb_nodes() ; ++pn )
   {
      ostringstream mesg ;
      mesg <<  TEST_NAME << " refinement equation (" << pn << ")" ;
      notify_one_test_result( mesg.str(), ok_parent( pn ) ) ;
   }

   pt->destroy() ;
   pt_ref_c->destroy() ;
}

//-----------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: display_error_coef(
                                            size_t ic,
                                            size_t child_n,
                                            size_t pn,
                                            double xx_1, double xx_2 ) const
//-----------------------------------------------------------------------
{
   std::ios_base::fmtflags original_flags = out().flags() ;
   out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;
   out() << "subcell " << ic
         << ", node " << child_n
         << " of parent " <<  pn << " : "
         << std::setw( 20 ) << xx_1
         << std::setw( 20 ) << xx_2
         << endl ;
   out().flags( original_flags ) ;
}

//-----------------------------------------------------------------------
void
PDE_ReferenceElementRefiner_TEST:: display_error( size_t pn,
                                            double xx_1, double xx_2 ) const
//-----------------------------------------------------------------------
{
   std::ios_base::fmtflags original_flags = out().flags() ;
   out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;
   out() << setw( 5 ) << pn
         << std::setw( 20 ) << xx_1
         << std::setw( 20 ) << xx_2
         << endl ;
   out().flags( original_flags ) ;
}
