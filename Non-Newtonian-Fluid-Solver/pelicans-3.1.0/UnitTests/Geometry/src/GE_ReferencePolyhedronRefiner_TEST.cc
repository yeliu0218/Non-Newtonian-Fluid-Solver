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

#include <GE_ReferencePolyhedronRefiner_TEST.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <boolVector.hh>
#include <stringVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_QuadratureRule.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>
#include <GE_SetOfPoints.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::ostringstream ;
using std::setw ;
using std::string ;

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner_TEST*
GE_ReferencePolyhedronRefiner_TEST:: REGISTRATOR
                                     = new GE_ReferencePolyhedronRefiner_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner_TEST:: GE_ReferencePolyhedronRefiner_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_ReferencePolyhedronRefiner", 
                     "GE_ReferencePolyhedronRefiner_TEST" )
{
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner_TEST:: ~GE_ReferencePolyhedronRefiner_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner_TEST:: process_one_test( 
                                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner_TEST:: process_one_test" ) ;

   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, 
                                 "list_of_GE_ReferencePolyhedronRefiner" ) ;
   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      PEL_ModuleExplorer const* se = ee->create_subexplorer( ee ) ;
      GE_ReferencePolyhedronRefiner const* cr = 
                               GE_ReferencePolyhedronRefiner::make( 0, se ) ;
      std::string const& name = se->string_data( "concrete_name" ) ;
      check_faces( name, cr ) ;

      GE_ReferencePolyhedron const* rp = cr->subcell_reference_polyhedron() ;

      size_t nb_dims = rp->dimension() ;

      GE_Point* pt_subcell = GE_Point::create( 0, nb_dims ) ;
      GE_Point* pt_ref_c = GE_Point::create( 0, nb_dims ) ;
      GE_Point* pt = GE_Point::create( 0, nb_dims ) ;

      GE_SetOfPoints* verts = GE_SetOfPoints:: create( 0, nb_dims ) ;
      for( size_t iv=0 ; iv<cr->nb_vertices() ; ++iv )
      {
         verts->append( cr->vertex( iv ) ) ;
      }

      GE_Mpolyhedron const* cpoly = create_coarse_polyhedron( verts, rp ) ;

      std::vector< GE_Mpolyhedron const* > rpolys ;
      build_new_polyhedra( cr, verts, rpolys ) ;

      // pour obtenir des points dans la maille de reference
      GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_5" ) ;
      GE_QuadratureRule const* qr = qrp->quadrature_rule( rp ) ;
      for( size_t ic=0 ; ic<cr->nb_subcells() ; ++ic )
      {
         bool ok = true ;
         GE_Mpolyhedron const* rpoly = rpolys[ic] ;
         for( size_t ip=0 ; ip<qr->nb_points() ; ++ip )
         {
            GE_Point const* pt_ref_r = qr->point( ip ) ;
            cr->compute_location_in_subcell( ic, pt_ref_r, pt_subcell ) ;

            cpoly->apply_mapping( pt_ref_r, pt ) ;
            rpoly->apply_inverse_mapping( pt, pt_ref_c ) ;

            ok = ok && pt_subcell->distance( pt_ref_c ) < 1.e-10 ;
         }
         ostringstream mesg ;
         mesg << name << " subcell " << ic ;
         notify_one_test_result( mesg.str(), ok ) ;
      }
      verts->destroy() ;
      pt_subcell->destroy() ;
      pt_ref_c->destroy() ;
      pt->destroy() ;
      cr->destroy() ;
   }
   ee->destroy() ;
}

//----------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner_TEST:: check_faces( 
                                      std::string const& classname,
                                      GE_ReferencePolyhedronRefiner const* cr )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner_TEST:: check_faces" ) ;

   boolVector traversed_face( cr->nb_subfaces() ) ;
   traversed_face.set( false ) ;
   
   GE_ReferencePolyhedron const* fpoly = cr->subface_reference_polyhedron() ;
   GE_ReferencePolyhedron const* cpoly = cr->subcell_reference_polyhedron() ;
   bool ok = true ;
   for( size_t i_cell=0 ; i_cell<cr->nb_subcells() ; ++i_cell )
   {
      for( size_t is=0 ; is<cpoly->nb_faces()  ; ++is )
      {
         size_t i_face = cr->subcell_face( i_cell, is ) ;
         std::set<size_t> face_1 ;
         std::set<size_t> face_2 ;
         for( size_t iv=0 ; iv<fpoly->nb_vertices() ; ++iv )
         {
            face_1.insert( cr->subface_vertex( i_face, iv ) ) ;
            face_2.insert( cr->subcell_vertex( i_cell, 
                                       cpoly->face_vertex( is, iv ) ) ) ;
         }
         ok &= ( face_1 == face_2 ) ;
         traversed_face( i_face ) = true ;
      }
   }
   
   GE_ReferencePolyhedron const* coarse_poly = cr->reference_polyhedron() ;
   for( size_t i_face=0 ; i_face<cr->nb_subfaces() ; ++i_face )
   {
      ok &= traversed_face( i_face ) ;
      size_t i_face_parent = PEL::bad_index() ;
      for( size_t jj=0 ; jj<coarse_poly->nb_faces() ; ++jj )
      {
         bool contains = true ;
         for( size_t iv=0 ; iv<fpoly->nb_vertices() ; ++iv )
         {
            size_t gv = cr->subface_vertex( i_face, iv ) ;
            GE_Point const* pt = cr->vertex( gv ) ;
            contains &= coarse_poly->face_contains( jj, pt ) ;
         }
         if( contains )
         {
            ok &= ( i_face_parent == PEL::bad_index() ) ;
            i_face_parent = jj ;
         }      
      }
      ok &= ( i_face_parent == cr->subface_parent( i_face ) ) ;      
   }
   
   notify_one_test_result( classname+" faces", ok ) ;
   
}

//----------------------------------------------------------------------------
GE_Mpolyhedron const*
GE_ReferencePolyhedronRefiner_TEST:: create_coarse_polyhedron( 
                                           GE_SetOfPoints* verts,
                                           GE_ReferencePolyhedron const* rp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner_TEST:: create_coarse_polyhedron" ) ;

   size_t nb_verts = rp->nb_vertices() ;
   size_t_vector vertices( nb_verts ) ;
   for( size_t iv=0 ; iv<nb_verts ; ++iv )
   {
      PEL_ASSERT( verts->has( rp->vertex( iv ) ) ) ;
      vertices( iv ) = verts->index( rp->vertex( iv ) ) ;
   }
   GE_Mpolyhedron const* result =
                         create_polyhedron( rp->name(), verts, vertices ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner_TEST:: build_new_polyhedra( 
                            GE_ReferencePolyhedronRefiner const* cr, 
                            GE_SetOfPoints* verts, 
                            std::vector< GE_Mpolyhedron const* >& rpolys )
//---------------------------------------------------------------------------
{
   GE_ReferencePolyhedron const* ref_poly = cr->subcell_reference_polyhedron() ;
   size_t nb_verts = ref_poly->nb_vertices() ;

   for( size_t ic=0 ; ic<cr->nb_subcells() ; ++ic )
   {
      size_t_vector vertices( nb_verts ) ;
      for( size_t iv=0 ; iv<nb_verts ; ++iv )
      {
         GE_Point const* pt = cr->vertex( cr->subcell_vertex( ic, iv ) ) ;
         PEL_ASSERT( verts->has( pt ) ) ;
         vertices( iv ) = verts->index( pt ) ;
      }
      GE_Mpolyhedron const* pp = 
                     create_polyhedron( ref_poly->name(), verts, vertices ) ;
      rpolys.push_back( pp ) ;
   }
}

//----------------------------------------------------------------------------
GE_Mpolyhedron const*
GE_ReferencePolyhedronRefiner_TEST:: create_polyhedron(
                                         std::string const& ref_poly_name,
                                         GE_SetOfPoints* verts,
                                         size_t_vector const& vertices )
//----------------------------------------------------------------------------
{
   std::string poly_name ;
   if( ref_poly_name == "GE_ReferenceSquare" )
   {
      poly_name = "GE_Rectangle" ;
   }
   else if( ref_poly_name == "GE_ReferenceTriangle" )
   {
      poly_name = "GE_Triangle" ;
   }
   else if( ref_poly_name == "GE_ReferenceCube" )
   {
      poly_name = "GE_Cuboid" ;
   }
   else
   {
      ostringstream mesg ;
      mesg << "*** GE_ReferencePolyhedronRefiner_TEST:" << endl ;
      mesg << "    the reference polyhedron \"" << ref_poly_name 
           << "\"" << endl ;
      mesg << "    is not handled" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;  
   }
   GE_Mpolyhedron const* result = 
                  GE_Mpolyhedron::create( poly_name, verts, vertices ) ;
   return( result ) ;
}
