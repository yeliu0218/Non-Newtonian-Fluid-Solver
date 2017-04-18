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

#include <GE_ReferencePolyhedronRefiner.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <iostream>
#include <map>
#include <set>

using std::string ;
using std::cout ; using std::endl ;

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner const*
GE_ReferencePolyhedronRefiner:: make( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string const& nn = exp->string_data( "concrete_name" ) ;
   GE_ReferencePolyhedronRefiner const* proto =
   static_cast<GE_ReferencePolyhedronRefiner const*>(plugins_map()->item(nn)) ;
   
   GE_ReferencePolyhedronRefiner const* result = 
                                 proto->create_replica( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner:: GE_ReferencePolyhedronRefiner(
                                PEL_Object* a_owner,
                                GE_ReferencePolyhedron const* ref_poly,
                                GE_ReferencePolyhedron const* subcell_ref_poly,
                                GE_ReferencePolyhedron const* subface_ref_poly )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , RPOLY( ref_poly )
   , SUBCELL_RPOLY( subcell_ref_poly )
   , SUBFACE_RPOLY( subface_ref_poly )
   , NB_VERTS( 0 )
   , NB_FACES( 0 )
   , NB_CELLS( 0 )
   , NB_SUB_PER_FACE( 0 )
   , VERTICES( 0 )
   , FACE_VERT( 0, 0 )
   , FACE_PARENT( 0 )
   , CELL_VERT( 0, 0 )
   , CELL_FACE( 0, 0 )
   , i_SUBFACE( 0 )
{
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner:: GE_ReferencePolyhedronRefiner( 
                                            std::string const& a_name )
//---------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , FACE_VERT( 0, 0 )
   , FACE_PARENT( 0 )
   , CELL_VERT( 0, 0 )
   , CELL_FACE( 0, 0 )
{
   PEL_LABEL( 
   "GE_ReferencePolyhedronRefiner:: GE_ReferencePolyhedronRefiner" ) ;
   
   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner:: ~GE_ReferencePolyhedronRefiner( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
bool
GE_ReferencePolyhedronRefiner:: is_a_prototype( void ) const
//---------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_ReferencePolyhedronRefiner:: reference_polyhedron( void ) const
//---------------------------------------------------------------------------
{
   return( RPOLY ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: nb_vertices( void ) const
//---------------------------------------------------------------------------
{
   return( NB_VERTS ) ;
}

//---------------------------------------------------------------------------
GE_Point const*
GE_ReferencePolyhedronRefiner:: vertex( size_t i ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: vertex" ) ;
   PEL_CHECK_PRE( i < nb_vertices() ) ;

   GE_Point const* result = static_cast<GE_Point*>( VERTICES->at( i ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_ReferencePolyhedronRefiner:: subface_reference_polyhedron( void ) const
//---------------------------------------------------------------------------
{
   return( SUBFACE_RPOLY ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: nb_subfaces_per_face( void ) const
//---------------------------------------------------------------------------
{
   return( NB_SUB_PER_FACE ) ;
}

//---------------------------------------------------------------------------
size_t 
GE_ReferencePolyhedronRefiner:: nb_subfaces( void ) const
//---------------------------------------------------------------------------
{
   return( NB_FACES ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: subface_vertex( size_t is, size_t iv ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: subface_vertex" ) ;
   PEL_CHECK_PRE( is < nb_subfaces() ) ;
   PEL_CHECK_PRE( iv < subface_reference_polyhedron()->nb_vertices() ) ;

   size_t result = FACE_VERT( is, iv ) ;

   PEL_CHECK_POST( result < nb_vertices() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: subface_parent( size_t is ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: subface_parent" ) ;
   PEL_CHECK_PRE( is < nb_subfaces() ) ;

   size_t result = FACE_PARENT( is ) ;

   PEL_CHECK_POST( result == PEL::bad_index() || 
                   result < reference_polyhedron()->nb_faces() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_ReferencePolyhedronRefiner:: subcell_reference_polyhedron( void ) const
//---------------------------------------------------------------------------
{
   return( SUBCELL_RPOLY ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: nb_subcells( void ) const
//---------------------------------------------------------------------------
{
   return( NB_CELLS ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: subcell_vertex( size_t ic, size_t iv ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: subcell_vertex" ) ;
   PEL_CHECK_PRE( ic < nb_subcells() ) ;
   PEL_CHECK_PRE( iv < subcell_reference_polyhedron()->nb_vertices() ) ;

   size_t result = CELL_VERT( ic, iv ) ;

   PEL_CHECK_POST( result < nb_vertices() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
GE_ReferencePolyhedronRefiner:: subcell_face( size_t ic, size_t is ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: subcell_face" ) ;
   PEL_CHECK_PRE( ic < nb_subcells() ) ;
   PEL_CHECK_PRE( is < subcell_reference_polyhedron()->nb_faces() ) ;

   size_t result = CELL_FACE( ic, is ) ;

   PEL_CHECK_POST( result < nb_subfaces() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner:: set_dimensions( size_t a_nb_vertices,
                                                size_t a_nb_subfaces,
                                                size_t a_nb_subcells,
                                                size_t a_nb_subfaces_per_face )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: set_dimensions" ) ;
   PEL_CHECK_PRE( nb_vertices() == 0 ) ;
   PEL_CHECK_PRE( nb_subfaces() == 0 ) ;
   PEL_CHECK_PRE( nb_subcells() == 0 ) ;
   PEL_CHECK_PRE( nb_subfaces_per_face() == 0 ) ;
   PEL_CHECK_PRE( a_nb_vertices != 0 ) ;
   PEL_CHECK_PRE( a_nb_subfaces != 0 ) ;
   PEL_CHECK_PRE( a_nb_subcells != 0 ) ;
   PEL_CHECK_PRE( a_nb_subfaces_per_face != 0 ) ;

   NB_VERTS = a_nb_vertices ;
   NB_FACES = a_nb_subfaces ;
   NB_CELLS = a_nb_subcells ;
   NB_SUB_PER_FACE = a_nb_subfaces_per_face ;
   VERTICES = PEL_Vector::create( this, a_nb_vertices ) ;
   FACE_VERT.re_initialize( a_nb_subfaces, SUBFACE_RPOLY->nb_vertices() ) ;
   FACE_VERT.set( PEL::bad_index() ) ;
   FACE_PARENT.re_initialize( a_nb_subfaces ) ;
   FACE_PARENT.set( PEL::bad_index() ) ;
   CELL_VERT.re_initialize( a_nb_subcells, SUBCELL_RPOLY->nb_vertices() ) ;
   CELL_VERT.set( PEL::bad_index() ) ;
   CELL_FACE.re_initialize( a_nb_subcells, SUBCELL_RPOLY->nb_faces() ) ;
   CELL_FACE.set( PEL::bad_index() ) ;

   PEL_CHECK_POST( nb_vertices() == a_nb_vertices ) ;
   PEL_CHECK_POST( nb_subfaces() == a_nb_subfaces ) ;
   PEL_CHECK_POST( nb_subcells() == a_nb_subcells ) ;
   PEL_CHECK_POST( nb_subfaces_per_face() == a_nb_subfaces_per_face ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner:: set_vertex( size_t i_vertex, GE_Point* pt )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: set_vertex" ) ;
   PEL_CHECK_PRE( i_vertex < nb_vertices() ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->owner() == this ) ;
   PEL_CHECK_PRE( pt->nb_coordinates() == 
                  subcell_reference_polyhedron()->dimension() ) ;

   VERTICES->set_at( i_vertex, pt ) ;

   PEL_CHECK_POST( vertex( i_vertex ) == pt ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedronRefiner:: set_subcell( size_t i_cell,
                                        size_t_vector const& vertex_indices )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedronRefiner:: set_subcell" ) ;
   PEL_CHECK_PRE( i_cell < nb_subcells() ) ;
   PEL_CHECK_PRE( vertex_indices.size() == 
                  subcell_reference_polyhedron()->nb_vertices() ) ;

   for( size_t iv=0 ; iv<vertex_indices.size() ; ++iv )
   {
      PEL_ASSERT( CELL_VERT( i_cell, iv ) == PEL::bad_index() ) ;
      CELL_VERT( i_cell, iv ) = vertex_indices( iv ) ;
   }

   for( size_t is=0 ; is<SUBCELL_RPOLY->nb_faces() ; ++is )
   {
      std::set<size_t> vertofside ;
      size_t nb_vvv = SUBCELL_RPOLY->nb_face_vertices( is ) ;
      size_t_vector vvv( nb_vvv ) ;
      for( size_t iv=0 ; iv<nb_vvv ; ++iv )
      {
         size_t gv = vertex_indices( SUBCELL_RPOLY->face_vertex( is, iv ) ) ;
         vertofside.insert( gv ) ;
         vvv( iv ) = gv ;
      }
      FACES_it = FACES.find( vertofside ) ;
      size_t current_face = PEL::bad_index() ;
      if( FACES_it == FACES.end() )
      {
         current_face = i_SUBFACE ;
         FACES[ vertofside ] = current_face ;
         for( size_t ii=0 ; ii<vvv.size() ; ++ii )
         {
            FACE_VERT( current_face, ii ) = vvv( ii ) ;
         }
         // search of the parent face
         for( size_t icf=0 ; icf<RPOLY->nb_faces() ; ++icf )
         {
            bool is_parent = true ;
            for( size_t iv=0 ; iv<nb_vvv && is_parent ; ++iv )
            {
               is_parent &= RPOLY->face_contains( icf, vertex( vvv( iv ) ) ) ;
            }
            if( is_parent )
            {
               PEL_ASSERT( FACE_PARENT( current_face ) == PEL::bad_index() ) ;
               FACE_PARENT( current_face ) = icf ;
            }
         }
         i_SUBFACE++ ;
      }
      else
      {
         current_face = FACES_it->second ;
      }
      PEL_ASSERT( current_face != PEL::bad_index() ) ;
      CELL_FACE( i_cell, is ) = current_face ;
   }

   PEL_CHECK_POST( 
      FORALL( 
      ( size_t iv=0 ; iv<subcell_reference_polyhedron()->nb_vertices() ; ++iv ),
         subcell_vertex( i_cell, iv ) == vertex_indices( iv ) ) ) ;
}

//---------------------------------------------------------------------------
bool
GE_ReferencePolyhedronRefiner:: compute_location_in_subcell_PRE( 
                                               size_t ic, 
                                               GE_Point const* pt_cell,
                                               GE_Point* pt_subcell ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( ic < nb_subcells() ) ;
   PEL_ASSERT( pt_cell != 0 ) ;
   PEL_ASSERT( pt_cell->nb_coordinates() == 
               subcell_reference_polyhedron()->dimension() ) ;
   PEL_ASSERT( pt_subcell != 0 ) ;
   PEL_ASSERT( pt_subcell->nb_coordinates() == 
               subcell_reference_polyhedron()->dimension() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
GE_ReferencePolyhedronRefiner:: create_replica_PRE( 
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
GE_ReferencePolyhedronRefiner:: create_replica_POST(  
                                  GE_ReferencePolyhedronRefiner const* result,
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
PEL_ObjectRegister*
GE_ReferencePolyhedronRefiner:: plugins_map( void )
//---------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "GE_ReferencePolyhedronRefiner descendant" ) ;
   return( result ) ;
}
