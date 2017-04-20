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

#include <GE_GambitMeshing.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <set>
#include <sstream>
#include <vector>

using std::cout ; using std::endl ;
using std::set ;
using std::string ;
using std::vector ;
using std::ostringstream ;

GE_GambitMeshing const*
GE_GambitMeshing::PROTOTYPE = new GE_GambitMeshing() ;

struct GE_GambitMeshing_ERROR
{
   static void n0( string const& fname ) ;
   static void n1( size_t gambit_dim, size_t pel_dim ) ;
   static void n2( size_t itype ) ;
   static void n3( string mesh, size_t ndp ) ;
   static void n4( string mesh ) ;
   static void n5( void ) ;
   static void n6( string gmesh, string cells, string sides ) ;
} ;

//-----------------------------------------------------------------------------
GE_GambitMeshing:: GE_GambitMeshing( void )
//-----------------------------------------------------------------------------
   : GE_Meshing( "GE_GambitMeshing" )
   , NTYPE( 0 )
   , COORD_OF_CURRENT_VERT( 0 )
   , SIDE_TREE(0)
   , SIDE_OF_CELL( 0, 0 )
   , VERT_OF_CURRENT_CELL( 0 )
   , SIDE_OF_CURRENT_CELL( 0 )
   , VERT_OF_SIDE( 0, 0 )
   , VERT_OF_CURRENT_SIDE( 0 )
{
}

//-----------------------------------------------------------------------------
GE_GambitMeshing*
GE_GambitMeshing:: create_replica( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp,
                                    size_t dim_space ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_GambitMeshing* result = new GE_GambitMeshing( a_owner, exp, dim_space ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_GambitMeshing:: GE_GambitMeshing( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp, 
                                     size_t dim_space )
//-----------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , INPUT( exp->string_data( "filename" ).c_str(), std::ios_base::binary )
   , NTYPE( PEL::bad_index() )
   , COLOR_OF_VERT( 0 )
   , COORD_OF_CURRENT_VERT( dim_space )
   , SIDE_TREE( PEL_BalancedBinaryTree:: create( this ) )
   , SIDE_OF_CELL( 0, 0 )
   , COLOR_OF_CELL( 0 )
   , VERT_OF_CURRENT_CELL( 0 )
   , SIDE_OF_CURRENT_CELL( 0 )
   , NB_SIDES( 0 )
   , VERT_OF_SIDE( 0, 0 )
   , COLOR_OF_SIDE( 0 )
   , VERT_OF_CURRENT_SIDE( 0 )
{
   PEL_LABEL( "GE_GambitMeshing:: GE_GambitMeshing" ) ;

   initialize_rounding_strategy( exp ) ;

   std::string const& fname = exp->string_data( "filename" ) ;
   if( !INPUT ) GE_GambitMeshing_ERROR::n0( fname ) ;

   read_mesh_polyhedron( exp, dim_space, SIDE_POLY_NAME, CELL_POLY_NAME ) ;

   string buffer_string ;
   size_t buffer_size_t = 0 ;

   go_i_th_next_line( 6 ) ;

   size_t ngrps, nbsets, ndfcd ;
   INPUT >> NB_VERTS >> NB_CELLS >> ngrps >> nbsets >> ndfcd ;
   if( ndfcd != nb_space_dimensions() ) 
      GE_GambitMeshing_ERROR::n1( ndfcd, nb_space_dimensions() ) ;

   // dimensionnement des tableaux internes
   GE_ReferencePolyhedron const* side_polyhedron
                     = GE_Mpolyhedron::reference_polyhedron( SIDE_POLY_NAME ) ;
   GE_ReferencePolyhedron const* cell_polyhedron
                     = GE_Mpolyhedron::reference_polyhedron( CELL_POLY_NAME ) ;
   NB_SIDES_PER_CELL = cell_polyhedron->nb_faces() ;
   VERT_OF_CURRENT_SIDE.re_initialize( side_polyhedron->nb_vertices() ) ;
   VERT_OF_SIDE.re_initialize( NB_SIDES_PER_CELL*NB_CELLS,
                               side_polyhedron->nb_vertices() ) ;
   SIDE_OF_CELL.re_initialize( NB_CELLS, NB_SIDES_PER_CELL ) ;
   COLOR_OF_VERT.resize( NB_VERTS, GE_Color::null_color() ) ;
   COLOR_OF_CELL.resize( NB_CELLS, GE_Color::null_color() ) ;

   go_i_th_next_line( 3 ) ; // next line... ENDOFSECTION... NODAL COORDINATES

   VERT_POS = INPUT.tellg() ;

   go_i_th_next_line( NB_VERTS ) ; // coordinates of all vertices
   go_i_th_next_line( 2 ) ;        // ENDOFSECTION... ELEMENTS/CELLS

   CELL_POS = INPUT.tellg() ;

   for( size_t igc=0 ; igc<NB_CELLS ; ++igc )
   {
      read_node_to_point_coordinate_data( NTYPE ) ;
      build_sides_connectivity( igc ) ;
   }
   check_meshes_consistency() ;

   go_i_th_next_line( 2 ) ; // next line... ELEMENT GROUP

   for( size_t igrp=0 ; igrp<ngrps ; ++igrp )
   {
      go_i_th_next_line( 1 ) ; // ELEMENT GROUP
      size_t nbc ;
      INPUT >> buffer_string >> buffer_size_t >> buffer_string >> nbc ;

      go_i_th_next_line( 1 ) ; // end of line GROUP

      INPUT >> buffer_string ;
      GE_Color::extend( buffer_string ) ;
      GE_Color const* color = GE_Color::object( buffer_string ) ;

      go_i_th_next_line( 2 ) ;

      for( size_t i=0 ; i<nbc ; ++i )
      {
         INPUT >> buffer_size_t ;
         COLOR_OF_CELL[ buffer_size_t-1 ] = color ;
      }

      go_i_th_next_line( 2 ) ; // next line... ENDOFSECTION
   }

   // SECTION : BOUNDARY CONDITIONS
   COLOR_OF_SIDE.resize( NB_SIDES, GE_Color::null_color() ) ;
   for( size_t iset=0 ; iset<nbsets ; ++iset )
   {
      go_i_th_next_line( 1 ) ; // BOUNDARY CONDITIONS
      INPUT >> buffer_string ; 
      GE_Color::extend( buffer_string ) ;
      GE_Color const* color = GE_Color::object( buffer_string ) ;
      INPUT >> buffer_size_t ;
      if( buffer_size_t != 1 )
         GE_GambitMeshing_ERROR::n2( buffer_size_t ) ;
      size_t nentry ;
      INPUT >> nentry ;
      go_i_th_next_line( 1 ) ;
      for( size_t ie=0 ; ie<nentry ; ++ie )
      {
         size_t igc, is ;
         INPUT >> igc >> buffer_size_t >> is ;
         COLOR_OF_SIDE[ SIDE_OF_CELL(igc-1,is-1) ] = color ;
      }
      go_i_th_next_line( 2 ) ; // ENDOFSECTION
   }

   build_colors_of_vertices() ;
   
   I_VERT = NB_VERTS+1 ;
   I_SIDE = NB_SIDES+1 ;
   I_CELL = NB_CELLS+1 ;
}

//-----------------------------------------------------------------------------
GE_GambitMeshing:: ~GE_GambitMeshing( void )
//-----------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//-----------------------------------------------------------------------------
size_t
GE_GambitMeshing:: nb_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: nb_vertices" ) ;
   return NB_VERTS ;
}

//-----------------------------------------------------------------------------
size_t
GE_GambitMeshing:: nb_cells( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: nb_cells" ) ;
   size_t result = NB_CELLS ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t
GE_GambitMeshing:: nb_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: nb_faces" ) ;
   size_t result = NB_SIDES ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: start_vertex_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   INPUT.seekg( VERT_POS ) ;
   I_VERT = 0 ;
   read_vertex() ;
}

//-----------------------------------------------------------------------------
bool
GE_GambitMeshing:: valid_vertex( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return I_VERT<NB_VERTS ;   
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: go_next_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++I_VERT ;
   if( valid_vertex() ) read_vertex() ;
}

//-----------------------------------------------------------------------------
doubleVector const& 
GE_GambitMeshing:: vertex_coordinates( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   doubleVector const& result = COORD_OF_CURRENT_VERT ;

   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: start_cell_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   I_CELL = 0 ;
   INPUT.seekg( CELL_POS ) ;
   read_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_GambitMeshing:: valid_cell( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: valid_cell" ) ;

   bool result = ( I_CELL< NB_CELLS )  ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: go_next_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   I_CELL++ ;
   if( valid_cell() ) read_cell() ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_GambitMeshing:: cell_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_GambitMeshing:: cell_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = VERT_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_GambitMeshing:: cell_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = SIDE_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: start_face_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   I_SIDE = 0 ;
   read_side() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_GambitMeshing:: valid_face( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: valid_face" ) ;

   bool result = ( I_SIDE<NB_SIDES ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: go_next_face( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   I_SIDE++ ;
   if( valid_face() )
   {
      read_side() ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_GambitMeshing:: face_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = SIDE_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_GambitMeshing:: face_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = VERT_OF_CURRENT_SIDE ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_GambitMeshing:: default_vertex_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = COLOR_OF_VERT[ I_VERT ] ;

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_GambitMeshing:: default_cell_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = COLOR_OF_CELL[ I_CELL ] ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_GambitMeshing:: default_face_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = COLOR_OF_SIDE[ I_SIDE ] ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: read_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: read_vertex" ) ;
   PEL_CHECK( valid_vertex() ) ;   
   PEL_CHECK_INV( invariant() ) ;
   PEL_ASSERT( !INPUT.eof() ) ;

   size_t ii ;
   INPUT >> ii ;
   PEL_ASSERT( ii == I_VERT+1 ) ;

   for( size_t i=0 ; i<nb_space_dimensions() ; ++i )
   {
      double coordinate = PEL::max_double() ;
      INPUT >> coordinate ;
      COORD_OF_CURRENT_VERT( i ) = roundoff( coordinate ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: read_node_to_point_coordinate_data( size_t& elm_geom_type )
//-----------------------------------------------------------------------------
{
   size_t buffer_size_t ;

   // ne    : global element number
   // ntype : element geometry type
   // ndp   : number of nodes that define the element
   size_t ne, ntype, ndp ;
   INPUT >> ne >> ntype >> ndp ;

   if( elm_geom_type == PEL::bad_index() )
   {
      elm_geom_type = ntype ;
   }
   else if( elm_geom_type != ntype )
   {
      if( elm_geom_type != ntype ) GE_GambitMeshing_ERROR::n5() ;
   }

   size_t_vector gcom( ndp ) ;   
   for( size_t iv=0 ; iv<ndp ; ++iv )
   {
      INPUT >> buffer_size_t ;
      gcom( iv ) = --buffer_size_t ;
   }

   VERT_OF_CURRENT_CELL.re_initialize( ndp ) ;
   if( ntype == 1 ) // Edge
   {
      GE_GambitMeshing_ERROR::n4( "Edge" ) ;
   }
   else if( ntype == 2 ) // Quadrilateral
   {
      if( ndp != 4 ) GE_GambitMeshing_ERROR::n3( "Quadrilateral", 4 ) ;
      VERT_OF_CURRENT_CELL( 0 ) = gcom( 0 ) ;
      VERT_OF_CURRENT_CELL( 1 ) = gcom( 1 ) ;
      VERT_OF_CURRENT_CELL( 2 ) = gcom( 2 ) ;
      VERT_OF_CURRENT_CELL( 3 ) = gcom( 3 ) ;
   }
   else if( ntype == 3 ) // Triangle
   {
      if( ndp != 3 ) GE_GambitMeshing_ERROR::n3( "Triangle", 3 ) ;
      VERT_OF_CURRENT_CELL( 0 ) = gcom( 0 ) ;
      VERT_OF_CURRENT_CELL( 1 ) = gcom( 1 ) ;
      VERT_OF_CURRENT_CELL( 2 ) = gcom( 2 ) ;
   }
   else if( ntype == 4 ) // Brick
   {
      if( ndp != 8 ) GE_GambitMeshing_ERROR::n3( "Brick", 8 ) ;
      VERT_OF_CURRENT_CELL( 0 ) = gcom( 0 ) ;
      VERT_OF_CURRENT_CELL( 1 ) = gcom( 4 ) ;
      VERT_OF_CURRENT_CELL( 2 ) = gcom( 5 ) ;
      VERT_OF_CURRENT_CELL( 3 ) = gcom( 1 ) ;
      VERT_OF_CURRENT_CELL( 4 ) = gcom( 2 ) ;
      VERT_OF_CURRENT_CELL( 5 ) = gcom( 6 ) ;
      VERT_OF_CURRENT_CELL( 6 ) = gcom( 7 ) ;
      VERT_OF_CURRENT_CELL( 7 ) = gcom( 3 ) ;
   }
   else if( ntype == 5 ) // Wedge (Prism)
   {
      GE_GambitMeshing_ERROR::n4( "Wedge (Prism)" ) ;
   }
   else if( ntype == 6 ) // Tetrahedron
   {
      if( ndp != 4 ) GE_GambitMeshing_ERROR::n3( "Tetrahedron", 4 ) ;
      VERT_OF_CURRENT_CELL( 0 ) = gcom( 3 ) ;
      VERT_OF_CURRENT_CELL( 1 ) = gcom( 0 ) ;
      VERT_OF_CURRENT_CELL( 2 ) = gcom( 2 ) ;
      VERT_OF_CURRENT_CELL( 3 ) = gcom( 1 ) ;
   }
   else if( ntype == 7 ) // Pyramid
   {
      GE_GambitMeshing_ERROR::n4( "Pyramid" ) ;
   }
   else
   {
      GE_GambitMeshing_ERROR::n4( "" ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: check_meshes_consistency( void )
//-----------------------------------------------------------------------------
{
   if( NTYPE == 2 )
   {
      if( CELL_POLY_NAME != "GE_Quadrilateral" && 
          CELL_POLY_NAME != "GE_Rectangle" &&
          CELL_POLY_NAME != "GE_Trapezoid" &&
          SIDE_POLY_NAME != "GE_Segment" )
         GE_GambitMeshing_ERROR::n6( "\"Quadrilateral\"",
            "\"GE_Quadrilateral\", \"GE_Rectangle\" or \"GE_Trapezoid\"",
            "\"GE_Segment\"" ) ;
   }
   else if( NTYPE == 3 )
   {
      if( CELL_POLY_NAME != "GE_Triangle" &&
          SIDE_POLY_NAME != "GE_Segment" )
         GE_GambitMeshing_ERROR::n6( "\"Triangle\"",
            "\"GE_Triangle\"", "\"GE_Segment\"" ) ;
   }
   else if( NTYPE == 4 )
   {
      if( CELL_POLY_NAME != "GE_Hexahedron" &&
          CELL_POLY_NAME != "GE_Cuboid" &&
          SIDE_POLY_NAME != "GE_Quadrilateral" &&
          SIDE_POLY_NAME != "GE_Trapezoid" &&
          SIDE_POLY_NAME != "GE_Rectangle" )
         GE_GambitMeshing_ERROR::n6( "\"Tetrahedron\"",
           "\"GE_Hexahedron\" or \"GE_Cuboid\"",
           "\"GE_Quadrilateral\", \"GE_Trapezoid\" or \"GE_Rectangle\"" ) ;
   }
   else if( NTYPE == 6 )
   {
      if( CELL_POLY_NAME != "GE_Tetrahedron" &&
          SIDE_POLY_NAME != "GE_Triangle" )
         GE_GambitMeshing_ERROR::n6( "\"Brick\"",
            "\"GE_Tetrahedron\"", "\"GE_Triangle\"" ) ;
   }
   else
   {
      PEL_Error::object()->raise_internal( "invalid NTYPE" ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: read_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: read_cell" ) ;
   PEL_CHECK( valid_cell() ) ;

   read_node_to_point_coordinate_data( NTYPE ) ;

   // Cell to its sides connectivity table
   SIDE_OF_CURRENT_CELL.re_initialize( NB_SIDES_PER_CELL ) ;
 
   // conversion : Gambit local face numbering to PELICANS local face numbering

   
   if( NTYPE == 1 || NTYPE == 2 || NTYPE == 3 )
   {
      // Edge          -> GE_ReferenceSegment
      // Quadrilateral -> GE_ReferenceSquare
      // Triangle      -> GE_ReferenceTriangle
      for( size_t i=0; i<NB_SIDES_PER_CELL ; i++ )
      {
         SIDE_OF_CURRENT_CELL( i ) = SIDE_OF_CELL( I_CELL, i ) ;
      }
   }
   else if( NTYPE == 4 )
   {
      // Brick -> GE_ReferenceCube
      SIDE_OF_CURRENT_CELL( 0 ) = SIDE_OF_CELL( I_CELL, 3 ) ;
      SIDE_OF_CURRENT_CELL( 1 ) = SIDE_OF_CELL( I_CELL, 5 ) ;
      SIDE_OF_CURRENT_CELL( 2 ) = SIDE_OF_CELL( I_CELL, 1 ) ;
      SIDE_OF_CURRENT_CELL( 3 ) = SIDE_OF_CELL( I_CELL, 4 ) ;
      SIDE_OF_CURRENT_CELL( 4 ) = SIDE_OF_CELL( I_CELL, 0 ) ;
      SIDE_OF_CURRENT_CELL( 5 ) = SIDE_OF_CELL( I_CELL, 2 ) ;
   }
   else if( NTYPE == 6 )
   {
      // Tetrahedron -> GE_ReferenceTetrahedron
      SIDE_OF_CURRENT_CELL( 0 ) = SIDE_OF_CELL( I_CELL, 3 ) ;
      SIDE_OF_CURRENT_CELL( 1 ) = SIDE_OF_CELL( I_CELL, 2 ) ;
      SIDE_OF_CURRENT_CELL( 2 ) = SIDE_OF_CELL( I_CELL, 1 ) ;
      SIDE_OF_CURRENT_CELL( 3 ) = SIDE_OF_CELL( I_CELL, 0 ) ;
   }
   else PEL_Error::object()->raise_plain( "bad NTYPE" ) ;

}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: read_side( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: read_side" ) ;
   PEL_CHECK( valid_face() ) ;

   for( size_t i=0; i<VERT_OF_CURRENT_SIDE.size() ; i++ )
   {
      VERT_OF_CURRENT_SIDE( i ) = VERT_OF_SIDE( I_SIDE, i ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: go_i_th_next_line( size_t i )
//-----------------------------------------------------------------------------
{
   for( size_t jj=0 ; jj<i ; ++jj )
   {
      while( INPUT.get() != '\n' )
      {
      }
   }
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: build_sides_connectivity( size_t cell_number )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: build_sides_connectivity" ) ;
   size_t side_number = 0 ;

   if( NTYPE==2 || NTYPE==3 )
   {
      size_t imax = VERT_OF_CURRENT_CELL.size()-1 ;
      for( size_t i=0; i<imax; i++ )
      {
         VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( i ) ;
         VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( i+1 ) ;
         side_number = extend_list_of_sides() ;
         SIDE_OF_CELL( cell_number, i ) = side_number ;
      }
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( imax ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, imax ) = side_number ;
   }
   else if( NTYPE==4 )
   {
      // First side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)0 ) = side_number ;
      // Second side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 7 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 6 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)1 ) = side_number ;
      // Third side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 7 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 4 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 5 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 6 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)2 ) = side_number ;
      // Fourth side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 4 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 5 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)3 ) = side_number ;
      // Fifth side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 4 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 7 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)4 ) = side_number ;
      // Sixth side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 6 ) ;
      VERT_OF_CURRENT_SIDE( 3 ) = VERT_OF_CURRENT_CELL( 5 ) ;
      side_number = extend_list_of_sides() ;
      SIDE_OF_CELL( cell_number, (size_t)5 ) = side_number ;
   }
   else if( NTYPE==6 )
   {
      // First side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)0 ) = side_number ;
      // Second side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)1 ) = side_number ;
      // Third side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 3 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)2 ) = side_number ;
      // Fourth side
      VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 2 ) ;
      VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 1 ) ;
      VERT_OF_CURRENT_SIDE( 2 ) = VERT_OF_CURRENT_CELL( 0 ) ;
      side_number = extend_list_of_sides( ) ;
      SIDE_OF_CELL( cell_number, (size_t)3 ) = side_number ;
   }
   else
   {
      PEL_Error::object()->raise_internal( "invalid NTYPE" ) ;
   }
}

//-----------------------------------------------------------------------------
size_t
GE_GambitMeshing:: extend_list_of_sides( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: extend_list_of_sides" ) ;

   size_t result = PEL::bad_index() ;

   PEL_IndexSet* indx = PEL_IndexSet::create( 0,
                                              VERT_OF_CURRENT_SIDE,
                                              SIDE_TREE->count() ) ;
   
   // Search if side has not already been registred
   bool const found = SIDE_TREE->has( indx ) ;
   if( found )
   {
      PEL_IndexSet const* l =
         static_cast<PEL_IndexSet const*>( SIDE_TREE->item( indx ) ) ;
      result = l->id() ;
      indx->destroy() ; indx = 0 ;
   }
   else
   {
      indx->set_owner( SIDE_TREE ) ;
      result = indx->id() ;
      SIDE_TREE->extend( indx ) ;
      PEL_ASSERT( result == NB_SIDES ) ;
      for( size_t j=0 ; j<VERT_OF_CURRENT_SIDE.size() ; ++j )
      {
         VERT_OF_SIDE( result, j ) = VERT_OF_CURRENT_SIDE( j ) ;
      }
      ++NB_SIDES ;
   }
   
   PEL_CHECK( result<NB_SIDES ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_GambitMeshing:: build_colors_of_vertices( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GambitMeshing:: build_colors_of_vertices" ) ;

   vector< set< string > > colors( NB_VERTS ) ;

   // nombre de sommets par face... supposé indépendant de la face.
   size_t nbv = VERT_OF_SIDE.index_bound( 1 ) ;

   for( size_t igs=0 ; igs<NB_SIDES ; ++igs )
   {
      for( size_t iv=0 ; iv<nbv ; ++iv )
      {
         set<string>& cc = colors[ VERT_OF_SIDE( igs, iv ) ] ;
         if( COLOR_OF_SIDE[ igs ] != GE_Color::null_color() )
         {
            cc.insert( COLOR_OF_SIDE[ igs ]->name() ) ;
         }
      }
   }

   for( size_t igv=0 ; igv<NB_VERTS ; ++igv )
   {
      if( !colors[igv].empty() )
      {
         set<string>::const_iterator it = colors[igv].begin() ;
         string color_name = *it ;
         for( ++it ; it != colors[igv].end() ; ++it )
         {
            color_name += "_" + (*it) ;
         }
         GE_Color::extend( color_name ) ;
         COLOR_OF_VERT[ igv ] = GE_Color::object( color_name ) ;
      }
   }
}

//-----------------------------------------------------------------------------
bool
GE_GambitMeshing:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;

   return true ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n0( string const& fname )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl<< "GE_GambitMeshing :" << endl ;
   mesg << "   unable to open the file : \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n1( size_t gambit_dim, size_t pel_dim )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   mesg << "   inconsistent number of space dimensions" << endl ;
   mesg << "   in the GAMBIT data file   : " << gambit_dim << endl ;
   mesg << "   in the PELICANS data file : " << pel_dim ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n2( size_t itype )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   mesg << "   invalid data type (" << itype << ") in the" << endl ;
   mesg << "   boundary condition control record" << endl ;
   mesg << "   (the only one supported is 1=element/cell)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n3( string mesh, size_t ndp )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   mesg << "   the element geometry type : " << mesh << endl ;
   mesg << "   should be defined with " << ndp << " nodes" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n4( string mesh )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   if( !mesh.empty() )
   {
      mesg << "   invalid element geometry type " ;
   }
   else
   {
      mesg << "   the element geometry type : " << mesh << endl ;
      mesg << "   is not supported" ;
   } 
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n5( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   mesg << "   all the cells should have the same geometry" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_GambitMeshing_ERROR:: n6( string gmesh, string cells, string sides )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_GambitMeshing :" << endl ;
   mesg << "   for the element geometry type : " << gmesh << endl ;
   mesg << "   the sides should be of type : " << endl ;
   mesg << "      " << sides << endl ;
   mesg << "   and the cells should be of type : " << endl ;
   mesg << "      " << cells ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

