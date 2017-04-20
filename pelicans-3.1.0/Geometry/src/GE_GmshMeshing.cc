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

#include <GE_GmshMeshing.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_BalancedBinaryTreeIterator.hh>
#include <PEL_IndexSet.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Data.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

#include <intVector.hh>
#include <stringVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

using std::string ;

struct GE_GmshMeshing_ERROR
{
   static void n0( std::string const& to_test, std::string const& sol ) ;
   static void n1( size_t elem_type, size_t face_type, size_t cell_type ) ;
   static void n2( size_t read_nb_verts, size_t nb_verts ) ;
   static void n3( void ) ;
   static void n4( size_t not_ascii ) ;
} ;

//-----------------------------------------------------------------------------
GE_GmshMeshing const* GE_GmshMeshing::prototype = new GE_GmshMeshing() ;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
GE_GmshMeshing:: GE_GmshMeshing( void )
//-----------------------------------------------------------------------------
   : GE_Meshing( "GE_GmshMeshing" )
   , FILENAME()
   , INPUT( )
   , INPUT_FILE_FORMAT( undefined )
   , N_C_KEYWORDS( 4 )
   , NBDIMS( PEL::bad_index() )
   , NB_VERTS( PEL::bad_index() )
   , I_VERT( PEL::bad_index() )
   , VERT_COORD( 0 )
   , VERTS_COLORS( 0 )
   , VERT_NUMB_COR( 0 )
   , NB_FACES( PEL::bad_index() )
   , FACE_POLY_NAME( "" )
   , FACE_TYPE( PEL::bad_index() )
   , I_FACE( PEL::bad_index() )
   , VERTS_OF_CURRENT_FACE( 0 )
   , VERTS_OF_FACES( 0, 0 )
   , FACE_POL( 0 )
   , FACE_TREE( 0 )
   , FACES_COLORS( 0 )
   , EDGE_POLY_NAME( "" )
   , EDGE_TYPE( PEL::bad_index() )
   , POINT_TYPE( PEL::bad_index() )
   , BOUNDS_TREE( 0 )
   , BOUNDS_COLORS( 0 )
   , NB_CELLS( PEL::bad_index() )
   , CELL_POLY_NAME( "" )
   , CELL_TYPE( PEL::bad_index() )
   , I_CELL( PEL::bad_index() )
   , FACES_OF_CELLS( 0, 0 )
   , VERTS_OF_CURRENT_CELL( 0 )
   , FACES_OF_CURRENT_CELL( 0 )
   , CELL_POL( 0 )
{
}

//-----------------------------------------------------------------------------
GE_GmshMeshing*
GE_GmshMeshing:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_GmshMeshing* result = new GE_GmshMeshing( a_owner, exp, dim_space ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_GmshMeshing:: GE_GmshMeshing( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space )
//-----------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , FILENAME( exp->string_data( "filename" ) )
   , INPUT( exp->string_data( "filename" ).c_str(), std::ios_base::binary )
   , INPUT_FILE_FORMAT( undefined )
   , N_C_KEYWORDS( 4 )
   , NBDIMS( dim_space )
   , NB_VERTS( PEL::bad_index() )
   , I_VERT( PEL::bad_index() )
   , VERT_COORD( dim_space )
   , VERTS_COLORS( PEL_Vector:: create( this, 0 ) )
   , VERT_NUMB_COR( 0 )
   , NB_FACES( PEL::bad_index() )
   , FACE_POLY_NAME( "" )
   , FACE_TYPE( PEL::bad_index() )
   , I_FACE( PEL::bad_index() )
   , VERTS_OF_CURRENT_FACE( 0 )
   , VERTS_OF_FACES( 0, 0 )
   , FACE_POL( 0 )
   , FACE_TREE( PEL_BalancedBinaryTree:: create( this ) )
   , FACES_COLORS( PEL_Vector:: create( this, 0 ) )
   , EDGE_POLY_NAME( "" )
   , EDGE_TYPE( PEL::bad_index() )
   , POINT_TYPE( PEL::bad_index() )
   , BOUNDS_TREE( PEL_BalancedBinaryTree:: create( this ) )
   , BOUNDS_COLORS( PEL_Vector:: create( this, 0 ) )
   , NB_CELLS( PEL::bad_index() )
   , CELL_POLY_NAME( "" )
   , CELL_TYPE( PEL::bad_index() )
   , I_CELL( PEL::bad_index() )
   , FACES_OF_CELLS( 0, 0 )
   , VERTS_OF_CURRENT_CELL( 0 )
   , FACES_OF_CURRENT_CELL( 0 )
   , CELL_POL( 0 )
{
   // For pattern:
   exp->test_file( "filename", "read" ) ;
   
   initialize_rounding_strategy( exp ) ;

   if( !INPUT )
   {
     PEL_Error::object()->raise_plain(
        "*** GE_GmshMeshing error:\n"
        "    Unable to open file \""+FILENAME+"\"" ) ;
   }
   read_mesh_polyhedron( exp, dim_space, FACE_POLY_NAME, CELL_POLY_NAME ) ;
   determine_mesh_types() ;
   read_input_file_header( exp ) ;   
   
   // Read vertices
   std::string trash_string( "" ) ;
   INPUT >> trash_string ; 
   GE_GmshMeshing_ERROR:: n0( trash_string, N_C_KEYWORDS( 0 ) ) ;
   go_i_th_next_line( 1 ) ;
   INPUT >> NB_VERTS ; 
   go_i_th_next_line( 1 ) ; 
   VERT_POS = INPUT.tellg() ;

   size_t trash_size_t ;
   size_t max_numb = 0 ;

   // First loop : determine the max number of vertices
   for( size_t i=0 ; i<NB_VERTS ; ++i )
   {
      INPUT >> trash_size_t ;
      if( trash_size_t>max_numb )
         max_numb = trash_size_t ;
      go_i_th_next_line( 1 ) ;
   }
   VERT_NUMB_COR.re_initialize( max_numb+1 ) ;
   
   // Second loop : PEL numbering to Gmsh numbering correspondance table.
   INPUT.seekg( VERT_POS ) ;
   for( size_t i=0 ; i<NB_VERTS ; ++i )
   {
      INPUT >> trash_size_t ;
      VERT_NUMB_COR( trash_size_t ) = i ;
      go_i_th_next_line( 1 ) ;
   }

   INPUT >> trash_string ; 
   GE_GmshMeshing_ERROR:: n0( trash_string, N_C_KEYWORDS( 1 ) ) ;
   go_i_th_next_line( 1 ) ;
   INPUT >> trash_string ;
   
   // Read elements : faces and cells are mixed
   GE_GmshMeshing_ERROR:: n0( trash_string, N_C_KEYWORDS( 2 ) ) ;
   go_i_th_next_line( 1 ) ;
   size_t nb_tot_elems = PEL::bad_index() ;
   INPUT >> nb_tot_elems ;
   go_i_th_next_line( 1 ) ;
   ELEM_POS = INPUT.tellg() ;

   initialize_tables( nb_tot_elems ) ;
  
   size_t_array2D vertices_colors( NB_VERTS, 
                                   2*CELL_POL->nb_vertices()+1 ) ;
   vertices_colors.set( 0 ) ;
   size_t elem_type = PEL::bad_index() ;
   size_t elem_nb = PEL::bad_index() ;
   size_t ic = 0 ;
   size_t ib = 0 ;
   for( size_t i=0 ; i<nb_tot_elems ; ++i )
   {
      PEL_ASSERT( INPUT ) ;
      INPUT >> elem_nb >> elem_type ;
      if( elem_type==POINT_TYPE )
      {
         read_one_point() ;
      }
      else if( elem_type==FACE_TYPE )
      {
         read_one_bound( ib, vertices_colors ) ;
         ++ib ;
      }
      else if( elem_type==CELL_TYPE )
      {
         preread_one_cell( ic, vertices_colors ) ;
         ++ic ;
      }
      else if( elem_type==EDGE_TYPE )
      {
         read_one_edge( ) ;
      }
      else
      {
         GE_GmshMeshing_ERROR:: n1( elem_type, FACE_TYPE, CELL_TYPE ) ;
      }
      go_i_th_next_line( 1 ) ;
   }
   PEL_ASSERT( ic == NB_CELLS ) ;
   INPUT >> trash_string ;
   GE_GmshMeshing_ERROR:: n0( trash_string, N_C_KEYWORDS( 3 ) ) ;
   go_i_th_next_line( 1 ) ;

   attribute_colors_to_vertices( vertices_colors ) ;
}

//------------------------------------------------------------------------
GE_GmshMeshing:: ~GE_GmshMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_GmshMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: nb_vertices" ) ;

   return NB_VERTS ;
}

//------------------------------------------------------------------------
size_t
GE_GmshMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: nb_cells" ) ;

   return NB_CELLS ;
}

//------------------------------------------------------------------------
size_t
GE_GmshMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: nb_faces" ) ;

   return NB_FACES ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   INPUT.seekg( VERT_POS ) ;
   I_VERT = 0 ;
   read_vertex() ;
}

//------------------------------------------------------------------------
bool
GE_GmshMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( I_VERT<NB_VERTS ) ;   
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   ++I_VERT ;
   if( valid_vertex() ) read_vertex() ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_GmshMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
   doubleVector const& result = VERT_COORD ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   I_CELL = 0 ;
   INPUT.seekg( ELEM_POS ) ;
   read_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_GmshMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: valid_cell" ) ;

   bool result = ( I_CELL<NB_CELLS )  ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   ++I_CELL ;
   if( valid_cell() )
   {
      read_cell() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_GmshMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_GmshMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = VERTS_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_GmshMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = FACES_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   I_FACE = 0 ;
   read_face() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_GmshMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: valid_face" ) ;

   bool result = ( I_FACE<NB_FACES ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   ++I_FACE ;
   if( valid_face() )
   {
      read_face() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_GmshMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = FACE_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_GmshMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = VERTS_OF_CURRENT_FACE ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   os << s << "filename: \"" << FILENAME << "\"" << std::endl ;
   os << s << "format: " ;
   if( INPUT_FILE_FORMAT==one )
   {
      os << "\"1.0\"" ;
   }
   else if( INPUT_FILE_FORMAT==two )
   {
      os << "\"2.0\"" ;
   }
   else
   {
      os << "???" ;
   }
   os << std::endl ;
   os << s << "nb_cells: " << NB_CELLS << std::endl ;
   os << s << "nb_faces: " << NB_FACES << std::endl ;
   os << s << "nb_vertices: " << NB_VERTS << std::endl ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_GmshMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = GE_Color::null_color() ; 
   if( VERTS_COLORS->at( I_VERT )!=0 )
   {
      result = static_cast<GE_Color const*>( VERTS_COLORS->at( I_VERT ) ) ;
   }

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_GmshMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = CELL_COLOR ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_GmshMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = FACE_COLOR ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_vertex( void )
//------------------------------------------------------------------------
{
   PEL_CHECK( valid_vertex() ) ;   
   PEL_CHECK_INV( invariant() ) ;
   PEL_ASSERT( !INPUT.eof() ) ;

   double x0, x1, x2 ;
   size_t trash_size_t ;
   INPUT >> trash_size_t >> x0 >> x1 >> x2 ;
   VERT_COORD( 0 ) = roundoff( x0 ) ;
   if( NBDIMS>1 )
   {
      VERT_COORD( 1 ) = roundoff( x1 ) ;
      if( NBDIMS>2 )
      {
         VERT_COORD( 2 ) = roundoff( x2 ) ;
      }
      PEL_CHECK( IMPLIES( NBDIMS==2, 
                          PEL::double_equality( x2, 0., PEL::epsilon_double(), 
                                                        PEL::epsilon_double() ) ) ) ;
   }
   else
   {
      PEL_CHECK( IMPLIES( NBDIMS==1, 
                          PEL::double_equality( x2, 0., PEL::epsilon_double(), 
                                                        PEL::epsilon_double() )
                       && PEL::double_equality( x1, 0., PEL::epsilon_double(), 
                                                        PEL::epsilon_double() ) ) ) ;
   }
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_cell( void )
//------------------------------------------------------------------------
{
   PEL_CHECK( valid_cell() ) ;

   size_t elem_type = PEL::bad_index() ;
   size_t trash_size_t = PEL::bad_index() ;
   while( !( INPUT >> trash_size_t >> elem_type && elem_type==CELL_TYPE ) )
   {
      go_i_th_next_line( 1 ) ;
   }
   
   // Colors and tags
   size_t c_col = PEL::bad_index() ;
   size_t const nb_vert = CELL_POL->nb_vertices() ;
   read_mesh_colors_and_tags( nb_vert, c_col ) ;

   // Connectivities
   for( size_t i=0; i<nb_vert; ++i )
   {
      INPUT >> trash_size_t ;
      VERTS_OF_CURRENT_CELL( i ) = VERT_NUMB_COR( trash_size_t ) ;
   }

   CELL_COLOR = provide_color( c_col ) ;

   // Cell to its faces connectivity table
   FACES_OF_CURRENT_CELL.re_initialize( CELL_POL->nb_faces() ) ;
   for( size_t i=0; i<CELL_POL->nb_faces() ; i++ )
   {
      FACES_OF_CURRENT_CELL( i ) = FACES_OF_CELLS( I_CELL, i ) ;
   }
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_face( void )
//------------------------------------------------------------------------
{
   PEL_CHECK( valid_face() ) ;
   for( size_t i=0; i<VERTS_OF_CURRENT_FACE.size() ; i++ )
   {
      VERTS_OF_CURRENT_FACE( i ) = VERTS_OF_FACES( I_FACE, i ) ;
   }
   FACE_COLOR =  static_cast<GE_Color const*>(FACES_COLORS->at( I_FACE ) ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: determine_mesh_types( void ) 
//------------------------------------------------------------------------
{
   CELL_POL = GE_Mpolyhedron::reference_polyhedron( CELL_POLY_NAME ) ;
   FACE_POL = GE_Mpolyhedron::reference_polyhedron( FACE_POLY_NAME ) ;
   if( CELL_POLY_NAME == "GE_Cuboid" || CELL_POLY_NAME == "GE_Hexahedron" )
   {
      CELL_TYPE = 5 ;
      FACE_TYPE = 3 ;
      EDGE_TYPE = 1 ;
      POINT_TYPE = 15 ;
      EDGE_POLY_NAME = "GE_Segment" ;
      EDGE_POL = GE_Mpolyhedron::reference_polyhedron( EDGE_POLY_NAME ) ;
      VERTS_OF_CURRENT_FACE.re_initialize( 4 ) ;
      VERTS_OF_CURRENT_CELL.re_initialize( 8 ) ;
      FACES_OF_CURRENT_CELL.re_initialize( 6 ) ;   
   }
   else if( CELL_POLY_NAME == "GE_Tetrahedron" )
   {
      CELL_TYPE = 4 ;
      FACE_TYPE = 2 ;
      EDGE_TYPE = 1 ;
      POINT_TYPE = 15 ;
      EDGE_POLY_NAME = "GE_Segment" ;
      EDGE_POL = GE_Mpolyhedron::reference_polyhedron( EDGE_POLY_NAME ) ;
      VERTS_OF_CURRENT_FACE.re_initialize( 3 ) ;
      VERTS_OF_CURRENT_CELL.re_initialize( 4 ) ;
      FACES_OF_CURRENT_CELL.re_initialize( 4 ) ;   
   }
   else if( CELL_POLY_NAME == "GE_Quadrilateral" || 
            CELL_POLY_NAME == "GE_Rectangle" || 
            CELL_POLY_NAME == "GE_Trapezoid" )
   {
      CELL_TYPE = 3 ;
      FACE_TYPE = 1 ;
      POINT_TYPE = 15 ;
      VERTS_OF_CURRENT_FACE.re_initialize( 2 ) ;
      VERTS_OF_CURRENT_CELL.re_initialize( 4 ) ;
      FACES_OF_CURRENT_CELL.re_initialize( 4 ) ;   
   }
   else if( CELL_POLY_NAME == "GE_Triangle" )
   {

      CELL_TYPE = 2 ;
      FACE_TYPE = 1 ;
      POINT_TYPE = 15 ;
      VERTS_OF_CURRENT_FACE.re_initialize( 2 ) ;
      VERTS_OF_CURRENT_CELL.re_initialize( 3 ) ;
      FACES_OF_CURRENT_CELL.re_initialize( 3 ) ;   
   }
   else
   {
      PEL_CHECK( CELL_POLY_NAME == "GE_Segment" ) ;
      CELL_TYPE = 1 ;
      FACE_TYPE = 15 ;
      VERTS_OF_CURRENT_FACE.re_initialize( 1 ) ;
      VERTS_OF_CURRENT_CELL.re_initialize( 2 ) ;
      FACES_OF_CURRENT_CELL.re_initialize( 2 ) ;   
   }   
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: initialize_tables( size_t nb_tot_elems )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: initialize_tables" ) ;

   // Loop on element to determine the number of cells and faces
   size_t elem_type = PEL::bad_index() ;
   size_t elem_nb = PEL::bad_index() ;
   NB_CELLS = 0 ;
   size_t nb_bounds = 0 ;
   for( size_t i=0 ; i<nb_tot_elems ; ++i )
   {
      PEL_ASSERT( INPUT ) ;
      INPUT >> elem_nb >> elem_type ;

      if( elem_type==CELL_TYPE )
      {
         ++NB_CELLS ;
      }
      else if( elem_type==FACE_TYPE )
      {
         ++nb_bounds ;
      }
      else if( !( elem_type==EDGE_TYPE || elem_type==POINT_TYPE ) )
      {
         GE_GmshMeshing_ERROR:: n1( elem_type, FACE_TYPE, CELL_TYPE ) ;
      }
      go_i_th_next_line( 1 ) ;
   }

   PEL_ASSERT( nb_tot_elems>=NB_CELLS ) ;
 
   // Initializations :
   VERTS_COLORS->re_initialize( NB_VERTS ) ;

   FACES_OF_CELLS.re_initialize( NB_CELLS, CELL_POL->nb_faces() ) ;
   VERTS_OF_FACES.re_initialize( NB_CELLS*CELL_POL->nb_faces(),
                                 FACE_POL->nb_vertices() ) ;
   FACES_COLORS->re_initialize( NB_CELLS*CELL_POL->nb_faces() ) ;
   BOUNDS_COLORS->re_initialize( nb_bounds ) ;
   NB_FACES = 0 ;

   // Back to first element description
   INPUT.seekg( ELEM_POS ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_one_bound( size_t bound_number,
                                 size_t_array2D& vertices_colors )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: read_one_bound" ) ;

   // Pure readings
   size_t trash_size_t ;
   size_t b_col ;
   size_t const nb_vert = FACE_POL->nb_vertices() ;
   read_mesh_colors_and_tags( nb_vert, b_col ) ;
   
   for( size_t i=0; i<nb_vert; ++i )
   {
      INPUT >> trash_size_t ;
      VERTS_OF_CURRENT_FACE( i ) = VERT_NUMB_COR( trash_size_t ) ;

      if( b_col!=0 )
      {
         // If vertex is only cell colored it turns to be bound colored
         if( vertices_colors( VERTS_OF_CURRENT_FACE( i ), 0 )==0 )
         {
             vertices_colors( VERTS_OF_CURRENT_FACE( i ), 0 ) = 1 ;
             for(  size_t d=1; d<vertices_colors.index_bound( 1 ); ++d )
             {
                vertices_colors( VERTS_OF_CURRENT_FACE( i ), d ) = 0 ;
             }
         }
         // Concatenate all bounds colors
         for( size_t d=1; d<vertices_colors.index_bound( 1 ); ++d )
         {
            if( vertices_colors( VERTS_OF_CURRENT_FACE( i ), d )==0 )
            {
               vertices_colors( VERTS_OF_CURRENT_FACE( i ), d ) = b_col ;
               break ;
            }
         }
      }
   }

   // Tree storage
   PEL_IndexSet* indx = PEL_IndexSet::create( BOUNDS_TREE,
                                              VERTS_OF_CURRENT_FACE,
                                              bound_number ) ;
   BOUNDS_TREE->extend( indx ) ;
   BOUNDS_COLORS->set_at( bound_number, provide_color( b_col ) ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_one_edge( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: read_one_edge" ) ;

   // Pure readings
   size_t e_col ;
   read_mesh_colors_and_tags( EDGE_POL->nb_vertices(), e_col ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: read_one_point( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: read_one_point" ) ;

   // Pure readings
   size_t p_col ;
   read_mesh_colors_and_tags( 1, p_col ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: preread_one_cell( size_t cell_number,
                                   size_t_array2D& vertices_colors )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: preread_one_cell" ) ;

   // Pure readings
   size_t trash_size_t ;
   size_t c_col ;
   size_t const nb_vert = CELL_POL->nb_vertices() ;
   read_mesh_colors_and_tags( nb_vert, c_col ) ;

   for( size_t i=0; i<nb_vert; ++i )
   {
      INPUT >> trash_size_t ;
      VERTS_OF_CURRENT_CELL( i ) = VERT_NUMB_COR( trash_size_t ) ;
      if( c_col!=0 )
      {
         // If vertex is only cell colored it turns to be bound colored
         if( VERTS_OF_CURRENT_CELL( i )>=vertices_colors.index_bound( 0 ) )
         {
            std::cout << " vertices_colors.index_bound( 0 ) " 
                      << vertices_colors.index_bound( 0 )
                      << " VERTS_OF_CURRENT_CELL( i ) " 
                      << VERTS_OF_CURRENT_CELL( i ) << std::endl ;
         }
         if( vertices_colors( VERTS_OF_CURRENT_CELL( i ), 0 )==0 )
         {


            // Concatenate all bounds colors
            for( size_t d=1; d<vertices_colors.index_bound( 1 ); ++d )
            {
               if( vertices_colors( VERTS_OF_CURRENT_CELL( i ), d )==0 )
               {
                  vertices_colors( VERTS_OF_CURRENT_CELL( i ), d ) = c_col ;
                  break ;
               }
            }
         }
      }
   }
   // Building of faces.
   faces_connectivity_building( cell_number ) ;
}

//-----------------------------------------------------------------------------
void
GE_GmshMeshing:: faces_connectivity_building( size_t cell_number )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: faces_connectivity_building" ) ;
   
   size_t face_number = 0 ;
   size_t NB_VERTS_PER_CELL = CELL_POL->nb_vertices() ;

   if( CELL_POLY_NAME=="GE_Segment" )
   {
      PEL_CHECK( FACE_POLY_NAME=="GE_Mpoint" ) ;
      for( size_t i=0; i<NB_VERTS_PER_CELL; i++ )
      {
         VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( i ) ;
         face_number = extend_list_of_faces() ;
         FACES_OF_CELLS( cell_number, i ) = face_number ;
      }
   }
   else if( CELL_POLY_NAME=="GE_Triangle"  || CELL_POLY_NAME=="GE_Quadrilateral" || 
       CELL_POLY_NAME=="GE_Rectangle" || CELL_POLY_NAME=="GE_Rectangle" )
   {
      PEL_CHECK( FACE_POLY_NAME=="GE_Segment" ) ;
      for( size_t i=0; i<NB_VERTS_PER_CELL-1; i++ )
      {
         VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( i ) ;
         VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( i+1 ) ;
         face_number = extend_list_of_faces() ;
         FACES_OF_CELLS( cell_number, i ) = face_number ;
      }
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( NB_VERTS_PER_CELL-1 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      face_number = extend_list_of_faces( ) ;
      FACES_OF_CELLS( cell_number, NB_VERTS_PER_CELL-1 ) = face_number ;
   }
   else if( CELL_POLY_NAME=="GE_Tetrahedron" )
   {
      PEL_CHECK( FACE_POLY_NAME=="GE_Triangle" ) ;
      // First face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      face_number = extend_list_of_faces( ) ;
      FACES_OF_CELLS( cell_number, 0 ) = face_number ;
      // Second face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 1 ) = face_number ;
      // Third face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 2 ) = face_number ;
      // Fourth face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 3 ) = face_number ;
   }
   else if( CELL_POLY_NAME=="GE_Hexahedron" || CELL_POLY_NAME=="GE_Cuboid" )
   {
      PEL_CHECK( FACE_POLY_NAME=="GE_Quadrilateral" || 
                 FACE_POLY_NAME=="GE_Trapezoid" ||
                 FACE_POLY_NAME=="GE_Rectangle" ) ;
      // First face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 5 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 4 ) ;
      face_number = extend_list_of_faces( ) ;
      FACES_OF_CELLS( cell_number, 0 ) = face_number ;
      // Second face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 6 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 5 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 1 ) = face_number ;
      // Third face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 7 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 6 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 2 ) = face_number ;
      // Fourth face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 4 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 7 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 3 ) = face_number ;
      // Fifth face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 0 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 3 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 2 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 1 ) ;
      face_number = extend_list_of_faces(  ) ;
      FACES_OF_CELLS( cell_number, 4 ) = face_number ;
      // Sixth face
      VERTS_OF_CURRENT_FACE( 0 ) = VERTS_OF_CURRENT_CELL( 4 ) ;
      VERTS_OF_CURRENT_FACE( 1 ) = VERTS_OF_CURRENT_CELL( 5 ) ;
      VERTS_OF_CURRENT_FACE( 2 ) = VERTS_OF_CURRENT_CELL( 6 ) ;
      VERTS_OF_CURRENT_FACE( 3 ) = VERTS_OF_CURRENT_CELL( 7 ) ;
      face_number = extend_list_of_faces( ) ;
      FACES_OF_CELLS( cell_number, 5 ) = face_number ;
   }
}

//-----------------------------------------------------------------------------
size_t
GE_GmshMeshing:: extend_list_of_faces( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: extend_list_of_faces" ) ;

   size_t result = PEL::bad_index() ;
   PEL_IndexSet* indx = PEL_IndexSet::create( 0,
                                              VERTS_OF_CURRENT_FACE,
                                              FACE_TREE->count() ) ;
   bool const found = FACE_TREE->has( indx ) ;
   if( found )
   {
      // Face already created
      PEL_IndexSet const* l =
         static_cast<PEL_IndexSet const*>( FACE_TREE->item( indx ) ) ;
      result = l->id() ;
      indx->destroy() ; indx = 0 ;
   }
   else
   {
      // Create new face
      indx->set_owner( FACE_TREE ) ;
      result = indx->id() ;
      FACE_TREE->extend( indx ) ;
      for( size_t d=0; d<VERTS_OF_CURRENT_FACE.size(); ++d )
      {
         VERTS_OF_FACES( NB_FACES, d ) = VERTS_OF_CURRENT_FACE( d ) ;
      }

      // Determine its color
      FACES_COLORS->set_at( result, 
                            const_cast<GE_Color*>( GE_Color::null_color() ) ) ;
      if( BOUNDS_TREE->has( indx ) )
      {
         PEL_IndexSet const* l =
            static_cast<PEL_IndexSet const*>( BOUNDS_TREE->item( indx ) ) ;
         FACES_COLORS->set_at( result, BOUNDS_COLORS->at( l->id() ) ) ;
      }
      ++NB_FACES ;
   }

   PEL_CHECK( result<NB_FACES ) ;

   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: attribute_colors_to_vertices( size_t_array2D const&
                                               vertices_colors )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_GmshMeshing:: attribute_colors_to_vertices" ) ;

   size_t_vector nb( vertices_colors.index_bound( 1 ) ) ;

   for( size_t i=0; i<NB_VERTS; ++i )
   {
      vertices_colors.extract_section( 0, i, nb ) ;
      VERTS_COLORS->set_at( i, provide_color( nb ) ) ;
   }
}

//------------------------------------------------------------------------
bool
GE_GmshMeshing:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;

   return true ;
}

//------------------------------------------------------------------------
void
GE_GmshMeshing:: go_i_th_next_line( size_t i )
//------------------------------------------------------------------------
{
   for( size_t jj=0 ; jj<i ; ++jj )
   {
      while( INPUT.get() != '\n' )
      {
      }
   }
}

//------------------------------------------------------------------------
GE_Color*
GE_GmshMeshing:: provide_color( size_t col ) const
//------------------------------------------------------------------------
{
   std::ostringstream oss ;
   oss << "r" ;
   oss << col ;
   std::string scol = oss.str() ;
   GE_Color::extend( scol ) ;

   return( const_cast<GE_Color*>( GE_Color::object( scol ) ) ) ;
}

//------------------------------------------------------------------------
GE_Color*
GE_GmshMeshing:: provide_color( size_t_vector const& cols  ) const
//------------------------------------------------------------------------
{
   std::ostringstream oss ;
   GE_Color* result = 0 ;

   size_t_vector x( cols ) ;
   x.remove_at( 0 ) ;
   x.sort_increasingly() ;
   for( size_t i=0; i<x.size(); ++i )
   {
      if( ( i==0 && x( i )!=0 ) || ( i>0 && x(i)>x(i-1) ) )
      {
         oss << "r" ;
         oss << x( i ) ;
      }
   }
   std::string scol = oss.str() ;
   if( scol.size()!= 0 )
   {
      GE_Color::extend( scol ) ;
      result = const_cast<GE_Color*>( GE_Color::object( scol ) ) ;
   }
   return( result ) ;
}

//------------------------------------------------------------------------
void 
GE_GmshMeshing:: read_input_file_header( PEL_ModuleExplorer const* exp )
//------------------------------------------------------------------------
{
   std::string const& ext_str = exp->string_data( "format" ) ;
   exp->test_data_in( "format", "1.0,2.0" ) ;
   if( ext_str=="1.0" )
   {
      INPUT_FILE_FORMAT = one ;
      N_C_KEYWORDS( 0 ) = "$NOD" ;
      N_C_KEYWORDS( 1 ) = "$ENDNOD" ;
      N_C_KEYWORDS( 2 ) = "$ELM" ;
      N_C_KEYWORDS( 3 ) = "$ENDELM" ;      
   }
   else if( ext_str=="2.0" )
   {
      INPUT_FILE_FORMAT = two ;
      N_C_KEYWORDS( 0 ) = "$Nodes" ;
      N_C_KEYWORDS( 1 ) = "$EndNodes" ;
      N_C_KEYWORDS( 2 ) = "$Elements" ;
      N_C_KEYWORDS( 3 ) = "$EndElements" ;
      std::string trash_string( "" ) ;  
      INPUT >> trash_string ;
      GE_GmshMeshing_ERROR:: n0( trash_string, "$MeshFormat" ) ;
      double trash_double ;
      INPUT >> trash_double ;
      size_t trash_size_t ;
      INPUT >> trash_size_t ;
      if( trash_size_t!=0 )
         GE_GmshMeshing_ERROR:: n4( trash_size_t ) ;
      INPUT >> trash_size_t ;
      INPUT >> trash_string ;
      GE_GmshMeshing_ERROR:: n0( trash_string, "$EndMeshFormat" ) ;
    }
   else 
   {
      GE_GmshMeshing_ERROR:: n3() ;
   }
}

//------------------------------------------------------------------------
void 
GE_GmshMeshing:: read_mesh_colors_and_tags( size_t nb_vert, size_t& c_col )
//------------------------------------------------------------------------
{
   size_t trash_size_t = PEL::bad_index() ;
   if( INPUT_FILE_FORMAT == one )
   {
      INPUT >> trash_size_t >> c_col >> trash_size_t ;
      if( trash_size_t!=nb_vert )
      {
         GE_GmshMeshing_ERROR:: n2( trash_size_t, nb_vert ) ;
      }
   }
   else
   {
      size_t nb_tags ;
      INPUT >> nb_tags ;
      if( nb_tags > 1 )
      {
         INPUT >> trash_size_t >> c_col ;
         for( size_t i=2 ; i<nb_tags ; ++i )
         {
            INPUT >> trash_size_t ;
         }
      }
      else
      {
         c_col = 0 ;
      }
   }
}

//Internal----------------------------------------------------------------
void 
GE_GmshMeshing_ERROR:: n0( std::string const& to_test,
                           std::string const& sol )
//Internal----------------------------------------------------------------
{
   if( to_test!=sol )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_GmshMeshing syntax error:\n"
         "    reading unexpected data: \""+to_test+"\"\n"
         "    (expected: \""+sol+"\")" ) ;
   }
}

   
//Internal--------------------------------------------------------------------
void 
GE_GmshMeshing_ERROR:: n1( size_t elem_type,
                           size_t face_type, size_t cell_type )
//Internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_GmshMeshing syntax error:\n" 
        << "    reading element type: \"" << elem_type << "\"\n"
        << "    while expecting \"" << face_type << "\" for faces\n"
        << "                and \"" << cell_type << "\" for cells" ;

   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//Internal--------------------------------------------------------------------
void 
GE_GmshMeshing_ERROR:: n2( size_t read_nb_verts, size_t nb_verts )
//Internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_GmshMeshing syntax error:" << std::endl
        << "    Error reading element vertices" << std::endl
        << "       number of vertices read: " << read_nb_verts << std::endl
        << "       expecting: " << nb_verts ;

   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//Internal--------------------------------------------------------------------
void 
GE_GmshMeshing_ERROR:: n3( void )
//Internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_GmshMeshing syntax error:\n" 
        << "    A format must be specified in the hierarchical data structure\n"
        << "    provided to GE_GmshMeshing objects by specifying a string data\n"
        << "    named \"format\" that can either take the values \"1.0\" or \"2.0\"\n "
        << "    see the Gmsh documentation section file formats." ;

   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//Internal--------------------------------------------------------------------
void 
GE_GmshMeshing_ERROR:: n4( size_t not_ascii )
//Internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_GmshMeshing syntax error:\n" 
        << "    only ascii format is supported,\n"
        << "    the first integer of the input file is: " << not_ascii  
        << std::endl
        << "    (a value of 0 is expected)." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
