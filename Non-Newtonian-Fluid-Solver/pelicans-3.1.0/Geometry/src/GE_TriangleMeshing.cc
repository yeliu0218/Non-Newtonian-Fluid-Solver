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

#include <GE_TriangleMeshing.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>

#include <size_t_vector.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::endl ; using std::string ; using std::ostringstream ;

GE_TriangleMeshing const* 
GE_TriangleMeshing::prototype = new GE_TriangleMeshing() ;

struct GE_TriangleMeshing_ERROR
{
   static void n0( string const& fname ) ;
   static void n1( size_t dim, size_t space_dim ) ;
   static void n3( size_t nb_pt ) ;
   static void n4( void ) ;
} ;

//-----------------------------------------------------------------------------
GE_TriangleMeshing:: GE_TriangleMeshing( void )
//-----------------------------------------------------------------------------
   : GE_Meshing( "GE_TriangleMeshing" )
   , NB_VERTS( PEL::bad_index() )
   , COLOR_OF_VERT( 0 )
   , COLOR_OF_VERT_INDEX( 0 )
   , I_VERT( PEL::bad_index() )
   , COORD_OF_CURRENT_VERT( 0 )
   , FROM_ZERO( false )
   , NB_AT( PEL::bad_index() )
   , HAS_BM( PEL::bad_index() )
   , GRID_VERTEX( 0 )
   , VERTEX_CORRECT_NB( 0 )
   , NB_CELLS( PEL::bad_index() )
   , SIDE_OF_CELL( 0, 0 )
   , COLOR_OF_CELL( 0 )
   , I_CELL( PEL::bad_index() )
   , VERT_OF_CURRENT_CELL( 0 )
   , SIDE_OF_CURRENT_CELL( 0 )
   , SIDE_TREE(0)
   , NB_SIDES( PEL::bad_index() )
   , VERT_OF_SIDE( 0, 0 )
   , COLOR_OF_SIDE( 0 )
   , I_SIDE( PEL::bad_index() )
   , VERT_OF_CURRENT_SIDE( 0 )
{
}

//-----------------------------------------------------------------------------
GE_TriangleMeshing*
GE_TriangleMeshing:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp,
                                     size_t dim_space ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_TriangleMeshing* result 
                         = new GE_TriangleMeshing( a_owner, exp, dim_space ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}
//-----------------------------------------------------------------------------
GE_TriangleMeshing:: GE_TriangleMeshing( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp, 
                                         size_t dim_space )
//-----------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , VERT_INPUT( exp->string_data( "filename_for_nodes" ).c_str(),
                 std::ios_base::binary )
   , CELL_INPUT( exp->string_data( "filename_for_cells" ).c_str(),
                 std::ios_base::binary )
   , NB_VERTS( 0 )
   , COLOR_OF_VERT( 0 )
   , COLOR_OF_VERT_INDEX( 0 )
   , I_VERT( 0 )
   , COORD_OF_CURRENT_VERT( dim_space )
   , FROM_ZERO( false )
   , NB_AT( PEL::bad_index() )
   , HAS_BM( PEL::bad_index() )
   , GRID_VERTEX( 0 )
   , VERTEX_CORRECT_NB( 0 )
   , NB_CELLS( 0 )
   , SIDE_OF_CELL( 0, 0 )
   , COLOR_OF_CELL( 0 )
   , I_CELL( 0 )
   , VERT_OF_CURRENT_CELL( 0 )
   , SIDE_OF_CURRENT_CELL( 0 )
   , SIDE_TREE( PEL_BalancedBinaryTree:: create( this ) )
   , NB_SIDES( 0 )
   , VERT_OF_SIDE( 0, 0 )
   , COLOR_OF_SIDE( 0 )
   , I_SIDE( 0 )
   , VERT_OF_CURRENT_SIDE( 0 )
{
   PEL_LABEL( "GE_TriangleMeshing:: GE_TriangleMeshing" ) ;

   initialize_rounding_strategy( exp ) ;

   size_t buffer_size_t = 0 ;

   // -- Prepare reading of vertices :
   std::string const& nname = exp->string_data( "filename_for_nodes" ) ;
   if( !VERT_INPUT ) GE_TriangleMeshing_ERROR::n0( nname ) ;

   // First line : nb_vertices/nb_dim/nb_attributes/nb_boundary_markers
   size_t dim ;
   VERT_INPUT >> NB_VERTS >> dim >> NB_AT >> buffer_size_t ;
   if( dim!=dim_space || dim!=2 )
      GE_TriangleMeshing_ERROR::n1(dim,dim_space) ;
   HAS_BM = ( buffer_size_t==1 ) ;
   VERT_POS = VERT_INPUT.tellg() ;
   COLOR_OF_VERT.resize( NB_VERTS, GE_Color::null_color() ) ;
   COLOR_OF_VERT_INDEX.re_initialize( NB_VERTS ) ;
   GRID_VERTEX.resize( NB_VERTS ) ;
   VERTEX_CORRECT_NB.resize( NB_VERTS ) ;
   VERT_INPUT >> buffer_size_t ;

   FROM_ZERO = ( buffer_size_t==0 ) ;
   if( !FROM_ZERO && buffer_size_t!=1 ) 
   {
      GE_TriangleMeshing_ERROR::n4() ;   
   }

   // -- Prepare reading of cells :
   std::string const& cname = exp->string_data( "filename_for_cells" ) ;
   if( !CELL_INPUT ) GE_TriangleMeshing_ERROR::n0( cname ) ;

   read_mesh_polyhedron( exp, dim_space, SIDE_POLY_NAME, CELL_POLY_NAME ) ;

   //      First line : nb_triangles/nb_points_per_triangle/nb_attributes
   size_t nb_pts ;
   CELL_INPUT >> NB_CELLS >> nb_pts >> buffer_size_t ;
   if( nb_pts!=3 ) GE_TriangleMeshing_ERROR::n3( nb_pts ) ; 
   CELL_POS = CELL_INPUT.tellg() ;
   remove_outofgrid_vertices() ;
   CELL_INPUT.clear() ;
   CELL_INPUT.seekg( CELL_POS ) ;
   PEL_ASSERT( !CELL_INPUT.fail() ) ;

   VERT_OF_CURRENT_SIDE.re_initialize( 2 ) ;
   VERT_OF_SIDE.re_initialize( 3*NB_CELLS, 2 ) ;
   SIDE_OF_CELL.re_initialize( NB_CELLS, 3 ) ;

   for( size_t igc=0 ; igc<NB_CELLS ; ++igc )
   {
      read_node_to_point_coordinate_data() ;
      build_sides_connectivity( igc ) ;
   }

   I_VERT = NB_VERTS+1 ;
   I_SIDE = NB_SIDES+1 ;
   I_CELL = NB_CELLS+1 ;


}

//-----------------------------------------------------------------------------
GE_TriangleMeshing:: ~GE_TriangleMeshing( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_TriangleMeshing:: nb_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: nb_vertices" ) ;
   return NB_VERTS ;
}

//-----------------------------------------------------------------------------
size_t
GE_TriangleMeshing:: nb_cells( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: nb_cells" ) ;
   return NB_CELLS ;
}

//-----------------------------------------------------------------------------
size_t
GE_TriangleMeshing:: nb_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: nb_faces" ) ;
   return NB_SIDES ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: start_vertex_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   VERT_INPUT.seekg( VERT_POS ) ;
   I_VERT = 0 ;
   read_vertex() ;
}

//-----------------------------------------------------------------------------
bool
GE_TriangleMeshing:: valid_vertex( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return ( I_VERT<GRID_VERTEX.size() ) ;   
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: go_next_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++I_VERT ;
   if( valid_vertex() ) read_vertex() ;
}

//-----------------------------------------------------------------------------
doubleVector const& 
GE_TriangleMeshing:: vertex_coordinates( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   doubleVector const& result = COORD_OF_CURRENT_VERT ;

   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: start_cell_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   I_CELL = 0 ;
   CELL_INPUT.seekg( CELL_POS ) ;
   read_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_TriangleMeshing:: valid_cell( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: valid_cell" ) ;

   bool result = ( I_CELL< NB_CELLS )  ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: go_next_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   ++I_CELL ;
   if( valid_cell() ) read_cell() ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_TriangleMeshing:: cell_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_TriangleMeshing:: cell_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = VERT_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_TriangleMeshing:: cell_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = SIDE_OF_CURRENT_CELL ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: start_face_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   I_SIDE = 0 ;
   read_side() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_TriangleMeshing:: valid_face( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: valid_face" ) ;

   bool result = ( I_SIDE<NB_SIDES ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: go_next_face( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   ++I_SIDE ;
   if( valid_face() )
   {
      read_side() ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_TriangleMeshing:: face_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = SIDE_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_TriangleMeshing:: face_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = VERT_OF_CURRENT_SIDE ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_TriangleMeshing:: default_vertex_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = COLOR_OF_VERT[ I_VERT ] ;

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_TriangleMeshing:: default_cell_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = GE_Color::null_color() ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_TriangleMeshing:: default_face_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = COLOR_OF_SIDE[ I_SIDE ] ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: read_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: read_vertex" ) ;
   PEL_CHECK( valid_vertex() ) ;   
   PEL_CHECK_INV( invariant() ) ;
   PEL_ASSERT( !VERT_INPUT.eof() ) ;

   size_t buffer_size_t ;
   VERT_INPUT >> buffer_size_t ;

   // vertices coordinates
   double x0, x1 ;
   VERT_INPUT >> x0 >> x1 ;

   // Attributes -> to trash
   double trash_double ;
   for( size_t i=0; i<NB_AT; ++i )
   {
      VERT_INPUT >> trash_double ;
   }

   size_t color = 0 ;
   if( HAS_BM )
   {
      VERT_INPUT >> color ;
   }

   if( GRID_VERTEX( I_VERT ) )
   {
      COORD_OF_CURRENT_VERT( 0 ) = roundoff( x0 ) ;
      COORD_OF_CURRENT_VERT( 1 ) = roundoff( x1 ) ;

      string str = "r" ;
      ostringstream oss_color;
      oss_color << color;
      str += oss_color.str();      
      GE_Color::extend( str ) ;

      COLOR_OF_VERT[ I_VERT ] = GE_Color::object( str ) ;
      COLOR_OF_VERT_INDEX( I_VERT ) = color ;
   }
   else
   {
      go_next_vertex() ;
   }
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: read_node_to_point_coordinate_data( void )
//-----------------------------------------------------------------------------
{
   // Each line : Cell nb / Vert nb / Vert nb / Vert nb 
    size_t buffer_size_t ;
    CELL_INPUT >> buffer_size_t ;

    VERT_OF_CURRENT_CELL.re_initialize( 3 ) ;

    CELL_INPUT >> buffer_size_t ;
    if( !FROM_ZERO )
    {
       if( buffer_size_t==0 )
          GE_TriangleMeshing_ERROR::n4() ;
       --buffer_size_t ;
    }
    VERT_OF_CURRENT_CELL( 0 ) = VERTEX_CORRECT_NB( buffer_size_t ) ;

    CELL_INPUT >> buffer_size_t ;
    if( !FROM_ZERO )
    {
       if( buffer_size_t==0 )
          GE_TriangleMeshing_ERROR::n4() ;
       --buffer_size_t ;       
    }
    VERT_OF_CURRENT_CELL( 1 ) = VERTEX_CORRECT_NB( buffer_size_t ) ;

    CELL_INPUT >> buffer_size_t ;
    if( !FROM_ZERO )
    {
       if( buffer_size_t==0 )
          GE_TriangleMeshing_ERROR::n4() ;
       --buffer_size_t ;       
    }
    VERT_OF_CURRENT_CELL( 2 ) = VERTEX_CORRECT_NB( buffer_size_t ) ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: read_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: read_cell" ) ;
   PEL_CHECK( valid_cell() ) ;

   read_node_to_point_coordinate_data() ;

   // Cell to its sides connectivity table
   SIDE_OF_CURRENT_CELL.re_initialize( 3 ) ;
   for( size_t i=0; i<3 ; ++i )
   {
      SIDE_OF_CURRENT_CELL( i ) = SIDE_OF_CELL( I_CELL, i ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: read_side( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: read_side" ) ;
   PEL_CHECK( valid_face() ) ;

   VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_SIDE( I_SIDE, 0 ) ;
   VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_SIDE( I_SIDE, 1 ) ;

   // Side color
   GE_Color const* scolor = GE_Color::null_color() ;
   std::string color = "" ;

   size_t cm = COLOR_OF_VERT_INDEX( VERT_OF_CURRENT_SIDE( 0 ) ) ;
   size_t cM = COLOR_OF_VERT_INDEX( VERT_OF_CURRENT_SIDE( 1 ) ) ;
   if( cm > cM )
   {
      size_t buffer_size_t = cm ; cm = cM ; cM = buffer_size_t ;
   }
   color += "r" ;
   ostringstream oss_color_cm;
   oss_color_cm << cm;
   color += oss_color_cm.str(); 
   
   if( cM > cm )
   {
      color += "r" ;
      ostringstream oss_color_cM;
      oss_color_cM << cM;
      color += oss_color_cM.str();
   }
   GE_Color::extend( color ) ;
   scolor = GE_Color::object( color ) ;

   COLOR_OF_SIDE.resize( COLOR_OF_SIDE.size()+1 ) ;
   COLOR_OF_SIDE[ I_SIDE ] = scolor ;
}

//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: build_sides_connectivity( size_t cell_number )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: build_sides_connectivity" ) ;

   size_t side_number = 0 ;

   // Side 0
   VERT_OF_CURRENT_SIDE.re_initialize( 2 ) ;
   VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 0 ) ;
   VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 1 ) ;
   side_number = extend_list_of_sides() ;
   SIDE_OF_CELL( cell_number, 0 ) = side_number ;

   // Side 1
   VERT_OF_CURRENT_SIDE.re_initialize( 2 ) ;
   VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 1 ) ;
   VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 2 ) ;
   side_number = extend_list_of_sides() ;
   SIDE_OF_CELL( cell_number, 1 ) = side_number ;

   // Side 2
   VERT_OF_CURRENT_SIDE.re_initialize( 2 ) ;
   VERT_OF_CURRENT_SIDE( 0 ) = VERT_OF_CURRENT_CELL( 2 ) ;
   VERT_OF_CURRENT_SIDE( 1 ) = VERT_OF_CURRENT_CELL( 0 ) ;
   side_number = extend_list_of_sides() ;
   SIDE_OF_CELL( cell_number, 2 ) = side_number ;
}

//-----------------------------------------------------------------------------
size_t
GE_TriangleMeshing:: extend_list_of_sides( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: extend_list_of_sides" ) ;

   size_t result = PEL::bad_index() ;
   PEL_IndexSet* indx = PEL_IndexSet::create( 0,
                                              VERT_OF_CURRENT_SIDE,
                                              SIDE_TREE->count() ) ;
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
      VERT_OF_SIDE( result, 0 ) = VERT_OF_CURRENT_SIDE( 0 ) ;
      VERT_OF_SIDE( result, 1 ) = VERT_OF_CURRENT_SIDE( 1 ) ;
      ++NB_SIDES ;
   }

   PEL_CHECK( result<NB_SIDES ) ;

   return( result ) ;
}


//-----------------------------------------------------------------------------
void
GE_TriangleMeshing:: remove_outofgrid_vertices( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: remove_outofgrid_vertices" ) ;

   size_t buffer_size_t ;
   for( size_t i=0; i<NB_CELLS; ++i )
   {
      CELL_INPUT >> buffer_size_t ;
      for( size_t d=0; d<3 ; ++d )
      {
         CELL_INPUT >> buffer_size_t ;
         if( !FROM_ZERO )
         {
            --buffer_size_t ;
         }
         GRID_VERTEX( buffer_size_t ) = true ;
      }
   }

   size_t switc = 0 ;
   NB_VERTS = 0 ;
   for( size_t i=0; i<GRID_VERTEX.size(); ++i )
   {
      if( !GRID_VERTEX( i ) )
      {
         ++switc ;
      }
      else
      {
         ++NB_VERTS ;
         VERTEX_CORRECT_NB( i ) = i - switc ; 
      }
   }

}

//-------------------------------------------------------------------------
void
GE_TriangleMeshing:: check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TriangleMeshing:: check_mesh_polyhedron" ) ;
   PEL_CHECK( check_mesh_polyhedron_PRE(
                            dim_space, face_poly_name, cell_poly_name ) ) ;

   // Check the cell reference polyhedron:
   GE_ReferencePolyhedron const* cell_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   if( cell_poly_ref->nb_vertices() != 3 )
   {
      raise_invalid_cell_polyhedron(
         cell_poly_name, "   inconsistent with triangular meshing" ) ;
   }

   // Check the side reference polyhedron:
   GE_ReferencePolyhedron const* side_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( face_poly_name ) ;
   if( side_poly_ref->nb_vertices() != 2 )
   {
      raise_invalid_face_polyhedron(
         face_poly_name, "   inconsistent with triangular meshing" ) ;
   }
}


//-----------------------------------------------------------------------------
bool
GE_TriangleMeshing:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;

   return true ;
}

//internal-------------------------------------------------------------------
void
GE_TriangleMeshing_ERROR:: n0( string const& fname )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl<< "GE_TriangleMeshing :" << endl ;
   mesg << "  Unable to open Triangle output file : \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_TriangleMeshing_ERROR:: n1( size_t dim, size_t space_dim )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_TriangleMeshing :" << endl ;
   mesg << "   inconsistent number of space dimensions" << endl ;
   mesg << "   in the Triangle data file : " << dim << endl ;
   mesg << "   in the PELICANS data file : " << space_dim ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_TriangleMeshing_ERROR:: n3( size_t nb_pts )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_TriangleMeshing :" << endl ;
   mesg << "   cells can only be triangles with 3 vertices while " << nb_pts
        << "   are declared in .ele file. " ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
GE_TriangleMeshing_ERROR:: n4( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "GE_TriangleMeshing : " << endl ;
   mesg << " Numbering problem while reading vertices " << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

