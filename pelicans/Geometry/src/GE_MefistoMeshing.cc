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

#include <GE_MefistoMeshing.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>

using std::string ;

class ColorGenerator
{
   public:
      
      ColorGenerator( void ) ;
     ~ColorGenerator( void ) ;
      void re_initialize( size_t dimension ) ;
      void add_new_color( std::string const& color, size_t number ) ; 
      GE_Color* build_color( size_t_vector const& my_colors ) const ;
      GE_Color* build_color( size_t my_color ) const ;
      void print( std::ostream& os, size_t indent_width ) const ;
      
   protected:

   private:
      
      size_t_vector CORRESPONDANCES ;
      stringVector COLORS ;
      size_t DIMENSION ;
      size_t CUMULATIVE_INDEX ;
} ;

struct GE_MefistoMeshing_ERROR
{
    static void n0( std::string const& message ) ;
} ;

GE_MefistoMeshing const*
GE_MefistoMeshing::PROTOTYPE = new GE_MefistoMeshing() ;

//-----------------------------------------------------------------------------
GE_MefistoMeshing:: GE_MefistoMeshing( void )
//-----------------------------------------------------------------------------
   : GE_Meshing( "GE_MefistoMeshing" )
   
   , NB_VERTS( PEL::bad_index() )
   , VERT_COORD( 0 )
   , VERT_COLORS( 0 )
   , I_VERT( PEL::bad_index() )
   , VERT_COLOR( 0 )

   , NB_VERTS_PER_CELL( PEL::bad_index() )
   , NB_SIDES_PER_CELL( PEL::bad_index() )
   , NB_LINES_PER_CELL( PEL::bad_index() )
   , NB_FACES_PER_CELL( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , CELL_POLY_NAME()
   , CELL_COLORS( 0 )
   , I_CELL( PEL::bad_index() )
   , CELL_COLOR( 0 )
   , CELLS_TO_SIDES( 0, 0 )
   , CELL_TO_VERTS( 0 )
   , CELL_TO_SIDES( 0 )

   , SIDE_TREE(0)
   , NB_SIDES( PEL::bad_index() )
   , SIDE_POLY_NAME()
   , SIDE_COLORS( 0 )
   , I_SIDE( PEL::bad_index() )
   , SIDE_COLOR( 0 )
   , SIDES_TO_VERTS( 0, 0 )
   , SIDE_TO_VERTS( 0 )
{
}

//-----------------------------------------------------------------------------
GE_MefistoMeshing*
GE_MefistoMeshing:: create_replica( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp,
                                    size_t dim_space ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   GE_MefistoMeshing* result = new GE_MefistoMeshing( a_owner,
                                                      exp, dim_space ) ;
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_MefistoMeshing:: GE_MefistoMeshing( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp, 
                                       size_t dim_space )
//-----------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , INPUT( exp->string_data( "filename" ).c_str(), std::ios_base::binary )
   
   , NB_VERTS( 0 )
   , VERT_COORD( dim_space )
   , VERT_COLORS( 0 )
   , I_VERT( PEL::bad_index() )
   , VERT_COLOR( 0 )
   
   , NB_VERTS_PER_CELL( PEL::bad_index() )
   , NB_SIDES_PER_CELL( PEL::bad_index() )
   , NB_LINES_PER_CELL( PEL::bad_index() )
   , NB_FACES_PER_CELL( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , CELL_POLY_NAME()
   , CELL_COLORS( 0 )
   , I_CELL( PEL::bad_index() )
   , CELL_COLOR( 0 )
   , CELLS_TO_SIDES( 0, 0 )
   , CELL_TO_VERTS( 0 )
   , CELL_TO_SIDES( 0 )

   , SIDE_TREE( PEL_BalancedBinaryTree:: create( this ) )
   , NB_SIDES( PEL::bad_index() )
   , SIDE_POLY_NAME()
   , SIDE_COLORS( 0 )
   , I_SIDE( PEL::bad_index() )
   , SIDE_COLOR( 0 )
   , SIDES_TO_VERTS( 0, 0 )
   , SIDE_TO_VERTS( 0 )     
{
   PEL_LABEL( "GE_MefistoMeshing:: GE_MefistoMeshing" ) ;
   
   initialize_rounding_strategy( exp ) ;

   std::string const& fname = exp->string_data( "filename" ) ;
   if( !INPUT )
   {
      GE_MefistoMeshing_ERROR:: n0( "Unable to open  file \""+fname+"\"" ) ;
   }
   read_mesh_polyhedron( exp, dim_space, SIDE_POLY_NAME, CELL_POLY_NAME ) ;

// First line (object_name + comment ...)
   string buffer_string ;
   getline( INPUT, buffer_string ) ;

// Second line ( space dimension + comment ...)
   size_t dimension = 0 ;
   INPUT >> dimension ;
   getline( INPUT, buffer_string ) ;
   if( dimension != dim_space )
   {
      GE_MefistoMeshing_ERROR:: n0( "invalid space dimension" ) ;
   }

// Third line ( nb of nodes + comment ...)
   INPUT >> NB_VERTS ;
   getline( INPUT, buffer_string ) ;

// Fourth line ( nb of tangent + comment ...)
   size_t buffer_size_t = 0 ;
   INPUT >> buffer_size_t ;
   getline( INPUT, buffer_string ) ;
   if( buffer_size_t != 0 )
   {
      GE_MefistoMeshing_ERROR:: n0( "No tangent should be defined" ) ;
   }

// NB_VERTS lines ( nodes coordinates )
   VERT_POS = INPUT.tellg() ;
   for( size_t iLine=0 ; iLine<NB_VERTS ; iLine++ )
   {
       getline( INPUT, buffer_string ) ;
   }
   size_t_vector vertices_color_order( NB_VERTS ) ;

// one more line + 2*nb_materials lines ( materials or sub-domains )
   size_t nb_materials = 0 ;
   ColorGenerator materials ;
   read_one_geometric_element_type( materials, nb_materials ) ;

// one more line ( number boundary objects )
   size_t nb_boundary_objects = 0 ;
   INPUT >> nb_boundary_objects ;
   getline( INPUT, buffer_string ) ;

// CASE 3D : one more line + nb_surfaces*2 lines (boundary surfaces)
   size_t nb_surfaces = 0 ;
   ColorGenerator surfaces ;
   if( dim_space == 3 )
   {
      read_one_geometric_element_type( surfaces, nb_surfaces ) ;
   }

// one more line + nb_lines*2 lines (boundary lines)
   size_t nb_lines = 0 ;
   ColorGenerator lines ;
   read_one_geometric_element_type( lines, nb_lines ) ;

// one more line + nb_points*2 lines (boundary points)
   size_t nb_points = 0 ;
   ColorGenerator points ;
   read_one_geometric_element_type( points, nb_points ) ;

// verification on the number of boundary objects
   PEL_CHECK( nb_boundary_objects == nb_surfaces + nb_lines + nb_points ) ;
   size_t max_nb_side_colors = PEL::max( nb_materials, 
                                         PEL::max( nb_surfaces, nb_lines ) ) ;
   size_t max_nb_vert_colors = PEL::max( nb_points, max_nb_side_colors ) ;

// seven more lines
   size_t nb_type_polyhedrons = 0 ;
   INPUT >> nb_type_polyhedrons ; 
   getline( INPUT, buffer_string ) ;
   
   if( nb_type_polyhedrons != 1 )
   {
      GE_MefistoMeshing_ERROR:: n0(
         "No more than one different type of polyhedron can now be accepted" ) ;
   }
   
   INPUT >> NB_CELLS ;
   getline( INPUT, buffer_string ) ;
   CELL_COLORS = PEL_Vector::create( this, NB_CELLS ) ;

   INPUT >> buffer_size_t ;
   getline( INPUT, buffer_string ) ;

   INPUT >> NB_VERTS_PER_CELL ;
   getline( INPUT, buffer_string ) ;
   if( buffer_size_t != NB_VERTS_PER_CELL )
   {
      GE_MefistoMeshing_ERROR:: n0(
         "No point except to polyhedron vertices can be defined" ) ;
   }
   if( NB_VERTS_PER_CELL == 6 )
   {
      GE_MefistoMeshing_ERROR:: n0( "No pentahedra can be defined" ) ;       
   }
   
   INPUT >> NB_LINES_PER_CELL ;
   getline( INPUT, buffer_string ) ;
   
   INPUT >> NB_FACES_PER_CELL ;
   getline( INPUT, buffer_string ) ;
   
   NB_SIDES_PER_CELL = NB_FACES_PER_CELL ;
   if( dim_space==2 )
   {
      NB_SIDES_PER_CELL = NB_LINES_PER_CELL ;
   }
   CELLS_TO_SIDES.re_initialize( NB_CELLS, NB_SIDES_PER_CELL ) ;

   GE_ReferencePolyhedron const* cell_polyhedron
                     = GE_Mpolyhedron::reference_polyhedron( CELL_POLY_NAME ) ;
   if( ( NB_VERTS_PER_CELL!=cell_polyhedron->nb_vertices() ) ||
       ( NB_SIDES_PER_CELL!=cell_polyhedron->nb_faces() ) )
   {
      GE_MefistoMeshing_ERROR:: n0( "Invalid definition of polyhedra" ) ;
   }

   size_t nb_volumes_this_polyhedron = 0 ;
   INPUT >> nb_volumes_this_polyhedron ;
   getline( INPUT, buffer_string ) ;
   if( nb_volumes_this_polyhedron != ( dim_space - 2 ) )
   {
      GE_MefistoMeshing_ERROR:: n0( "Invalid number of volumes" ) ;
   }
   
// Last lines minus one ( cell connectivities )
   GE_ReferencePolyhedron const* side_polyhedron
                     = GE_Mpolyhedron::reference_polyhedron( SIDE_POLY_NAME ) ;
   SIDE_TO_VERTS.re_initialize( side_polyhedron->nb_vertices() ) ;
   SIDES_TO_VERTS.re_initialize( NB_SIDES_PER_CELL*NB_CELLS,
                                     side_polyhedron->nb_vertices() ) ;
   size_t_array2D vertices_colors( NB_VERTS, max_nb_vert_colors ) ;
   size_t_array2D sides_colors( NB_SIDES_PER_CELL*NB_CELLS,
                                max_nb_side_colors ) ;
   size_t_vector sides_color_order( NB_SIDES_PER_CELL*NB_CELLS ) ;
   NB_SIDES = 0 ;

   CELL_POS = INPUT.tellg() ;
   for( size_t iCell=0 ; iCell<NB_CELLS ; iCell++ )
   {
      CELL_TO_VERTS.re_initialize( NB_VERTS_PER_CELL ) ;
      for( size_t iVertLoc=0; iVertLoc<NB_VERTS_PER_CELL; iVertLoc++ )
      {
         INPUT >> buffer_size_t ;
         CELL_TO_VERTS( iVertLoc ) = --buffer_size_t ;
         if( CELL_TO_VERTS( iVertLoc )>=NB_VERTS )
         {
            GE_MefistoMeshing_ERROR:: n0( "Invalid vertex number for cell" ) ;
         }
      }

      // Reference of points
      for( size_t iVertLoc=0; iVertLoc<NB_VERTS_PER_CELL; iVertLoc++ )
      {
         INPUT >> buffer_size_t ;
         size_t vert_global_nb = CELL_TO_VERTS( iVertLoc ) ;
         if( buffer_size_t!=0 )
	 {
	    // One reference of point.
            add_color_for_side_or_vertex( vertices_colors, vert_global_nb,
                                          buffer_size_t, 
                                          vertices_color_order( vert_global_nb )<4 ) ;
            vertices_color_order( vert_global_nb ) = (size_t)4 ;

	 }
      }
      // Building of sides.
      sides_connectivity_building( iCell, vertices_colors, 
                                   vertices_color_order, sides_colors,
                                   sides_color_order ) ;

      // Default color for vertices, the one of cell material
      size_t vol_surf_nb = 0 ;
      INPUT >> vol_surf_nb ;
      PEL_ASSERT( vol_surf_nb!=0 ) ;
      CELL_COLORS->set_at( iCell,  materials.build_color( vol_surf_nb ) ) ; 
      if( vol_surf_nb!=0 )
      {
         size_t lnn_order = 5 - dim_space ;
         for( size_t iVertLoc=0; iVertLoc<NB_VERTS_PER_CELL; iVertLoc++ )
         {
            size_t vert_global_nb = CELL_TO_VERTS( iVertLoc ) ;
            if( vertices_color_order( vert_global_nb  )< lnn_order )
	    {
	       // One reference of point.
               add_color_for_side_or_vertex( vertices_colors, vert_global_nb,
                                             vol_surf_nb, false ) ;
	       vertices_color_order( vert_global_nb ) = 1 ;
	    }
	 }
         for( size_t iSide=0; iSide<NB_SIDES_PER_CELL; iSide++ )
	 {
            size_t side_global_nb = CELLS_TO_SIDES( iCell, iSide ) ;
            if( sides_color_order( side_global_nb ) < lnn_order )
	    {
               add_color_for_side_or_vertex( sides_colors, side_global_nb,
                                             vol_surf_nb, false ) ;
               sides_color_order( side_global_nb ) = 1 ; 
	    }
	 }
      }
   }
   
   SIDE_COLORS = PEL_Vector::create( this, NB_SIDES ) ;
   VERT_COLORS = PEL_Vector::create( this, NB_VERTS ) ;

   colorize_vertices_and_sides( points, lines, surfaces, materials, 
                                vertices_colors, vertices_color_order, 
                                sides_colors, sides_color_order ) ;
}

//-----------------------------------------------------------------------------
GE_MefistoMeshing:: ~GE_MefistoMeshing( void )
//-----------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//-----------------------------------------------------------------------------
size_t
GE_MefistoMeshing:: nb_vertices( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_VERTS ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_MefistoMeshing:: nb_cells( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_CELLS ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_MefistoMeshing:: nb_faces( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_SIDES )  ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: start_vertex_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   INPUT.seekg( VERT_POS ) ;
   I_VERT = 0 ;
   read_vertex() ;
}

//-----------------------------------------------------------------------------
bool
GE_MefistoMeshing:: valid_vertex( void ) const
//-----------------------------------------------------------------------------
{
   return( I_VERT<NB_VERTS ) ;   
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: go_next_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   I_VERT++ ;
   if( valid_vertex() ) read_vertex() ;
}

//-----------------------------------------------------------------------------
doubleVector const& 
GE_MefistoMeshing:: vertex_coordinates( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
   doubleVector const& result = VERT_COORD ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: start_cell_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   I_CELL = 0 ;
   INPUT.seekg( CELL_POS ) ;
   read_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_MefistoMeshing:: valid_cell( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: valid_cell" ) ;

   bool result = ( I_CELL< NB_CELLS )  ;
   PEL_CHECK_POST( valid_cell_POST( result ) ) ;

   return result ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: go_next_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   I_CELL++ ;
   if( valid_cell() ) read_cell() ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_MefistoMeshing:: cell_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_MefistoMeshing:: cell_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = CELL_TO_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_MefistoMeshing:: cell_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = CELL_TO_SIDES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: start_face_iterator( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   I_SIDE = 0 ;
   read_side() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_MefistoMeshing:: valid_face( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: valid_face" ) ;

   bool result = ( I_SIDE<NB_SIDES ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: go_next_face( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   I_SIDE++ ;
   if( valid_face() )
   {
      read_side() ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_MefistoMeshing:: face_polyhedron_name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = SIDE_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
GE_MefistoMeshing:: face_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = SIDE_TO_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_MefistoMeshing:: default_vertex_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = VERT_COLOR ;

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_MefistoMeshing:: default_cell_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = CELL_COLOR ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_MefistoMeshing:: default_face_color( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = SIDE_COLOR ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: read_vertex( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: read_vertex" ) ;
   PEL_CHECK( valid_vertex() ) ;

   double coordinate = PEL::bad_double() ;
   INPUT >> coordinate ;
   VERT_COORD( 0 ) = roundoff( coordinate ) ;
   INPUT >> coordinate ;
   VERT_COORD( 1 ) = roundoff( coordinate ) ;
   INPUT >> coordinate ;
   if( nb_space_dimensions()==3 )
      VERT_COORD( 2 ) = roundoff( coordinate ) ;

   VERT_COLOR = static_cast<GE_Color const*>(VERT_COLORS->at( I_VERT )) ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: read_cell( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: read_cell" ) ;
   PEL_CHECK( valid_cell() ) ;

   // Mesh vertices
   CELL_TO_VERTS.re_initialize( NB_VERTS_PER_CELL ) ;
   for( size_t i=0 ; i<NB_VERTS_PER_CELL ; i++ )
   {
      INPUT >> CELL_TO_VERTS(i) ;
      CELL_TO_VERTS(i)-- ;
   }

   // Cell color the one of the surface in 2D, of the volume in 3D.
   size_t buffer_size_t ;
   for( size_t i=0; i<NB_VERTS_PER_CELL; i++ )
   {
      INPUT >> buffer_size_t ;
   }
   for( size_t i=0; i<NB_LINES_PER_CELL; i++ )
   {
      INPUT >> buffer_size_t ;
   }
   if( nb_space_dimensions()==3 )
   {
      for( size_t i=0; i<NB_FACES_PER_CELL; i++ )
      {
	 INPUT >> buffer_size_t ;
      }
   }
   INPUT >> buffer_size_t ;
   
   // Cell to its sides connectivity table
   CELL_TO_SIDES.re_initialize( NB_SIDES_PER_CELL ) ;

   if( CELL_POLY_NAME=="GE_Tetrahedron" ||
       CELL_POLY_NAME=="GE_Triangle"  || CELL_POLY_NAME=="GE_Quadrilateral" || 
       CELL_POLY_NAME=="GE_Rectangle" || CELL_POLY_NAME=="GE_Trapezoid" )
   {
      for( size_t i=0; i<NB_SIDES_PER_CELL ; i++ )
      {
         CELL_TO_SIDES( i ) = CELLS_TO_SIDES( I_CELL, i ) ;
      }
   }
   else if( CELL_POLY_NAME=="GE_Hexahedron" || CELL_POLY_NAME=="GE_Cuboid" )
   {
      CELL_TO_SIDES( 0 ) = CELLS_TO_SIDES( I_CELL, 2 ) ;
      CELL_TO_SIDES( 1 ) = CELLS_TO_SIDES( I_CELL, 4 ) ;
      CELL_TO_SIDES( 2 ) = CELLS_TO_SIDES( I_CELL, 5 ) ;
      CELL_TO_SIDES( 3 ) = CELLS_TO_SIDES( I_CELL, 1 ) ;
      CELL_TO_SIDES( 4 ) = CELLS_TO_SIDES( I_CELL, 0 ) ;
      CELL_TO_SIDES( 5 ) = CELLS_TO_SIDES( I_CELL, 3 ) ;
   }

   CELL_COLOR =  static_cast<GE_Color const*>(CELL_COLORS->at( I_CELL )) ;
   
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: read_side( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: read_side" ) ;
   PEL_CHECK( valid_face() ) ;

   for( size_t i=0; i<SIDE_TO_VERTS.size() ; i++ )
   {
      SIDE_TO_VERTS( i ) = SIDES_TO_VERTS( I_SIDE, i ) ;
   }
   SIDE_COLOR =  static_cast<GE_Color const*>(SIDE_COLORS->at( I_SIDE )) ;
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: read_one_geometric_element_type( 
                                                ColorGenerator& gene,
       		                                size_t& nb_geometric_elements )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: read_one_geometric_element_type" ) ;
   
   // First gets the number of boundary elements of that type
   INPUT >> nb_geometric_elements ;
   static string buffer_string ;
   getline( INPUT, buffer_string ) ;

   // Vector of colors
   gene.re_initialize( nb_geometric_elements ) ;
 
   // Loop on boundary elements of that type
   size_t current_number ;
   string current_name ;
   for( size_t i=0 ; i<nb_geometric_elements ; i++ )
   {
       // Number of current color
       INPUT >> current_number ;
       getline( INPUT, buffer_string ) ;

       // Corresponding color
       INPUT >> current_name ;
       getline( INPUT, buffer_string ) ;
       gene.add_new_color( current_name, current_number ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: sides_connectivity_building( 
                                          size_t cell_number,
                                          size_t_array2D& vertices_colors, 
                                          size_t_vector& vertices_color_order, 
                                          size_t_array2D&  sides_colors,
                                          size_t_vector& sides_color_order )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: sides_connectivity_building" ) ;
   
   size_t side_number = 0 ;

   if( CELL_POLY_NAME=="GE_Triangle"  || CELL_POLY_NAME=="GE_Quadrilateral" || 
       CELL_POLY_NAME=="GE_Rectangle" || CELL_POLY_NAME=="GE_Trapezoid" )
   {
      PEL_CHECK( SIDE_POLY_NAME=="GE_Segment" ) ;
      for( size_t i=0; i<NB_VERTS_PER_CELL-1; i++ )
      {
         SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( i ) ;
         SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( i+1 ) ;
         side_number = extend_list_of_sides( vertices_colors, 
					     vertices_color_order, 
					     sides_colors, 
					     sides_color_order ) ;
         CELLS_TO_SIDES( cell_number, i ) = side_number ;
      }
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( NB_VERTS_PER_CELL-1 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 0 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, NB_VERTS_PER_CELL-1 ) 
	 = side_number ;
   }
   else if( CELL_POLY_NAME=="GE_Tetrahedron" )
   {
      size_t buffer_size_t = 0 ;
      for( size_t i=0; i<NB_LINES_PER_CELL; i++ )
      {
         INPUT >> buffer_size_t ;
      }

      PEL_CHECK( SIDE_POLY_NAME=="GE_Triangle" ) ;
      // First side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 2 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 1 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 0 ) = side_number ;
      // Second side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 3 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 2 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 1 ) = side_number ;
      // Third side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 1 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 3 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 2 ) = side_number ;
      // Fourth side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 1 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 2 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 3 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, (size_t)3 ) = side_number ;
   }
   else if( CELL_POLY_NAME=="GE_Hexahedron" || CELL_POLY_NAME=="GE_Cuboid" )
   {
      PEL_CHECK( SIDE_POLY_NAME=="GE_Quadrilateral" || 
                 SIDE_POLY_NAME=="GE_Trapezoid" ||
                 SIDE_POLY_NAME=="GE_Rectangle" ) ;
      size_t buffer_size_t = 0 ;
      for( size_t i=0; i<NB_LINES_PER_CELL; i++ )
      {
         INPUT >> buffer_size_t ;
      }
      // First side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 3 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 2 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 1 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 0 ) = side_number ;
      // Second side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 4 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 7 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 3 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 1 ) = side_number ;
      // Third side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 0 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 1 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 5 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 4 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 2 ) = side_number ;
      // Fourth side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 4 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 5 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 6 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 7 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 3 ) = side_number ;
      // Fifth side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 1 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 2 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 6 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 5 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order ) ;
      CELLS_TO_SIDES( cell_number, 4 ) = side_number ;
      // Sixth side
      SIDE_TO_VERTS( 0 ) = CELL_TO_VERTS( 3 ) ;
      SIDE_TO_VERTS( 1 ) = CELL_TO_VERTS( 7 ) ;
      SIDE_TO_VERTS( 2 ) = CELL_TO_VERTS( 6 ) ;
      SIDE_TO_VERTS( 3 ) = CELL_TO_VERTS( 2 ) ;
      side_number = extend_list_of_sides( vertices_colors, 
                                          vertices_color_order, 
                                          sides_colors, 
                                          sides_color_order) ;
      CELLS_TO_SIDES( cell_number, 5 ) = side_number ;
   }
}

//-----------------------------------------------------------------------------
size_t
GE_MefistoMeshing:: extend_list_of_sides( size_t_array2D& vertices_colors, 
					  size_t_vector& vertices_color_order, 
					  size_t_array2D&  sides_colors,
					  size_t_vector& sides_color_order )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: extend_list_of_sides" ) ;

   size_t result = PEL::bad_index() ;
   
   // Search if side has not already been registred
   PEL_IndexSet* indx = PEL_IndexSet::create( 0,
                                              SIDE_TO_VERTS,
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
   }

   size_t current_side_color = 0 ;
   INPUT >> current_side_color ;
   
   size_t lnn_order = 2 ; // case 3D
   if( SIDE_TO_VERTS.size()==2 ) lnn_order = 3 ; // case 2D

   if( !found )
   {
      PEL_ASSERT( result == NB_SIDES ) ;
      for( size_t d=0 ; d<SIDE_TO_VERTS.size() ; ++d )
      {
         SIDES_TO_VERTS( NB_SIDES, d ) = SIDE_TO_VERTS( d ) ;

         if( current_side_color!=0 && 
             vertices_color_order( SIDE_TO_VERTS( d ) )<=lnn_order )
         {
            bool clean_prec
             = vertices_color_order( SIDE_TO_VERTS( d ) ) <lnn_order ;
            add_color_for_side_or_vertex( vertices_colors, 
                                          SIDE_TO_VERTS( d ),
                                          current_side_color, clean_prec ) ;

            vertices_color_order( SIDE_TO_VERTS( d ) ) = lnn_order ;
	 }
      }
      NB_SIDES++ ;
   }
   
   if( current_side_color!=0 )
   {
      add_color_for_side_or_vertex( sides_colors, result, 
				    current_side_color, 
                                    sides_color_order( result )<lnn_order ) ;
      sides_color_order( result ) = lnn_order ;
   }
   
   PEL_CHECK( result<NB_SIDES ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void 
GE_MefistoMeshing:: add_color_for_side_or_vertex( size_t_array2D& array, 
                                                  size_t posi,
                                                  size_t color,
                                                  bool clean_preceeding ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: add_color_for_side_or_vertex" ) ;
   
   size_t ind_bd = array.index_bound( 1 ) ;

   bool cont = true ;
   size_t index = 0 ;

   // One color of a geometrical element of higher level must be applied
   // clean preceeding colors.
   if( clean_preceeding )
   {
      for( size_t i=0; i<ind_bd && ( array( posi, i )!=0 );  i++ )
      {
	 array( posi, i ) = 0 ;
      }
   }
   else
   {
      for( size_t i=0; i<ind_bd && ( array( posi, i )!=0 ) && cont ; i++ )
      {
	 cont = ( array( posi, i )!=color ) ;
	 index++ ;
      }
   }

   if( cont ) // The color has not been found and should be addeed
   {
      PEL_CHECK( index<ind_bd ) ;
      array( posi, index ) = color ;
   }
   // Search for the first null color and test if color has not yet been inserted. 
}

//-----------------------------------------------------------------------------
void
GE_MefistoMeshing:: colorize_vertices_and_sides(
                            ColorGenerator Points,
                            ColorGenerator Lines,
                            ColorGenerator Surfaces,
                            ColorGenerator Materials,
                            size_t_array2D& vertices_colors,
                            size_t_vector& vertices_color_order,
                            size_t_array2D& sides_colors,
                            size_t_vector& sides_color_order )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_MefistoMeshing:: colorize_vertices_and_sides" ) ;
   
   static size_t_vector compose_color( 0 ) ;
   static string color_name ;
   
   compose_color.re_initialize( vertices_colors.index_bound( 1 ) ) ;
   for( size_t iVert=0 ; iVert<NB_VERTS ; ++iVert )
   {
      vertices_colors.extract_section( 0, iVert, compose_color ) ;

      if( vertices_color_order( iVert )==4 )
      {
	 VERT_COLORS->set_at( iVert, Points.build_color( compose_color ) ) ;
      }
      else if( vertices_color_order( iVert )==3 )
      {
	 VERT_COLORS->set_at( iVert, Lines.build_color( compose_color ) ) ;
      }
      else if( vertices_color_order( iVert )==2 )
      {
	 VERT_COLORS->set_at( iVert, Surfaces.build_color( compose_color ) ) ;
      }
      else if( vertices_color_order( iVert )==1 )
      {
	 VERT_COLORS->set_at( iVert, Materials.build_color( compose_color ) ) ;
      }
      else if( vertices_color_order( iVert )==0 )
	 PEL_Error::object()->raise_plain( "No color defined for one vertex") ;
   }

   compose_color.re_initialize( sides_colors.index_bound( (size_t)1 ) ) ;
   for( size_t iSide=0; iSide<NB_SIDES; iSide++ )
   {
      sides_colors.extract_section( 0, iSide, compose_color ) ;
      if( sides_color_order( iSide )==3 )
      {
	 SIDE_COLORS->set_at( iSide, Lines.build_color( compose_color ) ) ;
      }
      else if( sides_color_order( iSide )==2 )
      {
	 SIDE_COLORS->set_at( iSide, Surfaces.build_color( compose_color ) ) ;
      }
      else if( sides_color_order( iSide )==1 )
      {
	 SIDE_COLORS->set_at( iSide, Materials.build_color( compose_color ) ) ;
      }
      else if( sides_color_order( iSide )==0 )
	 PEL_Error::object()->raise_plain( "No color defined for one side" ) ;
   }
}

//-----------------------------------------------------------------------------
bool
GE_MefistoMeshing:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;
   return( true ) ;
}

//Internal--------------------------------------------------------------------
ColorGenerator:: ColorGenerator( void )
//Internal--------------------------------------------------------------------
   : CORRESPONDANCES( (size_t)0 )
   , COLORS( (size_t)0 )
   , DIMENSION( (size_t)0 )
   , CUMULATIVE_INDEX( (size_t)0 )
{
}

//Internal--------------------------------------------------------------------
ColorGenerator:: ~ColorGenerator( void )
//Internal--------------------------------------------------------------------
{
}

//Internal--------------------------------------------------------------------
void ColorGenerator:: re_initialize( size_t dimension )
//Internal--------------------------------------------------------------------
{
   DIMENSION = dimension ;
   CORRESPONDANCES.re_initialize( dimension ) ;
   COLORS.re_initialize( dimension ) ;
   CUMULATIVE_INDEX = (size_t)0 ;
}

//Internal--------------------------------------------------------------------
void 
ColorGenerator:: add_new_color( std::string const& color, size_t number )
//Internal--------------------------------------------------------------------
{
   PEL_CHECK( !color.empty() ) ;
   PEL_CHECK( CUMULATIVE_INDEX<DIMENSION ) ;

   size_t index = CUMULATIVE_INDEX ;
   for( size_t i=0; i<CUMULATIVE_INDEX; i++ )
   {
      if( COLORS(i)>color )
      {
         index = i ;
         break ;
      }
   }

   PEL_CHECK( index<DIMENSION ) ;

   if( color!=COLORS( index ) )
   {
      for( size_t i=CUMULATIVE_INDEX; i>index; i-- )
      {
	 COLORS( i ) = COLORS( i-1 ) ;
	 CORRESPONDANCES( i ) = CORRESPONDANCES( i-1 ) ; 
      }

      COLORS( index ) = color ;
      CORRESPONDANCES( index ) = number ;
      CUMULATIVE_INDEX++ ;
   }
// Add a new color sorting all colors by alphabetic order.

   PEL_CHECK( COLORS.has( color ) ) ;
   PEL_CHECK( CORRESPONDANCES.has( number ) ) ;
}

//Internal--------------------------------------------------------------------
GE_Color*
ColorGenerator:: build_color( size_t_vector const& my_colors ) const
//Internal--------------------------------------------------------------------
{
   PEL_CHECK(
      FORALL( ( size_t i=0 ; i<my_colors.size() && my_colors(i)!=0 ; i++ ),
              CORRESPONDANCES.has( my_colors(i) ) ) ) ;
   
   std::string namec ;
   bool first_name = true ;
   if( DIMENSION>0 )
   {
      for( size_t i=0; i<CUMULATIVE_INDEX; i++ )
      {
         for( size_t j=0 ; j<my_colors.size(); j++ )
	 {
	    if( my_colors( j )==CORRESPONDANCES( i ) )
	    {
               std::string const& a_col = COLORS( i ) ;
               if( !first_name ) namec.append( "_" ) ;
               namec.append( a_col ) ;
               
               first_name = false ;
               break ;
	    }
	 }
      }
   }
   GE_Color::extend( namec ) ;
   GE_Color const* result = GE_Color::object( namec ) ;
   
   return( const_cast<GE_Color*>( result ) ) ;
}

//Internal--------------------------------------------------------------------
GE_Color* 
ColorGenerator:: build_color( size_t my_color ) const
//Internal--------------------------------------------------------------------
{
   PEL_CHECK( CORRESPONDANCES.has( my_color ) ) ;
   size_t idx = CORRESPONDANCES.index_of( my_color ) ;
   
   
   std::string namec = COLORS( idx ) ;
   GE_Color::extend( namec ) ;
   return( const_cast<GE_Color*>( GE_Color::object( namec ) ) ) ;
}

//Internal--------------------------------------------------------------------
void 
ColorGenerator:: print( std::ostream& os, size_t indent_width ) const
//Internal--------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << " colors      |     number" << std::endl ;
   for( size_t i=0; i<CUMULATIVE_INDEX; i++ )
   {
      os << space << COLORS( i ) << " :" << CORRESPONDANCES( i ) << std::endl ;
   }
}

//Internal--------------------------------------------------------------------
void 
GE_MefistoMeshing_ERROR:: n0( std::string const& message )
//Internal--------------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** GE_MefistoMeshing error\n"
         "    Error reading the input file:\n"
         "       "+message ) ;
}
