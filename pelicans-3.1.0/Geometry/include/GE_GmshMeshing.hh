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

#ifndef GE_GMSH_MESHING_HH
#define GE_GMSH_MESHING_HH

#include <GE_Meshing.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL_BalancedBinaryTree.hh>
#include <PEL_IndexSet.hh>

#include <doubleArray2D.hh>
#include <stringVector.hh>
#include <size_t_array2D.hh>

#include <fstream>
#include <string>

/*
Meshings stored in a Gmsh output file with the .msh 
format, in version 1.0 or version 2.0 (this latter is
highly recommended)

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

   keyword                         data
   -------                         ----
filename                pathname of the Gmsh output file
mesh_polyhedron         polyhedron name for side and cell
format                  either "1.0" or "2.0"

Be careful, no prism and no pyramid can now be defined in 
PELICANS. 

Gmsh offers the capabilities to save boundary elements
and volumic elements associating them with a reference
number, named "reg-phys". This number corresponds to the
number of the physical entity the boundary or volumic
element belongs. This number is used to build colors
by adding the character "r" before, a.e. "r1", "r1r2"...
To 0 corresponds `GE_Color::null_color'.

Each cell is colored by the corresponding volumic element
color.

Each face (lines in 2D and facets in 3D) is colored 
by the color of the equivalent boundary element if any,
by `GE_Color::null_color' otherwise.

Each vertex is colored by the boundary elements it belongs
to, if any (with a composed color if it belongs to more than one
boundary element of different colors, a.e. "r1r2" if it is shared
by boundary elements of color "r1" and "r2" respectively). If
not, it is colored by the volumic elements it belongs to.
*/

class intVector ;
class PEL_ModuleExplorer ;
class PEL_KeywordDataIterator ;
class PEL_Vector ;
class GE_BoxWithBoxes ;

class PEL_EXPORT GE_GmshMeshing : public GE_Meshing
{
   public: //------------------------------------------------------------

   //-- Measurement

      virtual size_t nb_vertices( void ) const ;

      virtual size_t nb_cells( void ) const ;

      virtual size_t nb_faces( void ) const ;
      
   //-- Vertex-iterator movement

      virtual void start_vertex_iterator( void ) ;

      virtual bool valid_vertex( void ) const ;

      virtual void go_next_vertex( void ) ;

      virtual doubleVector const& vertex_coordinates( void ) const ;
      
   //-- Cell-iterator movement

      virtual void start_cell_iterator( void ) ;

      virtual bool valid_cell( void ) const ;

      virtual void go_next_cell( void ) ;

      virtual std::string const& cell_polyhedron_name( void ) const ;

      virtual size_t_vector const& cell_vertices( void ) const ;

      virtual size_t_vector const& cell_faces( void ) const ;

   //-- Face-iterator movement

      virtual void start_face_iterator( void )  ;

      virtual bool valid_face( void ) const ;

      virtual void go_next_face( void ) ;

      virtual std::string const& face_polyhedron_name( void ) const ;

      virtual size_t_vector const& face_vertices( void ) const ;
    
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_GmshMeshing( void ) ;
     ~GE_GmshMeshing( void ) ;
      GE_GmshMeshing( GE_GmshMeshing const& other ) ;
      GE_GmshMeshing const& operator=( GE_GmshMeshing const& other ) ;

      GE_GmshMeshing( PEL_Object* a_owner,
                      PEL_ModuleExplorer const* exp,
		      size_t dim_space ) ;

      virtual GE_GmshMeshing* create_replica( PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp,
                                              size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void read_vertex( void ) ;

      void read_cell( void ) ;

      void read_face( void ) ;
			 
      void go_i_th_next_line( size_t i ) ;

      void determine_mesh_types( void ) ;

      void initialize_tables( size_t nb_tot_elems ) ;

      void read_one_bound( size_t bound_number, 
                           size_t_array2D& vertices_colors ) ;
 
      void preread_one_cell( size_t cell_number, 
                             size_t_array2D& vertices_colors ) ;
      void read_one_edge( void ) ;

      void read_one_point( void ) ;

      void faces_connectivity_building( size_t cell_number ) ;

      size_t extend_list_of_faces( void ) ;

      void attribute_colors_to_vertices( size_t_array2D const& vertices_colors ) ;

      GE_Color* provide_color( size_t col ) const ;
      GE_Color* provide_color( size_t_vector const& cols ) const ;

      void read_input_file_header( PEL_ModuleExplorer const* exp ) ;

      void read_mesh_colors_and_tags( size_t nb_vert,  size_t& c_col ) ;

  //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_GmshMeshing const* prototype ;

   //-- Attributes

      // Input file
      enum FILE_FORMAT_TYPE{ undefined, one, two } ;
      std::string const FILENAME ;
      std::ifstream INPUT ;
      FILE_FORMAT_TYPE INPUT_FILE_FORMAT ;
      stringVector N_C_KEYWORDS ;

      // Generalities
      size_t NBDIMS ;

      // Vertices
      size_t NB_VERTS ;
      size_t I_VERT ;
      std::streampos VERT_POS ;
      doubleVector VERT_COORD ;
      PEL_Vector* VERTS_COLORS ;
      size_t_vector VERT_NUMB_COR ;

      // Elements
      std::streampos ELEM_POS ;

      // Faces
      size_t NB_FACES ;
      std::string FACE_POLY_NAME ;
      size_t FACE_TYPE ; // Gmsh dedicated index for faces
      size_t I_FACE ;
      mutable size_t_vector VERTS_OF_CURRENT_FACE ;
      size_t_array2D VERTS_OF_FACES ;
      GE_ReferencePolyhedron const* FACE_POL ;
      GE_ReferencePolyhedron const* EDGE_POL ;
      PEL_BalancedBinaryTree* FACE_TREE ;
      GE_Color const* FACE_COLOR ;
      PEL_Vector* FACES_COLORS ;
      std::string EDGE_POLY_NAME ;
      size_t EDGE_TYPE ; // Gmsh dedicated index for edges
      size_t POINT_TYPE ; // Gmsh dedicated index for points

      // Bounds
      PEL_BalancedBinaryTree* BOUNDS_TREE ;
      PEL_Vector* BOUNDS_COLORS ; 
      size_t I_BOUND ;     

      // Cells
      size_t NB_CELLS ;
      std::string CELL_POLY_NAME ;
      size_t CELL_TYPE ; // Gmsh dedicated index for cells
      size_t I_CELL ;
      GE_Color const* CELL_COLOR ;
      size_t_array2D FACES_OF_CELLS ;
      mutable size_t_vector VERTS_OF_CURRENT_CELL ;
      mutable size_t_vector FACES_OF_CURRENT_CELL ;
      GE_ReferencePolyhedron const* CELL_POL ;

} ;

#endif


