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

#ifndef GE_MEFISTO_MESHING_HH
#define GE_MEFISTO_MESHING_HH

#include <GE_Meshing.hh>

#include <PEL_Vector.hh>

#include <fstream>
#include <doubleVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <stringVector.hh>
#include <string>

/*
Meshings stored in an Mefisto output file under the xyznpef format.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

   keyword                         data
   -------                         ----
filename                pathname of the Mefisto output file
mesh_polyhedron         polyhedron name for sides and cells

One can found in the "MEFISTO-MAILLAGES, manuel de l'utilisateur"
documentation all necessary informations on the xyznpef format (pages
350-354 of the version of 03/07/2001).

To save a meshing, one has to use the "91;" option in the principal
menu of MEFISTO-MAILLAGES, gives the object name and then select "0"
to process without saving of any tangent.

Be careful, no pentahedra can be now defined in PELICANS. 

No 1D meshing can be obtained with Mefisto but `GE_LineWithSegments::'
provides all necessary functionalities.

One can give references to vertices, sides ( lines in 2D and facets in
3D ) and cells by associating to the saved object points, lines, surfaces
and volumes. A vertex take the reference of a points if they lie to the
same physical point, if not it takes the refence of the lines in 2D and
the surfaces in 3D it belongs to and as a default it takes the reference
of the surfaces in 2D and the volumes in 3D ( Mefisto says materials) it 
belongs to. A side takes the reference of the lines in 2D and the surface 
in 3D it belongs to, and as a default it takes the reference
of the surfaces in 2D and the volumes in 3D ( Mefisto says materials) it 
belongs to. A cell takes the reference of the surfaces in 2D and the 
volumes in 3D ( Mefisto says materials) it belongs to. If a vertex or a
side has to take the reference of several elements (boundary elements or
materials), the reference is compound of all reference names listed
alphabetically and separed by an underscore. As an exemple if a vertex
of the 3D space lies on surfaces EAST and NORTH, it will be associated
to reference EAST_NORTH, if a sides lies on volumes (materials) 
TOP and BOTTOM it will takes the reference BOTTOM_TOP.
*/

class PEL_BalancedBinaryTree ;
class PEL_ModuleExplorer ;
class ColorGenerator ;

using std::string ;

class PEL_EXPORT GE_MefistoMeshing : public GE_Meshing
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
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_MefistoMeshing( void ) ;
     ~GE_MefistoMeshing( void ) ;
      GE_MefistoMeshing( GE_MefistoMeshing const& other ) ;
      GE_MefistoMeshing const& operator=( GE_MefistoMeshing const& other ) ;

      GE_MefistoMeshing( PEL_Object* a_owner,
                         PEL_ModuleExplorer const* exp,
                         size_t dim_space ) ;

      virtual GE_MefistoMeshing* create_replica( PEL_Object* a_owner,
                                                 PEL_ModuleExplorer const* exp,
                                                 size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void read_vertex( void ) ;

      void read_cell( void ) ;

      void read_side( void ) ;

      void read_one_geometric_element_type( ColorGenerator& gene,
				            size_t& nb_geometric_elements) ;

      void sides_connectivity_building( size_t cell_number,
                                        size_t_array2D& vertices_colors, 
                                        size_t_vector& vertices_color_order, 
                                        size_t_array2D&  sides_colors,
                                        size_t_vector& sides_color_order ) ;
      
      size_t extend_list_of_sides( size_t_array2D& vertices_colors, 
                                   size_t_vector& vertices_color_order, 
                                   size_t_array2D&  sides_colors,
                                   size_t_vector& sides_color_order ) ;

      void add_color_for_side_or_vertex( size_t_array2D& array, size_t posi,
                                         size_t color, 
                                         bool clean_preceeding ) const ;

      void colorize_vertices_and_sides( ColorGenerator Points,
                                        ColorGenerator Lines,
                                        ColorGenerator Surfaces,
                                        ColorGenerator Materials,
                                        size_t_array2D& vertices_colors,
                                        size_t_vector& vertices_color_order,
                                        size_t_array2D&  sides_colors,
                                        size_t_vector& sides_color_order ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_MefistoMeshing const* PROTOTYPE ;

   //-- Attributes

      // General entry
      std::ifstream INPUT ;

      // Vertices ...
      size_t NB_VERTS ;
      std::streampos VERT_POS ;
      doubleVector VERT_COORD ;
      PEL_Vector* VERT_COLORS ;
      size_t I_VERT ;
      GE_Color const* VERT_COLOR ;

      // Cells ...
      size_t NB_VERTS_PER_CELL ;
      size_t NB_SIDES_PER_CELL ;
      size_t NB_LINES_PER_CELL ;
      size_t NB_FACES_PER_CELL ;
      size_t NB_CELLS ;
      std::streampos CELL_POS ;
      std::string CELL_POLY_NAME ;
      PEL_Vector* CELL_COLORS ;
      size_t I_CELL ;
      GE_Color const* CELL_COLOR ;
      size_t_array2D CELLS_TO_SIDES ;
      mutable size_t_vector CELL_TO_VERTS ;
      mutable size_t_vector CELL_TO_SIDES ;

      // Sides ...
      PEL_BalancedBinaryTree* SIDE_TREE ;
      size_t NB_SIDES ;
      std::string SIDE_POLY_NAME ;
      PEL_Vector* SIDE_COLORS ;
      size_t I_SIDE ;
      GE_Color const* SIDE_COLOR ;
      size_t_array2D SIDES_TO_VERTS ;
      mutable size_t_vector SIDE_TO_VERTS ;
      
} ;

#endif


