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

#ifndef GE_EXTRUDED_MESHING_HH
#define GE_EXTRUDED_MESHING_HH

#include <GE_Meshing.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <size_t_vector.hh>


class intVector ;
class PEL_BalancedBinaryTree ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class GE_Point ;
class GE_SetOfPoints ;

/*
3D meshings whose geometric domain and cells obtained by
by sweeping a planar meshing along a vertical straight line.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has 
   - a submodule of name "GE_Meshing" describing the associated 2D meshing
   - a submodule of name "vertices_coordinate_2" having the following 
     keyword-data pairs :
                keyword                          data
                -------                          ----
     name of a cell color     vector of the coordinates 2 of the 3D vertices 
     in the 2D meshing        obtained by sweeping the vertices of the
                              2D cells of that color
                ------                           ----
     default                  vector of the coordinates 2 of the 3D vertices 
     in the 2D meshing        obtained by sweeping the vertices of the
                              2D cells whose color has not be explicity quoted
   - the following keyword-data pair :
                keyword                          data
                -------                          ----
     mesh_polyhedron          vector of 2 strings : the side polyhedron name
                              and the cell polyhedron name
*/

class PEL_EXPORT GE_ExtrudedMeshing : public GE_Meshing
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

      virtual GE_ReferencePolyhedron const* cell_reference_polyhedron( 
                                                               void ) const ;

      virtual size_t_vector const& cell_vertices( void ) const ;

      virtual size_t_vector const& cell_faces( void ) const ;

   //-- Face-iterator movement

      virtual void start_face_iterator( void )  ;

      virtual bool valid_face( void ) const ;

      virtual void go_next_face( void ) ;

      virtual std::string const& face_polyhedron_name( void ) const ;

      virtual GE_ReferencePolyhedron const* face_reference_polyhedron( 
                                                               void ) const ;

      virtual size_t_vector const& face_vertices( void ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      virtual void display( std::ostream& os, size_t indent_width ) ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_ExtrudedMeshing( void ) ;
     ~GE_ExtrudedMeshing( void ) ;
      GE_ExtrudedMeshing( GE_ExtrudedMeshing const& other ) ;
      GE_ExtrudedMeshing const& operator=( GE_ExtrudedMeshing const& other ) ;

      GE_ExtrudedMeshing( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          size_t dim_space ) ;

      virtual GE_ExtrudedMeshing* create_replica(
                                                PEL_Object* a_owner,
                                                PEL_ModuleExplorer const* exp,
                                                size_t dim_space ) const ;
     
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void build_vertices( void ) ;

      static GE_Color const* merge_color( GE_Color const* first,
                                          GE_Color const* second ) ;

      void recover_initial_sides( intArray2D& old_side_vertices,
                                  PEL_Vector* initial_face_color ) ;

      void build_vertical_side( size_t_vector& loc,
                                GE_Color const*& s_color,
                                size_t old_side,
                                size_t iz,
                                size_t iMesh,
                                size_t_vector const& vert,
                                intArray2D const& old_side_vertices,
                                PEL_Vector const* initial_face_color ) const ;

      void build_horizontal_side( size_t_vector& loc,
                                  GE_Color const*& s_color,
                                  bool up_side,
                                  size_t idx_z,
                                  size_t iMesh,
                                  doubleVector const& elev,
                                  size_t_vector const& vert ) const ;

      size_t extend_side( size_t_vector const& loc,
                          GE_Color const* s_color ) ;
      
      void build_sides( void ) ;

      void build_meshes( void ) ;

      void build_local( void ) ;

      doubleVector const& elevation( GE_Color const* cell ) const ;

      void find_next_mesh( void ) ;

      void find_next_cell( void ) ;

      void find_next_side( void ) ;

      GE_Color const* cell_color_in_initial_meshing( 
                                              std::string const& str ) const ;

      virtual void check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const ;

   //-- Class attributes
      
      static GE_ExtrudedMeshing const* prototype ;
      
   //-- Attributes

      GE_Meshing * initial ;
      std::string the_mesh_polyhedron_name ;
      mutable size_t_vector local_mesh_vertices ;
      mutable size_t_vector local_mesh_connectivity ;
      std::string side_poly_name ;
      std::string  mesh_poly_name ;
      size_t mesh_iterator_dim ;
      bool connectivity ;
      bool transfer ;
      size_t i_mesh ;

      GE_SetOfPoints * sop ;
      intArray3D global_sides ;
      intArray3D global_vertices ;
      intArray2D side_connectivity ;

      PEL_Vector* SIDE_COLORS ;
      size_t_vector nb_meshes ;
      size_t idx_vertex ;
      size_t idx_mesh ;
      mutable doubleVector vertex_coord ;
      size_t local ;
      size_t last_side ;
      GE_ReferencePolyhedron const* side_ref ;
      GE_ReferencePolyhedron const* mesh_ref ;
      GE_ReferencePolyhedron const* ref_poly ;
      GE_Color const* behind_color ;
      GE_Color const* front_color ;

      PEL_Vector* IC_COLORS ;
      PEL_Vector* ELEVATIONS ;
      doubleVector DEFAULT_ELEVATION ;
      size_t max_elevation_size ;

      PEL_BalancedBinaryTree* const SIDE_TREE ;

      GE_Color const* m_color ;
} ;

#endif


