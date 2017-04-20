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

#ifndef GE_GAMBIT_MESHING_HH
#define GE_GAMBIT_MESHING_HH

#include <GE_Meshing.hh>

#include <fstream>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <string>
#include <vector>

/*
Meshings stored in an GAMBIT output file under the NEUTRAL format.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

   keyword                         data
   -------                         ----
filename                pathname of the GAMBIT output file
mesh_polyhedron         polyhedron name for sides and cells
*/

class PEL_BalancedBinaryTree ;
class PEL_ModuleExplorer ;

using std::string ;

class PEL_EXPORT GE_GambitMeshing : public GE_Meshing
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

      GE_GambitMeshing( void ) ;
     ~GE_GambitMeshing( void ) ;
      GE_GambitMeshing( GE_GambitMeshing const& other ) ;
      GE_GambitMeshing& operator=( GE_GambitMeshing const& other ) ;

      GE_GambitMeshing( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp,
                        size_t dim_space ) ;

      virtual GE_GambitMeshing* create_replica( PEL_Object* a_owner,
                                                PEL_ModuleExplorer const* exp,
                                                size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void read_vertex( void ) ;

      void read_node_to_point_coordinate_data( size_t& elm_geom_type ) ;

      void read_cell( void ) ;

      void read_side( void ) ;

      void go_i_th_next_line( size_t i ) ;

      void build_sides_connectivity( size_t cell_number ) ;

      void check_meshes_consistency( void ) ;
      
      size_t extend_list_of_sides( void ) ;

      void build_colors_of_vertices( void ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_GambitMeshing const* PROTOTYPE ;

   //-- Attributes

      // Input file stream
      std::ifstream INPUT ;
      size_t NTYPE ;
      std::streampos VERT_POS ;
      std::streampos CELL_POS ;

      // Vertices ...
      // COLOR_OF_VERT(igv) = color of the igv-th vertex
      size_t NB_VERTS ;
      std::vector< GE_Color const* > COLOR_OF_VERT ;
      size_t I_VERT ;
      doubleVector COORD_OF_CURRENT_VERT ;

      // Cells ...
      // SIDE_OF_CELL(igc,is) = global number of the is-th side
      //                        in the igc-th cell
      // COLOR_OF_CELL(igc) = color of the igc-th cell
      size_t NB_CELLS ;
      std::string CELL_POLY_NAME ;
      size_t NB_SIDES_PER_CELL ;
      PEL_BalancedBinaryTree* const SIDE_TREE ;
      size_t_array2D SIDE_OF_CELL ;
      std::vector< GE_Color const* > COLOR_OF_CELL ;
      size_t I_CELL ;
      mutable size_t_vector VERT_OF_CURRENT_CELL ;
      mutable size_t_vector SIDE_OF_CURRENT_CELL ;
      
      // Sides ...
      // VERT_OF_SIDE(igs,iv) = global number of the iv-th vertex
      //                        in the igs-th side
      size_t NB_SIDES ;
      std::string SIDE_POLY_NAME ;
      size_t_array2D VERT_OF_SIDE ;
      std::vector< GE_Color const* > COLOR_OF_SIDE ;
      size_t I_SIDE ;
      mutable size_t_vector VERT_OF_CURRENT_SIDE ;
} ;

#endif


