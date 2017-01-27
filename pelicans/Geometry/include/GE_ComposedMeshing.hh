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

#ifndef GE_COMPOSED_MESHING_HH
#define GE_COMPOSED_MESHING_HH

#include <GE_Meshing.hh>

#include <GE_SetOfPoints.hh>
#include <PEL_Vector.hh>

#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <stringVector.hh>

class GE_Point ;

class PEL_Vector ;

class size_t_vector ;

/*
Meshings obtained by merging a set of meshings.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that contains
   - a submodule of name "list_of_GE_Meshing" containing a set of modules
     that describe several not overlapping meshings
   - a keyword of name "check_meshing" that, if true, enables the checking
     of the created meshing (no overlappings). One must keep in mind that
     this operation is in O(`::nb_vertices()'*`::nb_cells()').

Example :

   MODULE GE_Meshing
      concrete_name = "GE_ComposedMeshing"
      check_meshing = true
      MODULE list_of_GE_Meshing
         MODULE GE_Meshing#1
            concrete_name = "GE_BoxWithBoxes"
            vertices_coordinate_0 = regular_vector( 0.0, 2, 1.0 )
            vertices_coordinate_1 = regular_vector( 0.0, 2, 1.0 )
            mesh_polyhedron = < "GE_Segment" "GE_Rectangle" > 
         END MODULE GE_Meshing#1
         MODULE GE_Meshing#2
            concrete_name = "GE_RectangleWithTriangles"
            MODULE GE_Meshing
               concrete_name = "GE_BoxWithBoxes"
               vertices_coordinate_0 = regular_vector( 0.5, 2, 1.5 )
               vertices_coordinate_1 = regular_vector( 1.0, 2, 2.0 )
               mesh_polyhedron = < "GE_Segment" "GE_Rectangle" >
            END MODULE GE_Meshing
            mesh_polyhedron = < "GE_Segment" "GE_Triangle" >
            MODULE refinement_strategy 
               default = "/" 
            END MODULE refinement_strategy
         END MODULE GE_Meshing#2
       END MODULE list_of_GE_Meshing
   END MODULE GE_Meshing

Remarks on element colors :
   - if a cell or an external bound has the color of name "colo" and is in
     the `i'-th meshing, its color is "m`i'_colo" in the merged meshing.
        examples: "m0_bottom_right"
                  "m1_top_right"
   - if a vertex or a side has the color of name "colo" and is in only
     one meshing (the `i'-th one),  its color is "m`i'_colo" in the merged
     meshing.
        examples: "m0_bottom_right"
                  "m1_top_right"
   - if a vertex or a side is shared by several meshings, its color is
     the concatenation of the differents colors described in the previous
     section.
        examples: "m0_top_m1_bottom"

   Meshings of default names "m0", "m1", ... can be rename:

        example:
        
           MODULE GE_Meshing
              concrete_name = "GE_ComposedMeshing"
              meshing_names = < "local1", "door", "local2" >
              MODULE list_of_GE_Meshing
                 ...
              END MODULE list_of_GE_Meshing
           END MODULE GE_Meshing

           colors are: "local1_right"
                       "door_bottom"
   

*/
class PEL_EXPORT GE_ComposedMeshing : public GE_Meshing
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

      GE_ComposedMeshing( void ) ;
     ~GE_ComposedMeshing( void ) ;
      GE_ComposedMeshing( GE_ComposedMeshing const& other ) ;
      GE_ComposedMeshing& operator=( GE_ComposedMeshing const& other ) ;

      GE_ComposedMeshing( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          size_t dim_space ) ;

      virtual GE_ComposedMeshing* create_replica(
                                                PEL_Object* a_owner,
                                                PEL_ModuleExplorer const* exp,
                                                size_t dim_space ) const ;

      virtual GE_Color const* default_vertex_color( void ) const ;
      virtual GE_Color const* default_cell_color( void ) const ;
      virtual GE_Color const* default_face_color( void ) const ;

  //-- Builder
      
      void build_vertices( void ) ;
      void build_cells( void ) ;
      void build_sides( void ) ;

      size_t_vector const& side_global_vertices(
                                size_t i_meshing ) const ;

   //-- Check

      void check_meshing( void ) ;
      
   //-- Class attributes

      static GE_ComposedMeshing const* PROTOTYPE ;
      
   //-- Attributes

      // Meshings :
      PEL_Vector* const MESHINGS ;
      stringVector MESHING_COLS ;
      
      // Vertices :
      GE_SetOfPoints* const VERTICES ;
      size_t_array2D VERT_LOC_TO_GLOB ;
      size_t VERT_INDEX ;

      // Cells :
      size_t NB_CELLS ;
      GE_Meshing* CELL_MESHING ;
      size_t CELL_MESHING_INDEX ;

      // Sides :
      size_t NB_SIDES ;
      size_t_array2D SIDE_LOC_TO_GLOB ;
      boolVector     SIDE_OK ;
      PEL_Vector*   SIDE_COLORS ;
      GE_Meshing* SIDE_MESHING ;
      size_t SIDE_MESHING_INDEX ;
      size_t SIDE_INDEX ;
} ;

#endif


