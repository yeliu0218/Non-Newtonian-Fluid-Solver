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

#ifndef GE_SPLIT_MESHING_HH
#define GE_SPLIT_MESHING_HH

#include <GE_Meshing.hh>

#include <GE_SplittingStrategy.hh>

#include <size_t_vector.hh>
#include <boolVector.hh>

/*
Meshings obtained by splitting another initial meshing.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that contains
   - a submodule of name "GE_Meshing" describing the entire meshing to be split
   - a submodule of name "splitting_strategy" describing how to split the global
     meshing.
     The implemented splitting strategy is a split along the coordinates :
       keyword                          data
       -------                          ----
     type                           "coordinate_splitting"
     coordinate_splitting_formula   formula that gives in function of its
                                    coordinates $DV_X, the rank which owns
                                    the different mesh centers of the entire 
                                    meshing                             
   - the following keyword-data pair :
       keyword                          data
       -------                          ----
     rank                           rank of the current split meshing. Modifing
                                    rank from 0 up to nb_ranks, one can build
                                    the several splits of the global meshing
     nb_ranks                       total number of ranks
     security_bandwith              number of meshes that are adjacent to the
                                    meshing defined with the "splitting_strategy"
                                    module that have to be added to the split
                                    meshing

Remarks on the colorization of the security bandwith geometrical elements :
        - vertices keep their initial colors ;
        - sides are colorized with `GE_Color::halo_color()' except for
          initial external bounds which keep their initial colors
          (used for boundary conditions, grid external hull, ...) ;
        - cells are colorized with `GE_Color::halo_color()'.

Example :

   MODULE GE_Meshing
      concrete_name = "GE_SplitMeshing"
      nb_ranks = 4
      rank = 2
      security_bandwith = 1
      MODULE GE_Meshing
         concrete_name = "GE_BoxWithBoxes"
         vertices_coordinate_0 = regular_vector( 0.0, 4, 4.0 )
         vertices_coordinate_1 = regular_vector( 0.0, 4, 4.0 )
         mesh_polyhedron = < "GE_Segment" "GE_Rectangle" > 
      END MODULE GE_Meshing
      MODULE splitting_strategy
         concrete_name = "GE_CoordinateSplitting"
         $DS_X = component( $DV_X, 0 )
         $DS_Y = component( $DV_X, 1 )
         coordinate_splitting_formula = 
            ( $DS_X < 1.0 && $DS_Y < 1.0 ? 2 :
              $DS_X < 1.0 && $DS_Y > 1.0 ? 0 :
              $DS_X > 1.0 && $DS_Y > 1.0 ? 3 :
              1 )
      END MODULE splitting_strategy
   END MODULE GE_Meshing


   Cells distribution (rank that owns the different cells) :
 
               0  3  3  3 
               0  3  3  3
               0  3  3  3
               2  1  1  1
         

   For rank=2, keeps the mesh at bottom left, and, a security bandwith
   including 3 meshes of color `GE_Color::halo_color()'.

PUBLISHED
*/

class PEL_EXPORT GE_SplitMeshing : public GE_Meshing
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
      
      virtual void display( std::ostream& os, size_t indent_width ) ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_SplitMeshing( void ) ;
     ~GE_SplitMeshing( void ) ;
      GE_SplitMeshing( GE_SplitMeshing const& other ) ;
      GE_SplitMeshing& operator=( GE_SplitMeshing const& other ) ;

      GE_SplitMeshing( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp,
                       size_t dim_space ) ;

      virtual GE_SplitMeshing* create_replica( PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp,
                                               size_t dim_space ) const ;

      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void find_next_vertex( void ) ;
      
      void find_next_cell( void ) ;

      void find_next_face( void ) ;
   
   //-- Meshing split
      
      void split_meshing( void ) ;
      
   //-- Class attributes

      static GE_SplitMeshing const* PROTOTYPE ;
      
   //-- Attributes

      GE_Meshing* INITIAL ;
      GE_SplittingStrategy* SPLIT_STRAT ;

      // Vertices:
      size_t NB_VERTS ;
      size_t_vector LOCAL_VERT ;
      size_t IDX_VERT ;

      // Faces:
      size_t NB_FACES ;
      size_t_vector LOCAL_FACE ;
      boolVector HALO_FACE ;
      size_t IDX_FACE ;

      // Cells:
      size_t NB_CELLS ;
      size_t_vector LOCAL_CELL ;
      size_t IDX_CELL ;
      
      // Connectivity:
      size_t_vector ELM_2_VERTS ;
      size_t_vector CELL_2_FACES ;
      
      // Halo zone:
      size_t const BAND ;
      GE_Color const* const HALO ;
} ;

#endif


