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

#ifndef GE_PerforatedMeshing_HH
#define GE_PerforatedMeshing_HH

#include <GE_Meshing.hh>

#include <size_t_vector.hh>
#include <boolVector.hh>

#include <vector>

class PEL_Data ;
class PEL_DoubleVector ;
class PEL_ModuleExplorer ;

/*
Meshings obtained by perforating another initial one i.e. by removing 
some defined zones (holes) in an original meshing.

The new mesh is obtained by selecting somes cells of the original
one using boolean formulas, by removing those selected cells and 
then, by building the new connectivities. A cell to be excluded
can not lie in more than one hole.

Internal sides of the original meshing that turn to be bounds of the new
meshing will take a color, the corresponding hole color, that have to be
precised while defining the selections of vertices to be excluded.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that contains
   - a submodule of name "GE_Meshing" defining the original meshing
   - a submodule of name "holes" containing a list of 
     keyword-data pairs, each one describes an hole:
       hole_color = boolean_formulae( $DV_X )
     ($DV_X beeing the center position of all initial meshing cells,
      boolean_formulae is true for cells to be excluded).

     Ex :
      MODULE GE_Meshing
         
         concrete_name = "GE_PerforatedMeshing"
         
        MODULE holes
           hole = in_box( $DV_X, <0.4 0.4 -0.1>, <0.6 0.6 1.0> )
        END MODULE holes

        MODULE GE_Meshing
            
            concrete_name = "GE_BoxWithBoxes"
            vertices_coordinate_0 = regular_vector( 0.0, 10, 1.0 )
            vertices_coordinate_1 = regular_vector( 0.0, 10, 2.0 )
            vertices_coordinate_2 = regular_vector( 0.0, 10, 3.0 )
            mesh_polyhedron = < "GE_Rectangle" "GE_Cuboid" >
        END MODULE GE_Meshing
        
      END MODULE GE_Meshing
*/


class PEL_EXPORT GE_PerforatedMeshing : public GE_Meshing
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

      GE_PerforatedMeshing( void ) ;
     ~GE_PerforatedMeshing( void ) ;
      GE_PerforatedMeshing( GE_PerforatedMeshing const& other ) ;
      GE_PerforatedMeshing const& operator=(
                            GE_PerforatedMeshing const& other ) ;

      GE_PerforatedMeshing( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            size_t dim_space ) ;

      virtual GE_PerforatedMeshing* create_replica(
                            PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            size_t dim_space ) const ;

      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      
      void find_next_vertex( void ) ;

      void find_next_cell( void ) ;

      void find_next_face( void ) ;
      
      void build_meshing( std::vector<PEL_Data*> const& hole_formula,
                          PEL_DoubleVector* coords ) ;
      
   //-- Class attributes

      static GE_PerforatedMeshing const* PROTOTYPE ;
      
   //-- Attributes

      // initial meshing:
      GE_Meshing* INITIAL ;

      // hole colors:
      std::vector<GE_Color const*> HOLE_COLORS ;

      // cells:
      size_t NB_CELLS ;
      boolVector REMOVED_CELLS ;
      
      // vertices:
      size_t NB_VERTS ;
      size_t_vector INITIAL_TO_LOCAL_VERTS ;
      size_t_vector VERTS_COL ;

      // faces:
      size_t NB_FACES ;
      size_t_vector INITIAL_TO_LOCAL_FACES ;
      size_t_vector FACES_COL ;

      // iterators:
      size_t IDX_VERT ;
      size_t IDX_FACE ;
      size_t IDX_CELL ;
      size_t_vector LOCAL_VERTS ;
      size_t_vector LOCAL_CONNS ;
} ;

#endif


