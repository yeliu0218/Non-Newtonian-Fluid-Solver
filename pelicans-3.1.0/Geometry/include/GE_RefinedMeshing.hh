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

#ifndef GE_REFINED_MESHING_HH
#define GE_REFINED_MESHING_HH

#include <GE_Meshing.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>

class PEL_ModuleExplorer ;
class PEL_Vector ;

class GE_Point ;
class GE_ReferencePolyhedronRefiner ;
class GE_SetOfPoints ;

/*
Meshings obtained by subdiving cells of another meshing.

Instances of subclasses are  delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has 
   - a submodule of name "GE_Meshing" describing the meshing to be refined
   - a submodule of name "list_of_GE_ReferencePolyhedronRefiner" that itself 
     contains a list of submodules describing the refinement pattern of each 
     of the reference polyhedra occuring in the meshing to be refined
   - the following keyword-data pair :
                keyword                          data
                -------                          ----
     mesh_polyhedron          vector of 2 strings : the side polyhedron name
                              and the cell polyhedron name
*/

class PEL_EXPORT GE_RefinedMeshing : public GE_Meshing
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

      GE_RefinedMeshing( void ) ;
     ~GE_RefinedMeshing( void ) ;
      GE_RefinedMeshing( GE_RefinedMeshing const& other ) ;
      GE_RefinedMeshing& operator=( GE_RefinedMeshing const& other ) ;

      GE_RefinedMeshing( PEL_Object* a_owner,
                         PEL_ModuleExplorer const* exp,
                         size_t dim_space ) ;

      virtual GE_RefinedMeshing* create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space ) const ;

      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void build_vertices( void ) ;

      void build_sides( void ) ;

      void build_cells( void ) ;

      GE_Color const* local_color( GE_Point const* pt ) const ;
      
      void find_next_cell( void ) ;

      void find_next_side( void ) ;
      
      void append_refiner( GE_Color const* col,
                           GE_ReferencePolyhedronRefiner const* refi ) ;
      
      GE_ReferencePolyhedronRefiner const* refiner( 
                                           GE_Meshing const* mm ) const ;
      
   //-- Class attributes
      
      static GE_RefinedMeshing const* PROTOTYPE ;
      
   //-- Attributes

      GE_Meshing* initial ;
      
      std::string the_mesh_polyhedron_name ;
      mutable size_t_vector MESH_VERTICES ;
      mutable size_t_vector CELL_SIDES ;
      std::string side_poly_name ;
      std::string  mesh_poly_name ;
      bool CELL_IT ;
      bool SIDE_IT ;
      bool connectivity ;
      bool transfer ;
      size_t i_mesh ;

      GE_SetOfPoints * sop ;
      intArray2D global_sides ;
      intArray2D global_vertices ;

      PEL_Vector* initial_side_poly ;
      PEL_Vector* initial_side_color ;
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
      bool verbose ;
            
      PEL_Vector* REFs ;
      PEL_Vector* REF_COLs ;
} ;

#endif


