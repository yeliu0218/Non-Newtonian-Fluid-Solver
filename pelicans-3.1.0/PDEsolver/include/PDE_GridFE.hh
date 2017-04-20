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

#ifndef PDE_GRID_FE_HH
#define PDE_GRID_FE_HH

#include <PEL_Object.hh>

#include <vector>

#include <size_t_vector.hh>

class PEL_BalancedBinaryTree ;
class PEL_Iterator ;
class PEL_IndexSet ;
class PEL_List ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class PEL_VectorIterator ;
class size_t_vector ;

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;
class GE_ReferencePolyhedronRefiner ;
class GE_SegmentSegment_INT ;
class GE_SetOfPoints ;
class GE_SimplePolygon2D ;

class PDE_BasisFunctionCell ;
class PDE_BoundFE ;
class PDE_CellFE ;
class PDE_DiscreteField ;
class PDE_FaceFE ;
class PDE_MeshFE ;
class PDE_ReferenceElement ;

class PEL_EXPORT PDE_GridFE : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance creation and initialization

      static PDE_GridFE* create( PEL_Object* a_owner,
                                 GE_SetOfPoints* a_vertex_manager,
                                 size_t a_nb_cells,
                                 size_t a_nb_faces ) ;
      
      // Update the specific attributes of `self' in case of moving of its
      // vertices.
      virtual void update( void ) ;
      
   //-- Discrete fields
         
      void extend_mesh_discretizations( PDE_DiscreteField const* ff,
                                        bool on_bounds,
                                        PEL_Vector const* elms ) ;
      
      void duplicate_discretization( PDE_DiscreteField const* model_f,
                                     PDE_DiscreteField* new_f,
                                     bool on_bounds ) ;
         
   //-- Building
      
      void build_face( size_t number,
                       std::string const& polyhedron_name,
                       size_t_vector const& face_vertices,
                       GE_Color const* color ) ;      
      
      void build_cell( size_t number,
                       std::string const& polyhedron_name,
                       size_t_vector const& cell_vertices,
                       size_t_vector const& cell_faces,
                       GE_Color const* color ) ;
      
      void set_periodicity_requests( PEL_ModuleExplorer const* exp ) ;
      
      void separate_internal_and_external_faces( void ) ;
      
      void set_decohesion_requests( PEL_ModuleExplorer const* exp ) ;
      
   //-- Clusters
      
      size_t nb_clusters( void ) const ;
      
      void define_mesh_as_cell_of_cluster( void ) ;
      
      void define_mesh_as_bound_of_cluster( void ) ;
      
      bool mesh_is_defined( void ) const ;
      
      size_t nb_meshes_for_all_clusters( void ) const ;
      
      bool mesh_belongs_to_cluster( size_t id_number_of_mesh, 
                                    size_t i_cluster ) const ; 
      
      void start_mesh_iterator( size_t i_cluster ) ;
      
      bool valid_mesh( void ) const ;
      
      void go_next_mesh( void ) ;
      
      PDE_MeshFE* current_mesh( void ) const ;
      
   //-- Node determination from a set of node points
      
      void start_nodepts_recording( void ) ;
      
      bool ongoing_nodepts_recording( void ) const ;
      
      void start_preexisting_nodepts_recording( void ) ;
      
      bool ongoing_preexisting_nodepts_recording( void ) const ;
      
      void record_preexisting_nodept( PDE_BasisFunctionCell const* bf, 
                                      size_t& bf_index ) ;
      
      void terminate_preexisting_nodepts_recording( void ) ;
      
      void record_nodept( GE_Mpolyhedron const* poly, 
                          GE_Point const* pt_ref,
                          size_t& nodept_index, bool& is_new ) ;
      
      void terminate_nodepts_recording( void ) ;
      
      size_t nb_recorded_nodepts( void ) const ;
      
      size_t nb_inferred_nodes( void ) const ;
      
      size_t inferred_node_of_recorded_nodept( size_t nodept_index ) const ;
      
      size_t bf_index_of_recorded_nodept( size_t nodept_index ) const ;
      
   //-- Vertices
            
      size_t nb_space_dimensions( void ) const ;

      GE_SetOfPoints* set_of_vertices( void ) const ;

   //-- Cells
      
      size_t nb_cells( void ) const ;

      PEL_Vector const* cells( void ) const ;
            
      void refine_cell( PDE_CellFE* ccell,
                        GE_ReferencePolyhedronRefiner const* crf ) ;

   //-- Sides (faces not located on the domain boundary)
      
      size_t nb_sides( void ) const ;

      PEL_Vector const* sides( void ) const ;
            
   //-- Bounds (faces located on the domain boundary)
      
      size_t nb_bounds( void ) const ;

      PEL_Vector const* bounds( void ) const ;
      
   //-- State
      
      void set_verbose_level( size_t a_verbose_level ) ;
      
      size_t geometric_state( void ) const ;
      
      void change_geometric_state( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
   
   //-- Hidden

      // mesh which contains `pt'
      PDE_CellFE* get_mesh(
                     GE_Point const* pt,
                     GE_SegmentSegment_INT const* seg_seg_intersector,
                     PDE_CellFE* guess_mesh=0 ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_GridFE( void ) ;
     ~PDE_GridFE( void ) ;
      PDE_GridFE( PDE_GridFE const& other ) ;
      PDE_GridFE& operator=( PDE_GridFE const& other ) ;

      PDE_GridFE( PEL_Object* a_owner,
                  GE_SetOfPoints* a_vertex_manager,
                  size_t a_nb_cells,
                  size_t a_nb_sides ) ;
      
   //-- Internals
      
      static void extend_mesh_discretizations( 
                                        bool elm_must_exist,
                                        PEL_Vector* discs,
                                        PDE_DiscreteField const* ff,
                                        PEL_Vector const* elms ) ;
      
      static void duplicate_discretization( PDE_DiscreteField const* model_f,
                                            PDE_DiscreteField* new_f,
                                            PEL_Vector const* discs,
                                            PEL_Vector const* meshes ) ;
      
      void add_cell( PDE_CellFE* cell ) ;
      
      void add_bound( PDE_BoundFE* bound ) ;

      void add_side( PDE_FaceFE* side ) ;

      void build_bound( size_t number, 
                        size_t side_number, 
                        GE_Color const* color ) ;
      
      void build_new_vertices( GE_ReferencePolyhedronRefiner const* crf,
                               PDE_CellFE const* ccell,
                               size_t_vector& idx_in_VERTS ) ;
      
      void record_already_existing_faces(
                               GE_ReferencePolyhedronRefiner const* crf,
                               PDE_CellFE const* ccell,
                               PEL_Vector* subfaces,
                               PEL_BalancedBinaryTree* tree ) ;
      
      void build_new_faces( GE_ReferencePolyhedronRefiner const* crf,
                            PDE_CellFE const* ccell,
                            size_t_vector const& idx_in_VERTS,
                            PEL_BalancedBinaryTree const* tree,
                            PEL_Vector* subfaces,
                            size_t_vector& idx_in_subfaces ) ;
      
      void build_new_cells( GE_ReferencePolyhedronRefiner const* crf,
                            PDE_CellFE* ccell,
                            size_t_vector const& idx_in_VERTS,
                            PEL_Vector const* subfaces,
                            size_t old_faces_limit,
                            size_t_vector const& idx_in_subfaces ) ;
      
      void replace_adjacent_cell( PDE_FaceFE* face,
                                  PDE_CellFE* old_cell,
                                  PDE_CellFE* new_cell ) ;

      void add_index_set( PEL_BalancedBinaryTree* tree,
                          PDE_MeshFE const* a_mesh,
                          size_t an_id ) ;
      
      void separate_sides_and_bounds( void ) ;

      void separate_sides_and_bounds( PEL_Iterator* it_tr ) ;
      
      void build_cluster_to_mesh_connectivity( PEL_Vector const* meshes ) ;
      
      static size_t associated_point( GE_SetOfPoints const* points, 
                                      PEL_ListIterator* tr_it,
                                      size_t ip,
                                      GE_Point* work_pt ) ;
      
      static void fill_pt_2_node( 
                        size_t ip, size_t_vector const& ass,
                        size_t_vector& pt_2_node, size_t& current_node ) ;
      
      bool rejected_mesh( void ) const ;
      
      static PDE_MeshFE* the_mesh( PEL_VectorIterator const* it ) ;
      
   //-- Hidden

      PDE_CellFE* search_mesh(
                       GE_Point const* pt,
                       GE_SegmentSegment_INT const* seg_seg_intersector,
                       PDE_CellFE* guess_mesh,
                       size_t recursive_max ) const ;
      
      
   //-- Attributes
      
      GE_SetOfPoints* const VERTS ;
      
      PEL_Vector* C_DISCS ;
      PEL_Vector* B_DISCS ;

      PEL_Vector* CELLS  ;  // list of PDE_CellFE*
      PEL_Vector* FACES  ;  // list of PDE_FaceFE*
      PEL_Vector* SIDES  ;  // list of PDE_FaceFE*
      PEL_Vector* BOUNDS ;  // list of PDE_BoundFE*
      
      PEL_VectorIterator* MESH_IT ;
      
      bool OK_MESH_IT ;
      size_t I_CURRENT_CLUSTER ;
      size_t NB_MESHES ;
      
      GE_SetOfPoints* NEW_NODEPTS ;
      GE_SetOfPoints* PRE_NODEPTS ;
      bool ON_PRE ;
      size_t NB_PRE_NODEPTS ;
      size_t NB_NEW_NODEPTS ;
      size_t_vector PT_2_NODE ;
      size_t NB_NODES ;
      
      size_t GEO_STATE ;
      
      PEL_List* TRANSFOS ;
      PEL_ListIterator* TRANSFOS_IT ;
      GE_Point* W_PT ;
      
      bool FULL_DECOHESION ;
      size_t NB_CLUSTERS ;
      std::vector< GE_Color const* > CLUSTER_COLOR ;
      size_t_vector CLUSTER_OF_MESH ;
            
      PEL_IndexSet* IDX_TEST ;
      size_t VERB ;
} ;

#endif
