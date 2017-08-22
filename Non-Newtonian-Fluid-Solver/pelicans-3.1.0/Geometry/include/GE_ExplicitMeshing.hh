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

#ifndef GE_EXPLICIT_TRIANGULATION_HH
#define GE_EXPLICIT_TRIANGULATION_HH

#include <GE_Meshing.hh>

class intVector ;
class PEL_Vector ;

/*
Arbitrary meshings.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that contains the
explicit description of the sets of vertices, cells and sides
with their connections.

PUBLISHED
*/

class PEL_EXPORT GE_ExplicitMeshing : public GE_Meshing
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

      GE_ExplicitMeshing( void ) ;
     ~GE_ExplicitMeshing( void ) ;
      GE_ExplicitMeshing( GE_ExplicitMeshing const& other ) ;
      GE_ExplicitMeshing& operator=( GE_ExplicitMeshing const& other ) ;

      GE_ExplicitMeshing( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          size_t dim_space ) ;

      virtual GE_ExplicitMeshing* create_replica(
                                               PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp,
                                               size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void recover_color_table( PEL_ModuleExplorer* color_identification,
                                PEL_ModuleExplorer* special_color,
                                size_t nb_elements,
                                PEL_Vector* color_table ) const ;

      void take_new_mesh_partition( void ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Class attributes
      
      static GE_ExplicitMeshing const* PROTOTYPE ;

   //-- Attributes

      PEL_ModuleExplorer const* const EXP ;
      PEL_ModuleExplorer* MESHES_EXP ;

      // Vertices:
      doubleVector const& VERT_COORDS ;
      PEL_Vector* const VERT_COLORS ;
      size_t IVERT ;
      mutable doubleVector A_VERT_COORDS ;

      // Cells/Faces:
      PEL_Vector* const MESH_COLORS ;
      std::string MESH_POLY_NAME ;
      bool CELL_IT ;
      bool FACE_IT ;
      size_t I_MESH ;
      size_t I_GLOB_MESH ;
      intVector const* MESH_IDX ;
      mutable size_t_vector MESH_2_VERTS ;
      bool HAS_MESH_2_FACES ;
      mutable size_t_vector MESH_2_FACES ;
      bool VALID_MESH ;
      size_t DECAL_MESH ;
} ;

#endif


