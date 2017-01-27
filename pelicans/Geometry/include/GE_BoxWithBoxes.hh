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

#ifndef GE_BOX_WITH_BOXES_HH
#define GE_BOX_WITH_BOXES_HH

#include <GE_Meshing.hh>

#include <doubleArray2D.hh>
#include <size_t_vector.hh>

class intVector ;
class PEL_Vector ;

/*
Meshings for which
   - the geometric domain is a box in 1D, 2D or 3D ;
   - the cells are boxes.

A box is an orthogonal parallelotope whose sides are parallel
to the coordinate axes.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

      keyword                          data
      -------                          ----
vertices_coordinate_0   vector of the coordinates 0 of all vertices
vertices_coordinate_1   vector of the coordinates 1 of all vertices
                        in the 2D or 3D case
vertices_coordinate_2   vector of the coordinates 2 of all vertices
                        in the 3D case
mesh_polyhedron         vector of 2 strings : the side polyhedron name
                        and the cell polyhedron name
*/

class PEL_EXPORT GE_BoxWithBoxes : public GE_Meshing
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

   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_BoxWithBoxes( void ) ;
     ~GE_BoxWithBoxes( void ) ;
      GE_BoxWithBoxes( GE_BoxWithBoxes const& other ) ;
      GE_BoxWithBoxes const& operator=( GE_BoxWithBoxes const& other ) ;

      GE_BoxWithBoxes( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp,
                       size_t dim_space ) ;

      virtual GE_BoxWithBoxes* create_replica( PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp,
                                               size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;
      
      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      GE_Color const* corner_color( size_t_vector const& idx ) const ;
      
      size_t vertex_index_ijk( size_t_vector const& idx ) const ;

      size_t side_index_ijk( size_t_vector const& idx, size_t normal ) const ;

      virtual void check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const ;

   //-- Class attributes

      static GE_BoxWithBoxes const* PROTOTYPE ;

   //-- Attributes

      // Reference polyhedra:
      std::string FACE_POLY_NAME ;
      GE_ReferencePolyhedron const* FACE_POLY ;
      std::string CELL_POLY_NAME ;
      GE_ReferencePolyhedron const* CELL_POLY ;
      
      // Meshing:
      doubleArray2D IJK_COORDS ;
      size_t_vector NB_VERTS ;
      size_t_vector NB_FACES ;
      size_t_vector IBD ;

      // Colors:
      PEL_Vector* COLOR_TABLE ;
      PEL_Vector* GE_COLOR_TABLE ;

      // Iterators:
      size_t_vector I_POLY ;
      size_t_vector IBD_POLY ;
      mutable size_t_vector POLY_TO_VERTS ;
      bool CELL_IT ; // in iterators: iteration on cells, or on faces ?

      // Iterator on vertices:
      size_t IVERT ;
      mutable doubleVector VCOORDS ;

      // Iterator on faces:
      size_t IFACE ;
      size_t FACE_NORMAL ;

      // Iterator on cells:
      mutable size_t_vector CELL_TO_FACES ;
} ;

#endif


