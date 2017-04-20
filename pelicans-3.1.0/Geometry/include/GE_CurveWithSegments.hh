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

#ifndef GE_CURVE_WITH_SEGMENTS_HH
#define GE_CURVE_WITH_SEGMENTS_HH

#include <GE_Meshing.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

/*
Meshings for which
   - the geometric domain is a curve made of connected straight lines 
   - the cells are segments.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

      keyword                          data
      -------                          ----
   bends                   DoubleArray2D representing the coordinates of
                           the points connecting the straight lines
   subdivisions            DoubleVector made of a concatenation of vectors
                           whose items grow from 0.0 to 1.0, each of these
                           vectors representing  the partitioning
                           of one straight line into segments
   mesh_polyhedron         vector of 2 strings : the side polyhedron name
                           and the cell polyhedron name
*/


class PEL_EXPORT GE_CurveWithSegments : public GE_Meshing
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

      GE_CurveWithSegments( void ) ;
     ~GE_CurveWithSegments( void ) ;
      GE_CurveWithSegments( GE_CurveWithSegments const& other ) ;
      GE_CurveWithSegments const& operator=( 
                           GE_CurveWithSegments const& other ) ;

      GE_CurveWithSegments( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp,
                           size_t dim_space ) ; 

      virtual GE_CurveWithSegments* create_replica(
                                                PEL_Object* a_owner,
                                                PEL_ModuleExplorer const* exp,
                                                size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      virtual void check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const ;

   //-- Class attributes
      
      static GE_CurveWithSegments const* prototype ;
      
   //-- Attributes

      size_t NB_SP_DIMS ;
      bool CLOSED_CURVE ;
      doubleArray2D BENDS ;
      size_t NB_VERTICES ;
      size_t NB_SEG ;
      size_t_vector NB_SEG_VERT ;
      doubleArray2D X_SEG_VERT ;
      size_t i_vert ;
      size_t i_seg ;
      size_t i_seg_vert ;
      size_t i_mesh ;
      bool CELL_IT ;
      bool SIDE_IT ;
      mutable doubleVector VERT_COORD ;
      mutable size_t_vector mesh2verts ;
      mutable size_t_vector mesh2sides ;
      std::string SIDE_POLY_NAME ;
      std::string CELL_POLY_NAME ;

} ;

#endif
