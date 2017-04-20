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

#ifndef GE_TRANSFORMED_MESHING_HH
#define GE_TRANSFORMED_MESHING_HH

#include <GE_Meshing.hh>

#include <doubleArray2D.hh>

class PEL_Bool ;
class PEL_ContextSimple ;
class PEL_Double ;
class PEL_DoubleVector ;
class PEL_Int ;

/*
Meshings obtained by applying a geometrical transformation to another meshing.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that contains 
   - a submodule of name "GE_Meshing" describing the meshing to be transformed
     (hereafter called initial meshing)
   - the following keyword-data pair :
                keyword                          data
                -------                          ----
     transformation           expression of type doubleVector, depending on the
                              variable $DV_X.
                              
Example :

      MODULE GE_Meshing
         concrete_name = "GE_TransformedMeshing"
         $DS_x = component( $DV_X, 0 )
         $DS_y = component( $DV_X, 1 )
         transformation = vector( 2.*$DS_x, -$DS_y )
         MODULE GE_Meshing
            concrete_name = "GE_GambitMeshing"
            mesh_polyhedron = < "GE_Segment" "GE_Quadrilateral" >
            filename = join( this_file_dir(), "cook.neu" )
         END MODULE GE_Meshing
      END MODULE GE_Meshing
      
The difference between `self' and the initial meshing relies only in the
coordinates of the vertices, which are calculated as follows.
   For each vertex of the original meshing
   -    assign the coordinates of this vertex to $DV_X ;
   -    evaluate the data of keyword "tranformation" ;
   -    the evaluation result gives the coordinates of the associated vertex 
        in `self'.
All other characteristics of `self' (connectivities, colors) are identical to
those of the initial meshing.

The data of keyword "tranformation" is evaluated in a context that may
contain the following variables:
   $DV_X, whose value is the coordinates of the current vertex ;
   $IS_VERT_IDX, whose value is the index of the current vertex ;
   $DS_VERT_MIN, whose value is the minimal distance between the 
                 current vertex and its neighbor vertices ;
   $BS_VERT_ON_BOUND, whose value `true' iff the current vertex
                      is located on the meshing boundary.
*/

class PEL_EXPORT GE_TransformedMeshing : public GE_Meshing
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

      GE_TransformedMeshing( void ) ;
     ~GE_TransformedMeshing( void ) ;
      GE_TransformedMeshing( GE_TransformedMeshing const& other ) ;
      GE_TransformedMeshing const& operator=( GE_TransformedMeshing const& other ) ;

      GE_TransformedMeshing( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t dim_space ) ;

      virtual GE_TransformedMeshing* create_replica( PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp,
                                              size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void set_meshing_data( void ) ;

   //-- Class attributes

      static GE_TransformedMeshing const* PROTOTYPE ;

   //-- Attributes
      
      GE_Meshing* INITIAL ;
      PEL_ContextSimple* CTX ;
      PEL_DoubleVector* COORDS ;
      PEL_Double* HMIN ;
      PEL_Bool* IS_BOUNDARY ;
      PEL_Int* VERT_IDX ;
      PEL_DataWithContext const* TRANSFO ;

      size_t IVERT ;
      doubleArray2D TRANSFORMED_COORDINATES ;
} ;

#endif


