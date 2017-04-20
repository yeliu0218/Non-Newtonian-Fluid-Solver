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

#ifndef GE_CACHED_MESHING_HH
#define GE_CACHED_MESHING_HH

#include <GE_Meshing.hh>

#include <PEL_String.hh>
#include <PEL_IntVector.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_StringVector.hh>
#include <PEL_Vector.hh>

#include <doubleArray2D.hh>
#include <size_t_array2D.hh>

class PEL_Vector ;

/*
  Cached Meshing.
  This class is provided to improve reading of file related meshing to
  avoid multiple readings of them.
  It assumes that only the same polyhedron is used for all cells and faces.

*/

class PEL_EXPORT GE_CachedMeshing : public GE_Meshing
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

      GE_CachedMeshing( void ) ;
     ~GE_CachedMeshing( void ) ;
      GE_CachedMeshing( GE_CachedMeshing const& other ) ;
      GE_CachedMeshing& operator=( GE_CachedMeshing const& other ) ;

      GE_CachedMeshing( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          size_t dim_space ) ;

      virtual GE_CachedMeshing* create_replica(
                                                PEL_Object* a_owner,
                                                PEL_ModuleExplorer const* exp,
                                                size_t dim_space ) const ;

      virtual GE_Color const* default_vertex_color( void ) const ;
      virtual GE_Color const* default_cell_color( void ) const ;
      virtual GE_Color const* default_face_color( void ) const ;

  //-- Builder
      
      void build_vertices( GE_Meshing* meshing ) ;
      void build_cells( GE_Meshing* meshing ) ;
      void build_faces( GE_Meshing* meshing ) ;

      void read_meshing( std::string const& filename ) ;
      void save_meshing( std::string const& filename ) ;
      
   //-- Class attributes

      static GE_CachedMeshing const* PROTOTYPE ;
      
   //-- Attributes
      
      // Vertices :
      doubleArray2D VERTICES ;
      PEL_Vector* VERTEX_COLORS ;
      size_t VERT_INDEX ;

      // Cells :
      std::string CELL_POLYHEDRON ;
      size_t_array2D CELL_FACES ;
      size_t_array2D CELL_VERTICES ;
      PEL_Vector* CELL_COLORS ;  
      size_t CELL_INDEX ;

      // Faces :
      std::string FACE_POLYHEDRON ;
      size_t_array2D FACE_VERTICES ;
      PEL_Vector* FACE_COLORS ;  
      size_t FACE_INDEX ;
} ;

#endif


