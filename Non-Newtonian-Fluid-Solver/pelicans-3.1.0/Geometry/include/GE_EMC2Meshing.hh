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

#ifndef GE_EMC2_MESHING_HH
#define GE_EMC2_MESHING_HH

#include <GE_Meshing.hh>

#include <fstream>
#include <doubleVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <stringVector.hh>

/*
Meshings stored in an EMC2 output file under one of the ftq or amdba format.

Instances are delivered by `GE_Meshing::create' whose second
argument is associated to a Module Hierarchy that has the
following keyword-data pairs :

   keyword                         data
   -------                         ----
filename                pathname of the EMC2 output file
format                  "ftq" or "ambda"

PUBLISHED
*/

class intVector ;
class PEL_ModuleExplorer ;
class PEL_KeywordDataIterator ;
class PEL_Vector ;
class GE_BoxWithBoxes ;

class PEL_EXPORT GE_EMC2Meshing : public GE_Meshing
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

      GE_EMC2Meshing( void ) ;
     ~GE_EMC2Meshing( void ) ;
      GE_EMC2Meshing( GE_EMC2Meshing const& other ) ;
      GE_EMC2Meshing const& operator=( GE_EMC2Meshing const& other ) ;

      GE_EMC2Meshing( PEL_Object* a_owner,
                      PEL_ModuleExplorer const* exp ) ;

      virtual GE_EMC2Meshing* create_replica( PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp,
                                              size_t dim_space ) const ;
      
      virtual GE_Color const* default_vertex_color( void ) const ;

      virtual GE_Color const* default_cell_color( void ) const ;

      virtual GE_Color const* default_face_color( void ) const ;

      void read_vertex( void ) ;

      void read_cell( void ) ;

      void read_face( void ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_EMC2Meshing const* PROTOTYPE ;

   //-- Attributes

      // Input file:
      std::ifstream INPUT_FILE ;
      enum  EMC2OutputFormat { ftq, amdba } ;
      EMC2OutputFormat const INPUT_FORMAT ;
      std::streampos VERT_POS ;
      std::streampos CELL_POS ;

      // Meshing:
      size_t NB_VERTS ;
      size_t NB_FACES ;
      size_t NB_TRIAS ;
      size_t NB_QUADS ;
      stringVector FACE_POLY_NAMES ;
      stringVector CELL_POLY_NAMES ;
      size_t_vector VERT_IDX ;
      size_t_array2D CELLS_2_FACES ;
      size_t_array2D FACES_2_VERTS ;

      // Iterator:
      size_t I_VERT ;
      size_t I_FACE ;
      size_t I_CELL ;
      GE_Color const* DEF_COLOR ;
      doubleVector VERT_COORDS ;
      size_t_vector MESH_2_VERTS ;
      size_t_vector CELL_2_FACES ;
} ;

#endif


