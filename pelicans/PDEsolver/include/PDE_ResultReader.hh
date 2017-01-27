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

#ifndef PDE_RESULT_READER_HH
#define PDE_RESULT_READER_HH

#include <PEL_Object.hh>
#include <boolVector.hh>
#include <doubleVector.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

class PEL_DataOnMeshingReader ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class doubleArray2D ;
class intArray2D ;

class GE_Mpolyhedron ;
class GE_Point ;

class PDE_ReferenceElement ;


class PEL_EXPORT PDE_ResultReader : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_ResultReader* create( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;

   //-- Geometry

      size_t nb_space_dimensions( void ) const ;

      bool is_in_grid( GE_Point const* pt ) const ;

      bool has_integration_domain( void ) const ;

      doubleArray2D const& inner_boundary( void ) const ;

   //-- Fields

      bool has_field( std::string const& field_name ) const ;

      doubleVector const&  field_value( std::string const& field_name,
                                        GE_Point const* pt ) const ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_ResultReader( void ) ;
     ~PDE_ResultReader( void ) ;
      PDE_ResultReader( PDE_ResultReader const& other ) ;
      PDE_ResultReader& operator=( PDE_ResultReader const& other ) ;

      PDE_ResultReader( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp,
                        PEL_DataOnMeshingReader const* data_reader ) ;

      size_t closest_vertex( GE_Point const* pt, bool& is_a_vertex ) const ;
      size_t closest_vertex_global_search( GE_Point const* pt, 
                                           bool& is_a_vertex ) const ;

      size_t cell_index( GE_Point const* pt, size_t i_closest_vertex,
                         GE_Mpolyhedron* cell_poly,
                         PEL_Vector* cell_poly_verts ) const ;

      void set_value_at_vertices(
                         GE_Point const* pt,
                         size_t i_cell, GE_Mpolyhedron const* cell_poly,
                         size_t i_comp, doubleArray2D const& value_table,
                         doubleVector& val ) const ;

      void set_value_at_cell_centers(
                         GE_Point const* pt,
                         size_t i_vert,
                         size_t i_comp,
                         GE_Mpolyhedron* dummy_cell_poly,
                         PEL_Vector*  dummy_cell_poly_verts,
                         doubleArray2D const& value_table,
                         doubleVector& val ) const ;

      double distance( GE_Point const* pt1, doubleVector const& pt2 ) const ;

      size_t index_of_neighbour( GE_Point const* pt, size_t const current_vert,
                                 double& d_min ) const ;
      
   //-- Attributes

      bool const DEBUG ;

      // Reader :
      PEL_DataOnMeshingReader const* const READER ;

      // Vertice table :
      size_t const NB_SP_DIMS ;
      size_t const NB_VERTS ;
      doubleArray2D const& VERTS ;

      // Cell connectivity table :
      size_t const NB_CELLS ;
      intArray2D const& CELL2VERT ;

      // Field names :
      stringVector FIELDS ;
      stringVector FNAMES ;
      size_t_vector COMP ; // Composed field

      // Cell search :
      size_t_array2D VERT2CELL ;
      size_t_vector  NB_CELL_PER_VERT ;
      double EPS_VERT ;
      mutable size_t CLOSEST_VERT ; // To speed up the algorithm :
                                    // we keep the last closest vertex found
      mutable boolVector VISITED_VERT ;
      mutable doubleVector VERT ;
      mutable GE_Mpolyhedron* POLY ;
      PEL_Vector* const POLY_VERTS ;
      PDE_ReferenceElement const* REF_ELM ;
      GE_Point* PT ;
      GE_Point* PT_REF ;
} ;

#endif
