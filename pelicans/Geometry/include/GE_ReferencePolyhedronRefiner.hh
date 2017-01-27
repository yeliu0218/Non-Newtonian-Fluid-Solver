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

#ifndef GE_REFERENCE_POLYHEDRON_REFINER_HH
#define GE_REFERENCE_POLYHEDRON_REFINER_HH

#include <PEL_Object.hh>

#include <size_t_array2D.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;
class PEL_Vector ;

class GE_Point ;
class GE_ReferencePolyhedron ;

#include <map>
#include <set>

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT GE_ReferencePolyhedronRefiner : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_ReferencePolyhedronRefiner const* make(
                                             PEL_Object* a_owner,
                                             PEL_ModuleExplorer const* exp ) ;
      
   //-- Coarse polyhedron
      
      GE_ReferencePolyhedron const* reference_polyhedron( void ) const ;

   //-- Vertices of the fine polyhedra

      size_t nb_vertices( void ) const ;

      GE_Point const* vertex( size_t i ) const ;

   //-- Faces of the fine polyhedra

      GE_ReferencePolyhedron const* subface_reference_polyhedron( void ) const ;
      
      size_t nb_subfaces_per_face( void ) const ;

      size_t nb_subfaces( void ) const ;

      size_t subface_vertex( size_t is, size_t iv ) const ;
      
      size_t subface_parent( size_t is ) const ;

   //-- Fine polyhedra

      GE_ReferencePolyhedron const* subcell_reference_polyhedron( void ) const ;

      size_t nb_subcells( void ) const ;

      size_t subcell_vertex( size_t ic, size_t iv ) const ;

      size_t subcell_face( size_t ic, size_t is ) const ;

      virtual void compute_location_in_subcell( 
                                       size_t ic, 
                                       GE_Point const* pt_cell,
                                       GE_Point* pt_subcell ) const = 0 ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~GE_ReferencePolyhedronRefiner( void ) ;

      GE_ReferencePolyhedronRefiner( std::string const& a_name ) ;
      
      GE_ReferencePolyhedronRefiner( 
                            PEL_Object* a_owner,
                            GE_ReferencePolyhedron const* ref_poly,
                            GE_ReferencePolyhedron const* subcell_ref_poly,
                            GE_ReferencePolyhedron const* subface_ref_poly ) ;

      virtual GE_ReferencePolyhedronRefiner const* create_replica(
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp) const = 0 ;
      
      bool is_a_prototype( void ) const ;
      
   //-- Internals building
      
      void set_dimensions( size_t a_nb_vertices,
                           size_t a_nb_subfaces,
                           size_t a_nb_subcells,
                           size_t a_nb_subfaces_per_face ) ;

      void set_vertex( size_t i_vertex, GE_Point* pt ) ;

      void set_subcell( size_t i_cell, size_t_vector const& vertex_indices ) ;

   //-- Preconditions, Postconditions, Invariant

      bool compute_location_in_subcell_PRE( size_t ic, 
                                            GE_Point const* pt_cell,
                                            GE_Point* pt_subcell ) const ;
      
      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      
      bool create_replica_POST( GE_ReferencePolyhedronRefiner const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;

   private: //----------------------------------------------------------

      GE_ReferencePolyhedronRefiner( void ) ;
      GE_ReferencePolyhedronRefiner( 
                            GE_ReferencePolyhedronRefiner const& other ) ;
      GE_ReferencePolyhedronRefiner& operator=( 
                            GE_ReferencePolyhedronRefiner const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
      
   //-- Attributes

      bool const IS_PROTO ;
      GE_ReferencePolyhedron const* RPOLY ;
      GE_ReferencePolyhedron const* SUBCELL_RPOLY ;
      GE_ReferencePolyhedron const* SUBFACE_RPOLY ;
      size_t NB_VERTS ;
      size_t NB_FACES ;
      size_t NB_CELLS ;
      size_t NB_SUB_PER_FACE ; 
      PEL_Vector* VERTICES ;
      size_t_array2D FACE_VERT ;
      size_t_vector FACE_PARENT ;
      size_t_array2D CELL_VERT ;
      size_t_array2D CELL_FACE ;

      std::map< std::set<size_t>, size_t > FACES ;
      std::map< std::set<size_t>, size_t >::iterator FACES_it ;
      size_t i_SUBFACE ;
} ;

#endif
