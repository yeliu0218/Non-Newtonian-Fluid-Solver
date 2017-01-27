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

#ifndef GE_M_POLYHEDRON_HH
#define GE_M_POLYHEDRON_HH

#include <PEL_Object.hh>
#include <iostream>
#include <string>
#include <size_t_vector.hh>

class doubleVector ;
class doubleArray2D ;
class doubleArray3D ;

class PEL_List ;
class PEL_ObjectRegister ;
class PEL_Vector ;

class GE_Matrix ;
class GE_Point ;
class GE_ReferencePolyhedron ;
class GE_SetOfPoints ;
class GE_Vector ;

/*

Meshes geometrical representation.

A GE_Mpolyhedron is a polyhedron i.e. a convex bounded intersection of half
spaces.

In addition to the polyhedron definition, its implies the definition
of an invertible mapping to a unique reference polyhedron that
belongs to an affine subspace of the canonical space E(n).

*/

class PEL_EXPORT GE_Mpolyhedron : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      // The vertices are given by the table of connectivities
      // `a_vertices_index_table' and the set of points `a_set_of_vertices'.
      static GE_Mpolyhedron* create(
                            std::string const& a_name,
                            GE_SetOfPoints* a_set_of_vertices,
                            size_t_vector const& a_vertices_index_table ) ;

      // Create and return an instance.
      // The vertices are given in the vector `vertices'.
      static GE_Mpolyhedron* create(
                            PEL_Object* a_owner,
                            std::string const& a_name,
                            PEL_Vector const* vertices ) ;

      // Update the specific attributes of `self' in case of moving of its
      // vertices.
      virtual void update( void ) ;

      static bool check_consistency( void ) ;
      static void unset_check_consistency( void ) ;
      static void set_check_consistency( void ) ;

   //-- Comparison

      virtual bool is_equal( PEL_Object const* other ) const ;

      virtual int three_way_comparison( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   //-- Status

      // characteristic name of `self'
      virtual std::string const& name( void ) const = 0 ;

   //-- Modification

      // Assuming that each vertex of `self' coincides with the vertex of
      // `other' (a fatal error is raised if not), modifie the index of
      // the vertices of `self' so that for all i, the i-th vertex of `self'
      // coincides with the i-th vertex of `other'.
      void reorder_vertices_according_to( GE_Mpolyhedron const* other ) ;

   //-- Geometrical characteristics

      // dimension of the space `self' is belongs to
      size_t nb_space_dimensions( void ) const ;

      // dimension of the subspace the reference polyhedron belongs to
      size_t dimension( void ) const ;

      // number of vertices
      size_t nb_vertices( void ) const ;

      // number of facets
      size_t nb_faces( void ) const ;

      // `i'-th vertex
      GE_Point const* vertex( size_t i ) const ;

      // measure of `self' : its length if `::dimension' is 1 ;
      // its area if `::dimension' is 2 ; its volume if `::dimension' is 3
      virtual double measure( void ) const = 0 ;

      // maximal distance between two vertices of `self'
      double inter_vertices_maximum_distance( void ) const ;

      // maximal distance between two vertices of `self' in the direction `dir'
      double inter_vertices_maximum_distance( size_t dir ) const ;

      // diameter of a ball which measure is the same than `self' ones
      double equivalent_ball_diameter( void ) const ;

      // distance between `pt1' and `pt2' in `self' characteristic length
      double reference_distance( GE_Point const* pt1,
                                 GE_Point const* pt2 ) const ;

      // geometrical center
      GE_Point const* center( void ) const ;

      // if it exists, geometrical point attached to `self' that is used in
      // finite volume discretizations, 0 otherwise
      // REMARK: the set of all cells equipped with their finite
      // volume centers must be admissible, meaning at least that the line
      // connecting the finite volume centers of two adjacent cells is
      // perpendicular to the adjacency face)
      virtual GE_Point const* finite_volume_center( void ) const = 0 ;

      // Is `pt' located in `self' ?
      virtual bool contains( GE_Point const* pt ) const = 0 ;

      // an unit normal of `self'
      // IMPLEMENTATION: a fatal error is raised, which means that the concrete
      // subclass did not implement this member function
      virtual GE_Vector const* unit_normal( void ) const ;

   //-- Reference polyhedron

      // reference polyhedron linked to the polyhedrons of name `a_name'
      // (fatal error raised if none)
      static GE_ReferencePolyhedron const* reference_polyhedron(
                                                 std::string const& a_name ) ;

      // reference polyhedron of `self'
      virtual GE_ReferencePolyhedron const* reference_polyhedron( void ) const = 0 ;

      // Projection from the reference polyhedron to `self' (mapping).
      virtual void apply_mapping( GE_Point const* pt_ref,
                                  GE_Point* pt ) const = 0 ;

      // Projection from `self' to its reference polyhedron (inverse mapping).
      virtual void apply_inverse_mapping( GE_Point const* pt,
                                          GE_Point* pt_ref ) const = 0 ;

      // Building the jacobian matrix of the mapping from the reference
      // polyhedron to `self'.
      virtual void build_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* jac ) const = 0 ;

      // Building the transpose of the jacobian matrix of the mapping from the
      // reference polyhedron to `self'.
      virtual void build_tr_mapping_derivative( GE_Point const* pt_ref,
                                                GE_Matrix* tjac ) const = 0 ;

      virtual void build_mapping_hessian( GE_Point const* pt_ref,
                                          doubleArray3D* hessian,
                                          bool& nonzero_hessian ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //---------------------------------------------------------------

      virtual ~GE_Mpolyhedron( void ) ;

      // Model registration constructor
      GE_Mpolyhedron( std::string const& a_name ) ;

      // Concrete polyhedron constructor
      GE_Mpolyhedron( PEL_Object* a_owner,
                      PEL_Vector* vertices ) ;

   //-- Vertices moving

      // Update internal specific attribute.
      virtual void update_internal( void ) = 0 ;

   //-- Internal status

      // is updating (lies in invariant calls)
      static bool is_updating( void ) ;

      // is `self' a prototype
      bool is_prototype( void ) const ;

      // Is `self' consistent
      // ( the criteria of consistence of `self' could be :
      //    - number of space dimension criteria
      //      ( in particular, `::dimension'<=`::nb_space_dimensions' )
      //    - non zero measure
      //    - vertices and facets criteria
      //      ( no cutting facets, congruent vertices, ...) ) ?
      virtual bool is_consistent( std::ostream& os=std::cout,
                                  bool verbose=false ) const = 0 ;

      void display_check( std::ostream& os,
                          std::string const& displayed_name,
                          bool success ) const ;

   //-- Utilities for the mapping

      void geobfs_2_mapping_hessian( doubleArray3D const& d2Ngeom,
                                     doubleArray3D* hessian ) const ;

   //-- Preconditions, Postconditions, Invariant

      bool finite_volume_center_PRE( void ) const ;
      virtual bool finite_volume_center_POST( GE_Point const* result ) const ;

      bool contains_PRE( GE_Point const* pt ) const ;

      bool measure_POST( double result ) const ;

      bool unit_normal_PRE( void ) const ;
      bool unit_normal_POST( GE_Vector const* result ) const ;

      bool apply_mapping_PRE( GE_Point const* pt_ref,
                              GE_Point* pt ) const ;
      bool apply_mapping_POST( GE_Point const* pt_ref,
                               GE_Point* pt ) const ;

      bool apply_inverse_mapping_PRE( GE_Point const* pt,
                                      GE_Point* pt_ref ) const ;
      bool apply_inverse_mapping_POST( GE_Point const* pt,
                                       GE_Point* pt_ref ) const ;

      bool build_mapping_derivative_PRE( GE_Point const* pt_ref,
                                         GE_Matrix* tjac ) const ;

      bool build_tr_mapping_derivative_PRE( GE_Point const* pt_ref,
                                            GE_Matrix* tjac ) const ;

      bool build_mapping_hessian_PRE( GE_Point const* pt_ref,
                                      doubleArray3D* hessian,
                                      bool& nonzero_hessian ) const ;

      bool update_internal_PRE( void ) const ;

      bool reference_polyhedron_PRE( void ) const ;
      virtual bool reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const ;

      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_Vector* vertices ) const ;
      bool create_replica_POST( PEL_Object* a_owner,
                                PEL_Vector* vertices,
                                GE_Mpolyhedron* result ) const ;

      virtual bool invariant( void  ) const ;

   private: //-----------------------------------------------------------------

      GE_Mpolyhedron( void ) ;
      GE_Mpolyhedron( GE_Mpolyhedron const& other ) ;
      GE_Mpolyhedron& operator=( GE_Mpolyhedron const& other ) ;

      virtual GE_Mpolyhedron* create_replica(
                                    PEL_Object* a_owner,
                                    PEL_Vector* vertices ) const = 0 ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Class attributes

      static bool CHECK_CONSISTENCY ;
      static bool UPDATING ; // internal used

   //-- Attributes

      // vertices vector
      PEL_Vector* VERTS ;

      // Center of `self' :
      mutable GE_Point* CENTER ;

} ;

#endif
