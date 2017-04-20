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

#ifndef GE_REFERENCE_POLYHEDRON_HH
#define GE_REFERENCE_POLYHEDRON_HH

#include <PEL_Object.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>

class doubleVector ;
class doubleArray2D ;
class doubleArray3D ;

class PEL_Vector ; 

class GE_Point ;
class GE_Vector ;

/*
Reference polyhedra.
*/

class PEL_EXPORT GE_ReferencePolyhedron : public PEL_Object
{
   public: //------------------------------------------------------------------
      
  //-- Status

      // characteristic name of `self'
      std::string const& name( void ) const ;

      // number, uniquely determining `self'
      size_t id_number( void ) const ;
      
      static size_t nb_objects( void ) ;

  //-- Geometrical properties

      // dimension
      size_t dimension( void ) const ;
      
      // number of vertices
      size_t nb_vertices( void ) const ;

      // `i'-th vertex
      GE_Point const* vertex( size_t i ) const ;

      // measure of `self' (its length if `::dimension' is 1 ; 
      // its area if `::dimension' is 2 ; its volume if `::dimension' is 3)
      double measure( void ) const ;

      // number of facets 
      size_t nb_faces( void ) const ;

      size_t nb_face_vertices( size_t i_face ) const ;

      size_t face_vertex( size_t i_face, size_t i_vert ) const ;
      
      GE_Vector const* face_outward_normal( size_t i_face ) const ;

      // geometrical center
      GE_Point const* center( void ) const ;

      // Does `self' contain `pt_ref' ?
      //   `tol'>=0.: Is `pt_ref' in the interior of `self' or located at
      //              a distance lower that `tol' from the  faces of `self' ?
      //   `tol'<0.:  Is `pt_ref' in the interior of `self' AND NOT located at
      //              a distance lower that -`tol' from the  faces of `self' ?
      virtual bool contains( GE_Point const* pt_ref,
                             double tol = epsilon() ) const = 0 ;

      // Is `pt_ref' located at a distance lower that `tol' from the 
      // `i_face'-th faces of `self' ?
      virtual bool face_contains( size_t i_face, 
                                  GE_Point const* pt_ref,
                                  double tol = epsilon() ) const = 0 ;

      // Project `pt_ref' on `self'.
      virtual void project( GE_Point* pt_ref ) const = 0 ;

      // Set `neighbor' as the point of the reference poyhedron that is at 
      // distance `delta' of `pt_ref' in direction `ic'.
      virtual void build_neighbor( GE_Point const* pt_ref,
                                   size_t ic,
                                   GE_Point* neighbor,
                                   double& delta ) const = 0 ;
      
      // tolerence allowed for the constaining of a point in `self'
      // IMPLEMENTATION: 1.E-5
      static double epsilon( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------------

      virtual ~GE_ReferencePolyhedron( void ) ;

      GE_ReferencePolyhedron( std::string const& a_name,
                              size_t const a_nb_vertices,
                              size_t const a_nb_faces,
                              size_t const a_dimension,
                              double const a_measure ) ;

      void set_vertex( size_t i, GE_Point* pt ) ;

      void append_face_vertex( size_t i_face, size_t vertex_id ) ;
      
      void set_face_normal( size_t i_face, GE_Vector* n ) ;
      
   //-- Preconditions, Postconditions, Invariant
      
      bool project_PRE( GE_Point const* pt_ref ) const ;
      bool project_POST( GE_Point const* pt_ref ) const ;

      bool contains_PRE( GE_Point const* pt_ref,
                         double tol ) const ;

      bool face_contains_PRE( size_t i_face,
                              GE_Point const* pt_ref,
                              double tol ) const ;

      bool build_neighbor_PRE( GE_Point const* pt_ref,
                               size_t ic,
                               GE_Point const* neighbor ) const ;
      bool build_neighbor_POST( GE_Point const* pt_ref,
                                size_t ic,
                                GE_Point const* neighbor,
                                double delta ) const ;

      virtual bool invariant( void ) const ;
      
   private: //-----------------------------------------------------------------

      GE_ReferencePolyhedron( void ) ;
      GE_ReferencePolyhedron( GE_ReferencePolyhedron const& other ) ;
      GE_ReferencePolyhedron& operator=( 
                              GE_ReferencePolyhedron const& other ) ;

   //-- Class attributes

      static size_t NB_INSTANCES ;

   //-- Attributes
      
      std::string const NAME ;
      size_t ID ;
      size_t const POLY_DIM ;

      size_t const NB_VERTS ;
      PEL_Vector* const VERTS ;

      size_t NB_FACES ;
      size_t_vector NB_FVERTS ;
      size_t_array2D IDX_FVERTS ;
      PEL_Vector* const F_NORMALS ;

      mutable GE_Point* CENTER ;
      double const MEASURE ;
} ;

#endif
