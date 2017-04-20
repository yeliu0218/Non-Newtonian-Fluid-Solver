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

#ifndef GE_SIMPLE_POLYGON_2D_HH
#define GE_SIMPLE_POLYGON_2D_HH

#include <GE_Polygon.hh>

class GE_Mpolyhedron ;
class GE_Polygon2D ;

class PEL_List ;
class PEL_Sequence ;
class PEL_Vector ;

/**

Polygons of the 2D space that are simple.

A polygon is simple if no two sides except for consecutive sides share
a point and if moreover two consecutive sides only share their common
endpoint. This definition excludes the so-called weakly simple polygons.

A simple polygon has a well-defined interior and it will be abusively 
associated to it, what gives a sense to methods `::area', `::has_in_interior' 
and `::create_triangulation'.

*/

class PEL_EXPORT GE_SimplePolygon2D : public GE_Polygon
{

   public: //-------------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual GE_SimplePolygon2D* create_clone( PEL_Object* a_owner ) const ;

      // Create and return an instance from the sequence of points `a_vertex_table'.
      static GE_SimplePolygon2D* create(
                        PEL_Object* a_owner,
                        PEL_Sequence const* a_vertex_table ) ;

      // Create and return an instance from the vector of double `coordinate_table'.
      static GE_SimplePolygon2D* create(
                        PEL_Object* a_owner,
                        doubleVector const& coordinate_table ) ;

      // The area of `self' has to be re-computed.
      virtual void update( void ) ;

   //-- Geometrical properties

      virtual double area( void ) const ;
      
      virtual bool has_in_interior( GE_Point const* pt ) const ;

      // orientation, result = 1 ( respectively -1 ) if `self' is counterclockwise 
      // ( respectively clockwise ) oriented 
      int orientation( void ) const ;

      // Create and return the triangulation of `self' without addition of any
      // new vertex by the ear cutting method (to avoid problems due to flat
      // angles, use the function `::suppress_flat_angles' before triangulate).
      PEL_Vector* create_triangulation( PEL_Object* seq_owner ) const ;
      
   //-- Modify the chain of vertices.

      // Suppress all vertices vi that are on the edge [vi-1,vi+1].
      void suppress_flat_angles( void ) ;
      
      // Suppress all vertices that define an angle of sine lower than 
      // `sine_of_flat_angle'.
      void suppress_flat_angles( double sine_of_flat_angle ) ;

   protected: //----------------------------------------------------------------

      

   private: //------------------------------------------------------------------

      GE_SimplePolygon2D( void ) ;
      GE_SimplePolygon2D( PEL_Object* a_owner ) ;
      GE_SimplePolygon2D( GE_SimplePolygon2D const& other ) ;
      GE_SimplePolygon2D const& operator=( GE_SimplePolygon2D const& other ) ;

      GE_SimplePolygon2D( PEL_Object* a_owner,
                          PEL_Sequence const* vertex_table ) ;
      GE_SimplePolygon2D( PEL_Object* a_owner,
                          doubleVector const& coordinate_table ) ;
      GE_SimplePolygon2D( PEL_Object* a_owner,
                          GE_SimplePolygon2D const* other ) ;
      ~GE_SimplePolygon2D( void ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool area_PRE( void ) const ;
      virtual bool has_in_interior_PRE( GE_Point const* pt ) const ;
      virtual bool invariant( void ) const ;
      
   //-- Internal status

      // IMPLEMENATION :
      //    - is_simple() ?
      virtual bool is_consistent( void ) const ;
      
      // Does `self' is simple ( two sides that share a point are compulsary  
      // consecutive sides ) ?
      bool is_simple( void ) const ;
      
      void compute_orientation( void ) const ;
      void compute_area( void ) const ;

   //-- Triangulation
      
      void set_ear_first_vertex( PEL_List* polygon ) const ;
      bool is_first_vertex_an_ear_vertex(
                                 PEL_List const* polygon ) const ;
      bool is_first_vertex_a_convex_vertex(
                                 PEL_List const* polygon ) const ;
      
   //-- Attributes

      // Area :
      mutable double SIGNED_AREA ;
} ;

#endif
