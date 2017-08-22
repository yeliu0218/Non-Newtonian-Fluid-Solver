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

#ifndef GE_POLYGON_HH
#define GE_POLYGON_HH

#include <PEL_Object.hh>

class GE_Point ;
class GE_PointIterator ;
class GE_PointPoint_INT ;
class GE_PointSegment_INT ;
class GE_SegmentSegment_INT ;

class PEL_List ;
class PEL_ListIterator ;
class PEL_Sequence ;

class doubleVector ;

/*
  
Closed curves of a plane composed of a finite number of straight 
line segments.

The space dimension can either be of 2 or of 3. 

Those segments are called the polygon sides and the endpoints where two
consecutive segments meet are called  the polygon vertices.
The number of vertices gives the polygon size.
The chain of vertices being circular each vertex is associated to a previous
one and to a next one.

*/

class PEL_EXPORT GE_Polygon : public PEL_Object
{

   public: //------------------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual void update( void ) = 0 ;

   //-- General properties

      // size of the polygon ( number of vertices )
      size_t size( void ) const ;

      // dimension of the space `self' belongs to
      size_t dimension( void ) const ;
      
   //-- Access to the chain of vertices

      // Create a chain of vertices iterator.
      GE_PointIterator* create_vertex_iterator( PEL_Object* owner ) const ;
      
      // location of the `i'-th vertex
      GE_Point const* vertex( size_t i ) const ;
      
      // location of the vertex that follows the `i'-th vertex
      GE_Point const* next_vertex( size_t i ) const ;

      // location of the vertex that precedes the `i'-th vertex
      GE_Point const* previous_vertex( size_t i ) const ;

      // index of the vertex that follows the `i'-th vertex
      size_t next_vertex_index( size_t i ) const ;

      // index of the vertex that precedes the `i'-th vertex
      size_t previous_vertex_index( size_t i ) const ;

   //-- Modify the chain of vertices.

      // Insert point `pt' has a new vertex of index `i' in the chain.
      // The function `::update' is called.
      void insert_vertex_at( size_t i, GE_Point const* pt ) ;

      // Insert point `pt' at the end of the chain of vertices.
      // The function `::update' is called.
      void append_vertex( GE_Point const* pt ) ;

      // Remove the vertex of index `i' in the chain (do NOT destroy it).
      // The function `::update' is called.
      void remove_vertex_at( size_t i ) ;

      // Replace the vertex of index `i' in the chain with `pt'.
      // The function `::update' is called.
      void replace_vertex( size_t i, GE_Point const* pt ) ;

      // Remove all the vertices in the chain (do NOT destroy them)
      // The function `::update' is called.
      void remove_vertices( void ) ;

   //-- Geometrical properties

      // area of `self'
      virtual double area( void ) const = 0 ;

      // Does `self' contains on its boundary the `GE_Point::' `pt' ?
      bool has_on_boundary( GE_Point const* pt ) const ;

      // Does `self' contains in its closed interior the `GE_Point::' `pt' ?
      virtual bool has_in_interior( GE_Point const* pt ) const = 0 ;

      // Project `pt' on the boundary of `self'.
      virtual void project_on_boundary( GE_Point* pt ) const ;

   //-- Intersector

      void set_point_point_intersector( GE_PointPoint_INT* a_pt_pt_intersector ) ;
      bool is_set_point_point_intersector( void ) const ;
      GE_PointPoint_INT* point_point_intersector( void ) const ;

      void set_point_segment_intersector( GE_PointSegment_INT* a_pt_seg_intersector ) ;
      bool is_set_point_segment_intersector( void ) const ;
      GE_PointSegment_INT* point_segment_intersector( void ) const ;

      void set_segment_segment_intersector( GE_SegmentSegment_INT* a_seg_seg_intersector ) ;
      bool is_set_segment_segment_intersector( void ) const ;
      GE_SegmentSegment_INT* segment_segment_intersector( void ) const ;

   //-- Comparison

      // Does `self' closed chain of vertices is equal to the one of `other'?
      virtual bool is_equal( PEL_Object const* other ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //----------------------------------------------------------------

      virtual ~GE_Polygon( void ) ;

      GE_Polygon( PEL_Object* a_owner,
                  size_t dim ) ;
      GE_Polygon( PEL_Object* a_owner,
                  size_t dim,
                  PEL_Sequence const* vertex_table ) ;
      GE_Polygon( PEL_Object* a_owner,
                  GE_Polygon const* other ) ;
      GE_Polygon( PEL_Object* a_owner,
                  size_t dim,
                  doubleVector const& coordinate_table ) ;
      
   //-- Internal status

      // Is `self' consistent
      // ( the criteria of consistence of `self' could be :
      //    - non zero measure 
      //    - vertices and facets criteria
      //      ( no cutting facets, congruent vertices, ...) ) ?
      // IMPLEMENATION :
      //    - if size()>=2 : for all i : !vertex(i)->is_close( next_vertex(i) )
      virtual bool is_consistent( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      virtual bool area_PRE( void ) const ;
      virtual bool has_in_interior_PRE( GE_Point const* pt ) const ;
      virtual bool project_on_boundary_PRE( GE_Point const* pt ) const ;
      virtual bool project_on_boundary_POST( GE_Point const* pt ) const ;
      
   private: //------------------------------------------------------------------

      GE_Polygon( void ) ;
      GE_Polygon( GE_Polygon const& other ) ;
      GE_Polygon const& operator=( GE_Polygon const& other) ;

   //-- Attributes

      PEL_List* const LIST_OF_POINTS ;
      size_t const DIM ;

      mutable PEL_ListIterator* POINTS_IT ;
      mutable bool OK_IT ;

      GE_PointPoint_INT* PT_PT_INTERSECTOR ;
      GE_PointSegment_INT* PT_SEG_INTERSECTOR ;
      GE_SegmentSegment_INT* SEG_SEG_INTERSECTOR ;
} ;

#endif
