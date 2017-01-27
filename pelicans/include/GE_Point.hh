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

#ifndef GE_POINT_HH
#define GE_POINT_HH

#include <PEL_Object.hh>

#include <PEL_assertions.hh>
#include <doubleVector.hh>

class PEL_DoubleComparator ;
class GE_Vector ;

/*
Geometrical points.
*/

class PEL_EXPORT GE_Point : public PEL_Object
{

   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create a point from its dimension.
      static GE_Point* create( PEL_Object* a_owner,
                               size_t a_dimension ) ;

      // Create a point from its table of coordinates.
      static GE_Point* create( PEL_Object* a_owner,
                               doubleVector const& a_coordinates_table ) ;

      // Create a 1D point from its coordinate.
      static GE_Point* create( PEL_Object* a_owner,
                               double x0 ) ;

      // Create a 2D point from its coordinates.
      static GE_Point* create( PEL_Object* a_owner,
                               double x0, double x1 ) ;

      // Create a 3D point from its coordinates.
      static GE_Point* create( PEL_Object* a_owner, 
                               double x0, double x1, double x2 ) ;

      // origin point
      static GE_Point const* origin( size_t dimension ) ;

      virtual GE_Point* create_clone( PEL_Object* a_owner ) const ;

   //-- Comparison

      static PEL_DoubleComparator const* coordinates_comparator( void ) ;
      
      // IMPLEMENATION: is three_way_comparison(other)==0 ?
      virtual bool is_equal( PEL_Object const* other ) const ;

      /*
      IMPLEMENTATION:
         use `::coordinates_comparator' in order to compare coordinates
         in lexicographical order
      */
      virtual int three_way_comparison( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   //-- Status

      // number of coordinates
      size_t nb_coordinates( void ) const ;

      // `ic'-th coordinate
      double coordinate( size_t ic ) const ;

      // vector of coordinate
      doubleVector const& coordinate_vector( void ) const ;

      // distance between `self' and `other'
      double distance( GE_Point const* other ) const ;

   //-- Modification

      // Reinitialize `self' with copy of `pt'.
      void copy( GE_Point const* pt ) ;
      
      // Set `self' with the coordinates of `pt'.
      void set( GE_Point const* pt ) ;

      // Set `self' with a table of coordinates.
      void set_coordinates( doubleVector const& a_coordinates_table ) ;

      // Set the `ic'-th coordinate of `self'.
      void set_coordinate( size_t ic, double x ) ;

      // Set `self' as the barycenter of (`ptA',1-`alpha') and 
      // (`ptB',`alpha').
      void set_as_barycenter( double alpha, 
                              GE_Point const* ptA, 
                              GE_Point const* ptB ) ;

      // Translate `self' of `vec'.
      void translate( double alpha, GE_Vector const* vec ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------- 

   private: //----------------------------------------------------------

      GE_Point( void ) ;
     ~GE_Point( void ) ;
      GE_Point( GE_Point const& other ) ;
      GE_Point const& operator=( GE_Point const& other ) ;

      GE_Point( PEL_Object* a_owner, size_t a_dimension ) ;
      GE_Point( PEL_Object* a_owner, doubleVector const& a_coordinates_table );
      GE_Point( PEL_Object* a_owner, double x0 ) ;
      GE_Point( PEL_Object* a_owner, double x0, double x1 ) ;
      GE_Point( PEL_Object* a_owner, double x0, double x1, double x2 ) ;
      GE_Point( PEL_Object* a_owner, GE_Point const* other ) ;
      

   //-- Overloaded to avoid implicit type conversion

      static GE_Point* create( PEL_Object* a_owner, int dimension ) ;
      GE_Point( PEL_Object* a_owner, int dimension ) ;
      
   //-- Attributes

      doubleVector xx ;
} ;

#ifndef OUTLINE
#include <GE_Point.icc>
#endif

#endif
