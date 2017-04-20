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

#ifndef GE_POLYGON_2D_HH
#define GE_POLYGON_2D_HH

#include <GE_Polygon.hh>

class GE_Mpolyhedron ;

class PEL_ListIdentity ;
class PEL_Vector ;

/**
   
Polygons of the 2D space.

*/

class PEL_EXPORT GE_Polygon2D : public GE_Polygon
{

   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual GE_Polygon2D* create_clone( PEL_Object* a_owner ) const ;

      // Create and return an instance.
      static GE_Polygon2D* create( PEL_Object* a_owner ) ;

      // Create and return an instance from the sequence of points `vertex_table'.
      static GE_Polygon2D* create( PEL_Object* a_owner,
                                   PEL_Sequence const* vertex_table ) ;

      // Create and return an instance from the vector of double `coordinate_table'.
      static GE_Polygon2D* create( PEL_Object* a_owner, 
				   doubleVector const& coordinate_table ) ;

      // The area of `self' and the simple polygon composing `self' has to be re-computed.
      virtual void update( void ) ;

      
   //-- Geometrical properties

      virtual double area( void ) const ;

      virtual bool has_in_interior( GE_Point const* pt ) const ;

      // Create and return a new split of simple polygons composing `self'.
      // Rem : the points potentially created at the intersections have the vector
      //       as owner.
      // Return the nil vector is the splitting has not succeeded.
      PEL_Vector* create_split_of_simple_polygons( PEL_Object* a_owner ) const ;

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      GE_Polygon2D( void ) ;
      GE_Polygon2D( GE_Polygon2D const& other ) ;
      GE_Polygon2D const& operator=( GE_Polygon2D const& other ) ;

      GE_Polygon2D( PEL_Object* a_owner ) ;
      GE_Polygon2D( PEL_Object* a_owner,
                    PEL_Sequence const* vertex_table ) ;
      GE_Polygon2D( PEL_Object* a_owner, GE_Polygon const* other ) ;
      GE_Polygon2D( PEL_Object* a_owner,
                    doubleVector const& coordinate_table ) ;
     ~GE_Polygon2D( void ) ;

   //-- Geometrical properties

      void compute_area( void ) const ;

      // split of simple polygons composing `self'
      // Rem : the points potentially created at the intersections is under
      //       ownership of `self'
      PEL_Vector const* split_of_simple_polygons( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;

   //-- Attributes
      
      // Split of simple polygons  :
      mutable PEL_Vector* SIMPLE_POLYGONS ;

      // Area :
      mutable double AREA ;      
} ;

#endif
