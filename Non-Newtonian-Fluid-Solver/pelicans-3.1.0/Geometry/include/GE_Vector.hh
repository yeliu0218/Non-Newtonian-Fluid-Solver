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

#ifndef GE_VECTOR_HH
#define GE_VECTOR_HH

#include <PEL_Object.hh>

#include <doubleVector.hh>

/*
Geometrical vectors.
*/

class GE_Point ;

class PEL_EXPORT GE_Vector : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create a vector from its dimension.
      static GE_Vector* create( PEL_Object* a_owner,
                                size_t a_dimension ) ;

      // Create a vector from its table of pomponents.
      static GE_Vector* create( PEL_Object* a_owner, 
                                doubleVector const& vector ) ;

      // Create a vector form `start' to 'end'.
      static GE_Vector* create( PEL_Object* a_owner, 
                                GE_Point const* end,
                                GE_Point const* start ) ;

      virtual GE_Vector* create_clone( PEL_Object* a_owner ) const ;

      // Reinitialize `self' to be a vector form `start' to 'end'.
      void re_initialize( GE_Point const* end, GE_Point const* start) ;

   //-- Comparison

      virtual bool is_equal( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   //-- Status

      // number of components
      size_t nb_components( void ) const ;

      // `ic'-th component
      double component( size_t ic ) const ;
      
      // vector of component
      doubleVector const& component_vector( void ) const ;
      
      // euclidian norm
      double norm( void ) const ;

      // euclidian dot between `self' and `vec'
      double dot_product( GE_Vector const* vec ) const ;

      // cosine between `self' and `vec'
      double cosine( GE_Vector const* vec ) const ;

   //-- Modification

      // Reinitialize `self' with copy of `vec'.
      void copy( GE_Vector const* vec ) ;

      // Set `self' with the coordinates of `vec'.
      void set( GE_Vector const* vec ) ;

      // Set the `ic'-th component of `self'.
      void set_component( size_t ic, double x ) ;

      // Set `self' as the cross product of `vec1' and `vec2'.
      void set_as_cross_product( GE_Vector const* vec1, 
                                 GE_Vector const* vec2 ) ;

      // Scale `self'.
      void scale( double alpha ) ;

      // Add `alpha'*`vec' to `self'.
      void sum( double alpha, GE_Vector const* vec ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------- 
   
   private: //----------------------------------------------------------

      GE_Vector( void ) ;
     ~GE_Vector( void ) ;
      GE_Vector( GE_Vector const& other ) ;
      GE_Vector const& operator=( GE_Vector const& other ) ;

      GE_Vector( PEL_Object* a_owner, GE_Vector const* other ) ;
      GE_Vector( PEL_Object* a_owner, size_t a_dimension ) ;
      GE_Vector( PEL_Object* a_owner,
                 GE_Point const* end, GE_Point const* start ) ;
      GE_Vector( PEL_Object* a_owner, doubleVector const& vector ) ;


   //-- Attributes

      doubleVector xx ;

} ;

#endif

