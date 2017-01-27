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

#ifndef GE_TRIANGLE_HH
#define GE_TRIANGLE_HH

#include <GE_Mpolyhedron.hh>

class GE_Vector ;

/*

Triangular meshes geometrical representations.

It can lies in subspaces of the two-dimensional canonical space or of the 
three-dimensional canonical space.

Its reference polyhedron is `GE_ReferenceTriangle::object()'.

*/

class PEL_EXPORT GE_Triangle : public GE_Mpolyhedron
{
   public: //----------------------------------------------------------

   //-- Status

      virtual std::string const& name( void ) const ;

   //-- Geometrical characteristics

      static double measure( GE_Point const* v0,
                             GE_Point const* v1,
                             GE_Point const* v2 ) ;
      
      virtual double measure( void ) const ;

      // IMPLEMENTATION: the circumcenter
      virtual GE_Point const* finite_volume_center( void ) const ;

      virtual bool contains( GE_Point const* pt ) const ;

      virtual GE_Vector const* unit_normal( void ) const ;
      
   //-- Reference polyhedron
      
      virtual GE_ReferencePolyhedron const* reference_polyhedron( void ) const ;
      
      virtual void apply_mapping( GE_Point const* pt_ref,
                                  GE_Point* pt ) const ;

      static void apply_inverse_mapping(
                                  GE_Point const* v0,
                                  GE_Point const* v1,
                                  GE_Point const* v2,
                                  double v0v1v2_measure,
                                  GE_Point const* pt,
                                  GE_Point* pt_ref ) ;
                             
      virtual void apply_inverse_mapping(
                                  GE_Point const* pt,
                                  GE_Point* pt_ref ) const ;
      
      virtual void build_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* jac ) const ;
      
      virtual void build_tr_mapping_derivative( GE_Point const* pt_ref,
                                                GE_Matrix* tjac ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      GE_Triangle( void ) ;
     ~GE_Triangle( void ) ;
      GE_Triangle( GE_Triangle const& other ) ;
      GE_Triangle& operator=( GE_Triangle const& other ) ;

      GE_Triangle( PEL_Object* a_owner, 
                   PEL_Vector* vertices ) ;

      virtual GE_Triangle* create_replica( PEL_Object* a_owner, 
                                           PEL_Vector* vertices ) const ;

   //-- Vertices moving

      virtual void update_internal( void ) ;

   //-- Internal status

      virtual bool is_consistent( std::ostream& os=std::cout, 
                                  bool verbose=false ) const ;

      static double signed_measure( GE_Point const* v0,
                                    GE_Point const* v1,
                                    GE_Point const* v2 ) ;
      
      static void compute_finite_volume_center( GE_Point const* a,
                                                GE_Vector const* ab,
                                                GE_Vector const* ac,
                                                GE_Point* result ) ;
 
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool finite_volume_center_POST( GE_Point const* result ) const ;
      
      virtual bool reference_polyhedron_POST(
                                GE_ReferencePolyhedron const* result ) const ;
      
      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_Triangle const* MODEL ;
      static GE_Vector* V0V1 ;
      static GE_Vector* V0V2 ;
      static GE_Vector* N ;

   //-- Attributes
      
      GE_Vector* UNIT_NORMAL ;
      double MEASURE ;
      bool ORIENT ;
      mutable GE_Point* FV_CENTER ;
} ;

#endif

