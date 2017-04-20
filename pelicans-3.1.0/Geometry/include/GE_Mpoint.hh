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

#ifndef GE_M_POINT_HH
#define GE_M_POINT_HH

#include <GE_Mpolyhedron.hh>

/*

Meshes reduced to a single point geometrical representations.

It can lies in subspaces of any canonical space of dimension greater
or equal to 1.

Its reference polyhedron is `GE_ReferencePoint::object()'.

*/

class PEL_EXPORT GE_Mpoint : public GE_Mpolyhedron
{
   public: //------------------------------------------------------------------
      
   //-- Status
      
      virtual std::string const& name( void ) const ;
      
   //-- Geometrical characteristics

      virtual double measure( void ) const ;
      
      // IMPLEMENTATION: 0
      virtual GE_Point const* finite_volume_center( void ) const ;

      virtual bool contains( GE_Point const* pt ) const ;
      
      virtual GE_Vector const* unit_normal( void ) const ;
      
   //-- Reference polyhedron
            
      virtual GE_ReferencePolyhedron const* reference_polyhedron( void ) const ;
      
      virtual void apply_mapping( GE_Point const* pt_ref,
                                  GE_Point* pt ) const ;
      
      virtual void apply_inverse_mapping(
                                  GE_Point const* pt,
                                  GE_Point* pt_ref ) const ;

      virtual void build_mapping_derivative( GE_Point const* pt_ref,
                                             GE_Matrix* jac ) const ;
      
      virtual void build_tr_mapping_derivative( GE_Point const* pt_ref,
                                                GE_Matrix* tjac ) const ;
      
   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      GE_Mpoint( void ) ;
     ~GE_Mpoint( void ) ;
      GE_Mpoint( GE_Mpoint const& other ) ;
      GE_Mpoint& operator=( GE_Mpoint const& other ) ;

      GE_Mpoint( PEL_Object* a_owner, PEL_Vector* vertices ) ;
      
      virtual GE_Mpoint* create_replica( PEL_Object* a_owner, 
                                         PEL_Vector* vertices ) const ;
      
   //-- Vertices moving

      virtual void update_internal( void ) ;

   //-- Internal status

      virtual bool is_consistent( std::ostream& os=std::cout, 
                                  bool verbose=false ) const ;      

   //-- Preconditions, Postconditions, Invariant

      virtual bool finite_volume_center_POST( GE_Point const* result ) const ;
      
      virtual bool reference_polyhedron_POST(
                                GE_ReferencePolyhedron const* result ) const ;
      
      virtual bool invariant( void ) const ;

   //-- Class Attributes

      static GE_Mpoint const* MODEL ;
      static GE_Vector const* UNIT_NORMAL_MODEL ;

   //-- Attributes

} ;

#endif
