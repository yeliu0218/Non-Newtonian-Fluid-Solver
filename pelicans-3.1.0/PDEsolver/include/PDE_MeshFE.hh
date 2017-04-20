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

#ifndef PDE_MESH_FE_HH
#define PDE_MESH_FE_HH

#include <PEL_Object.hh>

#include <vector>

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;

class PDE_BasisFunction ;
class PDE_DiscreteField ;
class PDE_DiscOnMeshFE ;
class PDE_ReferenceElement ;
class PDE_RefinementPatternProvider ;

class PEL_EXPORT PDE_MeshFE : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Comparison

      virtual bool is_equal( PEL_Object const* other ) const ;

      virtual int three_way_comparison( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   //-- Basic characteristics

      size_t id_number( void ) const ;

      GE_Mpolyhedron* polyhedron( void ) const ;

      GE_Color const* color( void ) const ;
      
   //-- Discrete fields

      void duplicate_discretization( PDE_DiscreteField const* model_f,
                                     PDE_DiscreteField* new_f ) ;

      bool has_discretization( PDE_DiscreteField const* ff ) const ;

      size_t index_of_reference_element( PDE_DiscreteField const* ff ) const ;

      virtual double value( PDE_DiscreteField const* ff,
                            size_t level,
                            GE_Point const* pt_ref,
                            size_t ic = 0 ) const ;

      size_t nb_discrete_fields( size_t ee ) const ;

      PDE_DiscreteField const* discrete_field( size_t ee, size_t i ) const ;

   //-- Reference Elements

      size_t index_of_element( PDE_ReferenceElement const* elm ) const ;

      size_t nb_reference_elements( void ) const ;

      PDE_ReferenceElement const* reference_element( size_t ee ) const ;

   //-- Basis functions
      
      size_t nb_basis_functions( size_t ee ) const ;
      
      // `ln'-th basis function of the `ee'-th reference element (may be 0)
      virtual PDE_BasisFunction* basis_function( size_t ee, 
                                                 size_t ln ) const = 0 ;

   //-- Adaptation

      static void set_refinement_pattern_provider( 
                                 PDE_RefinementPatternProvider const* r ) ;

      static PDE_RefinementPatternProvider const* 
                                  refinement_pattern_provider( void ) ;

      size_t refinement_level( void ) const ;

      virtual bool is_active( void ) const = 0 ;

      virtual PDE_MeshFE* parent( void ) const = 0 ;

   //-- Element change

      void set_color( GE_Color const* col ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_MeshFE( void ) ; 
     
      PDE_MeshFE( PEL_Object* a_owner,
                  size_t a_number,
                  GE_Mpolyhedron* a_polyhedron,
                  GE_Color const* a_color,
                  size_t a_refinement_level ) ;

      PDE_MeshFE( PEL_Object* a_owner,
                  size_t a_number,
                  GE_Mpolyhedron* a_polyhedron,
                  GE_Color const* a_color,
                  size_t a_refinement_level,
                  PDE_MeshFE const* a_pattern_mesh ) ;
      
      virtual PDE_DiscOnMeshFE const* disc( void ) const = 0 ;

      bool value_PRE( PDE_DiscreteField const* ff,
                      size_t level,
                      GE_Point const* pt_ref,
                      size_t ic ) const ;

      bool basis_function_PRE( size_t ee, size_t ln ) const ;
      
      bool basis_function_POST( PDE_BasisFunction* result,
                                size_t ee, size_t ln ) const ;

   private: //----------------------------------------------------------

      PDE_MeshFE( void ) ;
      PDE_MeshFE( PDE_MeshFE const& other ) ;
      PDE_MeshFE& operator=( PDE_MeshFE const& other ) ;

   //-- Class attributes

      static PDE_RefinementPatternProvider const* RPP ;

   //-- Attributes
      
      size_t const ID ;
      GE_Mpolyhedron* const POLY ;
      GE_Color const* COLOR ;
      size_t const RLEVEL ;
} ;

#endif
