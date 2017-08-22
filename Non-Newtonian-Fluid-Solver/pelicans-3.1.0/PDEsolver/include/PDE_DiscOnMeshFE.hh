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

#ifndef PDE_DISC_ON_MESH_FE_HH
#define PDE_DISC_ON_MESH_FE_HH

#include <PEL_Object.hh>

#include <vector>

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;
class GE_ReferencePolyhedron ;

class PDE_BasisFunction ;
class PDE_DiscreteField ;
class PDE_ReferenceElement ;
class PDE_RefinementPatternProvider ;

class PEL_EXPORT PDE_DiscOnMeshFE : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PDE_DiscOnMeshFE* create( PEL_Object* a_owner,
                                       PDE_DiscreteField const* ff,
                                       PDE_ReferenceElement const* elm ) ;
      
   //-- Discrete fields

      void add_discretization( PDE_DiscreteField const* ff,
                               PDE_ReferenceElement const* elm ) ;

      void duplicate_discretization( PDE_DiscreteField const* model_f,
                                     PDE_DiscreteField const* new_f ) ;

      bool has_discretization( PDE_DiscreteField const* ff ) const ;

      size_t index_of_reference_element( PDE_DiscreteField const* ff ) const ;

      size_t nb_discrete_fields( size_t ee ) const ;

      PDE_DiscreteField const* discrete_field( size_t ee, size_t i ) const ;

   //-- Reference Elements

      size_t index_of_element( PDE_ReferenceElement const* elm ) const ;

      size_t nb_reference_elements( void ) const ;

      PDE_ReferenceElement const* reference_element( size_t ee ) const ;

   //-- Basis functions
      
      size_t nb_basis_functions( size_t ee ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DiscOnMeshFE( void ) ;
     ~PDE_DiscOnMeshFE( void ) ; 
      PDE_DiscOnMeshFE( PDE_DiscOnMeshFE const& other ) ;
      PDE_DiscOnMeshFE& operator=( PDE_DiscOnMeshFE const& other ) ;

      PDE_DiscOnMeshFE( PEL_Object* a_owner ) ;
      
      PDE_DiscOnMeshFE( PEL_Object* a_owner,
                        PDE_DiscreteField const* ff,
                        PDE_ReferenceElement const* elm ) ;

      //-- Attributes
      
      std::vector< PDE_ReferenceElement const* > ELMS ;
      std::vector< size_t > FIELD_2_ELMS ;
      std::vector< std::vector< PDE_DiscreteField const* > > DFS ;         
} ;

#endif
