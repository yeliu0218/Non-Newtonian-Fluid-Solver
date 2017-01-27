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

#ifndef PDE_BOUND_FE_HH
#define PDE_BOUND_FE_HH

#include <PDE_MeshFE.hh>

#include <PDE_BasisFunctionBound.hh>

#include <string>

class PDE_BoundFE ;
class PDE_DiscreteField ;
class PDE_CellFE ;

class GE_Point ;
class GE_Vector ;

class PEL_List ;
class size_t_vector ;

class PEL_EXPORT PDE_BoundFE : public PDE_MeshFE
{
   public: //-----------------------------------------------------------

   //-- Instance creation and initialization

      static PDE_BoundFE* create( PEL_Object* a_owner,
                                  size_t a_number,
                                  GE_Mpolyhedron* a_polyhedron,
                                  GE_Color const* a_color,
                                  size_t a_refinement_level,
                                  PDE_DiscOnMeshFE const* a_disc ) ;

   //-- Basis functions

      void set_basis_function( size_t ee, 
                               size_t ln, 
                               PDE_BasisFunctionBound* bf ) ;

      virtual PDE_BasisFunctionBound* basis_function( size_t ee, 
                                                      size_t ln ) const ;

    //-- Discrete field reconstruction

      virtual double value( PDE_DiscreteField const* ff,
                            size_t level,
                            GE_Point const* pt_ref,
                            size_t ic ) const ;

   //-- Possible adjacent domain

      bool has_adjacent_bound( void ) const ;

      PDE_BoundFE const* adjacent_bound( void ) const ;

      void insert_adjacent_bound( PDE_BoundFE const* a_bound ) ;
    
   //-- Status

      // Test is the boundary is inflow :
      bool isInFlow( size_t level, PDE_DiscreteField const* velocity ) const ;
      bool isInFlow( size_t level,
                     PDE_DiscreteField const* velocity,
                     GE_Point const* ptRefCoord ) const ;

      // Outward normal (outward from the AdjacentMesh) :
      GE_Vector const* unit_outward_normal( void ) const ;

   //-- Adaptation

      void set_parent( PDE_BoundFE* a_parent ) ;

      virtual PDE_BoundFE* parent( void ) const ;

      virtual bool is_active( void ) const ;

   //-- Access

      bool has_adjacent_cell( void ) const  ;

      PDE_CellFE* adjacent_cell( void ) const ;

      void insert_adjacent_cell( PDE_CellFE* adjMesh ) ;

   //-- Input - Output

      void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_BoundFE( void ) ;
     ~PDE_BoundFE( void ) ;
      PDE_BoundFE( PDE_BoundFE const& other ) ;
      PDE_BoundFE& operator=( PDE_BoundFE const& other ) ;

      PDE_BoundFE( PEL_Object* a_owner,
                   size_t a_number,
                   GE_Mpolyhedron* a_polyhedron,
                   GE_Color const* a_color,
                   size_t a_refinement_level,
                   PDE_DiscOnMeshFE const* a_disc ) ;

      virtual PDE_DiscOnMeshFE const* disc( void ) const ;

   //--  Attributes

      PDE_DiscOnMeshFE const* DISC ;
      PDE_CellFE* adjacentMesh ;
      PDE_BoundFE const* ADJ_BOUND ;
      mutable GE_Vector* outwardNormal ;
      
      PDE_BoundFE* PARENT ;

      std::vector< std::vector< PDE_BasisFunctionBound* > > BFS ;
} ;

#endif
