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

#ifndef PDE_CELL_FE_HH
#define PDE_CELL_FE_HH

#include <PDE_MeshFE.hh>

#include <PDE_BasisFunctionCell.hh>

class PEL_Vector ;
class PDE_FaceFE ;
class PDE_FacesOfCellFE ;

class PEL_EXPORT PDE_CellFE : public PDE_MeshFE
{
   public: //-----------------------------------------------------------

   //-- Instance creation and initialization

      static PDE_CellFE* create( PEL_Object* a_owner,
                                 size_t a_number,
                                 GE_Mpolyhedron* a_polyhedron,
                                 GE_Color const* a_color,
                                 size_t a_refinement_level,
                                 PDE_DiscOnMeshFE const* a_disc ) ;

      static PDE_CellFE* create_with_discretizations_pattern( 
                                 PEL_Object* a_owner,
                                 size_t a_number,
                                 GE_Mpolyhedron* a_polyhedron,
                                 GE_Color const* a_color,
                                 size_t a_refinement_level,
                                 PDE_CellFE const* a_pattern_mesh ) ;

   //-- Basis functions

      bool all_basis_functions_are_dropped( void ) const ;

      void set_basis_function( size_t ee, 
                               size_t ln, 
                               PDE_BasisFunctionCell* bf ) ;

      void remove_basis_function( size_t ee, size_t ln ) ;

      virtual PDE_BasisFunctionCell* basis_function( size_t ee, 
                                                     size_t ln ) const ;

   //-- Adaptation

      void set_parent( PDE_CellFE* a_parent ) ;

      virtual PDE_CellFE* parent( void ) const ;

      void set_nb_childs( size_t a_nb ) ;

      size_t nb_childs( void ) const ;

      void set_child( size_t i, PDE_CellFE* a_child ) ;

      PDE_CellFE* child( size_t i ) const ;

      virtual bool is_active( void ) const ;

      void set_active( void ) ;

      void set_inactive( void ) ;

   //-- Element change

      void append_face( PDE_FaceFE* a_face ) ;

   //-- Access      

      size_t nb_faces( void ) const ;

      // list of the faces of `self' that have the same refinement level
      // than `self'
      PEL_Vector const* faces( void ) const ;
      
      // iterator over all the faces included in the boundary of `self'
      // (i.e. those accessible via `::faces()' and, recursively, all their
      // childs)
      PDE_FacesOfCellFE* create_faces_iterator( PEL_Object* a_owner ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------
      		
   private: //----------------------------------------------------------

      PDE_CellFE( void ) ;
     ~PDE_CellFE( void ) ;
      PDE_CellFE( PDE_CellFE const& other ) ;
      PDE_CellFE& operator =( PDE_CellFE const& other ) ;

      PDE_CellFE( PEL_Object* a_owner,
                  size_t a_number,
                  GE_Mpolyhedron* a_polyhedron,
                  GE_Color const* a_color,
                  size_t a_refinement_level,
                  PDE_DiscOnMeshFE const* a_disc ) ;

      PDE_CellFE( PEL_Object* a_owner,
                  size_t a_number,
                  GE_Mpolyhedron* a_polyhedron,
                  GE_Color const* a_color,
                  size_t a_refinement_level,
                  PDE_CellFE const* a_pattern_mesh ) ;

      virtual PDE_DiscOnMeshFE const* disc( void ) const ;

   //-- Attributes

      PDE_DiscOnMeshFE const* DISC ;
      PEL_Vector* FACES ;
      PDE_CellFE* PARENT ;
      PEL_Vector* CHILDS ;
      bool ACTIVE ;
      std::vector< std::vector< PDE_BasisFunctionCell* > > BFS ;
} ;

#endif
