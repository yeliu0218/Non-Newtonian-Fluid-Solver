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

#ifndef PDE_FACE_FE_HH
#define PDE_FACE_FE_HH

#include <PDE_CellFE.hh>

class PEL_List ;
class size_t_vector ;

class GE_Transform ;

class PDE_BoundFE ;
class PDE_CellFE ;

class PEL_EXPORT PDE_FaceFE : public PDE_MeshFE
{
   public: //--------------------------------------------------------------

   //-- Instance creation and initialization

      static PDE_FaceFE* create( PEL_Object* a_owner,
                                 size_t a_number,
                                 GE_Mpolyhedron* a_polyhedron,
                                 GE_Color const* a_color,
                                 size_t a_refinement_level ) ;

   //-- Basis functions

      virtual PDE_BasisFunction* basis_function( size_t ee, 
                                                 size_t ln ) const ;

   //-- Adaptation

      void set_parent( PDE_FaceFE* a_parent ) ;

      virtual PDE_FaceFE* parent( void ) const ;

      void set_nb_childs( size_t a_nb ) ;

      size_t nb_childs( void ) const ;

      void append_child( PDE_FaceFE* a_child ) ;

      PDE_FaceFE* child( size_t i ) const ;
      
      virtual bool is_active( void ) const ;

   //-- Element change

      void insert_adjacent_cell( PDE_CellFE* cell ) ;

      void replace_adjacent_cell( PDE_CellFE* old_cell, 
                                  PDE_CellFE* new_cell ) ;
      
      void set_adjacent_cells( PDE_CellFE* cell_a, PDE_CellFE* cell_b ) ;

      void insert_adjacent_bound( PDE_BoundFE* bound ) ;

      void set_side_id( size_t an_id ) ;
      
      void set_periodicity( PDE_FaceFE* side, GE_Transform const* tr ) ;

   //-- Access

      size_t side_id( void ) const ;

      size_t nb_adjacent_cells( void ) const ;
      
      bool has_adjacent_cell( PDE_CellFE const* m ) const ;

      PDE_CellFE* adjacent_cell( size_t i ) const ;

      PDE_CellFE* adjacent_cell_other_than( PDE_CellFE const* m ) const ;

      bool has_adjacent_bound( void ) const ;

      PDE_BoundFE* adjacent_bound( void ) const ;

      bool is_periodic( void ) const ;

      PDE_FaceFE* periodic_neighbour( void ) const ;
      
      GE_Transform const* periodic_transform( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------

      PDE_FaceFE( void ) ;
     ~PDE_FaceFE( void ) ;
      PDE_FaceFE( PDE_FaceFE const& other ) ;
      PDE_FaceFE& operator=( PDE_FaceFE const& other ) ;

      PDE_FaceFE( PEL_Object* a_owner,
                  size_t a_number,
                  GE_Mpolyhedron* a_polyhedron,
                  GE_Color const* a_color,
                  size_t a_refinement_level ) ;

      virtual PDE_DiscOnMeshFE const* disc( void ) const ;
      
   //-- Attributes

      PDE_CellFE* CELL_1 ;
      PDE_CellFE* CELL_2 ;

      PDE_BoundFE* bound ;
      size_t SIDE_ID ;
      PDE_FaceFE* PERIODIC_NEIGHBOUR ;
      GE_Transform const* PERIODIC_TRANSFORM ;

      PDE_FaceFE* PARENT ;
      PEL_Vector* CHILDS ;
} ;

#endif
