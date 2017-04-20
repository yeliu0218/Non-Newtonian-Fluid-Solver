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

#ifndef PDE_BASIS_FUNCTION_CELL_HH
#define PDE_BASIS_FUNCTION_CELL_HH

#include <PDE_BasisFunction.hh>
#include <PDE_SetOfBasisFunctions.hh>

class GE_Point ;

class PEL_EXPORT PDE_BasisFunctionCell : public PDE_BasisFunction
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_BasisFunctionCell* create( PEL_Object* a_owner ) ;

   //-- Activation and refinement

      virtual size_t refinement_level( void ) const ;

   //-- Parents

      size_t nb_parents( void ) const ;

      PDE_BasisFunctionCell* parent( size_t i ) const ;

      bool is_parent_of( PDE_BasisFunctionCell* a_child ) const ;

      bool parent_is_leading( size_t i ) const ;

   //-- Childs

      size_t nb_childs( void ) const ;

      PDE_BasisFunctionCell* child( size_t i ) const ;

      bool is_child_of( PDE_BasisFunctionCell const* a_parent ) const ;

      double refinement_coefficient( size_t i ) const ;

      // notify that `self' is a child of `a_parent' (`a_parent' is modified)
      // and that `a_parent' is a parent of `self'
      void set_child_parent_relationship( PDE_BasisFunctionCell* a_parent,
                                          double refinement_coef,
                                          bool is_leading ) ;

      void set_all_child_parent_relationship_on_cell(
                                      PDE_CellFE const* bf_cell, size_t ln,
                                      size_t ee, PDE_CellFE const* parent_cell,
                                      size_t ic, size_t verb_level ) ;

      bool check_child_parent_relationship( void ) const ;

   //-- Refinement deps

      size_t nb_ascendants( void ) const ;

      PDE_BasisFunctionCell* ascendant( size_t i ) const ;

      bool is_ascendant_of( PDE_BasisFunctionCell const* a_fine ) const ;

      size_t nb_descendants( void ) const ;

      PDE_BasisFunctionCell* descendant( size_t i ) const ;

      bool is_descendant_of( PDE_BasisFunctionCell const* a_coarse ) const ;

      // notify that `a_coarse' is a dep_down of `self'
      // and `self' is a dep_up of `a_coarse'
      void set_ascendant_relationship( PDE_BasisFunctionCell* a_coarse ) ;

      void set_all_ascendant_relationship_on_cell(
                                        PDE_CellFE const* bf_cell, size_t ln,
                                        size_t ee, PDE_CellFE const* parent_cell,
                                        size_t ic, size_t verb_level ) ;

   //-- Pieces

      void extend_pieces( PDE_CellFE* a_cell,
                          size_t elm_index,
                          size_t node_in_elm ) ;

      void remove_from_pieces( PDE_CellFE* a_cell,
                               size_t elm_index,
                               size_t node_in_elm ) ;

      size_t nb_cells( void ) const ;

      // support of the `i'-th piece
      PDE_CellFE* cell( size_t i ) const ;

      // reference element of the `i'-th piece
      size_t element_index_of_cell( size_t i ) const ;

      // local node index in the `i'-th piece
      size_t local_node_of_cell( size_t i ) const ;
      
      // associated geometrical node to `self'
      void geometrical_node( GE_Point* pt ) const ;
      
      // index of associated reference element group 
      // in `PDE_SetOfBasisFunctions'
      size_t ref_elts_grp_index( void ) const ;
      
   //-- Location

      virtual bool is_located_in_cell( PDE_CellFE const* cell ) const ;

      virtual bool is_located_on_bound( PDE_BoundFE const* bound ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      void print_2( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      void set_ref_elts_grp_index( size_t ind ) ;

      PDE_BasisFunctionCell( void ) ;
     ~PDE_BasisFunctionCell( void ) ;
      PDE_BasisFunctionCell( PDE_BasisFunctionCell const& other ) ;
      PDE_BasisFunctionCell& operator=( PDE_BasisFunctionCell const& other ) ;

      PDE_BasisFunctionCell( PEL_Object* a_owner ) ;

      bool check_location_in_cell( PDE_CellFE const* a_cell ) const ;

      bool check_location_on_bound( PDE_BoundFE const* a_bound ) const ;

      static size_t face_index( PDE_BoundFE const* bound,
                                PDE_CellFE const* cell ) ;

      static GE_Point* tmp_point( size_t dim ) ;
      
      friend void PDE_SetOfBasisFunctions:: add( PDE_BasisFunctionCell* bf,
                                                 size_t e ) ; 

   //-- Attributes

      struct Piece
      {
          Piece( PDE_CellFE* c, size_t i, size_t j )
             : cell( c ), elm_index( i ), node_in_elm( j ) {}
          PDE_CellFE* cell ;
          size_t elm_index ;
          size_t node_in_elm ;
      } ;
      std::vector< Piece > PIECES ;

      std::vector< PDE_BasisFunctionCell* > PARENTS ;
      std::vector< double > REFI_COEFS ;
      std::vector< PDE_BasisFunctionCell* > CHILDS ;

      std::vector<PDE_BasisFunctionCell*> ASCENDS ;
      std::vector<PDE_BasisFunctionCell*> DESCENDS ;

      size_t LEADING_PARENT ;
      
      size_t REF_ELT_INDEX ;
} ;

#endif
