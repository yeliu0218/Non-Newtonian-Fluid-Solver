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

#ifndef PDE_ADAPTATION_REQUEST_HH
#define PDE_ADAPTATION_REQUEST_HH

#include <PEL_Object.hh>

#include <vector>

class PEL_ListIdentity ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class PEL_Vector ;

class PDE_BasisFunctionCell ;
class PDE_CellFE ;
class PDE_CrossProcessBFNumbering ;
class PDE_ReferenceElement ;
class PDE_SetOfBasisFunctions ;

class PEL_EXPORT PDE_AdaptationRequest : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Internal configuration by applying adaptation rules

      void apply_criterion_to_infer_bfs_to_refine( PEL_Vector const* cells ) ;

      void compute_bfs_to_refine( PDE_SetOfBasisFunctions* bf_set = 0 ) ;

      bool something_to_refine( void ) const ;

      void build_info_on_cells_for_refinement ( void ) ;

      void extend_bfs_to_refine( PEL_ListIdentity* bfs_to_refine,
                                 PDE_SetOfBasisFunctions* bf_set = 0 ) ;

      void read_bfs_to_refine( PEL_ListIdentity const* bfs_to_refine ) ;

      void apply_criterion_to_infer_bfs_to_unrefine( PEL_Vector const* cells ) ;

      void compute_bfs_to_unrefine( PDE_SetOfBasisFunctions* bf_set = 0 ) ;

      bool something_to_unrefine( void ) const ;

      void build_info_on_cells_for_unrefinement( void ) ;

      void extend_bfs_to_unrefine( PEL_ListIdentity* bfs_to_unrefine,
                                   PDE_SetOfBasisFunctions* bf_set = 0 ) ;

      void read_bfs_to_unrefine( PEL_ListIdentity const* bfs_to_unrefine ) ;

   //-- Basis functions to refine : iterator movement

      void start_bf_to_refine( void ) ;

      bool valid_bf_to_refine( void ) const ;

      void go_next_bf_to_refine( void ) ;

      PDE_BasisFunctionCell* current_bf_to_refine( void ) const ;

   //-- Basis functions to unrefine : iterator movement

      void start_bf_to_unrefine( void ) ;

      bool valid_bf_to_unrefine( void ) const ;

      void go_next_bf_to_unrefine( void ) ;

      PDE_BasisFunctionCell* current_bf_to_unrefine( void ) const ;

   //-- Cells : iterator movement

      void start_cell_to_refine( void ) ;

      bool valid_cell_to_refine( void ) const ;

      void go_next_cell_to_refine( void ) ;

      PDE_CellFE* current_cell_to_refine( void ) const ;

      void start_cell_to_unrefine( void ) ;

      bool valid_cell_to_unrefine( void ) const ;

      void go_next_cell_to_unrefine( void ) ;

      PDE_CellFE* current_cell_to_unrefine( void ) const ;


   protected: //--------------------------------------------------------

      virtual ~PDE_AdaptationRequest( void ) ;

      PDE_AdaptationRequest( PEL_Object* a_owner,
                             size_t highest_refinement_level,
                             size_t verbose_level ) ;

   //-- Internal configuration by applying adaptation rules

      virtual bool to_be_refined( PDE_CellFE const* cell,
                                  PDE_BasisFunctionCell const* bf,
                                  PDE_ReferenceElement const* elm,
                                  size_t local_node ) const = 0 ;

      virtual bool to_be_unrefined( PDE_CellFE const* cell,
                                    PDE_BasisFunctionCell const* bf,
                                    PDE_ReferenceElement const* elm,
                                    size_t local_node ) const = 0 ;

   //-- Preconditions, Postconditions, Invariant

      bool to_be_refined_PRE( PDE_CellFE const* cell,
                              PDE_BasisFunctionCell const* bf,
                              PDE_ReferenceElement const* elm,
                              size_t local_node ) const ;

      bool to_be_unrefined_PRE( PDE_CellFE const* cell,
                                PDE_BasisFunctionCell const* bf,
                                PDE_ReferenceElement const* elm,
                                size_t local_node ) const ;

   private: //----------------------------------------------------------

      PDE_AdaptationRequest( void ) ;
      PDE_AdaptationRequest( PDE_AdaptationRequest const& other ) ;
      PDE_AdaptationRequest& operator=( PDE_AdaptationRequest const& other ) ;

   //-- Internals

      void extend_with_ascendants( PDE_BasisFunctionCell* bf,
                                   PEL_ListIdentity* bfs_to_refine,
                                   size_t indent ) const ;

      void consider_unrefining_bf( PDE_BasisFunctionCell* bf,
                                   PEL_ListIdentity* bfs_to_unrefine ) ;

      void synchronize_selected_bfs_to_refine( PEL_ListIdentity* bfs,
                                       PDE_SetOfBasisFunctions* bf_set ) ;

      void synchronize_selected_bfs_to_unrefine( PEL_ListIdentity* bfs,
                                      PDE_SetOfBasisFunctions* bf_set ) ;

   //-- Attributes

      std::vector< std::vector< PDE_CellFE* > > CELLS_TO_REFINE ;
      std::vector< std::vector< PDE_CellFE* > > CELLS_TO_UNREFINE ;

      size_t iM_R ;
      size_t iL_R ;

      size_t iM_U ;
      size_t iL_U ;

      size_t MAX_REF_LEVEL ;
      size_t VERB_LEVEL ;

      PEL_ListIdentity* BFS_DUM_R ;
      PEL_ListIdentity* BFS_DUM_U ;

      PEL_ListIterator* IT_R_CRITERIA ;
      PEL_ListIdentity* BFS_R_CRITERIA ;

      PEL_ListIterator* IT_U_CRITERIA ;
      PEL_ListIdentity* BFS_U_CRITERIA ;

      PEL_ListIterator* IT_R ;
      PEL_ListIdentity const* BFS_R ;

      PEL_ListIterator* IT_U ;
      PEL_ListIdentity const* BFS_U ;

      size_t SOMETHING_TO_REF ;
} ;

#endif
