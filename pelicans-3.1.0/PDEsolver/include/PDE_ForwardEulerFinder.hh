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

#ifndef PDE_FORWARD_EULER_FINDER_HH
#define PDE_FORWARD_EULER_FINDER_HH

#include <PDE_CFootFinder.hh>

#include <boolVector.hh>

class PEL_ModuleExplorer ;
class PEL_List ;

class GE_Point ;
class GE_Vector ;

class LA_Solver ;
class LA_DenseMatrix ;
class LA_SeqVector ;

class PDE_BoundFE ;
class PDE_CellFE ;
class PDE_DiscreteField ;
class PDE_GridFE ;
class PDE_MeshFE ;

class PEL_EXPORT PDE_ForwardEulerFinder : public PDE_CFootFinder
{
   public: //---------------------------------------------------------------

   //-- Instance creation and initialization

      static PDE_ForwardEulerFinder* create( PDE_LocalFEcell* fem,
                                             PEL_ModuleExplorer const* exp,
                                             PDE_GridFE const* grid ) ;

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------
      
      PDE_ForwardEulerFinder( void ) ;
     ~PDE_ForwardEulerFinder( void ) ;
      PDE_ForwardEulerFinder( PDE_ForwardEulerFinder const& other ) ;
      PDE_ForwardEulerFinder& operator=( 
                              PDE_ForwardEulerFinder const& other ) ;

      PDE_ForwardEulerFinder( PDE_LocalFEcell* fem,
                              PEL_ModuleExplorer const* exp,
                              PDE_GridFE const* grid ) ;

   //-- Characteristic foot

      virtual GE_Point const* foot( void ) const ;

      virtual GE_Point const* foot_ref( void ) const ;

   //-- Perform search

      virtual void search_foot( PDE_DiscreteField const* aa,
                                size_t level,
                                double time_step,
                                GE_Point const* M,
                                PDE_CellFE const* head_cell ) ;

      void search_in_cell( PDE_DiscreteField const* aa,
                           size_t level, 
                           double time_step,
                           GE_Point const* M,
                           PDE_CellFE const* trial_mesh ) ;

      void search_in_bound( PDE_DiscreteField const* aa,
                            size_t level, 
                            double time_step,
                            GE_Point const* M,
                            PDE_BoundFE const* trial_bound ) ;

      void initialize_trial_cells( PDE_CellFE const* head_cell ) ;

      void extend_trial_cells( PDE_CellFE const* a_cell,
                               PDE_CellFE const* head_cell ) ;

      PDE_CellFE const* pop_one_trial_cell( void ) ;

      bool more_trial_cells( void ) const ;

      void transport( PDE_DiscreteField const* aa,
                      size_t level,
                      double time_step,
                      PDE_CellFE const* cell_of_P,
                      GE_Point const* P_ref,
                      GE_Point const* P,
                      GE_Point* tr_P ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      PEL_List* MESHES ;
      boolVector ACCEPTABLE ;

      bool FOUND ;
      GE_Point* FOOT ;
      GE_Point* FOOT_REF ;

      double REF_DIST_MAX ;

      bool FOUND_IN_CELL ;
      GE_Point* CELL_FT ;
      GE_Point* CELL_FT_REF ;
      size_t CELL_ITER_MAX ;

      bool FOUND_IN_BOUND ;
      GE_Point* BOUND_FT ;
      GE_Point* BOUND_FT_REF ;
      size_t BOUND_ITER_MAX ;

      size_t DIM ;
      LA_Solver* SOLVER ;
      LA_DenseMatrix* A ;
      LA_SeqVector* b  ;
      LA_SeqVector* dX ;

      GE_Point* P ;
      GE_Point* Q ;
      GE_Point* tr_P ;
      GE_Point* tr_Q ;
      GE_Point* P_refC ;
      GE_Point* Q_refC ;
      GE_Point* P_refB ;
      GE_Point* Q_refB ;
      GE_Vector* PM ;
      GE_Vector* QM ;
      GE_Vector* aP ;
      GE_Point* PT_REF ;
} ;

#endif
