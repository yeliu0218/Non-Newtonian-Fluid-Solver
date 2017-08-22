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

#ifndef RR_REFINEMENT_FIC_DIFFUSION_HH
#define RR_REFINEMENT_FIC_DIFFUSION_HH

#include <FE_OneStepIteration.hh>

#include <boolVector.hh>

#include <vector>

class GE_Point ;
class GE_QRprovider ;
class GE_SetOfPoints ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_AdapterCHARMS ;
class PDE_CursorFEside ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class FE_Parameter ;

class FE_RefinementFICdiffusion : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_RefinementFICdiffusion( void ) ;
      FE_RefinementFICdiffusion( FE_RefinementFICdiffusion const& other ) ;
      FE_RefinementFICdiffusion& operator=( 
                                 FE_RefinementFICdiffusion const& other ) ;

      FE_RefinementFICdiffusion( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_RefinementFICdiffusion( void ) ;

      virtual FE_RefinementFICdiffusion* create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer* exp ) const ;
      
   //-- Internals
      
      void prepare_for_multigrid( size_t level, size_t level_max ) ;
      
      void build_and_solve( FE_TimeIterator const* t_it, size_t level ) ;
      
      void loop_on_sides( FE_TimeIterator const* t_it, 
                          size_t level, 
                          int mode ) ;
      
      void loop_on_bounds( FE_TimeIterator const* t_it, size_t level ) ;
      
      void loop_on_cells( FE_TimeIterator const* t_it, size_t level ) ;
      
      bool converged( void ) const ;

   //-- Class attributes

      static FE_RefinementFICdiffusion const* PROTOTYPE ;
  
    //-- Attributes

      PDE_DiscreteField* UU ;
      size_t L_UPDATE ;
      
      double ALPHA ;
      double KAPPA ;
      FE_Parameter* PI ;
      
      PDE_SetOfBCs const* BCs ;

      PDE_AdapterCHARMS* DA ;
      std::vector< boolVector > OBSERVED_NODES ;
      LA_SeqVector* RESID ;
      boolVector HAS_R ;
      
      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;
      
      GE_SetOfPoints const* VERTS ;

      GE_Point* PT ;
      
      GE_QRprovider const* QRP_PI ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;
      
      LA_SeqVector* d_U0 ;
      
      LA_Solver* SOLVER ;
      
      double TOL ;
} ;

#endif
