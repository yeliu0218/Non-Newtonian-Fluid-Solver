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

#ifndef AP_DIFFUSION_CG_HH
#define AP_DIFFUSION_CG_HH

#include <FE_Galerkin.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class AP_DiffusionCG : public FE_Galerkin
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void transfer_calculation_requirements_for_material_derivative( 
                                                     PDE_LocalFEcell* fe ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual void reset_discrete_problem( FE_TimeIterator const* t_it ) ;

      virtual GE_QRprovider const* QRprovider_for_material_derivative( 
                                                                 void ) const ;

      virtual void build_cell_contribution_to_material_derivative( 
              	                                FE_TimeIterator const* t_it,
                                                PDE_LocalFEcell* fe ) ;

      virtual void terminate_discrete_problem( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~AP_DiffusionCG( void ) ;
      AP_DiffusionCG( AP_DiffusionCG const& other ) ;
      AP_DiffusionCG& operator=( AP_DiffusionCG const& other ) ;
      
      AP_DiffusionCG( PEL_Object* a_owner, 
                      PDE_DomainAndFields const* dom,
                      FE_SetOfParameters const* prms,
                      PEL_ModuleExplorer const* exp ) ;

      void build_cell_contribution_to_creation( FE_TimeIterator const* t_it,
                                                PDE_LocalFEcell* fe ) ;

      void build_bound_contribution_to_creation( FE_TimeIterator const* t_it,
                                                 PDE_LocalFEbound* fe ) ;

   //-- Plug in

      AP_DiffusionCG( void ) ;

      virtual AP_DiffusionCG* create_replica( PEL_Object* a_owner, 
                                              PDE_DomainAndFields const* dom,
                                              FE_SetOfParameters const* prms,
                                              PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      static void add_row_col( PDE_LocalEquation* leq,
                               PDE_LocalFE const* fe,
                               double coef ) ;

      static void add_grad_row_dot_grad_col( PDE_LocalEquation* leq,
                                             PDE_LocalFEcell const* fe,
                                             double coef ) ;

      static void add_row( PDE_LocalEquation* leq,
                           PDE_LocalFE const* fe,
                           double coef ) ;

      void estimate_unknowns( void ) ;
      
      void update_fields( void ) ;
      
   //-- Class attributes

      static AP_DiffusionCG const* PROTOTYPE ;
      
   //-- Attributes

      PDE_DiscreteField* TT ;
      size_t L_UPDATE ;
      size_t L_EXPLICIT ;

      FE_Parameter* DENS ;
      FE_Parameter* COND ;
      FE_Parameter* CP ;
      FE_Parameter* POW ;

      PDE_SetOfBCs const* BCs ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* U ;
      LA_SeqVector* U_LOC ;
      
      LA_Solver* SOLVER ;
} ;

#endif
