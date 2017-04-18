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

#ifndef AP_STOKES_CG_HH
#define AP_STOKES_CG_HH

#include <FE_Galerkin.hh>

class doubleVector ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_TwoBlocksMethod ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;
class PDE_SystemNumbering ;

class AP_StokesCG : public FE_Galerkin
{ 
   public: //-----------------------------------------------------------

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

   protected: //-------------------------------------------------------

   private: //----------------------------------------------------------

     ~AP_StokesCG( void ) ;
      AP_StokesCG( AP_StokesCG const& other ) ;
      AP_StokesCG& operator=( AP_StokesCG const& other ) ;
      
      AP_StokesCG( PEL_Object* a_owner, 
                  PDE_DomainAndFields const* dom,
                  FE_SetOfParameters const* prms,
                  PEL_ModuleExplorer const* exp ) ;

      void build_cell_contribution_to_creation( FE_TimeIterator const* t_it,
                                                PDE_LocalFEcell* fe ) ;

   //-- Plug in

      AP_StokesCG( void ) ;

      virtual AP_StokesCG* create_replica( PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      static void add_D_row_dot_grad_col( PDE_LocalEquation* leq,
                                          PDE_LocalFEcell const* fe,
                                          double coef ) ;

      static void add_row_div_col( PDE_LocalEquation* leq,
                                   PDE_LocalFEcell const* fe,
                                   double coef ) ;

      static void add_row_dot_col( PDE_LocalEquation* leq,
                                   PDE_LocalFEcell const* fe,
                                   double coef ) ;

      static void add_row( PDE_LocalEquation* leq,
                           PDE_LocalFEcell const* fe,
                           doubleVector const& coef ) ;

      double density_for_body_force( PDE_LocalFEcell const* fe ) const ;
      
      void estimate_unknowns( void ) ;
      
      void update_fields( void ) ;
      
      void re_initialize_matrices_and_vectors( size_t nv_glob,
                                               size_t np_glob,
                                               size_t nv_loc,
                                               size_t np_loc ) ;

   //-- Class attributes

      static AP_StokesCG const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* UU ;
      PDE_DiscreteField* PP ;

      size_t L_UPDATE ;
      size_t L_EXPLICIT ;

      FE_Parameter* DENS ;
      FE_Parameter* VISC ;
      FE_Parameter* RHSU ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;

      PDE_SystemNumbering* NMB_U ;
      PDE_SystemNumbering* NMB_P ;
      
      LA_Matrix* A ;
      LA_Matrix* B ;
      LA_Matrix* C ;
      
      LA_Vector* F ; 
      LA_Vector* G ;
      LA_Vector* S ;
      
      LA_Vector* U ;
      LA_Vector* P ;
      
      LA_SeqVector* U_LOC ;
      LA_SeqVector* P_LOC ;
      
      LA_TwoBlocksMethod* SOLVER ;
      bool SOLUTION_ACHIEVED ;      
} ;

#endif
