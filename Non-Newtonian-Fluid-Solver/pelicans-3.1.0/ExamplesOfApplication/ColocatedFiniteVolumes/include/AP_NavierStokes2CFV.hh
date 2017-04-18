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

#ifndef AP_NAVIER_STOKES_2_CFV_HH
#define AP_NAVIER_STOKES_2_CFV_HH

#include <FE_OneStepIterationOpen.hh>

class PEL_ContextSimple ;
class PEL_DoubleVector ;
class PEL_Double ;

class GE_QRprovider ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_CursorFEside ;
class PDE_SetOfBCs ;

class FE_Parameter ;
class FE_SetOfParameters ;

/*
PUBLISHED
*/

class AP_NavierStokes2CFV : public FE_OneStepIterationOpen
{
   public: //-----------------------------------------------------------------

   //-- Jacobian

      virtual size_t nb_unknowns( void ) const ;

      virtual PDE_DiscreteField* field( size_t i_unk ) const ;

      virtual size_t level_of_field( size_t i_unk ) const ;

      virtual PDE_LinkDOF2Unknown const* link_DOF_2_unknown( 
                                                    size_t i_unk ) const ;

      virtual void build_function_and_jacobian( FE_TimeIterator const* t_it ) ;

      virtual LA_SeqVector const* create_function( PEL_Object* a_owner,
                                                   size_t i_unk ) const ;

      virtual LA_SeqMatrix const* create_jacobian( PEL_Object* a_owner,
                                                   size_t i_eq, 
                                                   size_t j_unk ) const ;

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_NavierStokes2CFV( void ) ;
      AP_NavierStokes2CFV( AP_NavierStokes2CFV const& other ) ;
      AP_NavierStokes2CFV& operator=( 
                                AP_NavierStokes2CFV const& other ) ;

      AP_NavierStokes2CFV( PEL_Object* a_owner, 
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_NavierStokes2CFV( void ) ;

      virtual AP_NavierStokes2CFV* create_replica( 
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_sides( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;

      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      bool convergence_achieved( size_t iter ) const ;

      void add_to_A_item( size_t i_row, size_t j_col, double xx ) ;
      void add_to_F_item( size_t i_row, double xx ) ;
      void add_to_B_item( size_t i_row, size_t j_col, double xx ) ;
      void add_to_D_item( size_t i_row, size_t j_col, double xx ) ;
      void add_to_G_item( size_t i_row, double xx ) ;
      void add_to_C_item( size_t i_row, size_t j_col, double xx ) ;

   //-- Class attributes

      static AP_NavierStokes2CFV const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* UU ;
      size_t L_UPDATE_UU ;
      size_t L_EXPLICIT_UU ;

      PDE_DiscreteField* PP ;
      size_t L_UPDATE_PP ;
      size_t L_EXPLICIT_PP ;

      bool IS_UNSTEADY ;
      double C_DT ;
      double C_VGRAD ;
      double MU ;
      FE_Parameter* PI ;

      size_t NB_ITER_MAX ;
      double RELAX ;
      double TOL ;

      double LAMBDA ;
      double h_EXPONENT ;
      GE_Color const* CLUSTER_BD ;

      PDE_SetOfBCs const* BCs ;

      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;

      GE_QRprovider const* QRP_PI ;

      size_t idx_UU ;
      size_t idx_PP ;
      
      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;

      LA_Solver* SOLVER ;
      
      double PRINTED_RES_UU ;
      double PRINTED_RES_PP ;

      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* XX ;
      PEL_Double* TT ;
} ;

#endif
