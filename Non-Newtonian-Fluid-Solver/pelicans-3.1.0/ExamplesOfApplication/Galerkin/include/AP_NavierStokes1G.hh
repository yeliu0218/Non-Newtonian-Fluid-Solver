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

#ifndef AP_NAVIER_STOKES_1_G_HH
#define AP_NAVIER_STOKES_1_G_HH

#include <FE_OneStepIteration.hh>

#include <vector>

class doubleVector ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_TwoBlocksMethod ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class FE_LocalBCsBuilder ;
class FE_Parameter ;

/*
PUBLISHED
*/

class AP_NavierStokes1G : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
            
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_NavierStokes1G( void ) ;
      AP_NavierStokes1G( AP_NavierStokes1G const& other ) ;
      AP_NavierStokes1G& operator=( AP_NavierStokes1G const& other ) ;

      AP_NavierStokes1G( PEL_Object* a_owner, 
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_NavierStokes1G( void ) ;

      virtual AP_NavierStokes1G* create_replica( 
                                        PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      enum TimeDisc{ Euler, BDF2, NoTime } ;
      
      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;
      
      void reset_discrete_problem( void ) ;

      void compute_coefs_at_IP( FE_TimeIterator const* t_it,
                                double& m_xx,
                                doubleVector& aa,
                                double& mu,
                                doubleVector& rhs ) const ;
                                
   //-- Class attributes

      static AP_NavierStokes1G const* PROTOTYPE ;

   //-- Attributes

      size_t NB_DIMS ;
      
      PDE_DiscreteField* UU ;
      size_t L_UU ;
      
      PDE_DiscreteField* UU_EXP ;
      size_t L_UU_EXP ;
      
      PDE_DiscreteField* UU_EXP_EXP ;
      size_t L_UU_EXP_EXP ;

      PDE_DiscreteField* PP ;
      size_t L_PP ;
      
      PDE_DiscreteField* PP_EXP ;
      size_t L_PP_EXP ;
      
      TimeDisc TDISC ;
      
      bool L2_STAB ;

      FE_Parameter* AA ;
      FE_Parameter* ALPHA ;
      FE_Parameter* ALPHA_EXP ;
      FE_Parameter* MU ;
      FE_Parameter* RHSU ;
      bool LAPL_UU ;

      PDE_SetOfBCs const* BCs ;

      FE_LocalBCsBuilder* LOCAL_BC ;

      PDE_LocalEquation* ELEMENT_EQ ;

      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      PDE_SystemNumbering* NMB_U ;
      PDE_SystemNumbering* NMB_P ;
      
      LA_Matrix* A ;
      LA_Matrix* B ;
      LA_Matrix* L ;
      LA_Matrix* M ;
      
      LA_Vector* F ; 
      LA_Vector* G ;
      LA_Vector* S ;
      LA_Vector* K ;
      
      LA_Vector* U ;
      LA_Vector* P ;
      
      LA_SeqVector* U_LOC ;
      LA_SeqVector* P_LOC ;
      
      bool INIT_DISCRETE_P ;
      bool INIT_DISCRETE_V ;
      
      LA_TwoBlocksMethod* SOLVER ;
} ;

#endif
