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

#ifndef AP_ADVECTION_DIFFUSION_1_CFV_HH
#define AP_ADVECTION_DIFFUSION_1_CFV_HH

#include <FE_OneStepIteration.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_CursorFEside ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class FE_Parameter ;
class FE_SetOfParameters ;

/*
PUBLISHED
*/

class AP_AdvectionDiffusion1CFV : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_AdvectionDiffusion1CFV( void ) ;
      AP_AdvectionDiffusion1CFV( AP_AdvectionDiffusion1CFV const& other ) ;
      AP_AdvectionDiffusion1CFV& operator=( 
                              AP_AdvectionDiffusion1CFV const& other ) ;

      AP_AdvectionDiffusion1CFV( PEL_Object* a_owner, 
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_AdvectionDiffusion1CFV( void ) ;

      virtual AP_AdvectionDiffusion1CFV* create_replica( 
                                         PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_sides( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;

      void loop_on_cells( FE_TimeIterator const* t_it ) ;
      
      double side_advective_flux( FE_TimeIterator const* t_it ) ;
      
      double bound_advective_flux( FE_TimeIterator const* t_it ) ;

   //-- Class attributes

      static AP_AdvectionDiffusion1CFV const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* UU ;
      size_t L_UPDATE ;
      size_t L_EXPLICIT ;

      FE_Parameter* ADV ;

      double ALPHA ;
      double KAPPA ;
      FE_Parameter* PI ;

      PDE_SetOfBCs const* BCs ;

      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;

      GE_QRprovider const* QRP_ADV ;
      GE_QRprovider const* QRP_PI ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;
      
      LA_Solver* SOLVER ;
} ;

#endif
