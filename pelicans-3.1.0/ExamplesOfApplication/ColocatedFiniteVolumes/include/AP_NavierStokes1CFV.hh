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

#ifndef AP_NAVIER_STOKES_1_CFV_HH
#define AP_NAVIER_STOKES_1_CFV_HH

#include <FE_OneStepIteration.hh>

class PEL_ContextSimple ;
class PEL_DoubleVector ;
class PEL_Double ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_CursorFEside ;
class PDE_SetOfBCs ;

class FE_Parameter ;
class FE_SetOfParameters ;

class AP_NavierStokes1System ;

/*
PUBLISHED
*/

class AP_NavierStokes1CFV : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_NavierStokes1CFV( void ) ;
      AP_NavierStokes1CFV( AP_NavierStokes1CFV const& other ) ;
      AP_NavierStokes1CFV& operator=( 
                              AP_NavierStokes1CFV const& other ) ;

      AP_NavierStokes1CFV( PEL_Object* a_owner, 
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_NavierStokes1CFV( void ) ;

      virtual AP_NavierStokes1CFV* create_replica( 
                                       PEL_Object* a_owner,
				       PDE_DomainAndFields const* dom,
				       FE_SetOfParameters const* prms,
				       PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_sides( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;

      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      double pressure_integral( size_t level ) ;

   //-- Class attributes

      static AP_NavierStokes1CFV const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* UU ;
      size_t L_UPDATE_UU ;
      size_t L_EXPLICIT_UU ;

      PDE_DiscreteField* PP ;
      size_t L_UPDATE_PP ;
      size_t L_EXPLICIT_PP ;

      double ALPHA ;
      double MU ;
      FE_Parameter* PI ;

      double LAMBDA ;
      double h_EXPONENT ;
      GE_Color const* CLUSTER_BD ;

      PDE_SetOfBCs const* BCs ;

      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;

      GE_QRprovider const* QRP_PI ;

      AP_NavierStokes1System* GLOBAL_EQ ;

      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* XX ;
      PEL_Double* TT ;
      double T_EXPLICIT ;
      double T_UPDATE ;
} ;

#endif
