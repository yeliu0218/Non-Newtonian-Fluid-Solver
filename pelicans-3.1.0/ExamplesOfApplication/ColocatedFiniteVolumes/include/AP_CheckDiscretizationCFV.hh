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

#ifndef AP_CHECK_DISCRETIZATION_CFV_HH
#define AP_CHECK_DISCRETIZATION_CFV_HH

#include <FE_OneStepIteration.hh>

#include <vector>

class GE_Mpolyhedron ;

class PDE_DiscreteField ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_CursorFEside ;

class FE_Parameter ;
class FE_SetOfParameters ;

class AP_CheckDiscretizationCFV : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_CheckDiscretizationCFV( void ) ;
      AP_CheckDiscretizationCFV( AP_CheckDiscretizationCFV const& other ) ;
      AP_CheckDiscretizationCFV& operator=( 
                               AP_CheckDiscretizationCFV const& other ) ;

      AP_CheckDiscretizationCFV( PEL_Object* a_owner, 
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_CheckDiscretizationCFV( void ) ;

      virtual AP_CheckDiscretizationCFV* create_replica( 
                                       PEL_Object* a_owner,
				       PDE_DomainAndFields const* dom,
				       FE_SetOfParameters const* prms,
				       PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      void display_side_pb( PDE_CursorFEside const* fe,
                            double dd0,
                            double dd1,
                            double dd,
                            bool ok_n ) const ;

      void display_bound_pb( PDE_LocalFEbound const* fe,
                             double dd ) const ;

      void display_not_in( GE_Mpolyhedron const* poly ) const ;

   //-- Class attributes

      static AP_CheckDiscretizationCFV const* PROTOTYPE ;

   //-- Attributes

      std::vector< PDE_DiscreteField const* > FIELDS ;
      std::vector< PDE_DiscreteField const* > BCFIELDS ;

      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;

      PDE_SetOfBCs const* BCs ;

      double MAX_NEG_DISTANCE_TO_SIDE ;
      double MAX_SCALAR_PRODUCT ;
      double MIN_DISTANCE_CENTERS ;
} ;

#endif
