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

#ifndef FE_GRID_MOTION_HH
#define FE_GRID_MOTION_HH

#include <FE_OneStepIteration.hh>

class GE_Color ;
class GE_QRprovider ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_GridMover ;
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

class PEL_EXPORT FE_GridMotion : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual void do_after_inner_iterations_stage( 
	                                          FE_TimeIterator const* t_it )  ;

   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
						                          FE_TimeIterator const* it,
						                          PDE_ResultSaver* rs ) ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~FE_GridMotion( void ) ;
      FE_GridMotion( FE_GridMotion const& other ) ;
      FE_GridMotion& operator=( FE_GridMotion const& other ) ;

      FE_GridMotion( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom,
                     FE_SetOfParameters const* prms,
                     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_GridMotion( void ) ;

      virtual FE_GridMotion* create_replica( PEL_Object* a_owner,
                                             PDE_DomainAndFields const* dom,
                                             FE_SetOfParameters const* prms,
                                             PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FE_GridMotion const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* CC ;
      PDE_LinkDOF2Unknown* CC_link ;
      size_t L_UPDATE ;

      FE_Parameter* ALPHA ;

      PDE_GridMover* GRID_MOVER ;

      PDE_SetOfBCs const* BCs ;

      FE_LocalBCsBuilder* LOCAL_BC ;

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
