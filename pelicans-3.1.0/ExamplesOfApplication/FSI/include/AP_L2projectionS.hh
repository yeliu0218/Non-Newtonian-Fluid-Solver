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

#ifndef AP_L2_PROJECTION_S_HH
#define AP_L2_PROJECTION_S_HH

#include <FE_OneStepIteration.hh>

class doubleVector ;
class doubleArray2D ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;

class AP_ConstitutiveLaw ;
class AP_KinematicState ;

/*
PUBLISHED
*/

class AP_L2projectionS : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
                                                  FE_TimeIterator const* t_it,
                                                  PDE_ResultSaver* rs ) ;

      void stress_to_smooth( doubleArray2D const& S,
                             doubleVector& s ) const ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_L2projectionS( void ) ;
      AP_L2projectionS( AP_L2projectionS const& other ) ;
      AP_L2projectionS& operator=( AP_L2projectionS const& other ) ;

      AP_L2projectionS( PEL_Object* a_owner, 
                        PDE_DomainAndFields const* dom,
                        FE_SetOfParameters const* prms,
                        PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_L2projectionS( void ) ;

      virtual AP_L2projectionS* create_replica( PEL_Object* a_owner,
                                                PDE_DomainAndFields const* dom,
                                                FE_SetOfParameters const* prms,
                                                PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static AP_L2projectionS const* PROTOTYPE ;

   //-- Attributes   

      PDE_DiscreteField* SS ;
      size_t L_SS ;

      PDE_DiscreteField* SD ;
      size_t L_SD ;

      size_t NB_COMP_S ;
      size_t NB_COMP_D ;

      AP_KinematicState* ST ;
      AP_ConstitutiveLaw* LAW ;

      PDE_LocalEquation* ELEMENT_EQ ;
     
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;

      LA_Solver* SOLVER ;
} ;

#endif
