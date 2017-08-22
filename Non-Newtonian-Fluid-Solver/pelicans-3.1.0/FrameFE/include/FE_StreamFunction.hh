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

#ifndef FE_STREAM_FUNCTION_HH
#define FE_STREAM_FUNCTION_HH

#include <FE_OneStepIteration.hh>

class GE_QRprovider ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DomainAndFields ;
class PDE_DiscreteField ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;
class PDE_LocalFEbound ;
class PDE_SystemNumbering ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_StreamFunction : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
                                                  FE_TimeIterator const* t_it,
                                                  PDE_ResultSaver* rs ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_StreamFunction( void ) ;
      FE_StreamFunction( FE_StreamFunction const& other ) ;
      FE_StreamFunction& operator=( FE_StreamFunction const& other ) ;

      FE_StreamFunction( PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_StreamFunction( void ) ;

      virtual FE_StreamFunction* create_replica(
                                        PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FE_StreamFunction const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* SF ;
      size_t L_UPDATE ;

      PDE_DiscreteField* VV ;
      size_t L_VV ;

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
