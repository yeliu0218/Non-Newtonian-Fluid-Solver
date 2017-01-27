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

#ifndef FE_MULTI_DOMAIN_SYSTEM_HH
#define FE_MULTI_DOMAIN_SYSTEM_HH

#include <FE_OneStepIteration.hh>

class PEL_Vector ;

class LA_Matrix ;
class LA_Solver ;
class LA_Vector ;

class PDE_SystemNumbering ;

class FE_MortarInterfaceDiscretizer ;
class FE_OneStepIterationOpen ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_MultiDomainSystem : public FE_OneStepIteration
{
   public: //------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

     ~FE_MultiDomainSystem( void ) ;
      FE_MultiDomainSystem( FE_MultiDomainSystem const& other ) ;
      FE_MultiDomainSystem& operator=( FE_MultiDomainSystem const& other ) ;

      FE_MultiDomainSystem( PEL_Object* a_owner,
                            PDE_SetOfDomains const* sdoms,
                            FE_SetOfParameters const* prms,
                            PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_MultiDomainSystem( void ) ;

      virtual FE_MultiDomainSystem* create_replica( 
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const ;
      
      virtual FE_MultiDomainSystem* create_replica( 
                                     PEL_Object* a_owner,
                                     PDE_SetOfDomains const* sdoms,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      FE_OneStepIterationOpen* subdomain_discretizer( size_t i ) const ;

      FE_MortarInterfaceDiscretizer* interface_discretizer( size_t i ) const ;

   //-- Class attributes

      static FE_MultiDomainSystem const* PROTOTYPE ;

   //-- Attributes

      PDE_SetOfDomains const* SDOMS ;

      PEL_Vector* D_PBS ;
      PEL_Vector* I_PBS ;
      
      PDE_SystemNumbering* NMB ;
      
      LA_Matrix* LHS ;
      LA_Vector* RHS ;
      LA_Vector* UNK ;
      
      LA_Solver* SOLVER ;
      
} ;

#endif
