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

#ifndef FE_POISSON_MORTAR_HH
#define FE_POISSON_MORTAR_HH

#include <FE_OneStepIterationOpen.hh>

class PEL_ModuleExplorer ;

class GE_Point ;
class GE_QRprovider ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;
class PDE_SystemNumbering ;

class FE_PoissonMortar : public FE_OneStepIterationOpen
{
   public: //-----------------------------------------------------------

   //-- Unknowns

      virtual size_t nb_unknowns( void ) const ;
      
      virtual PDE_DiscreteField* field( size_t i_unk ) const ;

      virtual size_t level_of_field( size_t i_unk ) const ;
      
      virtual PDE_LinkDOF2Unknown const* link_DOF_2_unknown( 
                                                    size_t i_unk ) const ;

   //-- Discretization

      virtual void assemble_contribution( FE_TimeIterator const* t_it,
                                          LA_Matrix* matrix,
                                          LA_Vector* vector,
                                          PDE_SystemNumbering const* nmb,
                                          size_t i_link ) const ;

      virtual void update_DOFs( LA_Vector* vector,
                                PDE_SystemNumbering const* nmb,
                                size_t i_link ) const ;

   //-- Substeps of the step by step progression

      void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_PoissonMortar( void ) ;
      FE_PoissonMortar( FE_PoissonMortar const& other ) ;
      FE_PoissonMortar& operator=( FE_PoissonMortar const& other ) ;

      FE_PoissonMortar( PEL_Object* a_owner,
                        PDE_DomainAndFields const* dom,
                        FE_SetOfParameters const* prms,
                        PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_PoissonMortar( void ) ;

      virtual FE_PoissonMortar* create_replica( 
                                    PEL_Object* a_owner,
                                    PDE_DomainAndFields const* dom,
                                    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer* exp ) const ;

   //-- Internals

      static void add_row_at_pt( PDE_LocalEquation* leq,
                                 PDE_LocalFEcell* fe,
                                 GE_Point const* pt ) ;

   //-- Class attributes
         
      static FE_PoissonMortar const* PROTOTYPE ;
     
   //-- Attributes

      PDE_DiscreteField* UU ;
      PDE_LinkDOF2Unknown* UU_link ;
      size_t L_UU ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;

      GE_Point const* ORIGIN ;
      
      PDE_SystemNumbering* NMB ;
      
      LA_Matrix* LHS ;
      LA_Vector* RHS ;
      LA_Vector* UNK ;
      LA_SeqVector* UNK_LOC ;
      
      LA_Solver* SOLVER ;
} ;

#endif
