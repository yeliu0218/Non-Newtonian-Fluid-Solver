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

#ifndef FE_REFINEMENT_APPLI_HH
#define FE_REFINEMENT_APPLI_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_Point ;
class GE_QRprovider ;

class PDE_AdapterCHARMS ;
class PDE_DomainAndFields ;
class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_ResultSaver ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class FE_Parameter ;
class FE_SetOfParameters ;
class FE_TimeIterator ;

class FE_RefinementAppli : public PEL_Application
{
   public: //-----------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_RefinementAppli( void ) ;
      FE_RefinementAppli( FE_RefinementAppli const& other ) ;
      FE_RefinementAppli& operator=( FE_RefinementAppli const& other ) ;

      FE_RefinementAppli( PEL_Object* a_owner, 
                          PEL_ModuleExplorer const* exp ) ;

    //-- Plug in

      FE_RefinementAppli( void ) ;

      virtual FE_RefinementAppli* create_replica( 
	                         PEL_Object* a_owner,
	                         PEL_ModuleExplorer const* exp ) const ;

    //-- Internals

      void do_saving( void ) ;

      void iterate( void ) ;

      void read_one_step_iteration_data( PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer const* exp ) ;

      void do_before_time_stepping( FE_TimeIterator const* t_it,
                                    PDE_AdapterCHARMS* da ) ;

      void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      bool check_consistency( FE_TimeIterator const* t_it ) const ;

      void check_faces_consistency( void ) ;
      void check_covering_of_cells_boundary( void ) ;
      
      void display_error( PDE_LocalFEcell const* fe,
                          double theo, double xx ) const ;

    //-- Class attributes

      static FE_RefinementAppli const* PROTOTYPE ;
  
    //-- Attributes

      PDE_DomainAndFields* DOM ;

      PDE_DiscreteField* UU ;
      size_t L_UU ;

      PDE_DiscreteField* UU_EXP ;
      size_t L_UU_EXP ;

      FE_TimeIterator* TIME_IT ;

      FE_SetOfParameters const* PRMS ;

      PDE_SetOfBCs const* BCs ;

      PDE_ResultSaver* SAVER ;

      size_t NB_ITER_MAX ;
      
      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      double ALPHA ;
      double LAMBDA ;

      FE_Parameter* FORCE_TERM ;
      FE_Parameter* AA ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;
      
      LA_Solver* SOLVER ;
      FE_Parameter* PRM_CHECK ;
      PDE_DiscreteField* UU_CHECK ;
      double EPS_CHECK ;
      double MIN_CHECK ;
} ;

#endif
