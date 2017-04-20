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

#ifndef FE_ADAPTATION_STEP_CHARMS_HH
#define FE_ADAPTATION_STEP_CHARMS_HH

#include <FE_OneStepIteration.hh>

class PEL_List ;
class PEL_ListIterator ;

class PDE_AdapterCHARMS ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_ResultSaver ;

class FE_AdaptationIndicator ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_AdaptationStepCHARMS : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_before_inner_iterations_stage( 
                                           FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual bool inner_iterations_are_completed( 
                                          FE_TimeIterator const* t_it ) const ;

      virtual void do_after_inner_iterations_stage( 
                                           FE_TimeIterator const* t_it ) ;

   //-- Savings for post-processing

      virtual void do_additional_savings( FE_TimeIterator const* t_it,
                                          PDE_ResultSaver* rs ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_AdaptationStepCHARMS( void ) ;
      FE_AdaptationStepCHARMS( FE_AdaptationStepCHARMS const& other ) ;
      FE_AdaptationStepCHARMS& operator=( 
                               FE_AdaptationStepCHARMS const& other ) ;

      FE_AdaptationStepCHARMS( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer* exp ) ;

   //-- Plug in

      FE_AdaptationStepCHARMS( void ) ;

      virtual FE_AdaptationStepCHARMS* create_replica( 
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const ;

   //-- Internals
      
      void give_time_to_indicators( FE_TimeIterator const* t_it ) ;

      void do_saving( FE_TimeIterator const* t_it ) ;
      
      void check_faces_consistency( void ) ;
      
      void check_covering_of_cells_boundary( void ) ;

   //-- Class attributes

      static FE_AdaptationStepCHARMS const* PROTOTYPE ;

   //-- Attributes

      PDE_DomainAndFields const* DOM ;
      PDE_SetOfDiscreteFields* DFS ;
      PEL_ModuleExplorer* EXP ;
      PDE_AdapterCHARMS* DA  ;
      PEL_List* INDICS ; // list of FE_AdaptationIndicator*
      PEL_ListIterator* INDICS_IT ;
      bool KEEP_ADAPTING ;
      PDE_ResultSaver* SAVER ;
      bool BUILD_INIT ;
      size_t NB_ITER_MAX_0 ;
      size_t NB_ITER_MAX_1 ;
      bool STOP_CV ;
      bool CHECK_FACES ;
      bool CHECK_BD_CELLS ;
} ;

#endif

