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

#ifndef FE_ITERATED_SYSTEM
#define FE_ITERATED_SYSTEM

#include <FE_OneStepIteration.hh>

class PEL_List ;
class PEL_ListIterator ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_IteratedSystem : public FE_OneStepIteration
{
   public: //------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_before_inner_iterations_stage( 
                                           FE_TimeIterator const* t_it ) ;

      virtual void do_inner_iterations_stage( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual bool inner_iterations_are_completed( 
                                      FE_TimeIterator const* t_it ) const ;

      virtual void do_after_inner_iterations_stage(  
                                      FE_TimeIterator const* t_it ) ;

      virtual void do_after_time_adaptation(  
                                      FE_TimeIterator const* t_it ) ;

      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;
      
   //-- Time iterator modification

      virtual void adapt_time_iterator( FE_TimeIteratorAdapter* t_adapter ) ;

      virtual void notify_inner_iterations_stage_failure( void ) ;
      
      virtual void reset_after_inner_iterations_stage_failure( void ) ;
      
   //-- Savings for post-processing

      virtual void do_additional_savings( FE_TimeIterator const* t_it,
                                          PDE_ResultSaver* rs ) ;
      
    //-- Persistence
      
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

     ~FE_IteratedSystem( void ) ;
      FE_IteratedSystem( FE_IteratedSystem const& other ) ;
      FE_IteratedSystem& operator=( FE_IteratedSystem const& other ) ;

      FE_IteratedSystem( PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp ) ;
      
      FE_IteratedSystem( PEL_Object* a_owner,
                         PDE_SetOfDomains const* sdoms,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      FE_IteratedSystem( void ) ;

      virtual FE_IteratedSystem* create_replica( 
                                        PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) const ;
      
      virtual FE_IteratedSystem* create_replica( 
                                        PEL_Object* a_owner,
                                        PDE_SetOfDomains const* sdoms,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FE_IteratedSystem const* PROTOTYPE ;

   //-- Attributes

      PEL_List* CMPS ;
      PEL_ListIterator* IT ;
      size_t NB_ITER_MAX ;
} ;

#endif
