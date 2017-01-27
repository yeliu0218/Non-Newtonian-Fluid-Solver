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

#ifndef FE_SPLIT_SYSTEM
#define FE_SPLIT_SYSTEM

#include <FE_OneStepIteration.hh>

class PEL_List ;
class PEL_ListIterator ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_SplitSystem : public FE_OneStepIteration
{
   public: //------------------------------------------------------------

   //-- Substeps of the step by step progression

      // IMPLEMENTATION : Call `FE_OneStepIteration::do_before_time_stepping'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_before_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_before_inner_iterations_stage( 
                                           FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_inner_iterations_stage( FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : raise a fatal error.
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // true if `FE_OneStepIteration::inner_iterations_are_completed' is true
      // for all handled `FE_OneStepIteration::' instances, false otherwise
      virtual bool inner_iterations_are_completed( 
                                      FE_TimeIterator const* t_it ) const ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_inner_iterations_stage(  
                                      FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_time_adaptation'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_time_adaptation(  
                                      FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_time_stepping'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Time iterator modification

      virtual void adapt_time_iterator( FE_TimeIteratorAdapter* t_adapter ) ;

      virtual void notify_inner_iterations_stage_failure( void ) ;
      
      virtual void reset_after_inner_iterations_stage_failure( void ) ;

   //-- Savings for post-processing

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_additional_savings'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_additional_savings( FE_TimeIterator const* t_it,
                                          PDE_ResultSaver* rs ) ;
      
    //-- Persistence
      
      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::add_storable_objects'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;
           
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

     ~FE_SplitSystem( void ) ;
      FE_SplitSystem( FE_SplitSystem const& other ) ;
      FE_SplitSystem& operator=( FE_SplitSystem const& other ) ;

      FE_SplitSystem( PEL_Object* a_owner,
                      PDE_DomainAndFields const* dom,
                      FE_SetOfParameters const* prms,
                      PEL_ModuleExplorer const* exp ) ;
      
      FE_SplitSystem( PEL_Object* a_owner,
                      PDE_SetOfDomains const* sdoms,
                      FE_SetOfParameters const* prms,
                      PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_SplitSystem( void ) ;

      virtual FE_SplitSystem* create_replica( 
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const ;
      
      virtual FE_SplitSystem* create_replica( 
                                     PEL_Object* a_owner,
                                     PDE_SetOfDomains const* sdoms,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FE_SplitSystem const* PROTOTYPE ;

   //-- Attributes

      PEL_List* CMPS ;
      PEL_ListIterator* IT ;
} ;

#endif
