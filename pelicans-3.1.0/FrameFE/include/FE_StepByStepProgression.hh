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

#ifndef FE_STEP_BY_STEP_PROGRESSION_HH
#define FE_STEP_BY_STEP_PROGRESSION_HH

#include <PEL_Application.hh>

#include <doubleVector.hh>

class PEL_ModuleExplorer ;
class PEL_Timer ;

class PDE_DomainAndFields ;
class PDE_ResultSaver ;
class PDE_SetOfDomains ;

class FE_OneStepIteration ;
class FE_SetOfParameters ;
class FE_TimeIteratorAdapter ;
class FE_TimeIterator ;


/*
Applications involving a time marching procedure in which unknowns are
computed at each step from their values at the previous steps.

Each instance is associated to
   1. an instance of `FE_TimeIterator::' that monitors the various steps
      of the time marching procedure ;
   2. an instance of `FE_OneStepIteration::' (constructed from the above 
      object) to which the progress of each step is delegated.

PUBLISHED 
*/

class PEL_EXPORT FE_StepByStepProgression : public PEL_Application
{
    public: //----------------------------------------------------------

    //-- Instance delivery and initialization

      static FE_StepByStepProgression* create( 
                                         PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) ;

    //-- Instance characteristics

      FE_TimeIterator const* time_iterator( void ) const ;

      FE_SetOfParameters const* set_of_parameters( void ) const ;

      PDE_DomainAndFields* domain_and_fields( void ) const ;

    //-- Program core execution

      /*
      Perform a complete time marching procedure, according to the following
      pseudo-code, where ONE_IT (resp. TIT) denotes the associated
      `FE_OneStepIteration::' (resp. `FE_TimeIterator::') instance.

      1. Call `FE_OneStepIteration::do_before_time_stepping' :
            ONE_IT->do_before_time_stepping()

      2. Perform the time loop, using `FE_TimeIterator::start',
         `FE_TimeIterator::is_finished', `FE_TimeIterator::go_next_time' :
 
         for( TIT->start() ; !TIT->is_finished() ; TIT->go_next_time() )
         {
            2.1 Call `FE_OneStepIteration::do_before_inner_iterations_stage':
                    ONE_IT->do_before_inner_iterations_stage( TIT ) ;

            2.2 Perform a stage of inner iterations, using 
                `FE_OneStepIteration::do_one_inner_iteration',
                `FE_OneStepIteration::inner_iterations_are_completed' :

                do
                   ONE_IT->do_one_inner_iteration( TIT ) ;
                while( !ONE_IT->inner_iterations_are_completed( TIT ) ) ;
      
            2.3 Call `FE_OneStepIteration::do_after_inner_iterations_stage':
                    ONE_IT->do_after_inner_iterations_stage( TIT ) ;

            2.4 If necessary and for subsequent postprocessing,
                save the current time (given by `FE_TimeIterator::time') 
                and the unknown fields, by calling
                   `PDE_ResultSaver::save_variable'
                   `PDE_ResultSaver::save_fields' 
                on behalf of ONE_IT->result_saver(). Then call 
                `FE_OneStepIteration::save_other_than_time_and_fields' :
                   ONE_IT->save_other_than_time_and_fields( TIT, 
                                                    ONE_IT->result_saver() ) ;

            2.5 Perform savings for persistence issues.
         }

      3. Call `FE_OneStepIteration::do_after_time_stepping' :
            ONE_IT->do_after_time_stepping()

      4. If necessary, perform final savings as in 2.4 .
      */
      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_StepByStepProgression( void ) ;
      FE_StepByStepProgression( FE_StepByStepProgression const& other ) ;
      FE_StepByStepProgression& operator=( 
                                FE_StepByStepProgression const& other ) ;

      FE_StepByStepProgression( PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_StepByStepProgression( void ) ;

      virtual FE_StepByStepProgression* create_replica( 
                                        PEL_Object* a_owner,
					PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      void save_for_restart( void ) ;

      void save_for_post_processing( bool force ) ;
      
      void display_timers( size_t indent_width ) const ;

      void do_save_for_post( PDE_ResultSaver* rs ) ;

      void check_times( std::string const& list_name,
                        doubleVector const& dates ) const ;
      
      static void display_memory_usage(
                                  std::ostream& os, size_t indent_width ) ;

   //-- Persistence
      
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static FE_StepByStepProgression const* PROTOTYPE ;

      static size_t graphics_level ;

   //-- Attributes

      FE_TimeIteratorAdapter* TIME_ADAPT ;
      FE_TimeIterator* TIME_IT ;
      
      PDE_DomainAndFields* DOM ;
      PDE_SetOfDomains* SDOMS ;
      
      FE_SetOfParameters* PRMS ;
      
      FE_OneStepIteration* ONE_IT ;

      doubleVector GRAPHICS_TIMES ;
      double GRAPHICS_NEXT_TIME ;

      doubleVector SAVER_TIMES ;
      double SAVER_NEXT_TIME ;

      PEL_Timer* overall ;
      PEL_Timer* POST_TIMER ;
      PEL_Timer* SAVE_TIMER ;
      bool SAVEFG ;
      size_t LAST_MEM_IT ;

      bool VERBOSE ;
} ;


#endif

