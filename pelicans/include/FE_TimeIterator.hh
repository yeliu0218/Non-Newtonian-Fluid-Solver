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

#ifndef FE_TIME_ITERATOR_HH
#define FE_TIME_ITERATOR_HH

#include <PEL_Object.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>

class PEL_Data ;
class PEL_Double ;
class PEL_ModuleExplorer ;

/*
Time steppers.

PUBLISHED
*/

class PEL_EXPORT FE_TimeIterator : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_TimeIterator* create( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) ;

   //-- Characteristics

      // initial iteration number
      size_t initial_iteration_number( void ) const ;

      // lower bound of the time interval
      double initial_time( void ) const ;

      // upper bound of the time interval
      double final_time( void ) const ;

      // time step storage depth
      size_t storage_depth( void ) const ;

   //-- Iteration

      // Start time stepping.
      void start( void ) ;

      // Advance one time step.
      void go_next_time( void ) ;

      // Is time stepping started ?
      bool is_started( void ) const ;

      // Is time stepping finished ?
      bool is_finished( void ) const ;

   //-- Access

      // current time
      double time( void ) const ;

      // time step of the past `level'-th iteration (eg current time step if 
      // `level' is zero, time step of the previous iteration if 
      // `level' is 1)
      double time_step( size_t level = 0 ) const ;

      // current iteration number
      size_t iteration_number( void ) const ;

      bool just_went_back_and_forth( void ) const ;
      
      enum FinishedReason
      {
         UserCheckpoint         = 1,
         FinishIterationsCalled = 2,
         FinalTimeReached       = 3,
         
         NotFinished = 0
      } ;
      
      FinishedReason finished_reason( void ) const ;

   //-- Adaptation(60.)
      
      void finish_iterations( void ) ;

      bool time_step_is_fixed( void ) const ;
      
      void set_next_time_step( double dt ) ;

      double next_time_step( void ) const ;
      
      void go_back( void ) ;

      bool just_went_back( void ) const ;

    //-- Utilities(70.)

      // Is `time1' greater than or equal to `time2' ?
      static bool greater_or_equal( double time1, double time2 ) ;

      /* 
      Is `times' a valid table of times, in the sense that:
          - items are sorted by increasing values,
          - all items are greater than (or equal to) `::initial_time()' and 
            lower than (or equal to) `::final_time()' ?
      */
      bool table_of_times_is_valid( doubleVector const& times ) const ;

      // Raise a fatal error due to the invalid table `times', this table
      // being the data of keyword `keyword' in a hierarchical data structure
      // read by the class of name `class_name'.
      void raise_invalid_table_of_times( std::string const& class_name,
                                         std::string const& keyword,
                                         doubleVector const& times ) const ;

      // first item of `times' that is `::greater_or_equal()' to `::time'()
      double next_time_in_table( doubleVector const& times ) const ;

      void set_time_offset( double offset ) ;
      // Advance one time step.
      void go_next_time( double t_end ) ;
    //-- Persistence

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FE_TimeIterator( void ) ;
     ~FE_TimeIterator( void ) ;
      FE_TimeIterator( FE_TimeIterator const& other ) ;
      FE_TimeIterator& operator=( FE_TimeIterator const& other ) ;

      FE_TimeIterator( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Internals
      
      void set_time_step( void ) ;
      
      void start_checkpointing( void ) ;
      
      void test_checkpoint( void ) ;

   //-- Attributes
      
      bool STARTED ;
      FinishedReason FINISHED_REASON ;

      size_t ITER_INIT ;
      size_t ITER ;

      double T_INIT ;
      double T_START ;
      double T_END  ;
      double T ;

      PEL_Data* VARYING_DT ;
      PEL_Double* TIME ;
      bool RESTART_IT ;
      bool BACK_AND_FORTH_IT ;
      doubleVector  DT ;
      double INI_DT ;
      double NEXT_DT ;
      
      std::string CHECKPOINT_FILE ;
      
       
} ;

#endif
