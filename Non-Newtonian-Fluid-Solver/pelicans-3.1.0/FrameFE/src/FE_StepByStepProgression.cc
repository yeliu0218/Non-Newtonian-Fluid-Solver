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

#include <FE_StepByStepProgression.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_MemoryTracer.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <PEL_Root.hh>
#include <PEL_Timer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDomains.hh>

#include <FE.hh>
#include <FE_OneStepIteration.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIteratorAdapter.hh>
#include <FE_TimeIterator.hh>

#include <intVector.hh>
#include <doubleVector.hh>

#include <iomanip>
#include <iostream>
#include <sstream>

using std::endl ;
using std::ostringstream ;
using std::string ;

size_t FE_StepByStepProgression::graphics_level = 0 ;

FE_StepByStepProgression const* 
FE_StepByStepProgression:: PROTOTYPE = new FE_StepByStepProgression() ;

//----------------------------------------------------------------------
FE_StepByStepProgression:: FE_StepByStepProgression( void )
//----------------------------------------------------------------------
   : PEL_Application( "FE_StepByStepProgression" )
   , TIME_ADAPT( 0 )
   , TIME_IT( 0 )
   , ONE_IT( 0 )
   , GRAPHICS_TIMES( 0 )
   , GRAPHICS_NEXT_TIME( PEL::bad_double() )
   , SAVER_TIMES( 0 )
   , SAVER_NEXT_TIME( PEL::bad_double() )
   , overall( 0 )
   , POST_TIMER( 0 )
   , SAVE_TIMER( 0 )
   , SAVEFG( true )
   , VERBOSE( true )
{
}

//----------------------------------------------------------------------
FE_StepByStepProgression*
FE_StepByStepProgression:: create( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   FE_StepByStepProgression* result = 
                              new FE_StepByStepProgression( a_owner, exp ) ;
   result->register_storable_objects() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_StepByStepProgression*
FE_StepByStepProgression:: create_replica( 
                                      PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   FE_StepByStepProgression* result = 
                              new FE_StepByStepProgression( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_StepByStepProgression:: FE_StepByStepProgression( 
                                            PEL_Object* a_owner,
                                            PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , TIME_ADAPT( 0 )
   , TIME_IT( 0 )
   , DOM( 0 )
   , SDOMS( 0 )
   , PRMS( 0 )
   , ONE_IT( 0 )
   , GRAPHICS_TIMES( 0 )
   , GRAPHICS_NEXT_TIME( PEL::max_double() )
   , SAVER_TIMES( 0 )
   , SAVER_NEXT_TIME( PEL::max_double() )
   , overall( PEL_Timer::create( this ) )
   , POST_TIMER( PEL_Timer::create( this ) )
   , SAVE_TIMER( PEL_Timer::create( this ) )
   , SAVEFG( true )
   , LAST_MEM_IT( PEL::bad_index() )
   , VERBOSE( true )
{
   PEL_LABEL( "FE_StepByStepProgression:: FE_StepByStepProgression" ) ;

   if( exp->has_entry( "verbose" ) ) VERBOSE = exp->bool_data( "verbose" ) ;
   
   if( exp->has_module( "memory_trace" ) )
   {
      PEL_MemoryTracer::object()->enable_memory_trace() ;
      PEL_ModuleExplorer* ee =
                       exp->create_subexplorer( 0, "memory_trace" ) ;
      std::string const& type = ee->string_data( "type" ) ;
      ee->test_data_in( "type", "trace_all_iterations,trace_first_iterations" ) ;
      if( type == "trace_first_iterations" )
      {
         LAST_MEM_IT = ee->int_data( "last_iteration_checked" ) ;
      }
      ee->destroy() ; ee = 0 ;
   }
   
   if( VERBOSE ) display_memory_usage( PEL::out(), 4 ) ;

   if( PEL_MemoryTracer::object()->memory_trace_enabled() )
   {
      PEL_MemoryTracer::object()->trace(
                   PEL_MemoryTracer::object()->message( "Start program" ) ) ;
   }
   
   std::string geom = "cartesian" ;
   if( exp->has_entry( "geometry" ) )
   {
      geom = exp->string_data( "geometry" ) ;
      exp->test_data_in( "geometry", "cartesian,axisymmetrical" ) ;
      exp->set_default( "geometry", "cartesian" ) ;
   }
   if( geom == "cartesian" )
   {
      FE::set_geometry( FE::cartesian ) ;
   }
   else if( geom == "axisymmetrical" )
   {
      FE::set_geometry( FE::axisymmetrical ) ;
   }
   if( VERBOSE )
   {   
      PEL::out() << "*** Building geometry" << std::endl ;
      PEL::out() << "    geometry : \"" << geom << "\"" << std::endl ;
      PEL::out() << std::endl ;
   }

   // PDE_DomainAndFields:
   if( exp->has_module( "PDE_DomainAndFields" ) )
   {
      PEL_ModuleExplorer* ee =
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
      PEL_MemoryTracer::object()->start_event( "Building PDE_DomainAndFields" ) ;
      DOM = PDE_DomainAndFields::create( this, ee, PEL_Exec::communicator() ) ;
      PEL_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;
      
      if( geom == "axisymmetrical" && DOM->nb_space_dimensions() != 2 )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "geometry",
            "\"axisymmetrical\" is only allowed in 2D geometry" ) ;
      }
   }

   // PDE_SetOfDomains:
   if( exp->has_module( "PDE_SetOfDomains" ) )
   {
      PEL_ModuleExplorer* ee =
                         exp->create_subexplorer( 0, "PDE_SetOfDomains" ) ;
      PEL_MemoryTracer::object()->start_event( "Building PDE_SetOfDomains" ) ;
      SDOMS = PDE_SetOfDomains::create( this, ee ) ;
      PEL_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;

      if( geom == "axisymmetrical" )
         for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
            if( SDOMS->domain( i )->nb_space_dimensions() != 2 )
               PEL_Error::object()->raise_bad_data_value(
                  exp, "geometry",
                  "\"axisymmetrical\" is only allowed in 2D geometry" ) ;
   }
   if( ( DOM == 0 && SDOMS == 0 ) || ( DOM != 0 && SDOMS != 0 ) )
   {
      PEL_Error::object()->raise_module_error(
         exp,
         "A sub-module \"PDE_DomainAndFields\"\n"
         "          or \"PDE_SetOfDomains\"\n"
         "is expected" ) ;
   }

   if( VERBOSE ) display_memory_usage( PEL::out(), 4 ) ;

   // FE_SetOfParameters:
   {
      PEL_ModuleExplorer* ee =
                     exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
      PEL_MemoryTracer::object()->start_event( "Building FE_SetOfParameters" ) ;
      PRMS = 
         ( DOM != 0 ? FE_SetOfParameters::create( this, DOM, ee )
                    : FE_SetOfParameters::create( this, SDOMS, ee ) ) ;
      PEL_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;
   }
   
   // FE_OneStepIteration:
   {
      if( VERBOSE ) PEL::out() << "*** building FE_OneStepIteration" << endl ;
      PEL_ModuleExplorer* ee =
                         exp->create_subexplorer( 0, "FE_OneStepIteration" ) ;
      PEL_MemoryTracer::object()->start_event( "Building FE_OneStepIteration" ) ;
      ONE_IT = 
         ( DOM != 0 ? FE_OneStepIteration::make( this, DOM, PRMS,  ee )
                    : FE_OneStepIteration::make( this, SDOMS, PRMS,  ee ) ) ;
      PEL_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;
      FE_OneStepIteration::reset_standard_times() ;
   }

   // FE_TimeIterator:
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
      TIME_IT = FE_TimeIterator::create( this, ee ) ;
      ee->destroy() ; ee = 0 ;
   
      if( exp->has_module( "FE_TimeIteratorAdapter" ) )
      {
         ee = exp->create_subexplorer( 0, "FE_TimeIteratorAdapter" ) ;
         TIME_ADAPT =
            ( DOM != 0 ? FE_TimeIteratorAdapter::make( TIME_IT, DOM, PRMS, ee )
                       : FE_TimeIteratorAdapter::make( TIME_IT, SDOMS, PRMS, ee ) ) ;
         ee->destroy() ; ee = 0 ;
      }
      else
      {
         TIME_ADAPT =
            ( DOM != 0 ? FE_TimeIteratorAdapter::make_default( TIME_IT,   DOM, PRMS )
                       : FE_TimeIteratorAdapter::make_default( TIME_IT, SDOMS, PRMS ) ) ;
      }
   }

   if( exp->has_entry( "graphics_output_times" ) )
   {
      GRAPHICS_TIMES = exp->doubleVector_data( "graphics_output_times" ) ;
   }

   if( exp->has_entry( "state_saving_times" ) )
   {
      SAVER_TIMES = exp->doubleVector_data( "state_saving_times" ) ;
   }
   if( exp->has_entry( "save_grid_and_fields_for_postprocessing" ) )
   {
      SAVEFG = exp->bool_data( "save_grid_and_fields_for_postprocessing" ) ;
   }

   if( VERBOSE && PEL_Exec::communicator()->rank() == 0 )
   {
      ONE_IT->print( PEL::out(), 4 ) ;
   }

   overall->start() ;
}

//----------------------------------------------------------------------
FE_StepByStepProgression:: ~FE_StepByStepProgression( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
FE_TimeIterator const* 
FE_StepByStepProgression:: time_iterator( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: time_iterator" ) ;

   FE_TimeIterator const* result = TIME_IT ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters const* 
FE_StepByStepProgression:: set_of_parameters( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: set_of_parameters" ) ;

   FE_SetOfParameters const* result = PRMS ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DomainAndFields* 
FE_StepByStepProgression:: domain_and_fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: domain_and_fields" ) ;

   PDE_DomainAndFields* result = DOM ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: display_memory_usage(
                                 std::ostream& os, size_t indent_width )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: display_memory_usage" ) ;
   PEL_CHECK_PRE( os ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "Memory usage: " ;
   PEL_MemoryTracer::display_memory( os, PEL_MemoryTracer::used_memory() ) ;
   os << endl
      << s << "Number of objects: "
      << PEL_Object::GetNumberOf_PEL_objects()
      << endl ;   
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: add_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: add_storable_objects" ) ;
   list->extend( TIME_ADAPT ) ;
   if( DOM != 0 )
   {
      list->extend( DOM ) ;
   }
   else
   {
      for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
      {
         list->extend( SDOMS->domain( i ) ) ;
      }
      for( size_t i=0 ; i<SDOMS->nb_interfaces() ; ++i )
      {
         list->extend( SDOMS->interface( i ) ) ;
      }
   }
   ONE_IT->register_storable_objects( list ) ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: run" ) ; 

   check_times( "graphics_output_times", GRAPHICS_TIMES ) ;
   check_times( "state_saving_times", SAVER_TIMES ) ;
   
   GRAPHICS_NEXT_TIME = TIME_IT->next_time_in_table( GRAPHICS_TIMES ) ;
   SAVER_NEXT_TIME = TIME_IT->next_time_in_table( SAVER_TIMES ) ;
   
   PEL::out() << endl << "++++++ TIME STEPPING " 
              << " *** INITIAL TIME = " << TIME_IT->initial_time()
              << " *** FINAL TIME = " << TIME_IT->final_time()
              << " ++++++" << endl << endl ;

   PEL_MemoryTracer::object()->start_event( "Time step initialization" ) ;
   ONE_IT->do_before_time_stepping( TIME_IT ) ;
   PEL_MemoryTracer::object()->stop_event() ;

   save_for_post_processing( true ) ;

   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      if( VERBOSE )
      {
         std::streamsize p = PEL::out().precision() ;
         PEL::out() << std::setprecision( 12 ) << endl ;;
         PEL::out() << "++++++ ITERATION = " << TIME_IT->iteration_number()
                    << " *** TIME = " << TIME_IT->time()
                    << " *** TIME STEP = " << TIME_IT->time_step()
                    << " ++++++" << endl << endl ;
         PEL::out() << std::setprecision( p ) ;
      }

      std::ostringstream mem_label ;
      mem_label << "Time iteration " << TIME_IT->iteration_number() ;
      PEL_MemoryTracer::object()->start_event( mem_label.str() ) ;

      TIME_ADAPT->initialize_time_step() ;
      ONE_IT->do_before_inner_iterations_stage( TIME_IT ) ;
      ONE_IT->do_inner_iterations_stage( TIME_IT ) ;
      if( ONE_IT->inner_iterations_stage_failed() )
      {
         TIME_ADAPT->set_time_iteration_failed() ;
         ONE_IT->reset_after_inner_iterations_stage_failure() ;
      }
      else
      {
         ONE_IT->do_after_inner_iterations_stage( TIME_IT ) ;
         ONE_IT->adapt_time_iterator( TIME_ADAPT ) ;
         TIME_ADAPT->adapt_time_iterator() ;
         if( !TIME_IT->just_went_back() )
         {
            ONE_IT->do_after_time_adaptation( TIME_IT ) ;
            save_for_post_processing( false ) ;
            save_for_restart() ;
         }
      }

      if( VERBOSE )
      {
         PEL::out() << std::endl ;
         display_memory_usage( PEL::out(), 3 ) ;
      }

      PEL_MemoryTracer::object()->stop_event() ;
      if( PEL_MemoryTracer::object()->memory_trace_enabled() &&
          TIME_IT->iteration_number() == LAST_MEM_IT )
      {
         PEL_MemoryTracer::object()->disable_memory_trace() ;
         PEL::out() << std::endl ;
      }
   }

   if( VERBOSE )
   {
      PEL::out() << endl ;
      if( TIME_IT->finished_reason() == FE_TimeIterator::UserCheckpoint )
         PEL::out() << "++++++ TIME STEPPING ENDED BY USER " ;
      else
         PEL::out() << "++++++ TIME STEPPING COMPLETED " ;
      PEL::out() << " ++++++" << endl << endl ;
   }

   ONE_IT->do_after_time_stepping() ;

   save_for_post_processing( true ) ;
   
   overall->stop() ;
   
   if( !VERBOSE ) 
   {
      PEL::out() << endl << "++++++ " << TIME_IT->iteration_number() 
                 << " ITERATIONS ++++++" << endl ;
   }
   
   display_timers( 0 ) ;
   PEL::out() << endl ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: save_for_post_processing( bool force )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: save_for_post_processing" ) ;
   POST_TIMER->start() ;
   
   static size_t last_iter_saved = PEL::bad_index() ;
   size_t const iter =
      TIME_IT->is_started() ? TIME_IT->iteration_number() : 0 ;
   bool const save_iteration =
      force || FE_TimeIterator::greater_or_equal( TIME_IT->time(),
                                                  GRAPHICS_NEXT_TIME ) ;
   
   if( save_iteration && iter!=last_iter_saved )
   {
      if( DOM != 0 )
      {
         do_save_for_post( DOM->result_saver() ) ;
         display_timers( 3 ) ;
      }
      else
      {
         for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
         {
            do_save_for_post( SDOMS->domain( i )->result_saver() )  ;
            display_timers( 3 ) ;
         }
      }
      last_iter_saved = iter ;
      GRAPHICS_NEXT_TIME = TIME_IT->next_time_in_table( GRAPHICS_TIMES ) ;
   }
   POST_TIMER->stop() ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: display_timers( size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: display_timers" ) ;
   
   std::string space( indent_width, ' ' ) ;
   PEL::out() << endl << space << "++++++ TIMERS  ++++++" << endl << endl ;
   PEL::out() << space << "total time : " ;
   overall->print( PEL::out(), 0 ) ; PEL::out() << endl ;
   PEL::out() << space << "   including \"save_for_post_processing\" : " ;
   POST_TIMER->print( PEL::out(), 0 ) ; PEL::out() << endl ;
   PEL::out() << space << "   including \"save_for_restart\" : " ;
   SAVE_TIMER->print( PEL::out(), 0 ) ; PEL::out() << endl ;
   FE_OneStepIteration::print_standard_times( PEL::out(), indent_width ) ;
   ONE_IT->print_additional_times( PEL::out(), indent_width ) ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: do_save_for_post( PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   rs->start_cycle() ;

   std::streamsize p = PEL::out().precision() ;
   PEL::out() << std::setprecision(12) ;
   PEL::out() << endl
              << "   +++ SAVE FOR POSTPROCESSING"
              << " *** CYCLE = " << rs->cycle_number()
              << " *** TIME = "  << TIME_IT->time()
              << " ++++++" << endl << endl ;
   PEL::out() << std::setprecision(p) ;

   if( SAVEFG ) if( !TIME_IT->is_started() ) rs->save_grid() ;
   ONE_IT->do_additional_savings( TIME_IT, rs ) ;
   if( SAVEFG ) rs->save_fields( graphics_level ) ;
   rs->save_time( TIME_IT->time() ) ;

   rs->terminate_cycle() ;
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: check_times( std::string const& list_name,
                                        doubleVector const& dates ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StepByStepProgression:: check_times" ) ;

   if( !TIME_IT->table_of_times_is_valid( dates ) )
   {
      TIME_IT->raise_invalid_table_of_times( "FE_StepByStepProgression",
                                             list_name,
                                             dates ) ;
   }
}

//----------------------------------------------------------------------
void
FE_StepByStepProgression:: save_for_restart( void )
//----------------------------------------------------------------------
{
   SAVE_TIMER->start() ;
   
   bool save_iteration =
      FE_TimeIterator::greater_or_equal( TIME_IT->time(), SAVER_NEXT_TIME ) ;
   if( save_iteration )
   {
      write_storable_objects() ;
      SAVER_NEXT_TIME = TIME_IT->next_time_in_table( SAVER_TIMES ) ;
   }
   SAVE_TIMER->stop() ;
}

//----------------------------------------------------------------------
bool
FE_StepByStepProgression:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Application::invariant() ) ;

   return( true ) ;
}



