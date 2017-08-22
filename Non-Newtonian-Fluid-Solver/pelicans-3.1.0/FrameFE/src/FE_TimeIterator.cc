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

#include <FE_TimeIterator.hh>

#include <PEL.hh>
#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>

#include <fstream>
#include <iostream>
#include <sstream>

using std::endl ;

//----------------------------------------------------------------------
FE_TimeIterator*
FE_TimeIterator:: create( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   FE_TimeIterator* result =  new FE_TimeIterator( a_owner,  exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ; 
   PEL_CHECK_POST( 
      result->initial_time() == exp->double_data( "time_initial" ) ) ;
   PEL_CHECK_POST( 
      result->final_time() == exp->double_data( "time_end" ) ) ;
   PEL_CHECK_POST(
      IMPLIES( !result->time_step_is_fixed(),
               result->time_step() == exp->double_data( "time_step" ) ) ) ;
   PEL_CHECK_POST(
      result->initial_iteration_number() == 1 ) ;
   PEL_CHECK_POST(
      IMPLIES( exp->has_entry( "storage_depth" ),
               result->storage_depth() == (size_t) exp->int_data( "storage_depth" ) ) ) ;
   PEL_CHECK_POST(
      IMPLIES( !exp->has_entry( "storage_depth" ),
	           result->storage_depth() == 5 ) ) ;
   PEL_CHECK_POST( 
      result->time() == result->initial_time() ) ;
   PEL_CHECK_POST( !result->is_started() ) ;
   PEL_CHECK_POST( result->time() == result->initial_time() ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
FE_TimeIterator:: FE_TimeIterator( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , STARTED( false )
   , FINISHED_REASON( NotFinished )
   , ITER_INIT( 1 )
   , ITER( PEL::bad_index() )
   , T_INIT( exp->double_data( "time_initial" ) )
   , T_START( exp->double_data( "time_initial" ) )
   , T_END( exp->double_data( "time_end" ) )
   , T( exp->double_data( "time_initial" ) )
   , VARYING_DT( 0 )
   , TIME( 0 )
   , RESTART_IT( false )
   , BACK_AND_FORTH_IT( false )
   , DT( doubleVector( 5 ) )
   , INI_DT( PEL::bad_double() )
   , NEXT_DT( PEL::bad_double() )
{
   PEL_LABEL( "FE_TimeIterator:: FE_TimeIterator" ) ;

   if( T_END < T_START )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "time_end",
         " a value greater than \"time_initial\" is expected" ) ;
   }
   if( exp->has_entry( "storage_depth" ) )
   {
      int s = exp->int_data( "storage_depth" ) ;
      if( s<=0 )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "storage_depth",
            "   a positive value is expected" ) ;
      }
      DT.re_initialize( (size_t) s ) ;
   }
   
   PEL_Data* cste_dt = exp->abstract_data( 0, "time_step" ) ;
   if( cste_dt->value_can_be_evaluated(0) )
   {
      double dt_ini = exp->double_data( "time_step" ) ;
      DT.set( dt_ini ) ;
   }
   else
   {
      PEL_ContextSimple* ctx = PEL_ContextSimple::create( this ) ;
      TIME = PEL_Double::create( ctx, T_INIT ) ;
      ctx->extend( PEL_Variable::object( "DS_T" ), TIME ) ;
      VARYING_DT = exp->abstract_data( this, "time_step", ctx ) ;
      if( !VARYING_DT->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
                    exp, "time_step", VARYING_DT->undefined_variables(0) ) ;
      }
      if( VARYING_DT->data_type() != PEL_Data::Double )
      {
         PEL_Error::object()->raise_bad_data_type(
                        exp, "time_step", PEL_Data::Double ) ;
      }
      DT.set( VARYING_DT->to_double() ) ;
   }
   cste_dt->destroy() ; cste_dt = 0 ;
   
   INI_DT = DT(0) ;
}

//----------------------------------------------------------------------
FE_TimeIterator:: ~FE_TimeIterator( void )
//----------------------------------------------------------------------
{
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   com->barrier() ;
   if( com->rank() == 0 )
   {
      PEL_System::erase( CHECKPOINT_FILE ) ;
   }
}

//----------------------------------------------------------------------
size_t
FE_TimeIterator:: initial_iteration_number( void ) const
//----------------------------------------------------------------------
{
   return( ITER_INIT ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: initial_time( void ) const
//----------------------------------------------------------------------
{
   return( T_INIT ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: final_time( void ) const
//----------------------------------------------------------------------
{
   return( T_END ) ;
}

//----------------------------------------------------------------------
size_t
FE_TimeIterator:: storage_depth( void ) const
//----------------------------------------------------------------------
{
   return( DT.size() ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: start" ) ;
   PEL_CHECK_PRE( !is_started() || is_finished() ) ;
   
   STARTED = true ;
   FINISHED_REASON = NotFinished ;
   ITER = ITER_INIT ;
   T = T_START ;
   set_time_step() ;
   start_checkpointing() ;

   PEL_CHECK_POST( is_started() ) ;
   PEL_CHECK_POST( !is_finished() ) ;
   PEL_CHECK_POST( iteration_number() == initial_iteration_number() ) ;
   PEL_CHECK_POST( FORMAL( time() == initial_time()+time_step() ) ) ;
   PEL_CHECK_POST( !just_went_back() ) ;
   PEL_CHECK_POST( !just_went_back_and_forth() ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: go_next_time( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: go_next_time" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   PEL_SAVEOLD( double, time, time() ) ;
   PEL_SAVEOLD( double, time_step, time_step() ) ;
   PEL_SAVEOLD( size_t, iteration_number, iteration_number() ) ;
   PEL_SAVEOLD( bool, just_went_back, just_went_back() ) ;
   
   test_checkpoint() ;
   if( ( FINISHED_REASON == NotFinished ) && ( T>T_END ) )
   {
      FINISHED_REASON = FinalTimeReached ;
   }
   
   if( FINISHED_REASON == NotFinished )
   {
      ITER++ ;
      set_time_step() ;
   }
   
   PEL_CHECK_POST( IMPLIES( is_finished(),
                            iteration_number() == OLD(iteration_number) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            iteration_number() == OLD(iteration_number)+1 ) ) ;
   PEL_CHECK_POST( IMPLIES( is_finished(), time() == OLD(time) ) ) ;
   PEL_CHECK_POST( IMPLIES( is_finished(),
                            time_step() == OLD(time_step) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            FORMAL( time() == OLD(time) + time_step() ) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            !just_went_back() ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            just_went_back_and_forth() == OLD(just_went_back) ) ) ;
   
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: is_started( void ) const
//----------------------------------------------------------------------
{
   return( STARTED ) ;
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: is_finished( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: is_finished" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   
   bool result = ( FINISHED_REASON != NotFinished ) ;
   
   PEL_CHECK_POST( EQUIVALENT( result, ( FINISHED_REASON != NotFinished ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: time( void ) const
//----------------------------------------------------------------------
{
   return( T ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: time_step( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: time_step" ) ;
   PEL_CHECK_PRE( level<storage_depth() ) ;
   return( DT(level) ) ;
}

//----------------------------------------------------------------------
size_t
FE_TimeIterator:: iteration_number( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: iteration_number" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   return( ITER ) ;
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: just_went_back_and_forth( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: just_went_back_and_forth" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   return( BACK_AND_FORTH_IT ) ;
}

//----------------------------------------------------------------------
FE_TimeIterator::FinishedReason
FE_TimeIterator:: finished_reason( void ) const
//----------------------------------------------------------------------
{
   return( FINISHED_REASON ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: finish_iterations( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: finish_iterations" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   PEL_CHECK_PRE( !is_finished() ) ;
   
   FINISHED_REASON = FinishIterationsCalled ;
   
   PEL_CHECK_POST( is_finished() ) ;
   PEL_CHECK_POST( finished_reason() == FinishIterationsCalled ) ;
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: time_step_is_fixed( void ) const
//----------------------------------------------------------------------
{
   return( VARYING_DT != 0 ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: set_next_time_step( double dt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: set_next_time_step" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   PEL_CHECK_PRE( !is_finished() ) ;
   PEL_CHECK_PRE( !time_step_is_fixed() ) ;
   PEL_CHECK_PRE( dt > 0. ) ;

   NEXT_DT = ( dt == DT(0) ? PEL::bad_double() : dt ) ;

   PEL_CHECK_POST( next_time_step() == dt ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: next_time_step( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: next_time_step" ) ;
   double result = ( NEXT_DT == PEL::bad_double() ? DT(0) : NEXT_DT ) ;
   PEL_CHECK_POST( result > 0. ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: go_back( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: go_back" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   PEL_CHECK_PRE( !just_went_back() ) ;
   PEL_CHECK_PRE( !is_finished() ) ;
   PEL_CHECK_PRE( !time_step_is_fixed() ) ;
   PEL_SAVEOLD( double, time_step, time_step() ) ;
   
   T = T-DT(0) ;
   RESTART_IT = true ;
 
   PEL_CHECK_POST( !is_finished() ) ;
   PEL_CHECK_POST( time_step() == OLD(time_step) ) ;
   PEL_CHECK_POST( FORMAL( time() == OLD(time) - time_step() ) ) ;
   PEL_CHECK_POST( just_went_back() ) ;
}


//----------------------------------------------------------------------
bool
FE_TimeIterator:: just_went_back( void ) const
//----------------------------------------------------------------------
{
   return( RESTART_IT ) ;
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: greater_or_equal( double time1, double time2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: greater_or_equal" ) ;

   double const epsilon = ( time2>0. ? 1.E-8 : -1.E-8 ) ;
   bool const result = ( time1 >= time2*( 1.0 - epsilon ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
FE_TimeIterator:: table_of_times_is_valid( doubleVector const& times ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: table_of_times_is_valid" ) ;

   bool result = true ;
   size_t const n = times.size() ;
   if( n>0 )
   {
      result = greater_or_equal( times(0), initial_time() )
            && greater_or_equal( final_time(), times(n-1) ) ;
      for( size_t i=0 ; result && i<n-1 ; ++i )
      {
         result = ( times(i+1)>times(i) ) ; // increasing values
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: raise_invalid_table_of_times( 
                                         std::string const& class_name,
                                         std::string const& keyword,
                                         doubleVector const& times ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: raise_invalid_table_of_times" ) ;
   PEL_CHECK_PRE( !table_of_times_is_valid( times ) ) ;

   std::ostringstream mesg ;
   mesg << "*** " << class_name << " error:" << endl ;
   mesg << "    invalid data of keyword \""<< keyword <<"\"" << endl ;
   mesg << "    the time sequence must " << endl ;
   mesg << "      - start after the initial time (which is: " 
        << initial_time() << ")" << endl ;
   mesg << "      - stop before the final time (which is: " 
        << final_time() << ")" << endl ;
   mesg << "      - be sorted by increasing values" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIterator:: next_time_in_table( doubleVector const& times ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: next_time_in_table" ) ;
   PEL_CHECK_PRE( table_of_times_is_valid( times ) ) ;
   
   double result = PEL::max_double() ;
   for( size_t i=0 ; i<times.size() ; i++ )
   {
      if( !greater_or_equal( time(), times(i) ) )
      {
         result = times(i) ;
         break ;
      }
   }
   
   PEL_CHECK_POST( !greater_or_equal( time(), result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: set_time_offset( double offset )
//----------------------------------------------------------------------
{
	T += offset ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: go_next_time( double t_local_end )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: go_next_time" ) ;
   PEL_CHECK_PRE( is_started() ) ;
   PEL_SAVEOLD( double, time, time() ) ;
   PEL_SAVEOLD( double, time_step, time_step() ) ;
   PEL_SAVEOLD( size_t, iteration_number, iteration_number() ) ;

   if( ( FINISHED_REASON == NotFinished ) && ( T>t_local_end ) )
   {
      FINISHED_REASON = FinalTimeReached ;
   }
   else
   {
      ITER++ ;
      set_time_step() ;
   }

   PEL_CHECK_POST( IMPLIES( is_finished(),
                            iteration_number() == OLD(iteration_number) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            iteration_number() == OLD(iteration_number)+1 ) ) ;
   PEL_CHECK_POST( IMPLIES( is_finished(), time() == OLD(time) ) ) ;
   PEL_CHECK_POST( IMPLIES( is_finished(),
                            time_step() == OLD(time_step) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(),
                            FORMAL( time() == OLD(time) + time_step() ) ) ) ;
   PEL_CHECK_POST( IMPLIES( !is_finished(), !just_went_back() ) ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;
   
   writer->start_new_object( "FE_TimeIterator" ) ;

   // Saving current time
   writer->add_entry( "t", PEL_Double::create( 0, T ) ) ;
   
   // Saving current iter
   writer->add_entry( "iter", PEL_Int::create( 0, ITER ) ) ;
   
   // Initial time step
   writer->add_entry( "initial_dt", PEL_Double::create( 0, INI_DT ) ) ;

   // Time step has been modified
   if( NEXT_DT != PEL::bad_double() )
   {
      writer->add_entry( "new_dt", PEL_Double::create( 0, NEXT_DT ) ) ;
   }
   
   // Restart
   if( RESTART_IT )
   {
     writer->add_entry( "just_went_back", PEL_Bool::create( 0, true ) ) ;
   }
   
   // Restart
   if( BACK_AND_FORTH_IT )
   {
     writer->add_entry( "just_went_back_and_forth", PEL_Bool::create( 0, true ) ) ;
   }
   
   // Saving current time step
   writer->add_entry( "dt", PEL_DoubleVector::create( 0, DT ) ) ;
   
   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FE_TimeIterator" ) ;
   
   T = reader->data_of_entry( "t" )->to_double() ;
   T_START = T ;
   ITER_INIT = reader->data_of_entry( "iter" )->to_int()+1 ;
   NEXT_DT = PEL::bad_double() ;
   if( !VARYING_DT )
   {
      double ini_dt = reader->data_of_entry( "initial_dt" )->to_double() ;
      if( ini_dt != INI_DT ) // time step is changed in the data deck
      {
         NEXT_DT = INI_DT ;
      }
      else if( reader->has_entry( "new_dt" ) )
      {
         NEXT_DT = reader->data_of_entry( "new_dt" )->to_double() ;
      }
   }
   RESTART_IT = reader->has_entry( "just_went_back" ) ;
   BACK_AND_FORTH_IT = reader->has_entry( "just_went_back_and_forth" ) ;
   DT = reader->data_of_entry( "dt" )->to_double_vector() ;
   ITER = PEL::bad_index() ;
   STARTED = false ;
   FINISHED_REASON = NotFinished ;
   reader->end_object_retrieval() ;
   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIterator:: set_time_step( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: set_time_step" ) ;

   if( !RESTART_IT )
   {
      for( size_t i=DT.size()-1 ; i>0 ; --i )
      {
         DT( i ) = DT( i-1 ) ;
      }
   }
   if( NEXT_DT != PEL::bad_double() )
   {
      DT(0) = NEXT_DT ;
      NEXT_DT =  PEL::bad_double() ;
   }
   else if( VARYING_DT != 0 )
   {
      TIME->set( T ) ;
      DT(0) = VARYING_DT->to_double() ;
   }
   T += DT(0) ;
   BACK_AND_FORTH_IT = RESTART_IT ;
   RESTART_IT = false ;
}

//------------------------------------------------------------------------
void
FE_TimeIterator:: start_checkpointing( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: start_checkpointing" ) ;
   
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   size_t rank = com->rank() ;
   int pid ;
   if( rank == 0 )
   {
      pid = PEL_System::process_id() ;
   }
   com->broadcast( pid, 0 ) ;
   std::ostringstream nstr ;
   nstr << "remove_to_stop_" << pid ;
   CHECKPOINT_FILE = nstr.str() ;
   if( rank == 0 )
   {
      std::ofstream file( CHECKPOINT_FILE.c_str(), std::ios::trunc ) ;
      file.close() ;
   }
}

//------------------------------------------------------------------------
void
FE_TimeIterator:: test_checkpoint( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator:: test_checkpoint" ) ;
   
   std::ifstream file( CHECKPOINT_FILE.c_str() ) ;
   if( !file ) FINISHED_REASON = UserCheckpoint ;
}
