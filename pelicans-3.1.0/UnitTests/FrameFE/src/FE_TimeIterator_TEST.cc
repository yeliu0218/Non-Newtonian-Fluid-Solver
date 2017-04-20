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

#include <FE_TimeIterator_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <intVector.hh>

#include <FE_TimeIterator.hh>

#include <iostream>

FE_TimeIterator_TEST const*
FE_TimeIterator_TEST::PROTOTYPE = new FE_TimeIterator_TEST() ;

//-------------------------------------------------------------------------
FE_TimeIterator_TEST:: FE_TimeIterator_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "FE_TimeIterator", "FE_TimeIterator_TEST" )
{
}

//-------------------------------------------------------------------------
FE_TimeIterator_TEST:: ~FE_TimeIterator_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_TimeIterator_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIterator_TEST:: process_one_test" ) ;

   out() << "| ... " << exp->name() << " : " << std::endl ;
   bool ok = false ;
   
   PEL_ModuleExplorer const* e =
                          exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   FE_TimeIterator* t_it = FE_TimeIterator::create( 0, e ) ;
   e->destroy() ; e = 0 ;

   PEL_ModuleExplorer const* r = exp->create_subexplorer( 0, "Result" ) ;
   double d_eps = r->double_data( "dbl_epsilon" ) ;
   double d_min = r->double_data( "dbl_minimum" ) ;

   // RESTART ONE ITERATION
   size_t restart_it = PEL::bad_index() ;
   double restart_dtime = PEL::bad_double() ;
   if( exp->has_module( "Restart_one_iteration" ) )
   {
      PEL_ModuleExplorer const* rest =
                    exp->create_subexplorer( 0, "Restart_one_iteration" ) ;
      restart_it = (size_t) rest->int_data( "restart_iteration" ) ;
      restart_dtime = rest->double_data( "restart_time_step" ) ;
      rest->destroy() ; rest = 0 ;
   }

   // FINISH A ONE ITERATION
   size_t finish_it = PEL::bad_index() ;
   if( exp->has_module( "Finish_at_iteration" ) )
   {
      PEL_ModuleExplorer const* rest =
                    exp->create_subexplorer( 0, "Finish_at_iteration" ) ;
      finish_it = (size_t) rest->int_data( "finish_iteration" ) ;
      rest->destroy() ; rest = 0 ;
   }

   // INITIAL TIME
   ok = ( PEL::double_equality( t_it->initial_time(),
                                r->double_data( "time_initial" ),
                                d_eps, d_min ) ) ;
   if( !ok )
   {
      out() << "initial time :   expected="
            << r->double_data( "time_initial" ) << std::endl ;
      out() << "                 computed="
            << t_it->initial_time() << std::endl ;
   }
   notify_one_test_result( "initial_time", ok ) ;

   // FINAL TIME
   ok = ( PEL::double_equality( t_it->final_time(),
                                r->double_data( "time_end" ),
                                d_eps, d_min ) ) ;
   if( !ok )
   {
      out() << "final time :   expected="
            << r->double_data( "time_end" ) << std::endl ;
      out() << "               computed="
            << t_it->final_time() << std::endl ;
   }
   notify_one_test_result( "final_time", ok ) ;

   // INITIAL ITERATION NUMBER
   ok = ( (int) t_it->initial_iteration_number() ==
                          r->int_data( "initial_iteration_number" ) ) ;
   if( !ok )
   {
      out() << "initial iteration number:  expected="
            << r->int_data( "initial_iteration_number" ) << std::endl ;
      out() << "                           computed="
            << t_it->initial_iteration_number() << std::endl ;
   }
   notify_one_test_result( "initial_iteration_number", ok ) ;

   // IS FIXED
   ok = ( t_it->time_step_is_fixed() == r->bool_data( "time_step_is_fixed" ) ) ;
   notify_one_test_result( "time_step_is_fixed", ok ) ;
   

   // TIME STEPPING
   doubleVector times(0) ;
   doubleVector dt(0) ;
   doubleVector dt1(0) ;
   doubleVector dt2(0) ;
   doubleVector dt3(0) ;
   intVector it(0) ;
   for( t_it->start() ; !t_it->is_finished() ; t_it->go_next_time() )
   {
      times.append( t_it->time() ) ;
      dt.append(  t_it->time_step(0) ) ;
      dt1.append( t_it->time_step(1) ) ;
      dt2.append( t_it->time_step(2) ) ;
      dt3.append( t_it->time_step(3) ) ;
      it.append( (int) t_it->iteration_number() ) ;
      if( t_it->iteration_number() == restart_it )
      {
         t_it->go_back() ;
         t_it->set_next_time_step( restart_dtime ) ;
      }
      if( t_it->iteration_number() == finish_it )
      {
         t_it->finish_iterations() ;
      }
   }
   ok = ( doubleVector_equality( times, r->doubleVector_data( "times" ),
                                 d_eps, d_min ) ) ;
   if( !ok )
   {
      out() << "times:" << std::endl ;
      out() << "   expected:" << r->doubleVector_data( "times" ) << std::endl ;
      out() << "   computed:" << times << std::endl ;
   }
   notify_one_test_result( "times", ok ) ;
   ok = ( doubleVector_equality( dt, r->doubleVector_data( "time_steps" ),
                                 d_eps, d_min ) ) ;
   if( !ok )
   {
      out() << "time_steps:" << std::endl ;
      out() << "   expected:" << r->doubleVector_data( "time_steps" ) << std::endl ;
      out() << "   computed:" << dt << std::endl ;
   }
   notify_one_test_result( "time_steps", ok ) ;
   if( r->has_entry( "time_steps_level_1" ) )
   {
      ok = ( doubleVector_equality(
                dt1, r->doubleVector_data( "time_steps_level_1" ),
                d_eps, d_min ) ) ;
      if( !ok )
      {
         out() << "time_steps_level_1:" << std::endl ;
         out() << "   expected:" << r->doubleVector_data( "time_steps_level_1" ) << std::endl ;
         out() << "   computed:" << dt1 << std::endl ;
      }
      notify_one_test_result( "time_steps_level_1", ok ) ;
   }
   if( r->has_entry( "time_steps_level_2" ) )
   {
      ok = ( doubleVector_equality(
                dt2, r->doubleVector_data( "time_steps_level_2" ),
                d_eps, d_min ) ) ;
      if( !ok )
      {
         out() << "time_steps_level_2:" << std::endl ;
         out() << "   expected:" << r->doubleVector_data( "time_steps_level_2" ) << std::endl ;
         out() << "   computed:" << dt2 << std::endl ;
      }
      notify_one_test_result( "time_steps_level_2", ok ) ;
   }
   if( r->has_entry( "time_steps_level_3" ) )
   {
      ok = ( doubleVector_equality(
                dt3, r->doubleVector_data( "time_steps_level_3" ),
                d_eps, d_min ) ) ;
      if( !ok )
      {
         out() << "time_steps_level_3:" << std::endl ;
         out() << "   expected:" << r->doubleVector_data( "time_steps_level_3" ) << std::endl ;
         out() << "   computed:" << dt3 << std::endl ;
      }
      notify_one_test_result( "time_steps_level_3", ok ) ;
   }
   
   ok = ( it == r->intVector_data( "iterations" ) ) ;
   if( !ok )
   {
      out() << "iterations:" << std::endl ;
      out() << "   expected:" << r->intVector_data( "iterations" ) << std::endl ;
      out() << "   computed:" << it << std::endl ;
   }
   notify_one_test_result( "iterations", ok ) ;

   r->destroy() ; r = 0 ;
   t_it->destroy() ; t_it = 0 ;
}

//-------------------------------------------------------------------------
bool
FE_TimeIterator_TEST:: doubleVector_equality(
                                  doubleVector const& v1,
                                  doubleVector const& v2,
                                  double d_eps, double d_min ) const
//-------------------------------------------------------------------------
{
   bool result = ( v1.size() == v2.size() ) ;
   for( size_t i=0 ; result && i<v1.size() ; ++i )
   {
      result = PEL::double_equality( v1(i), v2(i), d_eps, d_min ) ;
   }
   return( result ) ;
}
