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

#include <FE_TimeIteratorAdapter.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <FE_TimeIterator.hh>
#include <FE_SetOfParameters.hh>

FE_TimeIteratorAdapter const* FE_TimeIteratorAdapter::DEFAULT_PROTOTYPE = 0 ;

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: make( FE_TimeIterator* t_it,
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: make(single domain)" ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   PEL_CHECK_PRE( !t_it->is_started() ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string const& name = exp->string_data( "concrete_name" ) ;
   FE_TimeIteratorAdapter const* proto =
      static_cast<FE_TimeIteratorAdapter const*>(
                                            plugins_map()->item( name ) ) ;
      
   FE_TimeIteratorAdapter* result =
                            proto->create_replica( t_it, dom, prms, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == t_it ) ;
   PEL_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: make( FE_TimeIterator* t_it,
                         PDE_SetOfDomains const* sdoms,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: make(multi domain)" ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   PEL_CHECK_PRE( !t_it->is_started() ) ;
   PEL_CHECK_PRE( sdoms != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string const& name = exp->string_data( "concrete_name" ) ;
   FE_TimeIteratorAdapter const* proto =
      static_cast<FE_TimeIteratorAdapter const*>(
                                            plugins_map()->item( name ) ) ;
      
   FE_TimeIteratorAdapter* result =
                          proto->create_replica( t_it, sdoms, prms, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == t_it ) ;
   PEL_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: make_default( FE_TimeIterator* t_it,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: make_default(single domain)" ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   PEL_CHECK_PRE( !t_it->is_started() ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   if( DEFAULT_PROTOTYPE != 0 )
   {
      result = DEFAULT_PROTOTYPE->create_replica( t_it, dom, prms ) ;
   }
   else
   {
      result = new FE_TimeIteratorAdapter( t_it ) ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == t_it ) ;
   PEL_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: make_default( FE_TimeIterator* t_it,
                                 PDE_SetOfDomains const* sdoms,
                                 FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: make_default(multi domain)" ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   PEL_CHECK_PRE( !t_it->is_started() ) ;
   PEL_CHECK_PRE( sdoms != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   if( DEFAULT_PROTOTYPE != 0 )
   {
      result = DEFAULT_PROTOTYPE->create_replica( t_it, sdoms, prms ) ;
   }
   else
   {
      result = new FE_TimeIteratorAdapter( t_it ) ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == t_it ) ;
   PEL_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter( FE_TimeIterator* t_it )
//-------------------------------------------------------------------------
   : PEL_Object( t_it )
   , IS_PROTO( false )
   , T_IT( t_it )
   , DT( t_it->time_step() )
   , FINISHED( false )
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter" ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == t_it ) ;
   PEL_CHECK_POST( time_iterator() == t_it ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter(
                                       std::string const& a_concrete_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , T_IT( 0 )
   , DT( PEL::bad_double() )
   , FINISHED( false )
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter" ) ;

   plugins_map()->register_item( a_concrete_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter( void )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , T_IT( 0 )
   , DT( PEL::bad_double() )
   , FINISHED( false )
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: FE_TimeIteratorAdapter" ) ;

   static bool first = true ;
   if( !first )
   {
      PEL_Error::object()->raise_internal(
         "Several default FE_TimeIteratorAdapter objects are referenced" ) ;
   }
   first = false ;
   DEFAULT_PROTOTYPE = this ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter:: ~FE_TimeIteratorAdapter( void )
//-------------------------------------------------------------------------
{
   if( this == DEFAULT_PROTOTYPE )
   {
      DEFAULT_PROTOTYPE = 0 ;
   }
}

//-------------------------------------------------------------------------
FE_TimeIterator const*
FE_TimeIteratorAdapter:: time_iterator( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: time_iterator" ) ;
   FE_TimeIterator const* result = T_IT ;
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: initialize_time_step( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: initialize_time_step" ) ;
   PEL_CHECK_PRE( time_iterator()->is_started() ) ;
   PEL_CHECK_PRE( !time_iterator()->is_finished() ) ; 
   DT = T_IT->time_step() ;
   initialize_inner() ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: adapt_time_iterator( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: adapt_time_iterator" ) ;
   PEL_CHECK_PRE( time_iterator()->is_started() ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   FINISHED = com->boolean_and( FINISHED ) ;
   
   if( FINISHED )
   {
      T_IT->finish_iterations() ;
   }
   else
   {
      bool restart = false ;
      double next_dt = DT ;
      define_parameters_for_next_iteration( FINISHED, restart, next_dt ) ;
      if( restart )
      {
         restart_iteration_with_new_time_step( next_dt ) ;
      }
      else if( !T_IT->time_step_is_fixed() )
      {
         DT = next_dt ;
         T_IT->set_next_time_step( DT ) ;
      }
      FINISHED = com->boolean_and( FINISHED ) ;
      if( FINISHED )
      {
         T_IT->finish_iterations() ;
      }
   }
   
   PEL_CHECK_POST(
      EQUIVALENT( iterations_are_finished(),
                  time_iterator()->is_finished() ) ) ;
   PEL_CHECK_POST(
      IMPLIES( !iterations_are_finished() &&
                             !time_iterator()->time_step_is_fixed(),
               time_iterator()->next_time_step() == next_time_step() ) ) ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: propose_next_time_step( double dt )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: propose_next_time_step" ) ;
   PEL_CHECK_PRE( time_iterator()->is_started() ) ;
   PEL_CHECK_PRE( !time_iterator()->is_finished() ) ;
   PEL_CHECK_PRE( dt > 0. ) ;
   DT = next_time_step_from_proposed_one( dt ) ;
}

//-------------------------------------------------------------------------
double
FE_TimeIteratorAdapter:: next_time_step( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: next_time_step" ) ;
   double result = DT ;
   PEL_CHECK_POST( result > 0. ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: propose_to_finish_iterations( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: propose_to_finish_iterations" ) ;
   PEL_CHECK_PRE( time_iterator()->is_started() ) ;
   PEL_CHECK_PRE( !time_iterator()->is_finished() ) ;
   
   FINISHED = true ;
   
   PEL_CHECK_POST( iterations_are_finished() ) ;
}

//-------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: iterations_are_finished( void ) const
//-------------------------------------------------------------------------
{
   return( FINISHED ) ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: set_time_iteration_failed( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: set_time_iteration_failed" ) ;
   PEL_CHECK_PRE( set_time_iteration_failed_PRE() ) ;

   PEL_Error::object()->raise_plain(
      "*** Time step iteration error: iteration failed" ) ;
   
   PEL_CHECK_POST( set_time_iteration_failed_POST() ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   T_IT->save_state( writer ) ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   T_IT->restore_state( reader ) ;
   
   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: initialize_inner( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: initialize_inner" ) ;
   PEL_CHECK_PRE( initialize_inner_PRE() ) ;
   PEL_CHECK_POST( initialize_inner_POST() ) ;
}

//-------------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: define_parameters_for_next_iteration(
                           bool& finished, bool& restart, double& next_dt )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: define_parameters_for_next_iteration" ) ;
   PEL_CHECK_PRE( define_parameters_for_next_iteration_PRE( finished, restart, next_dt ) ) ;
   PEL_CHECK_POST( define_parameters_for_next_iteration_POST( finished, restart, next_dt ) ) ;
}

//----------------------------------------------------------------------
double
FE_TimeIteratorAdapter:: next_time_step_from_proposed_one( double dt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: next_time_step_from_proposed_one" ) ;
   PEL_CHECK_PRE( next_time_step_from_proposed_one_PRE( dt ) ) ;
   double result = dt ;
   PEL_CHECK_POST( next_time_step_from_proposed_one_POST( result, dt ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_TimeIteratorAdapter:: restart_iteration_with_new_time_step( double dt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: restart_iteration_with_new_time_step" ) ;
   PEL_CHECK_PRE( time_iterator()->is_started() ) ;
   PEL_CHECK_PRE( !time_iterator()->is_finished() ) ;
   PEL_CHECK_PRE( !time_iterator()->time_step_is_fixed() ) ;
   PEL_CHECK_PRE( !time_iterator()->just_went_back() ) ;
   PEL_CHECK_PRE( dt > 0. ) ;

   DT = dt ;
   T_IT->go_back() ;
   T_IT->set_next_time_step( DT ) ;

   PEL_CHECK_POST( time_iterator()->just_went_back() ) ;
   PEL_CHECK_POST( time_iterator()->next_time_step() == dt ) ;
   PEL_CHECK_POST( next_time_step() == dt ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: create_replica( FE_TimeIterator* t_it,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( t_it, dom, prms, exp ) ) ;

   PEL_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   
   PEL_CHECK_POST( create_replica_POST( result, t_it, dom, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: create_replica( FE_TimeIterator* t_it,
                                   PDE_SetOfDomains const* sdoms,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( t_it, sdoms, prms, exp ) ) ;

   PEL_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   
   PEL_CHECK_POST( create_replica_POST( result, t_it, sdoms, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: create_replica( FE_TimeIterator* t_it,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( t_it, dom, prms ) ) ;

   PEL_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   
   PEL_CHECK_POST( create_replica_POST( result, t_it, dom, prms ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_TimeIteratorAdapter:: create_replica( FE_TimeIterator* t_it,
                                   PDE_SetOfDomains const* sdoms,
                                   FE_SetOfParameters const* prms ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TimeIteratorAdapter:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( t_it, sdoms, prms ) ) ;

   PEL_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FE_TimeIteratorAdapter* result = 0 ;
   
   PEL_CHECK_POST( create_replica_POST( result, t_it, sdoms, prms ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: is_a_prototype( void ) const
//-------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//-------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: set_time_iteration_failed_PRE( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( time_iterator()->is_started() ) ;
   PEL_ASSERT( !time_iterator()->is_finished() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: set_time_iteration_failed_POST( void ) const
//--------------------------------------------------------------------------
{
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: initialize_inner_PRE( void ) const 
//--------------------------------------------------------------------------
{
   PEL_ASSERT( time_iterator()->is_started() ) ;
   PEL_ASSERT( !time_iterator()->is_finished() ) ;
   PEL_ASSERT( next_time_step() == time_iterator()->time_step() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: initialize_inner_POST( void ) const 
//--------------------------------------------------------------------------
{
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: define_parameters_for_next_iteration_PRE( 
                         bool finished, bool restart, double next_dt ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( time_iterator()->is_started() ) ;
   PEL_ASSERT( !time_iterator()->is_finished() ) ;
   PEL_ASSERT( finished == false ) ;
   PEL_ASSERT( restart == false ) ;
   PEL_ASSERT( next_dt == next_time_step() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: define_parameters_for_next_iteration_POST(
                         bool finished, bool restart, double next_dt ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( finished, !restart ) ) ;
   PEL_ASSERT( IMPLIES( restart, !finished ) ) ;
   PEL_ASSERT( IMPLIES( time_iterator()->time_step_is_fixed(), !restart ) ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: next_time_step_from_proposed_one_PRE(
                                                           double dt ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( time_iterator()->is_started() ) ;
   PEL_ASSERT( !time_iterator()->is_finished() ) ;
   PEL_ASSERT( dt > 0. ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: next_time_step_from_proposed_one_POST(
                                            double result, double dt ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result > 0. ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_PRE(
                                       FE_TimeIterator const* t_it,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( !t_it->is_started() ) ;
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_POST(
                                        FE_TimeIteratorAdapter const* result,
                                        FE_TimeIterator const* t_it,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == t_it ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_PRE(
                                       FE_TimeIterator const* t_it,
                                       PDE_SetOfDomains const* sdoms,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( !t_it->is_started() ) ;
   PEL_ASSERT( sdoms != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_POST(
                                        FE_TimeIteratorAdapter const* result,
                                        FE_TimeIterator const* t_it,
                                        PDE_SetOfDomains const* sdoms,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == t_it ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_PRE(
                                       FE_TimeIterator const* t_it,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( !t_it->is_started() ) ;
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_POST(
                                        FE_TimeIteratorAdapter const* result,
                                        FE_TimeIterator const* t_it,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == t_it ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_PRE(
                                       FE_TimeIterator const* t_it,
                                       PDE_SetOfDomains const* sdoms,
                                       FE_SetOfParameters const* prms ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( !t_it->is_started() ) ;
   PEL_ASSERT( sdoms != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_TimeIteratorAdapter:: create_replica_POST(
                                        FE_TimeIteratorAdapter const* result,
                                        FE_TimeIterator const* t_it,
                                        PDE_SetOfDomains const* sdoms,
                                        FE_SetOfParameters const* prms ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == t_it ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
FE_TimeIteratorAdapter:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "FE_TimeIteratorAdapter descendant" ) ;
   return( result ) ;
}
