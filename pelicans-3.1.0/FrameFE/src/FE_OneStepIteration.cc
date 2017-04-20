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

#include <FE_OneStepIteration.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Timer.hh>
#include <PEL_assertions.hh>

#include <GE_Color.hh>

#include <PDE_LocalFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GeometricMultilevel_PC.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDomains.hh>

#include <FE_Parameter.hh>
#include <FE_TimeIteratorAdapter.hh>
#include <FE_TimeIterator.hh>

#include <iomanip>
#include <iostream>
#include <sstream>

using std::ios_base ; using std::setprecision ; using std::setw ;
using std::string ;
using std::endl ;
using std::ostringstream ;
using std::set ;

std::string FE_OneStepIteration::INDENT ;

std::map< std::string, PEL_Timer* > FE_OneStepIteration::ASS_TIMER ;
std::map< std::string, PEL_Timer* > FE_OneStepIteration::SOL_TIMER ;
std::map< std::string, PEL_Timer* > FE_OneStepIteration::TOT_TIMER ;

//----------------------------------------------------------------------
FE_OneStepIteration*
FE_OneStepIteration:: make( PEL_Object* a_owner,
                            PDE_DomainAndFields const* dom,
                            FE_SetOfParameters const* prms,
                            PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: make" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   FE_OneStepIteration const* proto =
      static_cast<FE_OneStepIteration const*>(
                                        plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   FE_OneStepIteration* result = proto->create_replica( a_owner,
                                                        dom, prms, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
//??????????????????? beaucoup d'autres choses
   return( result ) ;
}

//----------------------------------------------------------------------
FE_OneStepIteration*
FE_OneStepIteration:: make( PEL_Object* a_owner,
                            PDE_SetOfDomains const* sdoms,
                            FE_SetOfParameters const* prms,
                            PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: make" ) ;
   PEL_CHECK_PRE( sdoms != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   FE_OneStepIteration const* proto =
      static_cast<FE_OneStepIteration const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   FE_OneStepIteration* result = proto->create_replica( a_owner,
                                                        sdoms, prms, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
//??????????????????? beaucoup d'autres choses
   return( result ) ;
}

//----------------------------------------------------------------------------
FE_OneStepIteration:: FE_OneStepIteration( PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , FAILURE( false )
   , VERBOSE_LEVEL( exp->has_entry( "verbose_level" ) ?
                                    exp->int_data( "verbose_level" ) : 2 )
{
   DOMS.insert( dom ) ;
}

//----------------------------------------------------------------------------
FE_OneStepIteration:: FE_OneStepIteration( PEL_Object* a_owner,
                                           PDE_SetOfDomains const* sdoms,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , FAILURE( false )
   , VERBOSE_LEVEL( exp->has_entry( "verbose_level" ) ?
                                    exp->int_data( "verbose_level" ) : 2 )
{
   for( size_t i=0 ; i<sdoms->nb_domains() ; ++i )
   {
      DOMS.insert( sdoms->domain( i ) ) ;
   }
}

//----------------------------------------------------------------------------
FE_OneStepIteration:: FE_OneStepIteration( std::string const& name )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , FAILURE( false )
   , VERBOSE_LEVEL( 0 )
{
   PEL_LABEL( "FE_OneStepIteration:: FE_OneStepIteration" ) ;

   plugins_map()->register_item( name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------------
FE_OneStepIteration*
FE_OneStepIteration:: create_replica( PEL_Object* a_owner,
                                      PDE_SetOfDomains const* sdoms,
                                      FE_SetOfParameters const* prms,
                                      PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, sdoms, prms, exp ) ) ;

   FE_OneStepIteration* result = 0 ;

   if( exp->has_entry( "domain" ) )
   {
      string const& dom_name = exp->string_data( "domain" ) ;
      PDE_DomainAndFields const* dom = sdoms->domain( dom_name ) ;
      result = create_replica( a_owner, dom, prms, exp ) ;
   }
   else
   {
      ostringstream mesg ;
      mesg << type_name() << " : " << endl ;
      mesg << "   an entry of keyword \"domain\" is required" << endl ;
      mesg << "   if there is more than one domain" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( create_replica_POST( result, a_owner, sdoms, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
FE_OneStepIteration:: ~FE_OneStepIteration( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_before_time_stepping( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_before_inner_iterations_stage(
                                               FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;
   FAILURE = false ;
   PEL_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_inner_iterations_stage(
                                               FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( do_inner_iterations_stage_PRE( t_it ) ) ;
   do_one_inner_iteration( t_it ) ;
}

//----------------------------------------------------------------------------
bool
FE_OneStepIteration:: inner_iterations_stage_failed( void ) const
//----------------------------------------------------------------------------
{
   return( FAILURE ) ;
}

//----------------------------------------------------------------------------
bool
FE_OneStepIteration:: inner_iterations_are_completed(
                                        FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: inner_iterations_are_completed" ) ;
   PEL_CHECK_PRE( inner_iterations_are_completed_PRE( t_it ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_after_inner_iterations_stage(
                                               FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: do_after_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_after_time_adaptation( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: do_after_time_adaptation" ) ;
   PEL_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_after_time_stepping( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: reset_standard_times( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: reset_standard_times" ) ;

   std::map< std::string, PEL_Timer* >::const_iterator it ;
   for( it = ASS_TIMER.begin() ; it != ASS_TIMER.end() ; ++it )
   {
      PEL_Timer* tt = it->second ;
      tt->reset() ;
   }
   for( it = SOL_TIMER.begin() ; it != SOL_TIMER.end() ; ++it )
   {
      PEL_Timer* tt = it->second ;
      tt->reset() ;
   }
   for( it = TOT_TIMER.begin() ; it != TOT_TIMER.end() ; ++it )
   {
      PEL_Timer* tt = it->second ;
      tt->reset() ;
   }
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: print_standard_times( std::ostream& os,
                                            size_t indent_width )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: print_standard_times" ) ;

   std::string space( indent_width, ' ' ) ;

   size_t l_name = 0 ;
   size_t l_ass = 0 ;
   size_t l_sol = 0 ;
   size_t l_tot = 0 ;
   std::map< std::string, PEL_Timer* >::const_iterator it ;
   std::set< std::string > names ;
   for( it = ASS_TIMER.begin() ; it != ASS_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = PEL::max( l_name, (size_t)nn.length() ) ;
      PEL_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_ass = PEL::max( l_ass, (size_t)time.str().length() ) ;
   }
   for( it = SOL_TIMER.begin() ; it != SOL_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = PEL::max( l_name, nn.length() ) ;
      PEL_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_sol = PEL::max( l_sol, (size_t)time.str().length() ) ;
   }
   for( it = TOT_TIMER.begin() ; it != TOT_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = PEL::max( l_name, nn.length() ) ;
      PEL_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_tot = PEL::max( l_tot, time.str().length() ) ;
   }

   l_ass = PEL::max( l_ass, (size_t)10 ) ; // 10 = length of "Assembling"
   l_sol = PEL::max( l_sol, (size_t) 7 ) ; //  7 = length of "Solving"
   l_tot = PEL::max( l_tot, (size_t) 5 ) ; //  5 = length of "Total"
   os << space
      << setw( l_ass+3 ) << "Assembling" << "  "
      << setw( l_sol+2 ) << "Solving" << "  "
      << setw( l_tot+2 ) << "Total" << std::endl ;
   std::set< std::string >::const_iterator itn = names.begin() ;
   for( ; itn != names.end() ; ++itn )
   {
      std::string const& nn = *itn ;

      std::ostringstream time_ass ;
      it = ASS_TIMER.find( nn ) ;
      if( it != ASS_TIMER.end() )
      {
         it->second->print( time_ass, 0 ) ;
      }
      std::ostringstream time_sol ;
      it = SOL_TIMER.find( nn ) ;
      if( it != SOL_TIMER.end() )
      {
         it->second->print( time_sol, 0 ) ;
      }
      std::ostringstream time_tot ;
      it = TOT_TIMER.find( nn ) ;
      if( it != TOT_TIMER.end() )
      {
         it->second->print( time_tot, 0 ) ;
      }

      os << space
         << setw( l_ass+3 ) << time_ass.str() << "  "
         << setw( l_sol+2 ) << time_sol.str() << "  "
         << setw( l_tot+2 ) << time_tot.str() << "  "
         << nn << std::endl ;
   }
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: print_additional_times( std::ostream& os,
                                              size_t indent_width ) const
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: adapt_time_iterator( FE_TimeIteratorAdapter* t_adapter )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: adapt_time_iterator" ) ;
   PEL_CHECK_PRE( adapt_time_iterator_PRE( t_adapter ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: notify_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: notify_inner_iterations_stage_failure" ) ;
   PEL_CHECK_PRE( notify_inner_iterations_stage_failure_PRE() ) ;
   FAILURE = true ;
   PEL_CHECK_POST( notify_inner_iterations_stage_failure_POST() ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: reset_after_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: reset_after_inner_iterations_stage_failure" ) ;
   PEL_CHECK_PRE( reset_after_inner_iterations_stage_failure_PRE() ) ;
   FAILURE = false ;
   PEL_CHECK_POST( reset_after_inner_iterations_stage_failure_POST() ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: do_additional_savings( FE_TimeIterator const* t_it,
                                             PDE_ResultSaver* rs )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: do_additional_savings" ) ;
   PEL_CHECK_PRE( do_additional_savings_PRE( t_it, rs ) ) ;

   set< PDE_DomainAndFields const* >::iterator it =
                                          DOMS.find( rs->attached_domain() ) ;
   if( it != DOMS.end() )
   {
      save_other_than_time_and_fields( t_it, rs ) ;
   }
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: register_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------------
{
   add_storable_objects( list ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: add_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( add_storable_objects_PRE( list ) ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: print" ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << type_name() << std::endl ;
}

//----------------------------------------------------------------------------
std::string const&
FE_OneStepIteration:: indent( void )
//----------------------------------------------------------------------------
{
   return( INDENT ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: increase_indent( void )
//----------------------------------------------------------------------------
{
   INDENT += "   " ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: decrease_indent( void )
//----------------------------------------------------------------------------
{
   if( INDENT.size() < 3 )
   {
      PEL_Error::object()->raise_plain( "attempt to decrease indentation "
                                        "below the zero limit" ) ;
   }
   INDENT.erase( INDENT.length()-3 ) ;
}

//----------------------------------------------------------------------------
bool
FE_OneStepIteration:: is_a_prototype( void ) const
//----------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: configure_multilevel_preconditioner(
                                PEL_ModuleExplorer const* exp,
                                std::string const& keyword,
                                PDE_DomainAndFields const* dom,
                                PDE_SystemNumbering const* nmb ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: configure_multilevel_preconditioner" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( !keyword.empty() ) ;
   PEL_CHECK_PRE( dom != 0 ) ;

   if( exp->has_entry( keyword ) )
   {
      PDE_GeometricMultilevel_PC* prec =
         PDE_GeometricMultilevel_PC::object( exp->string_data( keyword ) ) ;
      prec->set_discretization_scene( dom, nmb ) ;
   }
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: save_other_than_time_and_fields(
                                                FE_TimeIterator const* t_it,
                                                PDE_ResultSaver* rs )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;
}

//----------------------------------------------------------------------------
PEL_Communicator const*
FE_OneStepIteration:: communicator( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: communicator" ) ;

   static PEL_Communicator const* result = PEL_Exec::communicator() ;

   PEL_CHECK_POST( result == PEL_Exec::communicator() ) ;
   return result ;
}

//----------------------------------------------------------------------------
size_t
FE_OneStepIteration:: verbose_level( void ) const
//----------------------------------------------------------------------------
{
   return( VERBOSE_LEVEL ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: start_assembling_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: start_assembling_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   start_internal_timer( ASS_TIMER, "assemble... ", do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: stop_assembling_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: stop_assembling_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   stop_internal_timer( ASS_TIMER, do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: start_solving_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: start_solving_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   start_internal_timer( SOL_TIMER, "solve...", do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: stop_solving_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: stop_solving_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   stop_internal_timer( SOL_TIMER, do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: start_total_timer( std::string const& mesg,
                                         bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: start_total_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 1 ) ;
   start_internal_timer( TOT_TIMER, mesg+"...", do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: stop_total_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: stop_total_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 1 ) ;
   stop_internal_timer( TOT_TIMER, do_display ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: check_param_nb_components( FE_Parameter const* prm,
                                                 std::string const& keyword,
                                                 size_t required_nb_cmps )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: check_param_nb_components" ) ;
   PEL_CHECK_PRE( prm != 0 ) ;

   if( !( prm->nb_components() == required_nb_cmps ) )
   {
      ostringstream mesg ;
      mesg << "The parameter denoted by : \"" << keyword << "\"" << endl ;
      mesg << "should have " << required_nb_cmps << " components" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( prm->nb_components() == required_nb_cmps ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: check_field_nb_components(
                                          PDE_DiscreteField const* ff,
                                          size_t required_nb_cmps )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: check_field_nb_components" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   if( !( required_nb_cmps == ff->nb_components() ) )
   {
      ostringstream mesg ;
      mesg << "The field : \"" << ff->name() << "\"" << endl ;
      mesg << "should have " << required_nb_cmps << " components" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( ff->nb_components() == required_nb_cmps ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: check_field_storage_depth(
                                          PDE_DiscreteField const* ff,
                                          size_t requested_max_level ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: check_field_storage_depth" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   if( !( requested_max_level < ff->storage_depth() ) )
   {
      ostringstream mesg ;
      mesg << "field : \"" << ff->name() << "\"" << endl ;
      mesg << "   the level of storage number "
           << requested_max_level << " will be" << endl ;
      mesg << "   accessed in class " << type_name() << endl ;
      mesg << "   whereas the storage depth is only "
           << ff->storage_depth() ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( ff->storage_depth() > requested_max_level ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIteration:: check_nb_local_nodes( PDE_DiscreteField const* ff,
                                            PDE_LocalFE const* fe,
                                            size_t required_nb_nodes ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: check_nb_local_nodes" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->field_is_handled( ff ) ) ;

   if( fe->nb_local_nodes( ff ) != required_nb_nodes )
   {
      ostringstream mesg ;
      mesg << endl << "*** " << type_name() << endl << endl ;
      mesg << "   the field \"" << ff->name()
           << "\" should have only" << endl ;
      mesg << "   one node located in each cell" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( fe->nb_local_nodes( ff ) == required_nb_nodes ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIteration:: check_boundary_condition(
                                          PDE_SetOfBCs const* bcs,
                                          GE_Color const* color,
                                          PDE_DiscreteField const* field,
                                          size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: check_boundary_condition" ) ;
   PEL_CHECK_PRE( bcs != 0 ) ;
   PEL_CHECK_PRE( color != 0 ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   if( !bcs->has_BC( color, field, ic ) )
   {
      ostringstream mesg ;
      mesg << endl << "*** " << type_name() << endl << endl ;
      mesg << "   the field \"" << field->name()
           << "\" should have a boundary condition" << endl ;
      mesg << "   on all bounds" ;
      if( ic != PDE_SetOfBCs::all_components )
         mesg << "    for its component " << ic ;
      mesg << endl ;
      mesg << "   (such is not the case for at least one bound" << endl ;
      mesg << "    of color \"" << color->name() << "\")" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( bcs->has_BC( color, field, ic ) ) ;
}

//--------------------------------------------------------------------------
void
FE_OneStepIteration:: raise_bad_BC_type( std::string const& type,
                                         std::string const& allowed_types,
                                         PDE_DiscreteField const* field,
                                         size_t ic ) const
//--------------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "*** " << type_name() << ":" << endl << endl ;
   mesg << "   the boundary condition type \"" << type << "\"" << endl ;
   mesg << "   is not allowed for the field \"" << field->name() << "\"" ;
   if( ic != PDE_SetOfBCs::all_components )
      mesg << "(component " << ic << ")" ;
   mesg << endl << endl ;
   mesg << "   allowed types are:" << endl ;
   mesg << allowed_types ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_before_time_stepping_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( !t_it->is_started() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_before_inner_iterations_stage_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_before_inner_iterations_stage_POST(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_inner_iterations_stage_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_one_inner_iteration_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: inner_iterations_are_completed_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_after_inner_iterations_stage_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_after_time_adaptation_PRE(
                                         FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->just_went_back() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: do_additional_savings_PRE(
                                         FE_TimeIterator const* t_it,
                                         PDE_ResultSaver* rs ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( rs != 0 ) ;
   PEL_ASSERT( rs->has_an_opened_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: notify_inner_iterations_stage_failure_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: notify_inner_iterations_stage_failure_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: reset_after_inner_iterations_stage_failure_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: reset_after_inner_iterations_stage_failure_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: adapt_time_iterator_PRE(
                         FE_TimeIteratorAdapter const* t_adapter ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_adapter != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: save_other_than_time_and_fields_PRE(
                                         FE_TimeIterator const* t_it,
                                         PDE_ResultSaver* rs ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( rs != 0 ) ;
   PEL_ASSERT( rs->has_an_opened_cycle() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: add_storable_objects_PRE( PEL_ListIdentity* list ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( list != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: create_replica_PRE(
                                        PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: create_replica_POST(
                                    FE_OneStepIteration const* result,
                                    PEL_Object* a_owner,
                                    PDE_DomainAndFields const* dom,
                                    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: create_replica_PRE(
                                        PEL_Object* a_owner,
                                        PDE_SetOfDomains const* sdoms,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   PEL_ASSERT( sdoms != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_OneStepIteration:: create_replica_POST(
                                    FE_OneStepIteration const* result,
                                    PEL_Object* a_owner,
                                    PDE_SetOfDomains const* sdoms,
                                    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
FE_OneStepIteration:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
        PEL_ObjectRegister::create( PEL_Root::object(),
                                    "FE_OneStepIteration descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: start_internal_timer(
                            std::map< std::string, PEL_Timer* >& timers,
                            std::string const& display,
                            bool do_display ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: start_internal_timers" ) ;

   if( do_display )
   {
      increase_indent() ;
      PEL::out() << indent() << display << endl ;
   }

   PEL_Timer* tt = 0 ;
   std::string const& nn = type_name() ;
   std::map< std::string, PEL_Timer* >::const_iterator it = timers.find( nn ) ;
   if( it == timers.end() )
   {
      tt = PEL_Timer::create( PEL_Root::object() ) ;
      timers[ nn ] = tt ;
   }
   else
   {
      tt = it->second ;
   }
   tt->start() ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIteration:: stop_internal_timer(
                           std::map< std::string, PEL_Timer* > const& timers,
                           bool do_display ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIteration:: stop_internal_timers" ) ;

   if( do_display )
   {
      decrease_indent() ;
   }

   std::map< std::string, PEL_Timer* >::const_iterator it =
                                              timers.find( type_name() ) ;
   PEL_ASSERT( it != timers.end() ) ;

   PEL_Timer* tt = it->second ;
   tt->stop() ;
}

