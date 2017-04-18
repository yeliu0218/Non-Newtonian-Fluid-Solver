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

#include <FE_MorePostprocessing.hh>

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_ResultSaver.hh>

#include <FE_OneStepIteration.hh>
#include <FE_StepByStepProgression.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ; using std::ostringstream ;

struct FE_MorePostprocessing_ERROR
{
   static void n1( void ) ;
   static void n2( void ) ;
} ;

FE_MorePostprocessing const* 
FE_MorePostprocessing:: PROTOTYPE = new FE_MorePostprocessing() ;

//-------------------------------------------------------------------------
FE_MorePostprocessing:: FE_MorePostprocessing( void )
//-------------------------------------------------------------------------
   : PEL_Application( "FE_MorePostprocessing" )
   , EXP( 0 )
   , READER( 0 )
   , APPLI( 0 )
   , RS( 0 )
   , THE_CYCLE( PEL::bad_index() )
{
}
   
//---------------------------------------------------------------------------
FE_MorePostprocessing*
FE_MorePostprocessing:: create_replica( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MorePostprocessing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   FE_MorePostprocessing* result = new FE_MorePostprocessing( a_owner,
                                                                  exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_MorePostprocessing:: FE_MorePostprocessing( 
                                             PEL_Object* a_owner,
                                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , EXP( exp->create_clone( this ) )
   , READER( 0 )
   , APPLI( 0 )
   , RS( 0 )
   , THE_CYCLE( exp->has_entry( "cycle_number" ) 
                      ? exp->int_data( "cycle_number" ) : PEL::bad_index() )
{
   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "PEL_ObjectReader" ) ;
   READER = PEL_ObjectReader::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   PEL_Module* mod = create_modified_data_deck_module() ;
   ee = PEL_ModuleExplorer::create( 0, mod ) ;
   APPLI = FE_StepByStepProgression::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//---------------------------------------------------------------------------
FE_MorePostprocessing:: ~FE_MorePostprocessing( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_MorePostprocessing:: run( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MorePostprocessing:: run" ) ;

   FE_OneStepIteration::reset_standard_times() ;

   if( THE_CYCLE != PEL::bad_index() )
   {
      handle_one_cycle( THE_CYCLE ) ;
   }
   else
   {
      handle_one_cycle( 0 ) ; 
      for( size_t ic=1 ; ic<=READER->nb_cycles() ; ++ic )
      {
         handle_one_cycle( ic ) ;
      }
   }
}

//---------------------------------------------------------------------------
void
FE_MorePostprocessing:: handle_one_cycle( size_t cycle_number )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MorePostprocessing:: handle_one_cycle" ) ;

   if( cycle_number != 0 )
   {
      READER->seek_cycle( cycle_number ) ;
      APPLI->restore_registered_objects( READER ) ;
      READER->close_cycle() ;
   }

   FE_SetOfParameters const* prms = APPLI->set_of_parameters() ;
   FE_TimeIterator const* time_it = APPLI->time_iterator() ;

   PEL_ModuleExplorer const* ee = 
                          EXP->create_subexplorer( 0, "domain_and_fields" ) ;
   //                     -------------------------------------------------

   size_t verb = ee->int_data( "verbose_level" ) ;
   if( verb != 0 ) PEL::out() << endl ;

   PDE_DomainAndFields* dom = APPLI->domain_and_fields() ;
   if( dom == 0 ) FE_MorePostprocessing_ERROR::n2() ;

   if( ee->has_module( "interior_fields" ) || 
       ee->has_module( "boundary_fields" ) )
   {
      if( THE_CYCLE == PEL::bad_index() ) FE_MorePostprocessing_ERROR::n1() ;
      dom->append_fields( ee ) ;
   }

   PEL_ModuleExplorer* eee = 0 ;
   if( cycle_number == 0 || THE_CYCLE != PEL::bad_index() )
   {
      PEL_ASSERT( RS == 0 ) ;
      eee = ee->create_subexplorer( 0, "PDE_ResultSaver" ) ;
      RS = dom->novel_result_saver( eee ) ;
      eee->destroy() ; eee = 0 ;
   }

   FE_OneStepIteration* oit = 0 ;
   if( EXP->has_module( "FE_OneStepIteration" ) )
   {
      eee = EXP->create_subexplorer( 0, "FE_OneStepIteration" ) ;
      oit = FE_OneStepIteration::make( 0, dom, prms, eee ) ;
      eee->destroy() ; eee = 0 ;
   }

   ee->destroy() ; ee = 0 ;
   //              ------

   RS->start_cycle() ;

   PEL_Exec::out() << endl << "   +++ SAVE FOR POSTPROCESSING" 
        << " *** CYCLE = " << RS->cycle_number()
        << " *** TIME = "  << time_it->time() 
        << " ++++++" << endl << endl ;

   if( oit != 0 ) oit->do_additional_savings( time_it, RS ) ;

   if( cycle_number == 0 || THE_CYCLE != PEL::bad_index() )
   {
      RS->save_grid() ;
   }
   RS->save_fields( 0 ) ;   
   RS->save_variable( time_it->time(), "TIME" ) ;

   RS->terminate_cycle() ;

   if( oit != 0 ) oit->destroy() ;
}

//---------------------------------------------------------------------------
PEL_Module*
FE_MorePostprocessing:: create_modified_data_deck_module( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MorePostprocessing:: create_modified_data_deck_module" ) ;

   PEL_Module* m = 0 ;

   PEL_Module* header = READER->header_module() ;
   if( header->has_module( "PEL_Application" ) )
   {
      m = header->module( "PEL_Application" ) ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "invalid restart file" ) ; 
   }
   PEL_ASSERT( m != 0 ) ;
   PEL_Module* result = m->create_clone( this ) ;

   result->remove_module( "PEL_ObjectWriter" ) ;
   result->remove_module( "PDE_DomainAndFields/PDE_ResultSaver" ) ;

   change_owner( PEL_Root::object(), result ) ;
   
   PEL_CHECK( result != 0 ) ;
   PEL_CHECK( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
FE_MorePostprocessing_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** FE_MorePostprocessing" << endl << endl
        << "   A desired cycle number has to be specified when" << endl
        << "   new discrete fields are defined." << endl
        << "   --> Add an entry of keyword \"cycle_number\" " ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_MorePostprocessing_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** FE_MorePostprocessing" << endl << endl
        << "   multi-domain are not handled"  ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
