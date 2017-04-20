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

#include <PEL_ObjectTest.hh>

#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Timer.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

bool PEL_ObjectTest::FAILURE = false ;

//-------------------------------------------------------------------------
PEL_ObjectTest*
PEL_ObjectTest:: object( std::string const& a_name )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: object" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_ObjectTest* result =
      static_cast<PEL_ObjectTest*>( plugins_map()->item( a_name ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   PEL_CHECK_POST( result->registration_name() == a_name ) ; 
   PEL_CHECK_POST( !result->has_data_deck_explorer() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_ObjectTest:: ~PEL_ObjectTest( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
PEL_ObjectTest:: PEL_ObjectTest( std::string const& tested_class_name,
                                 std::string const& my_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , NAME( my_name )
   , TESTED_CLASS( tested_class_name )
   , EXP( 0 )
   , NB_ELEMENTARY_TESTS( 0 )
   , NB_ELEMENTARY_TESTS_OK( 0 )
{
   PEL_LABEL( "PEL_ObjectTest:: PEL_ObjectTest" ) ;
   
   plugins_map()->register_item( my_name, this ) ;

   PEL_CHECK_POST( owner() == PEL_Root::object() ) ;
   PEL_CHECK_POST( registration_name() == my_name ) ;
   PEL_CHECK_POST( !has_data_deck_explorer() ) ;
}

//-------------------------------------------------------------------------
std::string const&
PEL_ObjectTest:: registration_name( void ) const
//-------------------------------------------------------------------------
{
   return( NAME ) ;
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: set_data_deck_explorer( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: set_data_deck_explorer" ) ;
   PEL_CHECK_PRE( !has_data_deck_explorer() ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( exp->string_data("concrete_name") == registration_name() ) ;

   EXP = exp->create_clone( this ) ;
   
   PEL_CHECK_POST( has_data_deck_explorer() ) ;
}

//-------------------------------------------------------------------------
bool
PEL_ObjectTest:: has_data_deck_explorer( void ) const
//-------------------------------------------------------------------------
{
   return( EXP!=0 ) ;
}

//-------------------------------------------------------------------------
PEL_ModuleExplorer const*
PEL_ObjectTest:: data_deck_explorer( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: data_deck_explorer" ) ;
   PEL_CHECK_PRE( has_data_deck_explorer() ) ;

   PEL_ModuleExplorer const* result = EXP ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: run_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: run_all_tests" ) ;

   PEL_Timer* timer = PEL_Timer::create( 0 ) ;

   out() << endl ;
   out() << "----------------------------------------------------"    << endl ;
   out() << "|  Unit tests performed on class " << TESTED_CLASS << " :"
         << endl ;
   out() << "===================================================="    << endl ;
      
   reset_all_tests() ;
      
   timer->reset() ;
   timer->start() ;

   process_all_tests() ;

   timer->stop() ;
   out() << "===================================================="    << endl ;

   if( NB_ELEMENTARY_TESTS != NB_ELEMENTARY_TESTS_OK )
   {
      out() << "|  Class test FAILED !!!!!!!!!!!!!" << endl ;
   }
   out() << "|  End of " << NB_ELEMENTARY_TESTS << " test" ;
   if( NB_ELEMENTARY_TESTS> 1 )
   {
      out() << "s" ;
   }
   out() << " of class " <<  TESTED_CLASS 
         << " in " << timer->time() << " s " << endl ;
   out() << "----------------------------------------------------"    << endl ;

   timer->destroy() ;
}

//-------------------------------------------------------------------------
bool
PEL_ObjectTest:: tests_of_all_instances_are_successful( void )
//-------------------------------------------------------------------------
{
   return( !FAILURE ) ;
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: reset_all_tests( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: process_all_tests" ) ;

   if( EXP == 0 )
   {
      std::ostringstream m ;
      m << registration_name() << endl << endl ;
      m << "   The default implementation of \"process_all_tests\""  << endl ;
      m << "   provided by \"PEL_ObjectTest\" assumes that there is" << endl ;
      m << "   a specific data deck for \"" << registration_name() << "\"."
	<< endl << endl ;
      m << "   If not, the member function \"process_all_tests\""    << endl ;
      m << "   should be overridden." ;
      PEL_Error::object()->raise_plain( m.str() ) ;
   }
   
   EXP->start_module_iterator() ;
   for( ; EXP->is_valid_module() ; EXP->go_next_module() )
   {
      PEL_ModuleExplorer* ee = EXP->create_subexplorer( 0 ) ;
      process_one_test( ee ) ;
      ee->destroy() ;
   }
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: process_one_test" ) ;

   std::ostringstream m ;
   m << registration_name() << endl << endl ;
   m << "   member function \"process_one_test\" is not implemented"  << endl ;
   PEL_Error::object()->raise_plain( m.str() ) ;
   
}

//-------------------------------------------------------------------------
void
PEL_ObjectTest:: notify_one_test_result( std::string const& displayed_name, 
                                         bool success )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectTest:: notify_one_test_result" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   
   success = PEL_Exec::communicator()->boolean_and(success) ;
   
   out() << "| ... " << std::setw(40) << displayed_name << " :" ;
   if( success )
   {
      out() << "  OK" ;
      NB_ELEMENTARY_TESTS_OK++ ;
   }
   else
   {
      out() << " FAIL" ;
      FAILURE = true ;
   }
   out()  << endl ;
   NB_ELEMENTARY_TESTS++ ;
}

//----------------------------------------------------------------------------
void
PEL_ObjectTest:: print_time_result( std::string const& name, 
                                    double tt ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   PEL::out() << "| ... " << setw( 50 ) << name << " CPU "
              << setprecision( 6 ) << setw( 15 ) << tt << endl ;

   PEL::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------------
void
PEL_ObjectTest:: print_memory_result( std::string const& name, 
                                      size_t memory_size ) const
//----------------------------------------------------------------------------
{
   PEL::out() << "| ... " << setw( 50 ) << name << " MEM "
              << setw( 15 ) << memory_size << endl ;
}

//-------------------------------------------------------------------------
std::ostream &
PEL_ObjectTest:: out( void )
//-------------------------------------------------------------------------
{
   return( PEL::out() ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_ObjectTest:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
                PEL_ObjectRegister::create( PEL_Root::object(),
                                            "PEL_ObjectTest descendant" ) ;
   return( result ) ;
}

