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

#include <PEL_TestSuite.hh>

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Iterator.hh>
#include <PEL_List.hh>
#include <PEL_ObjectTest.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <iostream>
#include <sstream>

using std::endl ;

struct PEL_TestSuite_ERROR
{
   static void n0( std::string const& t_name ) ;
} ;

PEL_TestSuite const* PEL_TestSuite::PROTOTYPE = new PEL_TestSuite() ;

//-------------------------------------------------------------------------
PEL_TestSuite:: PEL_TestSuite( void )
//-------------------------------------------------------------------------
   : PEL_Application( "PEL_TestSuite" )
   , TESTS( 0 )
{
}

//-------------------------------------------------------------------------
PEL_TestSuite*
PEL_TestSuite:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TestSuite:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_TestSuite* result = new PEL_TestSuite( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_TestSuite:: PEL_TestSuite( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , TESTS( 0 )
{
   TESTS = PEL_List::create( this ) ;

   stringVector names( 0 ) ;

   if( exp->has_entry( "without_data_deck" ) )
   {
      names = exp->stringVector_data( "without_data_deck" ) ;
      for( size_t i=0 ; i<names.size() ; ++i )
      {
         TESTS->append( PEL_ObjectTest::object( names( i ) ) ) ;
      }
   }
   if( exp->has_module( "with_data_deck" ) )
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "with_data_deck" ) ;
   
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer* sexp = ee->create_subexplorer( 0 ) ;

         std::string const& t_name = sexp->string_data( "concrete_name" ) ;
         if( names.size()!=0 && names.has( t_name ) )
            PEL_TestSuite_ERROR::n0( t_name ) ;

         PEL_ObjectTest* ot = PEL_ObjectTest::object( t_name ) ;
         ot->set_data_deck_explorer( sexp ) ;
         TESTS->append( ot ) ;
         sexp->destroy() ;
      }
      ee->destroy() ; ee = 0 ;
   }
}

//-------------------------------------------------------------------------
PEL_TestSuite:: ~PEL_TestSuite( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_TestSuite:: run( void )
//-------------------------------------------------------------------------
{   
   PEL_Iterator* it = TESTS->create_iterator( this ) ;
   
   for( it->start(); it->is_valid() ; it->go_next() )
   {
      PEL_ObjectTest* ut = static_cast<PEL_ObjectTest*>( it->item() ) ;
      PEL_CHECK( dynamic_cast<PEL_ObjectTest*>( it->item() ) != 0 ) ;

      ut->run_all_tests() ;
   }

   int exit_code = 0 ;
   if( !PEL_ObjectTest::tests_of_all_instances_are_successful() )
   {
      PEL_Exec::out() << endl << "!!!! Failure of some unit tests !!!!" << endl << endl ;
      exit_code = 5 ;
   }
   else
   {
      PEL_Exec::out() << endl << "---- Success of all unit tests ---- " << endl << endl ;
   }
   PEL_Exec::set_exit_code( exit_code ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TestSuite_ERROR:: n0( std::string const& t_name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << endl << "*** PEL_TestSuite error:" << endl << endl ;
   mesg << "    \"" << t_name << "\" cannot appear both:" << endl ;
   mesg << "       - in a submodule of MODULE with_data_deck" << endl ;
   mesg << "       - in the entry of keyword without_data_deck";
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
