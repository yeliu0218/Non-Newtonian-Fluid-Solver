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

#include <Code1.hh>

#include <PEL.hh>
#include <PEL_CoupledApplications.hh>
#include <PEL_Double.hh>
#include <PEL_Exec.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_System.hh>
#include <PEL_assertions.hh>

#include <fstream>
#include <iostream>

Code1 const* Code1::PROTOTYPE = new Code1() ;

//----------------------------------------------------------------------
Code1:: Code1( void )
//----------------------------------------------------------------------
   : PEL_Application( "code1" )
{
}

//----------------------------------------------------------------------
Code1* 
Code1:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "Code1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   Code1* result = new Code1( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
Code1:: Code1( PEL_Object* a_owner,
               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , ERROR( PEL::bad_index() )
{
   PEL_LABEL( "Code1:: Code1" ) ;
   if( exp->has_entry( "fail_at_step" ) ) 
   {
      ERROR = exp->int_data( "fail_at_step" ) ;
   }

   if( ! PEL_CoupledApplications::has_coupled_applications() )
   {
      PEL_Error::object()->raise_plain(
         "*** Code1 error:\n"
         "    PEL_CoupledApplication is disabled." ) ;
   }
   if( ! PEL_CoupledApplications::has_application( "code1" ) )
   {
      PEL_Error::object()->raise_plain(
         "*** Code1 error:\n"
         "    \"code1\" is not defined in PEL_CoupledApplication." ) ;
   }
   if( ! PEL_CoupledApplications::has_application( "code2" ) )
   {
      PEL_Error::object()->raise_plain(
         "*** Code1 error:\n"
         "    \"code2\" is not defined in PEL_CoupledApplication." ) ;
   }
}

//----------------------------------------------------------------------
Code1:: ~Code1( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;  
   }
}

//----------------------------------------------------------------------
void
Code1:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "Code1:: run" ) ;

   PEL::out() << "code 1 " << "( pid : " 
              << PEL_System::process_id() << " ) " << std::endl ;
   
   for( size_t i=0 ; i<10 ; i++ ) 
   {
      // Retreive T value
      PEL_ModuleExplorer const* exp =
                  PEL_CoupledApplications::receive( 0, "code1", "code2" ) ;

      double T = exp->double_data( "T" ) ;

      // Compute P
      double P = 1.0/T ;

      PEL::out() << "Step " << i << " T= " << T << " P= " << P << std::endl ;
      
      // Send P value
      PEL_Module* mod = PEL_Module::create( 0, "SentByCode1" ) ;
      mod->add_entry( "P", PEL_Double::create( mod, P ) ) ;
      PEL_CoupledApplications::send(  mod, "code1", "code2" ) ;

      if( ERROR == i ) 
      {
         PEL_Error::object()->raise_internal( "I was said to crash now !" ) ;
      }
      
      mod->destroy() ;
      exp->destroy() ;
   }
}
