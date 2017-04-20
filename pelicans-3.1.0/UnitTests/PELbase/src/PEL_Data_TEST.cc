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

#include <PEL_Data_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Data.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>

#include <stringVector.hh>

#include <iostream>

PEL_Data_TEST*
PEL_Data_TEST::REGISTRATOR = new PEL_Data_TEST() ;

//-------------------------------------------------------------------------
PEL_Data_TEST:: PEL_Data_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_Data", "PEL_Data_TEST" )
   , MOD( PEL_Module::create( this, "empty" ) )
{
}

//-------------------------------------------------------------------------
PEL_Data_TEST:: ~PEL_Data_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_Data_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data_TEST:: process_one_test" ) ;
   
   bool request_exception =
      texp->has_entry( "exception" ) && texp->bool_data( "exception" ) ;
   bool catched = false ;
   try
   {
      process_one_comparison( texp ) ;
   }
   catch( PEL_Exceptions::UserError )
   {
      catched = true ;
   }
   std::string nn = texp->name() ;
   if( request_exception || catched )
      notify_one_test_result( nn+" exception", catched==request_exception ) ;
}

//-------------------------------------------------------------------------
void
PEL_Data_TEST:: process_one_comparison( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Data_TEST:: process_one_comparison" ) ;

   std::string const& type = texp->string_data( "type" ) ;
   
   if( type.find( "value_as_string" )==0 )
   {
      std::string nn = texp->name() ;
      std::string const& eval = "exp_to_eval" ;
      PEL_Data const* to_eval = texp->abstract_data( 0, eval ) ;
      std::string as_string = to_eval->value_as_string() ;
      
      PEL_Data const* parsed = MOD->create_evaluation( 0, as_string, 0 ) ;
      
      bool ok = parsed->data_type()==to_eval->data_type() &&
                parsed->value_as_string()==as_string ;
      
      notify_one_test_result( nn+" value_as_string", ok ) ;
      if( !ok ) 
      {
         out() << "to_eval=\""<<as_string<<"\" parsed= \"" << parsed->value_as_string() << std::endl ;
      }

      to_eval->destroy() ; to_eval = 0 ;
      parsed->destroy() ; parsed  = 0 ;
   }
}



