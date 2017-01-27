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

#include <PEL_ModulePattern_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <stringVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>


#include <string>

//-------------------------------------------------------------------------
PEL_ModulePattern_TEST*
PEL_ModulePattern_TEST::unique_instance = new PEL_ModulePattern_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
PEL_ModulePattern_TEST:: PEL_ModulePattern_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_ModulePattern", "PEL_ModulePattern_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_ModulePattern_TEST:: ~PEL_ModulePattern_TEST( void )
//-------------------------------------------------------------------------
{
}


//-------------------------------------------------------------------------
void
PEL_ModulePattern_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_ModuleExplorer * mand = texp->create_subexplorer( this, "Mandatory" ) ;

   double val = mand->double_data( "double_val" ) ;
   mand->test_data( "double_val", "double_val>0.0" ) ;
   mand->set_default( "double_val", "1.0" ) ;
   mand->set_help( "double_val", "A positive double value" ) ;
   notify_one_test_result(" val>0.0 ", val>0.0 ) ;
   
   mand->int_data( "int_val" ) ;
   mand->test_data( "int_val", "int_val>0" ) ;

   if( texp->has_module( "Optional" ) )
   {
      PEL_ModuleExplorer const* opt = texp->create_subexplorer( this, "Optional" ) ;
      if(opt->has_entry("double_val") )
         val = opt->double_data( "double_val" ) ;
   }

   std::string const& name = texp->string_data( "concrete_name" ) ;
   texp->set_default( "concrete_name", "first_instance" ) ;
   
   std::string some ;
   if( name=="first_instance" )
   {
      some = texp->string_data( "first_string" ) ;
      texp->test_data_in( "first_string", "choice1,choice2" ) ;
   }
   else
   {
      some = texp->string_data( "other_string" ) ;
   }

   
   std::string const& type = mand->string_data( "type" ) ;
   mand->set_default( "type", "first_type" ) ;
   
   std::string any ;
   if( type=="first_type" )
   {
     any = mand->string_data( "first_type" ) ;
   }
   else
   {
     any = mand->string_data( "other_type" ) ;
   }
   mand->doubleVector_data( "dv" ) ;
   mand->intVector_data( "iv" ) ;
   mand->test_data( "iv", "size(iv)=size(dv)" ) ;
   notify_one_test_result(" building pattern ", true ) ;

   if( mand->has_module( "LIST" ) )
   {
      PEL_ModuleExplorer * mlist = mand->create_subexplorer( this, "LIST" ) ;
      
      for( mlist->start_module_iterator() ;
           mlist->is_valid_module() ;
           mlist->go_next_module() ) 
      {
         PEL_ModuleExplorer * sexp = mlist->create_subexplorer(mlist) ;
         sexp->string_data( "name" ) ;
         sexp->test_data_as( "name", "/with_data_deck/*/*/Mandatory/ASSOCIATED/*/associated_name" ) ;
      }
      
   }
   if( mand->has_module( "ASSOCIATED" ) )
   {
      PEL_ModuleExplorer * massoc = mand->create_subexplorer( this, "ASSOCIATED" ) ;
      for( massoc->start_module_iterator() ;
           massoc->is_valid_module() ;
           massoc->go_next_module() ) 
      {
         PEL_ModuleExplorer * sexp = massoc->create_subexplorer(massoc) ;
         sexp->string_data( "associated_name" ) ;
         sexp->test_data_as( "associated_name", "../../LIST/*/name" ) ;
      }
      
   }
   if( texp->has_module( "TEST_FILE" ) )
   {
      PEL_ModuleExplorer * sexp = texp->create_subexplorer( this, "TEST_FILE" ) ;
      sexp->string_data( "existing_filename" ) ;
      sexp->test_file( "existing_filename", "read" ) ;
      
      sexp->string_data( "file_to_create" ) ;
      sexp->test_file( "file_to_create", "write" ) ;
   }
   
   
}



