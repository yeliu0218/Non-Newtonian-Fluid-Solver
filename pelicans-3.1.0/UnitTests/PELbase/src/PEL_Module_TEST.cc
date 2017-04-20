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

#include <PEL_Module_TEST.hh>

#include <PEL.hh>
#include <PEL_Double.hh>
#include <PEL_Module.hh>


#include <iostream>

using std::string ;
using std::cerr ;
using std::endl ;

//-------------------------------------------------------------------------
PEL_Module_TEST*
PEL_Module_TEST::unique_instance = new PEL_Module_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
PEL_Module_TEST:: PEL_Module_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_Module", "PEL_Module_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_Module_TEST:: ~PEL_Module_TEST( void )
//-------------------------------------------------------------------------
{
}


//-------------------------------------------------------------------------
void
PEL_Module_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_Module* root = PEL_Module::create( this, "root" ) ;
   
   notify_one_test_result("  name", root->name() == "root" ) ;
   
   PEL_Module* leaf = PEL_Module::create( root, "leaf" ) ;
   root->add_module( leaf ) ;
   
   notify_one_test_result(" absolute_path_name ",
                          leaf->absolute_path_name() == "/root/leaf" ) ;
   notify_one_test_result(" is_empty ",
                          !root->is_empty() && leaf->is_empty() ) ;
   notify_one_test_result(" has_module() ",
                          root->has_module() && !leaf->has_module() ) ;
   
   notify_one_test_result(" has_module(name) ", root->has_module("leaf") ) ;
   notify_one_test_result(" module(name) ", root->module("leaf")==leaf ) ;

   PEL_Double const* d = PEL_Double::create( leaf, 1.0 ) ;
   
   leaf->add_entry( "val", d ) ;
   notify_one_test_result(" has_entry() ",
                          !root->has_entry() && leaf->has_entry() ) ;
   notify_one_test_result(" has_entry(name) ",
                          leaf->has_entry("val") ) ;
   notify_one_test_result(" data_of_entry(name) ",
                          leaf->data_of_entry("val")==d ) ;
   PEL_Double const* r2d2 = PEL_Double::create( leaf, 2.0 ) ;
   root->replace_data_of_entry( "leaf/val", r2d2 ) ;
   notify_one_test_result(" replace_data_of_entry ",
                          leaf->data_of_entry("val")->to_double()==2.0 ) ;
   
   PEL_Module* leaf2 = PEL_Module::create( root, "leaf" ) ;
   PEL_Double const* d2 = PEL_Double::create( leaf2, 1.0 ) ;
   leaf2->add_entry( "val", d2 ) ;
   leaf->merge_module( leaf2 ) ;
   notify_one_test_result(" merge_module ",
                          root->data_of_entry("leaf/val")->to_double()==1.0 ) ;
   
   root->remove_entry( "leaf/val" ) ;
   notify_one_test_result(" remove_entry ",
                          !root->module("leaf")->has_entry() ) ;
   
   root->remove_module( "leaf" ) ;
   notify_one_test_result(" remove_module ",
                          !root->has_module() ) ;
   
   notify_one_test_result(" basename ",
                          PEL_Module::basename( "il/fau/drait/que/je/par/le/en/a/lex/an/drins" ) == "drins" ) ;
   notify_one_test_result(" dirname ",
                          PEL_Module::dirname( "il/fau/drait/que/je/par/le/en/a/lex/an/drins" ) == "il/fau/drait/que/je/par/le/en/a/lex/an" ) ;
   
}



