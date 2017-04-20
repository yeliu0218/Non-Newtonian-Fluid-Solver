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

#include <PEL_TICio_TEST.hh>

#include <PEL.hh>
#include <PEL_Data.hh>
#include <PEL_Module.hh>
#include <PEL_TICio.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

#include <iostream>

//-------------------------------------------------------------------------
PEL_TICio_TEST* PEL_TICio_TEST::registered_test =
new PEL_TICio_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void
PEL_TICio_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio_TEST:: process_all_tests" ) ;
   
   do_test("test.bgene",PEL_TICio::Binary ) ;
   do_test("test.gene",PEL_TICio::Text ) ;
}


//-------------------------------------------------------------------------
void
PEL_TICio_TEST:: do_test( std::string const& file,
                          PEL_TICio::TIC_FORMAT format )
//-------------------------------------------------------------------------
{
   size_t n = 10 ;
   size_t m = 5 ;
   
   PEL_TICio::create_file( file, format ) ;

   PEL_TICio::write_new_cycle( file, format, 2, 1 ) ;

   doubleArray2D da(n,m) ;
   
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
         da(i,j) = (double) i + j ;

   PEL_TICio::write_doubleArray2D_variable( file,
                                            format,
                                            "DA",
                                            da ) ;
   
   intArray2D ia(n,m) ;
   
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
         ia(i,j) = (int) i + j ;

   PEL_TICio::write_intArray2D_variable( file,
                                            format,
                                            "IA",
                                            ia ) ;
   
   PEL_Module* read = PEL_TICio::create_from_gene_file(
                                   this,
                                   "ROOT",
                                   file )   ;
   std::string path = "cycle_1/DA" ;
   
   notify_one_test_result( " read ",
                           read->has_entry( path ) ) ;

   PEL_Data const* data = read->data_of_entry( path ) ;
   notify_one_test_result( " type ",
                           data->data_type()==PEL_Data::DoubleArray2D ) ;
   doubleArray2D const& dar = data->to_double_array2D( ) ;
   bool ok = true ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
         ok = ok && PEL::equal( da(i,j), dar(i,j) ) ;
   
   notify_one_test_result( " equal ",
                           ok ) ;
   
   path = "cycle_1/IA" ;
   
   notify_one_test_result( " read ",
                           read->has_entry( path ) ) ;

   data = read->data_of_entry( path ) ;
   notify_one_test_result( " type ",
                           data->data_type()==PEL_Data::IntArray2D ) ;
   intArray2D const& iar = data->to_int_array2D( ) ;
   ok = true ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
         ok = ok && ia(i,j)== iar(i,j) ;
   
   notify_one_test_result( " equal ",
                           ok ) ;
   
}


//-------------------------------------------------------------------------
PEL_TICio_TEST:: PEL_TICio_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "PEL_TICio", "PEL_TICio_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_TICio_TEST:: ~PEL_TICio_TEST( void )
//-------------------------------------------------------------------------
{
}



