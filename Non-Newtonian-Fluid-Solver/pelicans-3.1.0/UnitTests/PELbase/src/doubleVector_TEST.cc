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

#include <doubleVector_TEST.hh>

#include <PEL.hh>
#include <PEL_System.hh>
#include <doubleVector.hh>

doubleVector_TEST const* 
doubleVector_TEST:: REGISTRATOR = new doubleVector_TEST() ;

//-------------------------------------------------------------------------
doubleVector_TEST:: doubleVector_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "doubleVector", "doubleVector_TEST" )
{
}

//-------------------------------------------------------------------------
doubleVector_TEST:: ~doubleVector_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   do_tests_1() ;
   do_tests_2() ;
   do_tests_3() ;
   do_tests_4() ;
   do_tests_5() ;
   do_tests_6() ;
   do_tests_7() ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_1( void )
//-------------------------------------------------------------------------
{
   doubleVector initial_table( 3 ) ;
   initial_table( 0 ) = 0. ;
   initial_table( 1 ) = 1.0 ;
   initial_table( 2 ) = 2.0 ;
   
   notify_one_test_result( "  index_of ",
        ( initial_table.index_of( 3.0 )==PEL::bad_index() ) &&
        ( initial_table.index_of( 1.0 )==1 ) &&
        ( initial_table.index_of( 2.0, 1.1, 0.0 )==1 )
      ) ;
   notify_one_test_result( "  has ",
        ( initial_table.has( 3.0 )==false ) &&
        ( initial_table.has( 3.0, 1.1 )==true ) &&
        ( initial_table.has( 1.0 )==true )  
      ) ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_2( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_2" ) ;

   bool ok = true ;

   doubleVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   notify_one_test_result( "size", xxx.size()==10 ) ;
   notify_one_test_result( "capacity", xxx.capacity() == 10 ) ;
   
   xxx.resize( 0 ) ;
   ok  = ( xxx.size() == 0 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   notify_one_test_result( "(nullify size via) resize", ok ) ;
   
   ok = true ;
   for( size_t i=0 ; i<10 ; ++i )
   {
      xxx.append( 3.0*((double)i) ) ;
      ok &= ( xxx.size() == i+1 ) ;
      ok &= ( xxx.capacity() == 10 ) ;
   }
   notify_one_test_result( "append", ok ) ;
   
   xxx.resize( 15 ) ;
   size_t cap = PEL_System::new_block_size( 10, 15 ) ;
   ok  = ( xxx.size() == 15 ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=10 ; i<15 ; ++i )
   {
      ok = ok && ( xxx( i )==0.0 ) ;
      xxx( i ) = 5.0*((double)i) ;
   }
   notify_one_test_result( "resize", ok ) ;
   
   xxx.resize( cap-1 ) ;
   ok  = ( xxx.size() == cap-1 ) ;
   ok &= ( xxx.capacity() == cap ) ; 
   
   xxx.append( 321984330.0 ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == 3.0*((double)i) ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == 5.0*((double)i) ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == 0.0   ) ;
   ok &= ( xxx( cap-1 ) == 321984330.0 ) ;
   notify_one_test_result( "(unchanged capacity) append", ok ) ;
   
   xxx.append( 789.0 ) ;
   size_t cap2 = PEL_System::new_block_size( cap, cap+1 ) ;
   ok  = ( cap2 > cap ) ;
   ok &= ( xxx.size() == cap+1 ) ;
   ok &= ( xxx( cap ) == 789.0 ) ;
   ok &= ( xxx.capacity() == cap2  ) ; 
   notify_one_test_result( "append", ok ) ;
   
   xxx.remove_at( cap-1 ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == 3.0*((double)i) ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == 5.0*((double)i) ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == 0.0   ) ;
   ok &= ( xxx( cap-1 ) == 789.0 ) ;
   notify_one_test_result( "remove_at", ok ) ;
}


//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_3( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_3" ) ;

   bool ok = true ;
   
   doubleVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   
   ok  = ( xxx.size() == 10 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == ((double)i)*2.0 ) ;
   notify_one_test_result( "(already existing item) extend", ok ) ;
   
   xxx.append( 3.0 ) ;
   size_t cap = PEL_System::new_block_size( 10, 11 ) ;
   ok  = ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == ((double)i)*2.0 ) ;
   ok &= ( xxx( 10 ) == 3.0 ) ;
   notify_one_test_result( "(new item, with reallocation) extend", ok ) ;
   
   xxx.resize( 20 ) ;
   size_t cap2 = cap ;
   if( 20 > cap )
   {
      cap2 = PEL_System::new_block_size( 11, 20 ) ;
   }
   ok  = ( xxx.size() == 20 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10 ; ++i ) ok &= ( xxx( i ) == ((double)i)*2.0 ) ;
   ok &= ( xxx( 10 ) == 3.0 ) ;
   for( size_t i=11 ; i<20 ; ++i ) ok &= ( xxx( i ) == 0.0 ) ;
   xxx.resize( 11 ) ;
   ok &= ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   notify_one_test_result( "resize", ok ) ;   
   
   xxx.append( 5.0 ) ;
   xxx.append( 7.0 ) ;
   ok  = ( xxx.size() == 13 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == ((double)i)*2.0 ) ;
   ok &= ( xxx( 10 ) == 3.0 ) ;
   ok &= ( xxx( 11 ) == 5.0 ) ;
   ok &= ( xxx( 12 ) == 7.0 ) ;
   notify_one_test_result( "(new item, without reallocation) extend", ok ) ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_4( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_4" ) ;
   
   bool ok = true ;
   
   doubleVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   xxx.resize( 50 ) ;
   xxx.resize( 10 ) ;
   
   doubleVector yyy( 10 ) ;
   yyy = xxx ;
   
   ok  = ( xxx.capacity() >= 50 ) ; 
   ok &= ( yyy.capacity() == 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( yyy( i ) == xxx( i ) ) ;
   notify_one_test_result( "(without reallocation) operator=", ok ) ;
   
   size_t cap = xxx.capacity() ;
   xxx.resize( 12 ) ;
   yyy = xxx ;
   ok  = ( xxx.capacity() == cap ) ; 
   ok &= ( yyy.capacity() == 12 ) ;
   ok &= ( yyy.size() == 12 ) ;
   for( size_t i=0 ; i<12 ; ++i ) ok &= ( yyy( i ) == xxx( i ) ) ;
   notify_one_test_result( "(with reallocation) operator=", ok ) ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_5( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_5" ) ;
   
   bool ok = true ;
   doubleVector xxx( 8 ) ;
   doubleVector yyy( 8 ) ;
   
   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   xxx.resize( 28 ) ;
   xxx.re_initialize( 8 ) ;
   
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == 0 ) ;
   notify_one_test_result( "(with reallocation) re_initialize", ok ) ;
   
   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   xxx.re_initialize( 8 ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == 0.0 ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;

   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   xxx.resize( 5 ) ;   
   xxx.re_initialize( 8 ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == 0.0 ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_6( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_6" ) ;
   
   bool ok ;
   
   doubleVector xxx( 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = ((double)i)*2.0 ;
   ok  = ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   
   xxx.resize( 10 ) ;
   xxx.resize( 8 ) ;
   ok  = ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() >= 10 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == ((double)i)*2.0 ) ;
   
   doubleVector yyy( xxx ) ;
   ok &= ( yyy.size() == xxx.size() ) ;
   ok &= ( yyy.capacity() == xxx.capacity() ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( yyy( i ) == xxx( i ) ) ;
   
   notify_one_test_result( "copy constructor", ok ) ;
}

//-------------------------------------------------------------------------
void
doubleVector_TEST:: do_tests_7( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector_TEST:: do_tests_7" ) ;
   
   bool ok ;
   
   {
      doubleVector xxx( 8, 10 ) ;
      ok = ( xxx.size() == 8 ) ;
      for( size_t i=0 ; ok && i<8 ; ++i )
      {
         ok &= ( xxx(i) == 10 ) ;
      }
      notify_one_test_result( "doubleVector( dim, val)", ok ) ;

      xxx.resize( 15, 2 ) ;
      ok = ( xxx.size() == 15 ) ;
      for( size_t i=0 ; ok && i<8 ; ++i )
      {
         ok &= ( xxx(i) == 10 ) ;
      }
      for( size_t i=8 ; ok && i<15 ; ++i )
      {
         ok &= ( xxx(i) == 2 ) ;
      }
      notify_one_test_result( "resize( dim>size, val)", ok ) ;

      xxx.resize( 6, 7 ) ;
      ok = ( xxx.size() == 6 ) ;
      for( size_t i=0 ; ok && i<6 ; ++i )
      {
         ok &= ( xxx(i) == 10 ) ;
      }
      notify_one_test_result( "resize( dim<size, val)", ok ) ;
      
      xxx.resize( 0, 7 ) ;
      ok = ( xxx.size() == 0 ) ;
      notify_one_test_result( "resize( 0, val)", ok ) ;

      xxx.re_initialize( 20, 1 ) ;
      ok = ( xxx.size() == 20 ) ;
      for( size_t i=0 ; ok && i<20 ; ++i )
      {
         ok &= ( xxx(i) == 1 ) ;
      }
      notify_one_test_result( "reinitialize( dim, val)", ok ) ;
      
      xxx.re_initialize( 0, 10 ) ;
      ok = ( xxx.size() == 0 ) ;
      notify_one_test_result( "reinitialize( 0, val)", ok ) ;

      doubleVector yyy( 0, 10 ) ;
      ok = ( yyy.size() == 0 ) ;
      notify_one_test_result( "doubleVector( 0, val)", ok ) ;
   }
}
