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

#include <boolVector_TEST.hh>

#include <PEL.hh>
#include <PEL_System.hh>
#include <boolVector.hh>

#include <iostream>

using std::cout ;
using std::endl ;

boolVector_TEST const* 
boolVector_TEST:: REGISTRATOR = new boolVector_TEST() ;

//-------------------------------------------------------------------------
boolVector_TEST:: boolVector_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "boolVector", "boolVector_TEST" )
{
}

//-------------------------------------------------------------------------
boolVector_TEST:: ~boolVector_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   do_tests_2() ;
   do_tests_3() ;
   do_tests_4() ;
   do_tests_5() ;
   do_tests_6() ;
   do_tests_7() ;
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: do_tests_2( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_2" ) ;

   bool ok = true ;

   boolVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = true ;
   notify_one_test_result( "size", xxx.size()==10 ) ;
   notify_one_test_result( "capacity", xxx.capacity() == 10 ) ;
   
   xxx.resize( 0 ) ;
   ok  = ( xxx.size() == 0 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   notify_one_test_result( "(nullify size via) resize", ok ) ;
   
   ok = true ;
   for( size_t i=0 ; i<10 ; ++i )
   {
      xxx.append( false ) ;
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
      ok = ok && ( xxx( i )==false ) ;
      xxx( i ) = true ;
   }
   notify_one_test_result( "resize", ok ) ;
   
   xxx.resize( cap-1 ) ;
   ok  = ( xxx.size() == cap-1 ) ;
   ok &= ( xxx.capacity() == cap ) ; 
   
   xxx.append( true ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == false ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == true ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == false   ) ;
   ok &= ( xxx( cap-1 ) == true ) ;
   notify_one_test_result( "(unchanged capacity) append", ok ) ;
   
   xxx.append( false ) ;
   size_t cap2 = PEL_System::new_block_size( cap, cap+1 ) ;
   ok  = ( cap2 > cap ) ;
   ok &= ( xxx.size() == cap+1 ) ;
   ok &= ( xxx( cap ) == false ) ;
   ok &= ( xxx.capacity() == cap2  ) ; 
   notify_one_test_result( "append", ok ) ;
   
   xxx.remove_at( cap-1 ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == false ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == true ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == false ) ;
   ok &= ( xxx( cap-1 ) == false ) ;
   notify_one_test_result( "remove_at", ok ) ;
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: do_tests_3( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_3" ) ;

   bool ok = true ;
   
   boolVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = true ;
   
   ok  = ( xxx.size() == 10 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == true ) ;
   notify_one_test_result( "initialization", ok ) ;
   
   xxx.append( false ) ;
   size_t cap = PEL_System::new_block_size( 10, 11 ) ;
   ok  = ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == true ) ;
   ok &= ( xxx( 10 ) == false ) ;
   notify_one_test_result( "(new item, with reallocation) extend", ok ) ;
   
   xxx.resize( 20 ) ;
   size_t cap2 = cap ;
   if( 20 > cap )
   {
      cap2 = PEL_System::new_block_size( 11, 20 ) ;
   }
   ok  = ( xxx.size() == 20 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10 ; ++i ) ok &= ( xxx( i ) == true ) ;
   ok &= ( xxx( 10 ) == false ) ;
   for( size_t i=11 ; i<20 ; ++i ) ok &= ( xxx( i ) == false ) ;
   xxx.resize( 11 ) ;
   ok &= ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   notify_one_test_result( "resize", ok ) ;   
   
   xxx.append( true ) ;
   xxx.append( false ) ;
   ok  = ( xxx.size() == 13 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == true ) ;
   ok &= ( xxx( 10 ) == false ) ;
   ok &= ( xxx( 11 ) == true ) ;
   ok &= ( xxx( 12 ) == false ) ;
   notify_one_test_result( "(new item, without reallocation) extend", ok ) ;
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: do_tests_4( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_4" ) ;
   
   bool ok = true ;
   
   boolVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = i*2 ;
   xxx.resize( 50 ) ;
   xxx.resize( 10 ) ;
   
   boolVector yyy( 10 ) ;
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
boolVector_TEST:: do_tests_5( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_5" ) ;
   
   bool ok = true ;
   boolVector xxx( 8 ) ;
   boolVector yyy( 8 ) ;
   
   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = true ;
   xxx.resize( 28 ) ;
   xxx.re_initialize( 8 ) ;
   
   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == false ) ;
   notify_one_test_result( "(with reallocation) re_initialize", ok ) ;
   
   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = true ;
   xxx.re_initialize( 8 ) ;
   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == false ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;

   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = true ;
   xxx.resize( 5 ) ;   
   xxx.re_initialize( 8 ) ;
   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == false ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: do_tests_6( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_6" ) ;
   
   bool ok ;
   
   {
      boolVector xxx( 8 ) ;
      for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = true ;
      ok  = ( xxx.size() == 8 ) ;
      ok &= ( xxx.capacity() == 8 ) ;

      xxx.resize( 10 ) ;
      xxx.resize( 8 ) ;
      ok  = ( xxx.size() == 8 ) ;
      ok &= ( xxx.capacity() >= 10 ) ;
      for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == true ) ;

      boolVector yyy( xxx ) ;
      ok &= ( yyy.size() == xxx.size() ) ;
      ok &= ( yyy.capacity() == xxx.capacity() ) ;
      for( size_t i=0 ; i<8 ; ++i ) ok &= ( yyy( i ) == xxx( i ) ) ;

      notify_one_test_result( "boolVector(boolVector const&)", ok ) ;
   }
}

//-------------------------------------------------------------------------
void
boolVector_TEST:: do_tests_7( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "boolVector_TEST:: do_tests_7" ) ;
   
   bool ok ;
   
   {
      boolVector xxx( 8, true ) ;
      ok = ( xxx.size() == 8 ) ;
      for( size_t i=0 ; ok && i<8 ; ++i )
      {
         ok &= ( xxx(i) == true ) ;
      }
      notify_one_test_result( "boolVector( dim, val)", ok ) ;

      xxx.resize( 15, false ) ;
      ok = ( xxx.size() == 15 ) ;
      for( size_t i=0 ; ok && i<8 ; ++i )
      {
         ok &= ( xxx(i) == true ) ;
      }
      for( size_t i=8 ; ok && i<15 ; ++i )
      {
         ok &= ( xxx(i) == false ) ;
      }
      notify_one_test_result( "resize( dim>size, val)", ok ) ;

      xxx.resize( 6, false ) ;
      ok = ( xxx.size() == 6 ) ;
      for( size_t i=0 ; ok && i<6 ; ++i )
      {
         ok &= ( xxx(i) == true ) ;
      }
      notify_one_test_result( "resize( dim<size, val)", ok ) ;

      xxx.resize( 0, false ) ;
      ok = ( xxx.size() == 0 ) ;
      notify_one_test_result( "resize( 0, val)", ok ) ;

      xxx.re_initialize( 20, true ) ;
      ok = ( xxx.size() == 20 ) ;
      for( size_t i=0 ; ok && i<20 ; ++i )
      {
         ok &= ( xxx(i) == true ) ;
      }
      notify_one_test_result( "reinitialize( dim, val)", ok ) ;

      xxx.re_initialize( 0, true ) ;
      ok = ( xxx.size() == 0 ) ;
      notify_one_test_result( "reinitialize( 0, val)", ok ) ;

      boolVector yyy( 0, true ) ;
      ok = ( yyy.size() == 0 ) ;
      notify_one_test_result( "boolVector( 0, val)", ok ) ;
   }
}
