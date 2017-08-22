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

#include <stringVector_TEST.hh>

#include <PEL.hh>
#include <PEL_System.hh>
#include <stringVector.hh>
#include <iostream>
#include <sstream>

stringVector_TEST const*
stringVector_TEST:: REGISTRATOR = new stringVector_TEST() ;

//-------------------------------------------------------------------------
stringVector_TEST:: stringVector_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "stringVector", "stringVector_TEST" )
{
}

//-------------------------------------------------------------------------
stringVector_TEST:: ~stringVector_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   do_tests_0() ;
   do_tests_2() ;
   do_tests_3() ;
   do_tests_4() ;
   do_tests_5() ;
   do_tests_6() ;
   do_tests_7() ;
   do_tests_8() ;
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_0( void )
//-------------------------------------------------------------------------
{
   // Sort:
   {
      stringVector initial_table( 6 ) ;
      initial_table( 0 ) = "b" ;
      initial_table( 1 ) = "aa" ;
      initial_table( 2 ) = "bb" ;
      initial_table( 3 ) = "a" ;
      initial_table( 4 ) = "a" ;

      stringVector sorted = initial_table ;
      sorted.sort() ;

      stringVector result( 6 ) ;
      result( 1 ) = "a" ;
      result( 2 ) = "a" ;
      result( 3 ) = "aa" ;
      result( 4 ) = "b" ;
      result( 5 ) = "bb" ;

      bool ok = ( result == sorted ) ;
      if( !ok ) print_diff( sorted, result ) ;
      notify_one_test_result("sort", ok ) ;
   }

   // String constructor
   {
      stringVector initial_table( 5 ) ;
      initial_table( 0 ) = "b" ;
      initial_table( 1 ) = "aa" ;
      initial_table( 2 ) = "bb" ;
      initial_table( 3 ) = "a" ;
      initial_table( 4 ) = "a" ;

      std::string const s1 = "b,aa,bb,a,a" ;
      stringVector values( s1 ) ;
      bool ok = ( values == initial_table ) ;
      if( !ok ) print_diff( values, initial_table ) ;
      notify_one_test_result("stringVector(string,',')", ok ) ;

      stringVector values1(0) ;
      values1 = s1 ;
      bool ok1 = ( values1 == initial_table ) ;
      if( !ok1 ) print_diff( values1, initial_table ) ;
      notify_one_test_result("operator=(string)", ok1 ) ;

      std::string const s2 = "b aa bb a a" ;
      stringVector values2( s2, ' ' ) ;
      bool ok2 = ( values2 == initial_table ) ;
      if( !ok2 ) print_diff( values2, initial_table ) ;
      notify_one_test_result("stringVector(string,' ')", ok2 ) ;

      std::string const s3 = "b,   aa bb ,bb,a,a" ;
      stringVector values3( s3 ) ;
      initial_table( 1 ) = "   aa bb " ;
      bool ok3 = ( values3 == initial_table ) ;
      if( !ok3 ) print_diff( values3, initial_table ) ;
      notify_one_test_result("stringVector(string_with_blanks)", ok3 ) ;

      stringVector values4( "b,   aa bb ,bb,a,a" ) ;
      initial_table( 1 ) = "   aa bb " ;
      bool ok4 = ( values4 == initial_table ) ;
      if( !ok4 ) print_diff( values4, initial_table ) ;
      notify_one_test_result("stringVector(const char*)", ok4 ) ;
   }
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_2( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_2" ) ;

   bool ok = true ;

   stringVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = sss(i*2) ;
   notify_one_test_result( "size", xxx.size()==10 ) ;
   notify_one_test_result( "capacity", xxx.capacity() == 10 ) ;

   xxx.resize( 0 ) ;
   ok  = ( xxx.size() == 0 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   notify_one_test_result( "(nullify size via) resize", ok ) ;

   ok = true ;
   for( size_t i=0 ; i<10 ; ++i )
   {
      xxx.append( sss(3*i) ) ;
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
      ok = ok && ( xxx( i )=="" ) ;
      xxx( i ) = sss( 5*i ) ;
   }
   notify_one_test_result( "resize", ok ) ;

   xxx.resize( cap-1 ) ;
   ok  = ( xxx.size() == cap-1 ) ;
   ok &= ( xxx.capacity() == cap ) ;

   xxx.append( sss(321984330) ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == sss(3*i) ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == sss(5*i) ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == ""   ) ;
   ok &= ( xxx( cap-1 ) == "321984330" ) ;
   notify_one_test_result( "(unchanged capacity) append", ok ) ;

   xxx.append( sss(789) ) ;
   size_t cap2 = PEL_System::new_block_size( cap, cap+1 ) ;
   ok  = ( cap2 > cap ) ;
   ok &= ( xxx.size() == cap+1 ) ;
   ok &= ( xxx( cap ) == "789" ) ;
   ok &= ( xxx.capacity() == cap2  ) ;
   notify_one_test_result( "append", ok ) ;

   xxx.remove_at( cap-1 ) ;
   ok  = ( xxx.size() == cap ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10      ; ++i )  ok &= ( xxx( i ) == sss(3*i) ) ;
   for( size_t i=10 ; i<15      ; ++i )  ok &= ( xxx( i ) == sss(5*i) ) ;
   for( size_t i=16 ; i<(cap-2) ; ++i )  ok &= ( xxx( i ) == ""   ) ;
   ok &= ( xxx( cap-1 ) == "789" ) ;
   notify_one_test_result( "remove_at", ok ) ;
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_3( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_3" ) ;

   bool ok = true ;

   stringVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = sss(i*2) ;

   xxx.extend( sss(2) ) ;
   ok  = ( xxx.size() == 10 ) ;
   ok &= ( xxx.capacity() == 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == sss(i*2) ) ;
   notify_one_test_result( "(already existing item) extend", ok ) ;

   xxx.extend( sss(3) ) ;
   size_t cap = PEL_System::new_block_size( 10, 11 ) ;
   ok  = ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == sss(i*2) ) ;
   ok &= ( xxx( 10 ) == "3" ) ;
   notify_one_test_result( "(new item, with reallocation) extend", ok ) ;

   xxx.resize( 20 ) ;
   size_t cap2 = cap ;
   if( 20 > cap )
   {
      cap2 = PEL_System::new_block_size( 11, 20 ) ;
   }
   ok  = ( xxx.size() == 20 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0  ; i<10 ; ++i ) ok &= ( xxx( i ) == sss(i*2) ) ;
   ok &= ( xxx( 10 ) == "3" ) ;
   for( size_t i=11 ; i<20 ; ++i ) ok &= ( xxx( i ) == "" ) ;
   xxx.resize( 11 ) ;
   ok &= ( xxx.size() == 11 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   notify_one_test_result( "resize", ok ) ;

   xxx.extend( sss(5) ) ;
   xxx.extend( sss(7) ) ;
   ok  = ( xxx.size() == 13 ) ;
   ok &= ( xxx.capacity() == cap2 ) ;
   for( size_t i=0 ; i<10 ; ++i ) ok &= ( xxx( i ) == sss(i*2) ) ;
   ok &= ( xxx( 10 ) == "3" ) ;
   ok &= ( xxx( 11 ) == "5" ) ;
   ok &= ( xxx( 12 ) == "7" ) ;
   notify_one_test_result( "(new item, without reallocation) extend", ok ) ;
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_4( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_4" ) ;

   bool ok = true ;

   stringVector xxx( 10 ) ;
   for( size_t i=0 ; i<10 ; ++i ) xxx( i ) = sss(i*2) ;
   xxx.resize( 50 ) ;
   xxx.resize( 10 ) ;

   stringVector yyy( 10 ) ;
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
stringVector_TEST:: do_tests_5( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_5" ) ;

   bool ok = true ;
   stringVector xxx( 8 ) ;
   stringVector yyy( 8 ) ;

   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = sss(i*2) ;
   xxx.resize( 28 ) ;
   xxx.re_initialize( 8 ) ;

   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == "" ) ;
   notify_one_test_result( "(with reallocation) re_initialize", ok ) ;

   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = sss(i*2) ;
   xxx.re_initialize( 8 ) ;
   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == "" ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;

   for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = sss(i*2) ;
   xxx.resize( 5 ) ;
   xxx.re_initialize( 8 ) ;
   ok  = ( xxx == yyy ) ;
   ok &= ( xxx.size() == 8 ) ;
   ok &= ( xxx.capacity() == 8 ) ;
   for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == "" ) ;
   notify_one_test_result( "(without reallocation) re_initialize", ok ) ;
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_6( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_6" ) ;

   bool ok ;

   {
      stringVector xxx( 8 ) ;
      for( size_t i=0 ; i<8 ; ++i ) xxx( i ) = sss(i*2) ;
      ok  = ( xxx.size() == 8 ) ;
      ok &= ( xxx.capacity() == 8 ) ;

      xxx.resize( 10 ) ;
      xxx.resize( 8 ) ;
      ok  = ( xxx.size() == 8 ) ;
      ok &= ( xxx.capacity() >= 10 ) ;
      for( size_t i=0 ; i<8 ; ++i ) ok &= ( xxx( i ) == sss(i*2) ) ;

      stringVector yyy( xxx ) ;
      ok &= ( yyy.size() == xxx.size() ) ;
      ok &= ( yyy.capacity() == xxx.capacity() ) ;
      for( size_t i=0 ; i<8 ; ++i ) ok &= ( yyy( i ) == xxx( i ) ) ;

      notify_one_test_result( "stringVector(stringVector const&)", ok ) ;
   }
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_7( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "stringVector_TEST:: do_tests_7" ) ;

   bool ok ;

   {
      stringVector xxx( "a,a,a,a" ) ;

      xxx.resize( 15, "b" ) ;
      ok = ( xxx.size() == 15 ) ;
      for( size_t i=0 ; ok && i<4 ; ++i )
      {
         ok &= ( xxx(i) == "a" ) ;
      }
      for( size_t i=4 ; ok && i<15 ; ++i )
      {
         ok &= ( xxx(i) == "b" ) ;
      }
      notify_one_test_result( "resize( dim>size, val)", ok ) ;

      xxx.resize( 3, "c" ) ;
      ok = ( xxx.size() == 3 ) ;
      for( size_t i=0 ; ok && i<3 ; ++i )
      {
         ok &= ( xxx(i) == "a" ) ;
      }
      notify_one_test_result( "resize( dim<size, val)", ok ) ;

      xxx.resize( 0, "t" ) ;
      ok = ( xxx.size() == 0 ) ;
      notify_one_test_result( "resize( 0, val)", ok ) ;
   }
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: do_tests_8( void )
//-------------------------------------------------------------------------
{
   // Comparisons: operator==, operator!=,
   //              three_way_comparison, operator<, operator>
   {
      std::string const s1 = "b,aa,bb,a,a" ;
      stringVector v1( s1 ) ;
      std::string const s2 = "b,aa,bb,a,a" ;
      stringVector v2( s2 ) ;
      std::string const s3 = "b,aa,bb,a" ;
      stringVector v3( s3 ) ;
      std::string const s4 = "b,aa,b,a,a" ;
      stringVector v4( s4 ) ;

      bool ok = true ;
      ok &=  ( v1 == v2 ) &&  ( v2 == v1 ) ;
      ok &= !( v1 == v3 ) && !( v3 == v1 ) ;
      ok &= !( v1 == v4 ) && !( v4 == v1 ) ;
      ok &= !( v3 == v4 ) && !( v4 == v3 ) ;
      notify_one_test_result( "operator==", ok ) ;

      ok = true ;
      ok &= !( v1 != v2 ) && !( v2 != v1 ) ;
      ok &=  ( v1 != v3 ) &&  ( v3 != v1 ) ;
      ok &=  ( v1 != v4 ) &&  ( v4 != v1 ) ;
      ok &=  ( v3 != v4 ) &&  ( v4 != v3 ) ;
      notify_one_test_result( "operator!=", ok ) ;

      ok = true ;
      ok &= !( v1 > v2 ) && !( v2 > v1 ) ;
      ok &=  ( v1 > v3 ) && !( v3 > v1 ) ;
      ok &=  ( v1 > v4 ) && !( v4 > v1 ) ;
      ok &=  ( v3 > v4 ) && !( v4 > v3 ) ;
      notify_one_test_result( "operator>", ok ) ;

      ok = true ;
      ok &= !( v1 < v2 ) && !( v2 < v1 ) ;
      ok &= !( v1 < v3 ) &&  ( v3 < v1 ) ;
      ok &= !( v1 < v4 ) &&  ( v4 < v1 ) ;
      ok &= !( v3 < v4 ) &&  ( v4 < v3 ) ;
      notify_one_test_result( "operator<", ok ) ;

      ok = true ;
      int c12 = v1.three_way_comparison( v2 ) ;
      int c21 = v2.three_way_comparison( v1 ) ;
      ok &= ( c12 == 0 ) && ( c21 == 0 ) ;

      int c13 = v1.three_way_comparison( v3 ) ;
      int c31 = v3.three_way_comparison( v1 ) ;
      ok &= ( c13 > 0 ) && ( c31 < 0 ) ;

      int c14 = v1.three_way_comparison( v4 ) ;
      int c41 = v4.three_way_comparison( v1 ) ;
      ok &= ( c14 > 0 ) && ( c41 < 0 ) ;

      int c34 = v3.three_way_comparison( v4 ) ;
      int c43 = v4.three_way_comparison( v3 ) ;
      ok &= ( c34 > 0 ) && ( c43 < 0 ) ;

      notify_one_test_result( "three_way_comparison", ok ) ;
   }
}

//-------------------------------------------------------------------------
void
stringVector_TEST:: print_diff( stringVector const& result,
                                stringVector const& expected ) const
//-------------------------------------------------------------------------
{
   out() << "Result:" << std::endl ;
   for( size_t i=0 ; i<result.size() ; ++i )
   {
      out() << "   - \"" << result(i) << "\"" << std::endl ;
   }
   out() << "Expected:" << std::endl ;
   for( size_t i=0 ; i<expected.size() ; ++i )
   {
      out() << "   - \"" << expected(i) << "\"" << std::endl ;
   }
}

//-------------------------------------------------------------------------
std::string
stringVector_TEST:: sss( size_t i )
//-------------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << i ;
   return( mesg.str() ) ;
}

