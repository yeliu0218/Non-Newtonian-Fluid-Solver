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

#include <PEL_assertions_TEST.hh>

#include <PEL_assertions.hh>
#include <PEL_Exceptions.hh>
#include <PEL.hh>
#include <doubleVector.hh>

#include <iostream>

//-------------------------------------------------------------------------
PEL_assertions_TEST* PEL_assertions_TEST::registered_test = new PEL_assertions_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void
PEL_assertions_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   bool ok = false ;

   PEL_ASSERT( FORALL( ( size_t i=0 ; i<10 ; i++ ),
                       EXISTS( ( size_t j=0 ; j<10 ; j++ ), (j==i) ) ) ) ;
   notify_one_test_result( "PEL_ASSERT( FORALL EXIST)", true ) ;

   try
   {
      PEL_ASSERT( FORALL( ( size_t i=0 ; i<10 ; i++ ),
                          ! EXISTS( ( size_t j=0 ; j<10 ; j++ ), (j==i) ) ) ) ;
   }
   catch( PEL_Exceptions::InternalError )
   {
      std::cout << std::endl << std::endl ;
      ok = true ;
   }
   notify_one_test_result( "PEL_ASSERT( FORALL !EXIST)", ok ) ;

   PEL_ASSERT( ! EXISTS( ( size_t j=0 ; j<10 ; j++ ), (j==10) ) ) ;
   PEL_ASSERT( ! FORALL( ( size_t j=0 ; j<10 ; j++ ), (j==8) ) ) ;
   PEL_ASSERT( FORALL( ( size_t j=0 ; j<10 ; j++ ), true ) ) ;
   PEL_ASSERT( ! FORALL( ( size_t j=0 ; j<10 ; j++ ), false ) ) ;

   PEL_ASSERT(   FORALL( ( size_t j=0 ; j<10 ; j++ ), true ) &&
               ! FORALL( ( size_t j=0 ; j<10 ; j++ ), (j==8) ) ) ;

   try
   {
      PEL_ASSERT( ! FORALL( ( size_t j=0 ; j<10 ; j++ ), true ) ||
                    FORALL( ( size_t j=0 ; j<10 ; j++ ), (j==8) ) ) ;
   }
   catch( PEL_Exceptions::InternalError )
   {
      std::cout << std::endl << std::endl ;
      ok = true ;
   }
   notify_one_test_result( "PEL_ASSERT( !FORALL || FORALL)", ok ) ;

   FORALL( ( size_t j=0 ; j<10 ; j++ ), true ) ;
   notify_one_test_result("  FORALL true", PEL_Assertion::result ) ;

   FORALL( ( size_t j=0 ; j<10 ; j++ ), j<10 ) ;
   notify_one_test_result("  FORALL j<10", PEL_Assertion::result ) ;

   !FORALL( ( size_t j=0 ; j<10 ; j++ ), false ) ;
   notify_one_test_result("  FORALL false", PEL_Assertion::result ) ;

   !FORALL( ( size_t j=0 ; j<10 ; j++ ), j<9 ) ;
   notify_one_test_result("  FORALL j<9", PEL_Assertion::result ) ;

   EXISTS( ( size_t j=0 ; j<10 ; j++ ), true ) ;
   notify_one_test_result("  EXISTS true", PEL_Assertion::result ) ;

   EXISTS( ( size_t j=0 ; j<10 ; j++ ), j>8 ) ;
   notify_one_test_result("  EXISTS j>8", PEL_Assertion::result ) ;

   !EXISTS( ( size_t j=0 ; j<10 ; j++ ), false ) ;
   notify_one_test_result("  EXISTS false", PEL_Assertion::result ) ;

   !EXISTS( ( size_t j=0 ; j<10 ; j++ ), j>9 ) ;
   notify_one_test_result("  EXISTS j>9", PEL_Assertion::result ) ;

   notify_one_test_result("  1=>1", ( IMPLIES( true, true ) ) ) ;
   notify_one_test_result("  0=>1", ( IMPLIES( false, true ) ) ) ;
   notify_one_test_result("  0=>0", ( IMPLIES( false, false ) ) ) ;
   notify_one_test_result("  !1=>0", ( !IMPLIES( true, false ) ) ) ;

   notify_one_test_result("  1<=>1", ( EQUIVALENT( true, true ) ) ) ;
   notify_one_test_result("  !0<=>1", !( EQUIVALENT( false, true ) ) ) ;
   notify_one_test_result("  !1<=>0", !( EQUIVALENT( true, false ) ) ) ;
   notify_one_test_result("  0<=>0", ( EQUIVALENT( false, false ) ) ) ;

   try
   {
      PEL_ASSERT( false ) ;
   }
   catch( PEL_Exceptions::InternalError )
   {
      std::cout << std::endl << std::endl ;
      ok = true ;
   }
   notify_one_test_result( "PEL_ASSERT(false)", ok ) ;

   PEL_ASSERT( true && true ) ;
   PEL_ASSERT(
      FORALL( (size_t i=0 ; i<10 ; ++i ),
              ( IMPLIES( i!=3, true ) ) &&
              ( IMPLIES( i==3, true ) ) ) ) ;

   ok = false ;
   try
   {
      PEL_ASSERT(
                 FORALL( (size_t i=0 ; i<10 ; ++i ),
                         FORALL( ( size_t j=0 ; j<3 ; ++j ),
                                 ( IMPLIES( i!=3, true ) ) &&
                                 ( IMPLIES( i==3, false ) ) ) ) ) ;
   }
   catch( PEL_Exceptions::InternalError )
   {
      std::cout << std::endl << std::endl ;
      ok = true ;
   }
   notify_one_test_result( "PEL_ASSERT FORALL FORALL", ok ) ;


   PEL_ASSERT( FORALL( ( size_t j=0 ; j<10 ; j++ ), true ) ) ;
   PEL_ASSERT( EXISTS( ( size_t j=0 ; j<10 ; j++ ), true ) ) ;
   PEL_ASSERT( !FORALL( ( size_t j=0 ; j<10 ; j++ ), false ) ) ;
   PEL_ASSERT( !EXISTS( ( size_t j=0 ; j<10 ; j++ ), false ) ) ;

   doubleVector xx( 10 ) ;
   PEL_ASSERT( FORALL( ( size_t j=0 ; j<10 ; j++ ),
                       xx(j) <= 1.0 && xx(j) >= -1.0 ) ) ;

   ok = false ;
   try
   {
      PEL_ASSERT( FORALL( ( size_t j=0 ; j<10 ; j++ ),
                          xx(j) <= -1.0 ) ) ;
   }
   catch( PEL_Exceptions::InternalError )
   {
      std::cout << std::endl << std::endl ;
      ok = true ;
   }
   notify_one_test_result( "PEL_ASSERT FORALL", ok ) ;
}



//-------------------------------------------------------------------------
PEL_assertions_TEST:: PEL_assertions_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "PEL_assertions", "PEL_assertions_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_assertions_TEST:: ~PEL_assertions_TEST( void )
//-------------------------------------------------------------------------
{
}



