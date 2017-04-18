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

#include <PEL_IndexSet_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_IndexSet.hh>

#include <iostream>

PEL_IndexSet_TEST const*
PEL_IndexSet_TEST::PROTOTYPE = new PEL_IndexSet_TEST() ;

//-------------------------------------------------------------------------
PEL_IndexSet_TEST:: PEL_IndexSet_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_IndexSet", "PEL_IndexSet_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_IndexSet_TEST:: ~PEL_IndexSet_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_IndexSet_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet_TEST:: process_all_tests" ) ;

   bool ok = true ;
   size_t_vector vec(4) ;
   vec(0)=3 ;
   vec(1)=1 ;
   vec(2)=2 ;
   vec(3)=4 ;
   
   size_t_vector sorted(vec) ;
   sorted.sort_increasingly() ;
   
   PEL_IndexSet* A = PEL_IndexSet::create( 0, vec, 10 ) ;
   
   ok = ( A->elements() == sorted ) ;
   if( !ok )
   {
      A->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "elements", ok ) ;
   
   ok = ( A->id() == 10 ) ;
   if( !ok )
   {
      A->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "id", ok ) ;
   
   // Test8 :
   vec(3) = 5 ;
   PEL_IndexSet* B = PEL_IndexSet::create( A, vec, 12 ) ;
   ok = ( !A->is_equal( B ) && !B->is_equal( A )
              && A->three_way_comparison( B )==-1
              && B->three_way_comparison( A )==1 ) ;
   if( !ok )
   {
      A->print( out(), 3 ) ;
      out() << std::endl ;
      B->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "compare#0", ok ) ;
   

   // Test9 :
   vec(3) = 3 ;
   PEL_IndexSet* C = PEL_IndexSet::create( A, vec, 12 ) ;
   ok = ( !C->is_equal( B ) && !B->is_equal( C )
              && C->three_way_comparison( B )==-1
              && B->three_way_comparison( C )==1 ) ;
   if( !ok )
   {
      B->print( out(), 3 ) ;
      out() << std::endl ;
      C->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "compare#1", ok ) ; 
   
   vec(3) = 4 ;
   B = PEL_IndexSet::create( A, vec, 15 ) ;
   ok = ( B->is_equal(A) && A->is_equal( B )
              && A->three_way_comparison( B )==0
              && B->three_way_comparison( A )==0 ) ;
   if( !ok )
   {
      B->print( out(), 3 ) ;
      out() << std::endl ;
      C->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "compare#2", ok ) ;
   
   ok = ( B->is_equal(B) && B->three_way_comparison( B )==0 ) ;
   if( !ok )
   {
      B->print( out(), 3 ) ;
      out() << std::endl ;
   }
   notify_one_test_result( "compare#3", ok ) ;

   A->destroy() ; A = 0 ;
}
