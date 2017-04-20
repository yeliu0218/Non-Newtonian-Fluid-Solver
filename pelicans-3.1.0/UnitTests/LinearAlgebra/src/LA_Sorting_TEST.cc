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

#include <LA_Sorting_TEST.hh>

#include <LA_SeqVector.hh>
#include <LA_Sorting.hh>
#include <PEL_Root.hh>
#include <PEL_Timer.hh>

#include <iostream>

//-------------------------------------------------------------------------
LA_Sorting_TEST*
LA_Sorting_TEST:: registered = new LA_Sorting_TEST() ;
//-------------------------------------------------------------------------




//-------------------------------------------------------------------------
LA_Sorting_TEST:: LA_Sorting_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "LA_Sorting", "LA_Sorting_TEST" )
{
}



//-------------------------------------------------------------------------
LA_Sorting_TEST:: ~LA_Sorting_TEST( void )
//-------------------------------------------------------------------------
{
}



//-------------------------------------------------------------------------
void
LA_Sorting_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Sorting_TEST:: process_all_tests" ) ;
   
   const int n = 2000 ;
   
   PEL_Timer* tim = PEL_Timer::create( this ) ;
   
   LA_SeqVector * r_rand = LA_SeqVector::create( this, n ) ;
   for( int i=0 ; i<n ; i++ )
   {
      r_rand->set_item( i, 1.0*PEL::rand() ) ;
   }
   LA_SeqVector * r_sorted = LA_SeqVector::create( this, n ) ;
   size_t_vector r_rank(  n ) ;
   tim->start() ;
   r_rand->synchronize() ;
   
   LA_Sorting::sort( r_rand, r_sorted, r_rank ) ;
   tim->stop() ;
   double insert_t = tim->time() ;
   
   bool ok = true ;
   for( int i=1 ; i<n ; i++ )
   {
      bool okt1 = r_sorted->item( i ) >= r_sorted->item( i-1 ) ;
      if( !okt1 )
      {
         out() << "Bad order : " << r_sorted->item( i )
               << " <= " <<r_sorted->item( i-1 ) << std::endl ;
      }
      bool okt2 = r_rand->item( i ) == r_sorted->item( r_rank( i ) ) ;
      if( !okt2 )
      {
         out() << "Rank def : " << r_rand->item( i )
               << " = " << r_sorted->item( r_rank( i ) ) << std::endl ;
      }
      ok = ok && okt1 && okt2 ;
   }
   notify_one_test_result( "insert sorting", ok ) ;
   
   r_sorted->set( r_rand ) ;
   tim->reset() ;
   tim->start() ;
   LA_Sorting::quick_sort( r_sorted, 0, n-1 ) ;
   tim->stop() ;
   double quick_t = tim->time() ;
   ok = true ;
   for( int i=1 ; i<n ; i++ )
   {
      bool okt1 = r_sorted->item( i ) >= r_sorted->item( i-1 ) ;
      if( !okt1 )
      {
         out() << "Bad order : " << r_sorted->item( i )
               << " <= " <<r_sorted->item( i-1 ) << std::endl ;
      }
      ok = okt1 ;
   }
   notify_one_test_result( "quick sorting", ok ) ;
   out() << "insert : " << insert_t << "s quick : "<< quick_t<<"s"<<std::endl ;
}   


