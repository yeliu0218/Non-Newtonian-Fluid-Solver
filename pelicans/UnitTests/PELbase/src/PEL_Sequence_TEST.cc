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

#include <PEL_Sequence_TEST.hh>

#include <iostream>
#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Iterator.hh>
#include <PEL_List.hh>

#include <PEL_Double.hh>


//-------------------------------------------------------------------------
void
PEL_Sequence_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   // Base ancestor tests
    PEL_Collection_TEST::process_all_tests() ;

   // Test of PEL_Sequence only fonctionalities
   PEL_Sequence * list = emptyList( this ) ;
   for( size_t i = 0 ; i<10 ; i++ )
   {
      list->append( PEL_Double::create( this, 2.0 ) ) ;
   }
   notify_one_test_result( "append", list->count()==10 ) ;
   list->clear() ;
   
   const size_t n = 10 ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      list->prepend( PEL_Double::create( this, (double) i ) ) ;
   }

   bool ok1 = true ;
   bool ok2 = list->count()==n ;
   for( int i = n-1 ; i>=0 ; i-- )
   {
      PEL_Double * d = dynamic_cast<PEL_Double*>( list->at(i) ) ;
      if( d==0 )
      {
         ok1 = false ;
      }
      else
      {
         ok2 = ok2 && d->to_double() == (double)(n-1-i) ;
      }
   }
   notify_one_test_result( "at", ok1 ) ;
   notify_one_test_result( "prepend", ok2 ) ;
   
   PEL_Double * dLast = dynamic_cast<PEL_Double*>( list->at(n-1) ) ;
   notify_one_test_result( "last", ( dLast!=0 && dLast->to_double()==0.0 ) ) ;
  
   PEL_Double * dFirst = dynamic_cast<PEL_Double*>( list->at(0) ) ;
   notify_one_test_result( "first", ( dFirst!=0 && dFirst->to_double()==(double)(n-1) ) ) ;

   dLast = dynamic_cast<PEL_Double*>( list->at(n-1) ) ;
   list->remove_at( n-1 ) ;
   notify_one_test_result( "removeLast",
               ( dLast!=0 && dLast->to_double()==0.0 && list->count()==(n-1) ) ) ;
   
   dFirst = dynamic_cast<PEL_Double*>( list->at(0) ) ;
   list->remove_at( 0 ) ;
   notify_one_test_result( "removeFirst",
               ( dFirst!=0 && dFirst->to_double()==(double)(n-1) && list->count()==n-2 ) ) ;
  
   PEL_Double * dMiddle = dynamic_cast<PEL_Double*>( list->at(0) ) ;
   list->remove_at(0) ;
   bool ok = ( dMiddle !=0 && dMiddle->to_double()==(double)(n-2) && list->count()==n-3 ) ;
   dMiddle = dynamic_cast<PEL_Double*>( list->at(n-4) ) ;
   list->remove_at(n-4) ;
   ok = ok && ( dMiddle !=0 && dMiddle->to_double()==1.0 && list->count()==n-4 ) ;
   
   size_t m = 2 ;
   PEL_ASSERT( n>7 ) ;
   dMiddle = dynamic_cast<PEL_Double*>( list->at(m) ) ;
   list->remove_at(m) ;
   ok = ok && ( dMiddle !=0 && dMiddle->to_double()==(double)(n-5) && list->count()==n-5 ) ;

   notify_one_test_result( "remove_at", ok ) ;
   
   list->insert_at(0, PEL_Double::create( this,-1) ) ;
   PEL_Double * dIns = dynamic_cast<PEL_Double*>( list->at(0) ) ;
   ok = dIns->to_double() == -1.0 && list->count()==n-4 ;
   
   list->insert_at(n-5, PEL_Double::create( this,-2) ) ;
   dIns = dynamic_cast<PEL_Double*>( list->at(list->count()-2) ) ;
   ok = ok && dIns->to_double() == -2.0 ;
   notify_one_test_result( "insert_at", ok ) ;

   list->clear() ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      list->extend( PEL_Double::create( this, (double) i ) ) ;
   }
   PEL_Sequence * list2 = list->create_clone( this ) ;
   PEL_Iterator * it = list2->create_iterator( this ) ;
   size_t cpt = 0 ;
   ok = true ;
   bool oki = true ;
   
   for( ; it->is_valid() ; it->go_next() )
   {
      PEL_Double * d = dynamic_cast<PEL_Double*>( it->item() ) ;
      ok = ok && ( d->to_double() == (double)cpt ) ;
      PEL_Double * cmp = PEL_Double::create( this, cpt) ;      
      oki = oki && list->index_of( cmp ) == cpt ;
      cpt++ ;
   }
   notify_one_test_result( "create_iterator", cpt==n ) ;
   notify_one_test_result( "clone", ok ) ;
   notify_one_test_result( "index", oki ) ;
   
   PEL_Sequence * vector = emptyList( this ) ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      vector->append( PEL_Double::create( this, 0 ) ) ;
   }
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      vector->set_at( i, PEL_Double::create( this, i ) ) ;
   }
   ok = true ;
   for( int i = n-1 ; i>=0 ; i-- )
   {
      PEL_Object* o = vector->at( i ) ;
      PEL_Double * d = dynamic_cast<PEL_Double * >(o) ;
      ok = ok && d->to_double()==(double)i ;
   }
   
   notify_one_test_result( "at, set_at", ok ) ;
   
   list->clear() ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      PEL_Object * a_owner = this ;
      if( i==4 || i==5 )
      {
         a_owner = 0 ;
      }
      list->extend( PEL_Double::create( a_owner, (double) i ) ) ;
   }
   list->remove_section( 2, 2 ) ;
   PEL_Double * d = dynamic_cast<PEL_Double * >( list->at( 2 ) ) ;
   ok = list->count()==n-2 && d->to_double()==(double)4 ;
   notify_one_test_result( "remove_section", ok ) ;
   
   list->destroy_items_and_remove_section( 2, 2 ) ;
   d = dynamic_cast<PEL_Double * >( list->at( 2 ) ) ;
   ok = list->count()==n-4 && d->to_double()==(double)6 ;
   notify_one_test_result( "destroy_items_and_remove_section", ok ) ;
   
   list->remove_section( 1, list->count()-1 ) ;
   ok = list->count()==1 ;
   notify_one_test_result( "remove_section to end", ok ) ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      vector->set_at( i, PEL_Double::create( this, n-i-1 ) ) ;
   }
   vector->sort() ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      PEL_Double* dbl = static_cast<PEL_Double*>( vector->at(i) ) ;
      ok=ok&&PEL::equal(dbl->to_double(),(double)i) ;
      if(!ok)
      {
         out()<<i<<" "<<dbl->to_double()<<" <> "<<dbl->to_double()-i<<std::endl ;
      }
   }
   notify_one_test_result( "sort", ok ) ;
   
}



//-------------------------------------------------------------------------
PEL_Sequence_TEST:: PEL_Sequence_TEST( 
                      std::string const& tested_class,
                      std::string const& registration ) 
//-------------------------------------------------------------------------
   :   PEL_Collection_TEST( tested_class, registration )
{
}



//-------------------------------------------------------------------------
PEL_Sequence_TEST:: ~PEL_Sequence_TEST( void )
//-------------------------------------------------------------------------
{
}



