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

#include <PEL_Collection_TEST.hh>

#include <iostream>

#include <PEL_Double.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>

#include <iostream>

using std::endl ;

//-------------------------------------------------------------------------
PEL_Collection_TEST:: PEL_Collection_TEST( 
                      std::string const& tested_class,
                      std::string const& registration )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( tested_class, registration )
{
}



//-------------------------------------------------------------------------
PEL_Collection_TEST:: ~PEL_Collection_TEST( void )
//-------------------------------------------------------------------------
{
}



//-------------------------------------------------------------------------
void
PEL_Collection_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_Collection * list = emptyList( this ) ;
   PEL_Iterator * it = list->create_iterator( this ) ;
   notify_one_test_result( "Empty list", list->count()==0 ) ;

   list->clear() ;
   notify_one_test_result( "Clear empty list", list->count()==0 ) ;
   
   list->extend( PEL_Double::create( this, 0.0) ) ;
   notify_one_test_result( "extend", list->count()==1 ) ;
   list->extend( PEL_Double::create( this, 1.0) ) ;
   PEL_Double * d2 = PEL_Double::create( this, 2.0) ;
   list->extend( d2 ) ;
   notify_one_test_result( "Size", list->count()==3 ) ;
   
   PEL_Double * doubled2 = PEL_Double::create( this, 2.0) ;
   PEL_Double * d2Find = dynamic_cast<PEL_Double * >( list->item( doubled2 ) ) ;
   
   notify_one_test_result( "find", d2Find==d2 ) ;
   
   notify_one_test_result( "add", list->count()==3 ) ;
   size_t nbElem = 0 ;
   
   bool ok = true ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {      
      PEL_Double* d = dynamic_cast<PEL_Double*>( it->item() ) ;
      bool okt = ( d != 0 ) ;
      // PEL_Double* dd = PEL_Double::create( this,  *d ) ;
      PEL_Double* dd = d->create_clone( this ) ;
      
      okt = okt && list->item( dd ) == d ;
      if( !okt )
      {
         out() << "Error : Iterator elem " << nbElem << " = " ;
         d->print( out(), 0 ) ;
      }
       ok = ok && okt ;
      nbElem++ ;
   }
   bool okt = nbElem==list->count() ;
   if( !okt )
   {
      out() << "Iterator nb elem : " << nbElem
            << " different from list nb elem : " << list->count() << endl ;
   }
   ok = ok && okt ;
   notify_one_test_result( "iterator", ok ) ;
   
   list->clear() ;
   notify_one_test_result( "clear", list->count()==0 ) ;

   for( int i=0 ; i< 100 ; i++ )
   {
      list->extend( PEL_Double::create( this, 1.0*i) ) ;
   }
   
   while( list->count()>0  )
   {      
      PEL_Iterator* itd = list->create_iterator( this ) ;
      PEL_Double* d = dynamic_cast<PEL_Double*>( itd->item() ) ;
      bool bon = d != 0 ;
      PEL_Double* dd = d->create_clone( this ) ;

      bon = bon && list->has( d ) ;
      if( bon ) list->remove( dd ) ;
      
      if( !bon )
      {
         out() << "Error : Iterator elem " << nbElem << " = " ;
         d->print( out(), 0 ) ;
      }
      
      ok = ok && bon ;
   }

   notify_one_test_result( "remove", ok ) ;

   for( int i=0 ; i< 100 ; i++ )
   {
      list->extend( PEL_Double::create( this, 1.0*i) ) ;
   }
   
   PEL_Collection* list2 = list->create_clone( this ) ;
   PEL_Iterator* it2 = list2->create_iterator( this ) ;
   ok = true ;
   for( it->start() ; it->is_valid() && ok ; it->go_next() ) 
   {
      ok = ok && it2->is_valid() ;
      if( ok )
      {
         ok = ok && it->item()->is_equal( it2->item() ) ;
         it2->go_next() ;
      }
   }
   ok = ok && !it2->is_valid() ;
   notify_one_test_result( "clone", ok ) ;

//   notify_one_test_result( "is_equal", list->is_equal( list2 ) && list2->is_equal( list ) ) ;
//   notify_one_test_result( "hash_code", list->hash_code() == list2->hash_code() ) ;
   
//   list->extend( PEL_Double::create( this, -1) ) ;
//   int ret1 = list->three_way_comparison( list2 ) ;
//   int ret2 = list2->three_way_comparison( list ) ;
   
//   notify_one_test_result( "three_way_comparison", ret1*ret2==-1  ) ;
   
}



//-------------------------------------------------------------------------
PEL_Collection_TEST
const&
PEL_Collection_TEST:: operator=( PEL_Collection_TEST const& thePEL_Collection_TEST )
//-------------------------------------------------------------------------
{
   return( *this ) ;
}


