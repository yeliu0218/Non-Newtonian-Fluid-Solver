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

#include <PEL_NumberedDoubleVectors_TEST.hh>

#include <PEL_DoubleComparatorExact.hh>
#include <PEL_NumberedDoubleVectors.hh>
#include <PEL_Randomizer.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

#include <iostream>

PEL_NumberedDoubleVectors_TEST* 
PEL_NumberedDoubleVectors_TEST:: UNIQUE_INSTANCE =
                                       new PEL_NumberedDoubleVectors_TEST() ;

//-------------------------------------------------------------------------
PEL_NumberedDoubleVectors_TEST:: PEL_NumberedDoubleVectors_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_NumberedDoubleVectors", 
                     "PEL_NumberedDoubleVectors_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_NumberedDoubleVectors_TEST:: ~PEL_NumberedDoubleVectors_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_NumberedDoubleVectors_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors_TEST:: process_all_tests" ) ;

   PEL_DoubleComparator const* dbl_comp = PEL_DoubleComparatorExact::object() ;
   
   size_t const& nb_items = 10 ;
   PEL_Randomizer* ran = PEL_Randomizer::create( this, 1 ) ;
   
   ran->start() ;
   
   for( size_t d=1 ; d<5 ; d++ )
   {
      PEL_NumberedDoubleVectors* num = 
         PEL_NumberedDoubleVectors::create( this, dbl_comp, d ) ;
      doubleArray2D array( d, nb_items ) ;
      for( size_t j=0 ; j<nb_items ; j++ )
         for( size_t i=0 ; i<d ; i++ )
         {
            array(i,j) = ran->item() ;
            ran->go_next() ;
         }
      
      fill( num, array ) ;
      doubleArray2D ordered = num->ordered_items() ;
      size_t n_ordered = num->nb_items() ;
      
      destroy_possession( num ) ;
      num = PEL_NumberedDoubleVectors::create( this, dbl_comp, d ) ;
      
      size_t_vector perm( nb_items ) ;
      ran->build_permutation( perm ) ;
      
      doubleArray2D permuted( d, nb_items ) ;
      for( size_t j=0 ; j<nb_items ; j++ )
      {
         size_t jp = perm(j) ;
            
         for( size_t i=0 ; i<d ; i++ )
            permuted(i,j) = array(i,jp) ;
      }
      
      fill(num, permuted) ;
      doubleArray2D const& ordered2 = num->ordered_items() ;
      bool ok_order = n_ordered == num->nb_items() ;
      notify_one_test_result( "nb_items", ok_order ) ;
      
      if( ok_order )
         for( size_t i=0 ; i<n_ordered ; i++ )
            for( size_t j=0 ; j<d ; j++ )
               ok_order = ok_order && ordered(j,i) == ordered2(j,i) ;
      
      
      notify_one_test_result( "independant order", ok_order ) ;
   }
}

//-------------------------------------------------------------------------
void
PEL_NumberedDoubleVectors_TEST:: fill( PEL_NumberedDoubleVectors* num, 
                                       doubleArray2D const& array )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors_TEST:: fill" ) ;
   
   size_t d = array.index_bound(0) ;
   size_t n = array.index_bound(1) ;
   size_t_vector index(n)  ;
   
   doubleVector vec( d ) ;
   bool ok_fill = true ;
   bool ok_order = true ;
   
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<d ; j++ ) vec(j) = array(j,i) ;
      num->extend(vec) ;
      
      index(i) = num->index(vec) ;
   }
   ok_fill = ok_fill && num->nb_items() <= n ;
   size_t old_nb = num->nb_items() ;
   
   size_t_vector const& perm = num->order() ;
   doubleArray2D const& ordered = num->ordered_items() ;
   doubleVector tmp( d ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<d ; j++ ) vec(j) = array(j,i) ;
      ok_fill = ok_fill && index(i) == num->index(vec) ;
      size_t p = perm( index(i) ) ;
      
      for( size_t j=0 ; j<d ; j++ )
         ok_order = ok_order && ordered( j, p )==vec(j) ;
   }
   ok_fill = ok_fill && old_nb==num->nb_items() ;
   notify_one_test_result( "filling", ok_fill ) ;
   notify_one_test_result( "ordering", ok_order ) ;
   if( !( ok_fill && ok_order ) )
   {
      std::cout << " array "  << array   << std::endl ;
      std::cout << " ordered "<< ordered << std::endl ;
      std::cout << " order "  << perm    << std::endl ;
   }
}
