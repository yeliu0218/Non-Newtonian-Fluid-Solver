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

#include <LA_SeqVector_TEST.hh>

#include <LA_SeqVector.hh>

#include <PEL_Root.hh>
#include <PEL.hh>

#include <size_t_vector.hh>

//-------------------------------------------------------------------------
LA_SeqVector_TEST const* 
LA_SeqVector_TEST::PROTOTYPE = new LA_SeqVector_TEST()  ;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
LA_SeqVector_TEST:: LA_SeqVector_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest("LA_SeqVector", "LA_SeqVector_TEST" )
{
}

//-------------------------------------------------------------------------
LA_SeqVector_TEST:: ~LA_SeqVector_TEST( void )
//-------------------------------------------------------------------------
{
   PROTOTYPE = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqVector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   doBasicTests() ;
   doBlockTests() ;
   doBlas1Tests() ;
   doStatisticalTests() ;
}

//-------------------------------------------------------------------------
void
LA_SeqVector_TEST:: doBasicTests( void )
//-------------------------------------------------------------------------
{
   size_t const dim = 20 ;
   
   LA_SeqVector* vec = createVector( dim ) ;
   notify_one_test_result( "Constructor size", vec->nb_rows()==dim );
   
   double sum = 0.0 ;
   double sum1 = 0.0 ;
   for( size_t i=0 ; i<dim ; i++ )
   {
      sum += PEL::abs( vec->item(i) ) ;
      vec->set_item(i , 1.0*i ) ;
      sum1 += vec->item(i) ;
   }
   notify_one_test_result( "Constructor initialisation", PEL::toler(sum) );
   notify_one_test_result( "Accessors-Modifiers",
                        PEL::equal( sum1, dim*(dim-1)/2.0 ) ) ;

   size_t newDim = dim ;
   newDim++ ;
   vec->re_initialize( newDim ) ;
   notify_one_test_result( "Reshape size", vec->nb_rows()==newDim );
   
   sum = 0.0 ;
   vec->synchronize() ;
   for( size_t i=0 ; i<newDim ; i++ )
   {
      sum += PEL::abs( vec->item(i) ) ;
      vec->add_to_item(i , 1.0*i ) ;
   }
   notify_one_test_result( "Constructor reInit initialisation", PEL::toler(sum) );
   vec->synchronize() ;
   LA_SeqVector* vecBis = createVector( newDim ) ;
   vecBis->set( vec ) ;
   bool ok = true ;
   for( size_t i=0 ; i<newDim ; i++ )
   {
      ok = ok && PEL::equal( vec->item(i), vecBis->item(i) ) ;
   }
   sum/=(dim+1) ;
   notify_one_test_result( "Vector copy", ok );

   vec->nullify() ;
   ok = true ;
   for( size_t i=0 ; i<newDim ; i++ )
   {
      ok = ok && PEL::toler( vec->item(i) ) ;
   }
   notify_one_test_result( "Nullify", ok );
   
}

//-------------------------------------------------------------------------
void
LA_SeqVector_TEST:: doBlockTests( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_SeqVector_TEST:: doBlas1Tests( void )
//-------------------------------------------------------------------------
{
   const int m = 100 ;

   LA_SeqVector* vec1 = createVector(  m ) ;
   LA_SeqVector* vec2 = createVector(  m ) ;

   for( int i=0 ; i<m ; i++ )
   {
      vec1->set_item( i, 0.01*i ) ;
   }
   vec1->synchronize() ;

   vec2->set( vec1 ) ;
   
   double alpha = 0.314 ;
   vec1->scale( alpha ) ;
   bool ok = true ;
   for( int i=0 ; i<m ; i++ )
   {
      ok &= PEL::equal( vec1->item( i ),
                        i*0.01*alpha );
   }
   notify_one_test_result( "scale(a)", ok );
   
   vec2->sum( vec1, alpha ) ;
   ok = true ;
   for( int i=0 ; i<m ; i++ )
   {
      ok &= PEL::equal( vec2->item( i ),
                        i*0.01*(1.0+alpha*alpha) );
   }
   notify_one_test_result( "sum", ok );

   double dot_prod = 0.0 ;
   for( int i=0 ; i<m ; i++ )
   {
      dot_prod += vec1->item( i ) * vec2->item( i ) ;
   }
   notify_one_test_result( "dot product", PEL::equal( dot_prod, vec1->dot( vec2 ) ) );

   notify_one_test_result( "2 norm", PEL::equal( vec1->two_norm(), PEL::sqrt( vec1->dot( vec1 ) ) ) );

   double grand = 1e6 ;
   vec1->set_item( 10, grand ) ;
   vec1->synchronize() ;
   
   notify_one_test_result( "max norm", PEL::equal( vec1->max_norm(), grand ) ) ;   

}

//-------------------------------------------------------------------------
void
LA_SeqVector_TEST:: doStatisticalTests( void )
//-------------------------------------------------------------------------
{
   const int m = 100 ;

   LA_SeqVector* vec = createVector(  m ) ;
   double mean = 0 ;
   for( int i=0 ; i<m ; i++ )
   {
      vec->set_item( i, 1e-3*i ) ;
      mean += vec->item( i ) ;
   }
   mean /= m ;
   
   notify_one_test_result( "Mean value", PEL::equal( vec->mean(), mean ) ) ;
   double deviation = 0 ;
   for( int i=0 ; i<m ; i++ )
   {
      deviation += PEL::pow( mean-vec->item( i ), 2 ) ;
   }
   deviation = PEL::sqrt( deviation/m );
   notify_one_test_result( "Standard deviation",
               PEL::equal( vec->standard_deviation(), deviation ) ) ;

}

//-------------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector_TEST:: createVector( int m )
//-------------------------------------------------------------------------
{
   LA_SeqVector* vector = LA_SeqVector::create( this, m ) ;
   return( vector ) ;
}
