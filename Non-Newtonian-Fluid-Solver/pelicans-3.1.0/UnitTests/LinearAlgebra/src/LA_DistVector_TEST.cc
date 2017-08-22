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

#include <LA_DistVector_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>
#include <doubleVector.hh>
#include <intVector.hh>

#include <LA_DistVector.hh>
#include <LA_DistScatter.hh>

#include <iostream>

//-------------------------------------------------------------------------
LA_DistVector_TEST* 
LA_DistVector_TEST:: REGISTRATOR = new LA_DistVector_TEST() ;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
LA_DistVector_TEST:: LA_DistVector_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "LA_DistVector", "LA_DistVector_TEST" )
{
}

//-------------------------------------------------------------------------
LA_DistVector_TEST:: ~LA_DistVector_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_DistVector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector_TEST:: process_all_tests" ) ;

   com = PEL_Exec::communicator() ;

   out() << "| Communicator is " << com->name() << std::endl ;

   rank = com->rank() ;
   size = com->nb_ranks() ;
   size_t const N = 10 ;

   
   LA_DistVector* local_sized_vector =
      LA_DistVector::create( this, PEL::bad_index(), N, 
                             LA::FromLocalSize ) ;
   
   test_vector( local_sized_vector, "local sized", N, N*size ) ;

   local_sized_vector->re_initialize( PEL::bad_index(), 2*N ) ;
   
   test_vector( local_sized_vector, "resized local sized", 2*N, 2*N*size ) ;
   
   local_sized_vector->re_initialize( PEL::bad_index(), N*rank ) ;
   
   test_vector( local_sized_vector, "non-uniform local sized", N*rank, N*size*(size-1)/2 ) ;  

   LA_DistVector* global_sized_vector =
      LA_DistVector::create( this, N*size, PEL::bad_index(), 
                             LA::FromGlobalSize ) ;
   test_vector( global_sized_vector, "global sized", N, N*size ) ;
   
   global_sized_vector->re_initialize( 2*N*size, PEL::bad_index() ) ; 
   test_vector( global_sized_vector, "resized global sized", 2*N, 2*N*size ) ;

   
   global_sized_vector->re_initialize( size-1, PEL::bad_index() ) ; 
   test_vector( global_sized_vector, "under sized", (rank==0?size-1:0), size-1 ) ;

   global_sized_vector = LA_DistVector::create( this, 0, PEL::bad_index(), 
                                                LA::FromGlobalSize ) ;
   test_vector( global_sized_vector, "null sized", 0, 0 ) ;
   
   global_sized_vector->re_initialize( N*size, PEL::bad_index() ) ; 
   test_vector( global_sized_vector, "resized null sized", N, N*size ) ;
   
   global_sized_vector = LA_DistVector::create( this, PEL::bad_index(), 0, 
                                                LA::FromLocalSize ) ;
   test_vector( global_sized_vector, "local null sized", 0, 0 ) ;
}

//-------------------------------------------------------------------------
void
LA_DistVector_TEST:: test_vector( LA_DistVector* vec,
                                  std::string const& test,
                                  size_t local_size,
                                  size_t global_size )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector_TEST:: test_vector" ) ;

   out() << "| Testing " << test << std::endl ;
   notify_one_test_result(" local_nb_rows ",
                          vec->row_distribution()->local_number()==local_size ) ;
   notify_one_test_result(" nb_rows ",
                          vec->nb_rows()==global_size ) ;
   notify_one_test_result(" synchronized inital ",
                          vec->is_synchronized() ) ;
   
   notify_one_test_result("first_local_index <= local_index_limit  ",
                          vec->row_distribution()->first_local_index() <= vec->row_distribution()->local_index_limit() ) ;
   
   notify_one_test_result("index_limit - first = local size ",
                          vec->row_distribution()->local_index_limit() - vec->row_distribution()->first_local_index()
                          == local_size ) ;
   
   LA_DistVector* vec2 = vec->create_vector( vec ) ;
   notify_one_test_result(" create_vector nb_rows ",
                          vec->nb_rows()==vec2->nb_rows() ) ;
   notify_one_test_result(" create_vector local_nb_rows",
                          vec->row_distribution()->local_number()==vec2->row_distribution()->local_number() ) ;

   
   for( size_t i=0 ; i<global_size ; i++ )
   {
      vec->add_to_item( i, (1.0)*rank ) ;
   }
   vec->synchronize() ;
   notify_one_test_result(" synchronized",
                          (vec->is_synchronized()) ) ;
   bool ok = true ;
   double S = (size*(size-1))/2.0 ;
   
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == S ;
   }
   if( !ok ) vec->print( out(), 0 ) ;
   
   notify_one_test_result(" add_to_item", ok ) ;
   
   for( size_t i=vec->row_distribution()->first_local_index() ; i<vec->row_distribution()->local_index_limit() ; i++ )
   {
      double dd = 0.5 ;
      vec->set_item( i, dd ) ;
   }
   vec->synchronize() ;
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == 0.5 ;
   }
   notify_one_test_result(" set local items", ok ) ;
   
   vec->nullify() ;

   if( global_size > 0 )
   {
      vec->set_item( 0, 1.0 ) ;
      vec->synchronize() ;
      ok = ( vec->row_distribution()->first_local_index()==0 && local_size > 0 ? vec->item(0)==1.0 : true ) ;
      com->boolean_and( ok ) ;
      
      notify_one_test_result(" add one item", ok ) ;
     
   }
   
   for( size_t i=0 ; i<global_size ; i++ )
   {
      vec->set_item( i, (1.0)*i ) ;
   }
   vec->synchronize() ;
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == i ;
   }
   notify_one_test_result(" set_item", ok ) ;

   vec2->nullify() ;
   notify_one_test_result(" nullify is synchronized", vec2->is_synchronized() ) ;
   
   for( size_t i=0; i<global_size ; i++ )
   {
      vec->set_item( i, 2.0 ) ;
      vec2->set_item( i , ( i%2==0 ? 1.0 : -1.0 ) ) ;
   }
   vec->synchronize() ;
   vec2->synchronize() ;

   if( global_size > 0 )
   {
      ok = vec->dot( vec )==global_size*4.0  &&
         vec->dot( vec2 )==( global_size % 2==0 ? 0 : 2.0 ) ;
      notify_one_test_result(" dot", ok ) ;
      double waited_2_norm = PEL::sqrt( global_size * 4.0 ) ;

      ok = PEL::double_equality(vec->two_norm(),waited_2_norm,1.0e-15,0.0) ;
      if( !ok )
         out() << "Calculated 2-norm : " << vec->two_norm() << " waited " << waited_2_norm 
               << " eps= " << vec->two_norm()-waited_2_norm << std::endl ;
   
      notify_one_test_result(" two_norm", ok ) ;
   
      notify_one_test_result(" max_norm", vec2->max_norm()==1.0 ) ;
   }
   
   vec->scale( 0.5 ) ;
   
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == 1.0 ;
   }
   notify_one_test_result(" scale ", ok ) ;
   
   vec->set( vec2 ) ;
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == vec2->item(i) ;
   }

   notify_one_test_result(" set ", ok ) ;
   
   vec->set_as_v_product( vec, vec2 ) ;
   ok = true ;
   
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == vec2->item(i)*vec2->item(i) ;
   }
  
   notify_one_test_result(" set_as_v_product ", ok ) ;

   vec->set( vec2 ) ;
   vec->sum( vec2 ) ;
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == 2.0*vec2->item(i) ;
   }
   notify_one_test_result(" sum ", ok ) ;

   vec->sum( vec2, 3.0 ) ;
   ok = true ;
   for( size_t i=vec->row_distribution()->first_local_index() ;
        i<vec->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok &= vec->item( i ) == 5.0*vec2->item(i) ;
   }
   notify_one_test_result(" sum (alpha) ", ok ) ;

   if( global_size>7 )
   {
      size_t_vector sub( 0 ) ;
      size_t_vector local( 0 ) ;
      for( size_t i=0 ; i*7<global_size ; i ++ )
      {
         sub.append( 7*i ) ;
         local.append(i) ;
      }

      for( size_t i=0; i<global_size ; i++ )
      {
         vec->set_item( i, 1.0*i ) ;
      }
      vec->synchronize() ;
      
      LA_SeqVector* lavec = LA_SeqVector::create( vec, sub.size() ) ;
      LA_SeqVector* global = LA_SeqVector::create( vec, global_size ) ;
      vec->recover_global_vector( global ) ;
      ok = true ;
      for( size_t i=0 ; i<global->nb_rows() ; i++ )
      {
         ok &= global->item( i ) == 1.0*i ;
      }
      notify_one_test_result(" recover_global_vector ", ok ) ;
      
      vec->recover_global_vector( global ) ;
      LA_DistScatter* scatter = LA_DistScatter::create(
         this, vec->row_distribution() ) ;
      scatter->set_unsorted( sub, local ) ;
      scatter->get( vec, lavec ) ;
      
      // vec->get( lavec, sub ) ;
      ok = true ;
      for( size_t i=0 ; i<sub.size() ; i++ )
      {
         ok &= lavec->item( i ) == global->item( sub(i) ) ;
         lavec->set_item( i, 0.0 ) ;
      }
      notify_one_test_result(" get ", ok ) ;

      vec->set( lavec, sub ) ;
      vec->synchronize() ;

      vec->recover_global_vector( global ) ;
      ok = true ;
      for( size_t i=0 ; i<global_size ; i++ )
      {
         ok &= global->item( i ) == ( i%7 == 0 ? 0.0 : i ) ;
      }
      notify_one_test_result(" set ", ok ) ;
   }
}
