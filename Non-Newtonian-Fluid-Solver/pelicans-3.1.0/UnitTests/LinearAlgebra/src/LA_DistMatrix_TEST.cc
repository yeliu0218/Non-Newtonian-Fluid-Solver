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

#include <LA_DistMatrix_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>
#include <PEL_System.hh>
#include <doubleVector.hh>
#include <intVector.hh>

#include <LA_DistMatrix.hh>
#include <LA_DistVector.hh>

#include <iostream>
#include <iomanip>

LA_DistMatrix_TEST const*
LA_DistMatrix_TEST:: PROTOTYPE = new LA_DistMatrix_TEST() ;

//-------------------------------------------------------------------------
LA_DistMatrix_TEST:: LA_DistMatrix_TEST( void )
//-------------------------------------------------------------------------
   : LA_Matrix_TEST( "LA_DistMatrix", "LA_DistMatrix_TEST" )
{
}

//-------------------------------------------------------------------------
LA_DistMatrix_TEST:: ~LA_DistMatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_DistMatrix_TEST:: test( PEL_ModuleExplorer const* exp,
                           LA_Matrix const* matrix,
                           std::string const& test_matrix_name )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistMatrix_TEST:: process_all_tests" ) ;

   LA_Matrix_TEST::test( exp, matrix, test_matrix_name ) ;

   LA_DistMatrix const* dmat = dynamic_cast<LA_DistMatrix const*>( matrix) ;
   notify_one_test_result("type", dmat ) ;

   com = PEL_Exec::communicator() ;

   out() << "| Communicator is " << com->name() << std::endl ;

   rank = com->rank() ;
   size = com->nb_ranks() ;
   size_t const N = 10 ;


   LA_DistMatrix* local_sized_matrix =
      LA_DistMatrix::create( this, PEL::bad_index(), PEL::bad_index(), 
                             N, N, LA::FromLocalSize ) ;
   test_matrix( local_sized_matrix, "local sized", N, N, N*size, N*size ) ;

   local_sized_matrix->re_initialize( PEL::bad_index(), PEL::bad_index(),
                                      rank, rank ) ;
   test_matrix( local_sized_matrix, "resized local sized", rank, rank,
                (size*(size-1))/2, (size*(size-1))/2 ) ;

   LA_DistMatrix* global_sized_matrix =
      LA_DistMatrix::create( this, N*size, N*size,
                                   PEL::bad_index(), PEL::bad_index(),
                                   LA::FromGlobalSize ) ;
   test_matrix( global_sized_matrix, "local sized", N, N, N*size, N*size ) ;

   global_sized_matrix->re_initialize( 2*N*size, 2*N*size,
                                       PEL::bad_index(), PEL::bad_index() ) ;
   test_matrix( global_sized_matrix, "resized global sized", 2*N, 2*N,
                2*N*size, 2*N*size ) ;

   global_sized_matrix->re_initialize( size-1, size-1,
                                       PEL::bad_index(), PEL::bad_index() ) ;
   test_matrix( global_sized_matrix, "under sized", (rank==0?size-1:0),
                (rank==0?size-1:0), size-1, size-1 ) ;

   local_sized_matrix->re_initialize( PEL::bad_index(), PEL::bad_index(),
                                      N, 2*N ) ;
   test_matrix( local_sized_matrix, "larger", N, 2*N, N*size, 2*N*size ) ;

   local_sized_matrix->re_initialize( PEL::bad_index(), PEL::bad_index(),
                                      2*N, N ) ;
   test_matrix( local_sized_matrix, "higher", 2*N, N, 2*N*size, N*size ) ;

   if( size>1 ) process_parallel_tests() ;

}

//-------------------------------------------------------------------------
void
LA_DistMatrix_TEST:: test_matrix( LA_DistMatrix* mat,
                                         std::string const& test_name,
                                         size_t local_row_size,
                                         size_t local_col_size,
                                         size_t global_row_size,
                                         size_t global_col_size )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistMatrix_TEST:: test_matrix" ) ;

   out() << "| Testing " << test_name << std::endl ;
   notify_one_test_result(" nb_rows ",
                          mat->nb_rows()==global_row_size ) ;
   notify_one_test_result(" nb_cols ",
                          mat->nb_cols()==global_col_size ) ;

   notify_one_test_result("first_local_index <= local_index_limit  ",
                          mat->row_distribution()->first_local_index()
                          <= mat->row_distribution()->local_index_limit() ) ;

   notify_one_test_result("local row size ",
                          mat->row_distribution()->local_index_limit()
                        - mat->row_distribution()->first_local_index()
                          == local_row_size ) ;

   notify_one_test_result("local col size ",
                          mat->col_distribution()->local_index_limit()
                        - mat->col_distribution()->first_local_index()
                          == local_col_size ) ;
   mat->synchronize() ;

   LA_DistMatrix* mat2 = mat->create_matrix( mat ) ;
   bool ok = mat->nb_rows() == mat2->nb_rows()
          && mat->nb_cols() == mat2->nb_cols()
          && mat->row_distribution()->is_compatible( mat2->row_distribution() )
          && mat->col_distribution()->is_compatible( mat2->col_distribution() )
          && !mat2->is_synchronized() ;
   notify_one_test_result(" create_matrix  ", ok ) ;

   for( size_t i=0 ; i<global_row_size ; i++ )
   {
      mat->add_to_item( i, 0, (1.0)*rank ) ;
   }
   if( global_row_size > 0 )
   {
      notify_one_test_result( " !synchronized", mat->state() != LA::Sync ) ;
      mat->synchronize() ;
   }
   notify_one_test_result( " synchronized", mat->is_synchronized() ) ;
   ok = true ;
   double S = (size*(size-1))/2.0 ;

   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok = ok &&
         mat->item( i, 0 ) == S ;
   }
   notify_one_test_result(" add_item", ok ) ;
   mat->nullify() ;
   size_t local_size = PEL::min( local_row_size, local_col_size ) ;

   if( global_row_size > 0 )
   {
      mat->set_item( 0, 0, 1.0 ) ;
      mat->synchronize() ;
      ok = ( mat->row_distribution()->first_local_index() == 0 &&
             local_size > 0 ? mat->item(0,0) == 1.0 : true ) ;
      com->boolean_and( ok ) ;

      notify_one_test_result(" add one item", ok ) ;
   }

   for( size_t i=0 ; i<global_row_size ; i++ )
   {
      mat->set_item( i, 0, (1.0)*i ) ;
   }
   if( global_row_size > 0 )
   {
      mat->synchronize() ;
   }
   ok = true ;
   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() ;
        i++ )
   {
      ok = ok &&
         mat->item( i, 0 ) == i ;
   }
   notify_one_test_result(" set_item", ok ) ;

   size_t global_size = PEL::min( global_row_size, global_col_size ) ;

   for( size_t i=0; i<global_size ; i++ )
   {
      mat->set_item( i, i, 2.0 ) ;
      mat2->set_item( i , i, ( i%2==0 ? 1.0 : -1.0 ) ) ;
   }
   mat->synchronize() ;
   mat2->synchronize() ;

   mat->scale( 0.5 ) ;

   ok = true ;
   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() && i<global_size ;
        i++ )
   {
      ok = ok &&
         mat->item( i, i ) == 1.0 ;
   }
   notify_one_test_result(" scale ", ok ) ;

   mat->set( mat2 ) ;
   mat->synchronize() ;
   ok = true ;
   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() && i<global_size ;
        i++ )
   {
      ok = ok &&
         mat->item( i, i ) == mat2->item(i,i) ;
   }

   notify_one_test_result(" set ", ok ) ;


   mat->set( mat2 ) ;
   mat->add_Mat( mat2, 1.5 ) ;
   mat->synchronize() ;

   ok = true ;
   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() && i<global_size;
        i++ )
   {
      ok = ok &&
         mat->item( i, i ) == 2.5*mat2->item(i, i) ;
   }
   notify_one_test_result(" add_Mat ", ok ) ;

   LA_DistVector* vec_x =
      LA_DistVector::create( mat, PEL::bad_index(), local_row_size, 
                             LA::FromLocalSize ) ;
   LA_DistVector* vec_y =
      LA_DistVector::create( mat, PEL::bad_index(), local_col_size, 
                             LA::FromLocalSize ) ;
   mat->nullify() ;

   size_t hd = global_size/2 ;
   double d = 1.0 ;
   double dhd = 0.5 ;


   for( size_t i=0; i<global_size ; i++ )
   {
      mat->set_item( i, hd, dhd ) ;
      mat->set_item( i, i, d*i ) ;
   }
   for( size_t i=0; i<global_row_size ; i++ )
   {
      vec_x->set_item( i, 1.0*(i) ) ;
   }
   for( size_t i=0; i<global_col_size ; i++ )
   {
      vec_y->set_item( i, 1.0*(i) ) ;
   }
   mat->synchronize() ;
   vec_x->synchronize() ;
   vec_y->synchronize() ;
   double alpha = 1.5 ;
   double beta = 2.3 ;

   mat->multiply_vec_then_add( vec_y, vec_x, alpha, beta ) ;
   ok = true ;
   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() ;
        i++ )
   {
      double aw = beta*i ;
      if( i<global_size )
      {
         aw += alpha*i*d*i ;
         if( i!=hd ) aw += alpha*dhd*hd ;
      }


      bool okt = PEL::double_equality(vec_x->item(i),aw,1.0e-15,1.0e-15)  ;
      if( !okt ) out() << " |ERR res=" << vec_x->item( i )
                       << " waited=" << aw <<std::endl ;

      ok = ok && okt ;
   }
   notify_one_test_result(" multiply_vec_then_add ", ok ) ;

   for( size_t i=0; i<global_col_size ; i++ )
   {
      vec_y->set_item( i, 1.0*(i) ) ;
   }
   for( size_t i=0; i<global_row_size ; i++ )
   {
      vec_x->set_item( i, 1.0*(i) ) ;
   }
   vec_x->synchronize() ;
   vec_y->synchronize() ;

//    out() << " vec_y avant tr_multiply_vec_then_add : " << std::endl ;
//    vec_y->print( out(), 1 ) ;

   mat->tr_multiply_vec_then_add( vec_x, vec_y, alpha, beta ) ;

   ok = true ;
   for( size_t i=mat->col_distribution()->first_local_index() ;
        i<mat->col_distribution()->local_index_limit() ;
        i++ )
   {
      double aw = beta*i ;
      if( i<global_size ) aw += alpha*i*d*i ;
      if( i==hd ) aw += alpha*(
         dhd*( global_size*(global_size-1)/2.0 - hd ) ) ;

      bool okt = PEL::double_equality(vec_y->item(i),aw,1.0e-15,1.0e-15)  ;
      if( !okt ) out() << " |ERR res=" << vec_y->item( i )
                       << " waited=" << aw <<std::endl ;

      ok = ok && okt ;
   }
   if( !ok )
   {
      out() << " vec_x : " << std::endl ;
      vec_x->print( out(), 1 ) ;

      out() << " vec_y apres tr_multiply_vec_then_add : " << std::endl ;
      vec_y->print( out(), 1 ) ;

      out() << " mat : " << std::endl ;
      mat->print( out(), 1 ) ;
   }

   notify_one_test_result(" tr_multiply_vec_then_add ", ok ) ;

   mat2->re_initialize( global_row_size, global_col_size,
                        local_col_size, local_row_size ) ;

   mat2->add_tMat( mat ) ;
   mat2->synchronize() ;

   ok = true ;

   for( size_t i=mat->row_distribution()->first_local_index() ;
        i<mat->row_distribution()->local_index_limit() ;
        i++ )
   {
      for( size_t j=mat->col_distribution()->first_local_index() ;
           j<mat->col_distribution()->local_index_limit() ;
           j++ )
      {
         bool okt = mat->item(i,j)==mat2->item(j,i)  ;
         if( !okt )
            out() << " |ERR mat=" << mat->item(i,j)
                  << " mat2=" << mat2->item(j,i) <<std::endl ;

         ok = ok && okt ;
      }
   }
   notify_one_test_result("  add_tMat ", ok ) ;

}



//-------------------------------------------------------------------------
void
LA_DistMatrix_TEST:: process_parallel_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistMatrix_TEST:: process_parallel_tests" ) ;

   com = PEL_Exec::communicator() ;

   out() << "| Communicator is " << com->name() << std::endl ;

   rank = com->rank() ;
   size = com->nb_ranks() ;
   size_t tot = 51 ;
   double const P = 1.0e-20 ;

   LA_DistMatrix* mat =
      LA_DistMatrix::create( this, tot, tot,
                             PEL::bad_index(), PEL::bad_index(),
                             LA::FromGlobalSize ) ;
   LA_DistMatrix* mat1 =
      LA_DistMatrix::create( this, tot, tot,
                             PEL::bad_index(), PEL::bad_index(),
                             LA::FromGlobalSize ) ;
   LA_DistVector* vec_x =
      LA_DistVector::create( mat, mat->nb_rows(), PEL::bad_index(),
                             LA::FromGlobalSize ) ;
   LA_DistVector* vec_y = vec_x->create_vector( vec_x ) ;
   size_t deb = mat->row_distribution()->first_local_index() ;
   size_t end = mat->row_distribution()->local_index_limit() ;

   for( size_t attempt=0 ; attempt<2 ; attempt++ )
   {
      mat->nullify() ;
      mat1->nullify() ;

      for( size_t i=0 ; i<tot ; i++ )
      {
         for( size_t j=0 ; j<tot ; j++ )
         {
            if( j>i) mat->add_to_item(i,j,P/(j+1)) ;
            else mat1->add_to_item(i,j,P/(j+1)) ;
         }
         vec_x->set_item(i, 1.0*(i+1) ) ;
      }
      vec_x->synchronize() ;
      mat->synchronize() ;
      mat1->synchronize() ;
   }
   mat->add_Mat( mat1 ) ; // Scattering should be adapted

   mat->multiply_vec_then_add( vec_x, vec_y ) ;
   double R = P*size*(tot*1.0);
   bool ok  = true ;
   for( size_t i=deb ; i<end ; i++ )
   {
      for( size_t j=0 ; j<tot ; j++ )
      {
         double M = P*(size*1.0)/(j+1) ;
         bool okm = PEL::abs( mat->item(i,j) - M ) < M*1.0e-14 ;
         if( !okm ) out() << std::setprecision(16) << " mat (" << i << ","
                          << j <<") = " << mat->item(i,j)
                          << " Expected = " << M << std::endl ;
         ok = ok && okm ;
      }
      double V = 1.0*(i+1) ;
      bool okv = PEL::abs( vec_x->item(i) - V ) < V*1.0e-14 ;
      if( !okv ) out() << std::setprecision(16) << " vec_x (" << i << ") = "
                       << vec_x->item(i) << " Expected = " << V << std::endl ;

      bool okt = PEL::abs( vec_y->item(i) - R ) < R*1.0e-14 ;
      if( !okt ) out() << std::setprecision(16) << " vec_y (" << i << ") = "
                       << vec_y->item(i) << " Expected = " << R << std::endl ;
      ok = ok && okt ;
   }
   notify_one_test_result(" multiply_vec_then_add  ", ok ) ;

}

