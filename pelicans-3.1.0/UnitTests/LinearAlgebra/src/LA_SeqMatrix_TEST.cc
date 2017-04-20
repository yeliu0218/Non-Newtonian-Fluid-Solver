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

#include <LA_SeqMatrix_TEST.hh>
#include <LA_Matrix_TEST.hh>

#include <LA_DenseMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <doubleArray2D.hh>

#include <iostream>

LA_SeqMatrix_TEST const* 
LA_SeqMatrix_TEST::PROTOTYPE = new LA_SeqMatrix_TEST()  ;

//-------------------------------------------------------------------------
LA_SeqMatrix_TEST:: LA_SeqMatrix_TEST(void )
//-------------------------------------------------------------------------
   : LA_Matrix_TEST( "LA_SeqMatrix", "LA_SeqMatrix_TEST" )
   , SYM(false)
{
}

//-------------------------------------------------------------------------
LA_SeqMatrix_TEST:: LA_SeqMatrix_TEST( std::string const& matrixClass,
                                 std::string const& registration )
//-------------------------------------------------------------------------
   : LA_Matrix_TEST( matrixClass, registration )
   , SYM(false)
{
}

//-------------------------------------------------------------------------
LA_SeqMatrix_TEST:: ~LA_SeqMatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: test( PEL_ModuleExplorer const* exp,
                          LA_Matrix const* matrix,
                          std::string const& test_matrix )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: test" ) ;
   LA_Matrix_TEST::test(exp,matrix,test_matrix);
   LA_SeqMatrix const* smat = dynamic_cast<LA_SeqMatrix const* >( matrix ) ;
   PEL_ASSERT( smat!=0 ) ;

   SYM = matrix->is_symmetric() ;
   
   doBasicTests() ;
   doBlockVectorTests() ;
   doBlas2Tests() ;
   doBlas3Tests() ;
   
   doBasicSparseTests() ;
   doSolveLUTest(smat) ;
   PEL_Module* mod = PEL_Module::create(this,"a_module") ;
   smat->save(mod) ;
   PEL_ModuleExplorer const* a_exp = PEL_ModuleExplorer::create(this,mod);
   LA_SeqMatrix* cmat = smat->create_matrix(this) ;
   cmat->restore(a_exp) ;      
   notify_one_test_result(  " save and restore",
                            matrix_equality( smat, cmat ) ) ;
      
   
   
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doBasicTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doBasicTests" ) ;
   
   size_t const n = 100 ;
   size_t const m = 100 ;
   double const petit = 1.0e-10 ;
   
   LA_SeqMatrix* mat = create_matrix( 0, n, m ) ;

   notify_one_test_result( "nb_rows", mat->nb_rows()==n ) ;
   notify_one_test_result( "nb_cols", mat->nb_cols()==m ) ;

   
   for( size_t i = 0 ; i<n ; i++ )
   {
      for( size_t j = ( SYM ? i : 0 ) ; j<m ; j++ )
      {     
         mat->set_item( i, j, petit*(i+1)*(j+1) ) ;
      }
   }
   
   for( size_t i = 0 ; i< PEL::min( n, m ) ; i++ )
   {
     mat->set_item( i, i, i ) ;
   }
   for( size_t i = 0 ; i< PEL::min( n, m ) ; i++ )
   {
     mat->add_to_item( i, i, i ) ;
   }
   mat->synchronize() ;
   bool ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      for( size_t i = 0 ; i<n ; i++ )
      {
         double v = mat->item( i, j ) ;
         double v0 ;        
         if( i==j )
         {
            v0 = 2.0*i ;
         }
         else
         {
            v0 = petit*(i+1)*(j+1) ;
         }
         ok &= LA_Matrix_TEST::double_equality( v, v0 ) ;
      }
   }
   notify_one_test_result( "set_item,add_to_item,() ", ok ) ;
   mat->set_item( 0, 0, 1. ) ;
   mat->set_item( 0, 0, 0. ) ;
   mat->synchronize() ;
   notify_one_test_result( "set_item(0) ", mat->item(0,0)==0. ) ;
   size_t nb = mat->nb_stored_items() ;
   mat->nullify() ;
   ok = mat->nb_stored_items()==nb ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      for( size_t i = 0 ; i<n ; i++ )
      {
         double v = mat->item( i, j ) ;
         ok = ok && v==0.0 ;
      }
   }
   notify_one_test_result( "nullify ", ok ) ;

   mat->destroy() ; mat = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doBlockVectorTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doBlockVectorTests" ) ;
   
   size_t const n = 10 ;
   size_t const m = ( SYM ? n : 10 ) ;
   size_t const iRow = n/2 ;
   size_t const iCol = m/2 ;

   LA_SeqMatrix* mat = create_matrix( 0, n, m ) ;
   LA_SeqVector* row = LA_SeqVector::create( mat, m ) ;
   LA_SeqVector* col = LA_SeqVector::create( mat, n ) ;
   for( size_t k=0 ; k<m*n/2 ; k++ )
   {
      size_t i = PEL::rand()%n ;
      size_t j = PEL::rand()%m ;

      if( !SYM || i<=j )
      {
         double v = 1.0*PEL::rand() ;
         mat->set_item( i, j, v ) ;
      
         if( i==iRow )
         {
            row->set_item( j, v ) ;
         }
         if( j==iCol )
         {
            col->set_item( i, v ) ;
         }
         if( SYM )
         {
            if( j==iRow )
            {
               row->set_item( i, v ) ;
            }
            if( i==iCol )
            {
               col->set_item( j, v ) ;
            }
         }
      }
   }
   mat->synchronize() ;
   
   LA_SeqVector* rowExtracted = LA_SeqVector::create( mat, m ) ;
   mat->extract_row( iRow, rowExtracted ) ;
   bool ok = true ;
   for( size_t j=0 ; j<m ; j++ )
   {
      ok &= LA_Matrix_TEST::double_equality( rowExtracted->item(j), row->item(j) ) ;
   }
   notify_one_test_result( "extract_row", ok ) ;

   ok = true ;
   LA_SeqVector* colExtracted = LA_SeqVector::create( mat, n ) ;
   mat->extract_col( iCol, colExtracted ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      ok &= LA_Matrix_TEST::double_equality( colExtracted->item(i), col->item(i) ) ;
   }
   notify_one_test_result( "extract_col", ok ) ;
   
   mat->destroy() ; mat = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doBlas2Tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doBlas2Tests" ) ;
   
   size_t const n = 100 ;
   size_t const m = ( SYM ? n : 200 ) ;
   double const alpha = 0.5 ;
   double const beta = 0.11 ;
   
   PEL::out() << "LA_SeqMatrix_TEST:: doBlas2Tests" << std::endl ;
   
   LA_SeqMatrix * mat = create_matrix( 0, n, m ) ;
   LA_SeqVector * y = LA_SeqVector::create( mat, n ) ;
   LA_SeqVector * x = LA_SeqVector::create( mat, m ) ;
   double res[ n ] ;
   double * res_t = new double [ m ] ;

   for( size_t i = 0 ; i<n ; i++ )
   {
      y->set_item( i, i*1.0/n ) ;
      res[ i ] = (beta * i)/n ;
   }
   for( size_t j = 0 ; j<m ; j++ )
   {
      x->set_item( j, j*1.1/m ) ;
      res_t[ j ] = (beta * 1.1*j)/m ;
   }
   
   for( size_t k=0 ; k<5*(m+n) ; k++ )
   {
      size_t i = PEL::rand()%n ;
      size_t j = PEL::rand()%m ;

      if( !SYM || i<=j )
      {
         double v = PEL::rand() ;
         if( v > 0 )
         {
            v = 1e-3/v ;
         }

         mat->add_to_item( i, j, v ) ;
         
         double r = (v*alpha*1.1*j)/m ;
         res[ i ] += r ;
         
         double r_t = (v*alpha*i)/n ;
         res_t[ j ] += r_t ;

         if( SYM && i!=j ) 
         {
            res[ j ] += (v*alpha*1.1*i)/m ;
            res_t[ i ] += (v*alpha*j)/n ;
         }
      }
      
   }
   mat->synchronize() ;
   x->synchronize() ;
   y->synchronize() ;
   PEL_ASSERT( mat->is_synchronized() ) ;
   
   mat->multiply_vec_then_add( x, y, alpha, beta ) ;
   bool ok = true ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      ok &= LA_Matrix_TEST::double_equality( res[ i ], y->item( i ) ) ;
   }
   notify_one_test_result( "multiply_vec_then_add", ok ) ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      y->set_item( i, i*1.0/n ) ;
   }
   for( size_t j = 0 ; j<m ; j++ )
   {
      x->set_item( j, j*1.1/m ) ;
   }
   x->synchronize() ;
   y->synchronize() ;
   mat->tr_multiply_vec_then_add( y, x, alpha, beta ) ;
   ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      ok &= LA_Matrix_TEST::double_equality( res_t[ j ], x->item( j ) ) ;
   }
   notify_one_test_result( "tr_multiply_vec_then_add", ok ) ;

   mat->destroy() ; mat = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doBlas3Tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doBlas3Tests" ) ;
   
   size_t const n = 10 ;
   size_t const m = ( SYM ? n : 20 ) ;

   double alpha = 0.5 ;
   
   LA_SeqMatrix* matA = create_matrix( 0, n, m ) ;
   doubleArray2D mA(n,m) ;
   fill_sparse_matrix( matA, mA ) ;
   
   LA_SeqMatrix* matB = create_matrix( 0, m, n ) ;
   doubleArray2D mB(m,n) ;
   fill_sparse_matrix( matB, mB ) ;
   
   LA_SeqMatrix* matC = create_matrix( 0, m, n ) ;
   doubleArray2D mC(m,n) ;
   fill_sparse_matrix( matC, mC ) ;
   
//   test_clone( matA ) ;
   
   LA_SeqMatrix* mat = create_matrix( 0, n, m ) ;
   doubleArray2D mm(n,m) ;
   LA_SeqMatrix* mat2 = create_matrix( 0, n, n ) ;
   doubleArray2D mm2(n,n) ;
   LA_SeqMatrix* mat3 = create_matrix( 0, m, m ) ;
   doubleArray2D mm3(m,m) ;

   // add_Mat(sparse)
   {
      fill_sparse_matrix( mat, mm ) ;
      mat->synchronize() ;
      matA->synchronize() ;
      mat->synchronize() ;
      
      mat->add_Mat( matA, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<n ; i++ )
      {
         for( size_t j=0 ; j<m ; j++ )
         {
            double r = mm(i,j)+alpha*mA(i,j)  ;
            ok &= LA_Matrix_TEST::double_equality( r, mat->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_Mat(sparse)", ok ) ;
   }

   // add_tMat(sparse,sparse)
   {
      fill_sparse_matrix( mat, mm ) ;
      mat->synchronize() ;
      matB->synchronize() ;
      mat->add_tMat( matB, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<n ; i++ )
      {
         for( size_t j=0 ; j<m ; j++ )
         {
            double r = mm(i,j)+alpha*mB(j,i)  ;
            ok &= LA_Matrix_TEST::double_equality( r, mat->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_tMat(sparse)", ok ) ;
   }
   
   if( !SYM )
   // add_Mat_Mat(sparse,sparse)
   {
      fill_sparse_matrix( mat2, mm2 ) ;
      mat2->synchronize() ;
      mat2->add_Mat_Mat( matA, matB, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<n ; i++ )
      {
         for( size_t j=0 ; j<n ; j++ )
         {
            double r = mm2(i,j) ;
            for( size_t k=0 ; k<m ; k++ )
            {
               r += alpha * mA(i,k) * mB(k,j) ;
            }
            ok &= LA_Matrix_TEST::double_equality( r, mat2->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_Mat_Mat(sparse,sparse)", ok ) ;
   }

   if( !SYM )
   // add_Mat_tMat(sparse,sparse)
   {
      fill_sparse_matrix( mat3, mm3 ) ;
      mat3->synchronize() ;
      matC->synchronize() ;
      mat3->add_Mat_tMat( matC, matB, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<m ; i++ )
      {
         for( size_t j=0 ; j<m ; j++ )
         {
            double r = mm3(i,j) ;
            for( size_t k=0 ; k<n ; k++ )
            {
               r += alpha * mC(i,k) * mB(j,k) ;
            }
            ok &= LA_Matrix_TEST::double_equality( r, mat3->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_Mat_tMat(sparse,sparse)", ok ) ;
   }
   
   if( !SYM )
   // add_tMat_Mat(sparse,sparse)
   {
      fill_sparse_matrix( mat2, mm2 ) ;
      mat2->synchronize() ;
      mat2->add_tMat_Mat( matC, matB, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<n ; i++ )
      {
         for( size_t j=0 ; j<n ; j++ )
         {
            double r = mm2(i,j) ;
            for( size_t k=0 ; k<m ; k++ )
            {
               r += alpha * mC(k,i) * mB(k,j) ;
            }
            ok &= LA_Matrix_TEST::double_equality( r, mat2->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_tMat_Mat(sparse,sparse)", ok ) ;
   }
   
   if( !SYM )
   // add_tMat_tMat(sparse,sparse)
   {
      fill_sparse_matrix( mat3, mm3 ) ;
      mat3->synchronize() ;
      mat3->add_tMat_tMat( matA, matB, alpha ) ;
      bool ok = true ;
      for( size_t i=0 ; i<m ; i++ )
      {
         for( size_t j=0 ; j<m ; j++ )
         {
            double r = mm3(i,j) ;
            for( size_t k=0 ; k<n ; k++ )
            {
               r += alpha * mA(k,i) * mB(j,k) ;
            }
            ok &= LA_Matrix_TEST::double_equality( r, mat3->item( i, j )  ) ;
         }
      }
      notify_one_test_result( "add_tMat_tMat(sparse,sparse)", ok ) ;
   }
   
   matA->destroy() ; matA = 0 ;
   matB->destroy() ; matB = 0 ;
   matC->destroy() ; matC = 0 ;
   
   mat->destroy() ; mat = 0 ;
   mat2->destroy() ; mat2 = 0 ;
   mat3->destroy() ; mat3 = 0 ;
}


//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doBasicSparseTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doBasicSparseTests" ) ;
   
   const size_t n = 100 ;
   size_t const m = ( SYM ? n : 200 ) ;
   double petit = 1e-10 ;
   
   LA_SeqMatrix* mat = create_matrix( 0, n, m ) ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      for( size_t j = (SYM ? i : 0) ; j<m ; j++ )
      {     
         mat->set_item( i, j, petit*(i+1)*(j+1) ) ;
      }
   }
   for( size_t i = 0 ; i< PEL::min( n, m ) ; i++ )
   {
      mat->set_item( i, i, i ) ;
   }
   mat->synchronize() ;
   for( size_t i = 0 ; i< PEL::min( n, m ) ; i++ )
   {
      mat->add_to_item( i, i, i ) ;
   }
   bool ok = true ;

   mat->synchronize() ;
   LA_MatrixIterator * it = mat->create_stored_item_iterator( mat ) ;

   it->start_all_items() ;
   notify_one_test_result( "Iterator nb_rows ", mat->nb_rows() == n ) ;
   ok = true ;
   
   for( ; it->is_valid() ; it->go_next() )
   {
      size_t i = it->row() ;
      size_t j = it->col() ;
      double v = it->item() ;
      double v0 ;        
      if( i==j )
      {
         v0 = 2.0*i ;
      }
      else
      {
         v0 = petit*(i+1)*(j+1) ;
      }
      ok &= LA_Matrix_TEST::double_equality( v, v0 ) ;      
   }
   notify_one_test_result( "Iterator", ok ) ;
   
   LA_SeqMatrix * mat_set = create_matrix( 0, n, m ) ;
   mat_set->set( mat ) ;
   bool ok_set = true ;
   PEL_ASSERT( mat_set->state() == LA::Sync ) ;
   
   LA_MatrixIterator* it_set = mat_set->create_stored_item_iterator( mat_set ) ;
   for( it_set->start_all_items() ; it_set->is_valid() ; it_set->go_next() )
   {
      size_t i = it_set->row() ;
      size_t j = it_set->col() ;
      if( !SYM || i<=j )
      {
         double v = it_set->item() ;
         double v0 = mat->item(i,j) ;        
         ok_set &= LA_Matrix_TEST::double_equality( v, v0 ) ;
      }
   }
   mat_set->synchronize() ;
   notify_one_test_result( "set", ok_set ) ;

   ok_set = true ;
   for( it_set->start_all_items() ; it_set->is_valid() ; it_set->go_next() )
   {
      size_t i = it_set->row() ;
      size_t j = it_set->col() ;
      if( !SYM || i<=j )
         it_set->set_item( i+j ) ;      
   }
   mat_set->synchronize() ;
   for( it_set->start_all_items() ; it_set->is_valid() ; it_set->go_next() )
   {
      size_t i = it_set->row() ;
      size_t j = it_set->col() ;
      if( !SYM || i<=j )
         ok_set &= LA_Matrix_TEST::double_equality( it_set->item(), i+j ) ;
   }
   notify_one_test_result( "set_item", ok_set ) ;
   for( size_t k = 0 ; k<n+m ; k++ )
   {
      size_t i = PEL::rand()%n ;
      size_t j = PEL::rand()%m ;
      if( !SYM || i<=j )
         mat->set_item( i, j, 1.0 ) ;
   }
   bool ok_fill=true ;
   mat->set_stored_items(PEL::bad_double()) ;
   mat->synchronize() ;
   
   it->start_all_items() ;
   for( it_set->start_all_items() ; it_set->is_valid() ; it_set->go_next() )
   {
      ok_fill = ok_fill && it->item()==PEL::bad_double() ;
   }
   notify_one_test_result( "set_stored_items", ok ) ;
   
   mat->nullify() ;
   
   ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      for( size_t i = 0 ; i<n ; i++ )
      {
         double v = mat->item( i, j ) ;
         ok &= LA_Matrix_TEST::double_equality( 0, v ) ;
      }
   }
   notify_one_test_result( "nullify", ok ) ;

   mat->destroy() ; mat = 0 ;
   mat_set->destroy() ; mat_set = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: doSolveLUTest( LA_SeqMatrix const* matA )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: doSolveLUTest" ) ;
   
   size_t n = matA->nb_rows() ;

   LA_SeqVector * rhs = LA_SeqVector::create( 0, n ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      rhs->set_item( i, 1.0/(i+1) ) ;
   }
   rhs->synchronize() ;

   LA_SeqVector* sol = LA_SeqVector::create( rhs, n ) ;

   matA->solve_LU( rhs, sol ) ;

   LA_SeqVector * x = LA_SeqVector::create( rhs, n ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=i ; j<n ; j++ )
      {
         x->add_to_item( i, matA->item(i,j)*sol->item(j) ) ;
      }
   }
   
   for( size_t i=0 ; i<n ; i++ )
   {
      double xx = x->item(i) ;      
      for( size_t j=0 ; j<i ; j++ )
      {
        xx += matA->item(i,j)*x->item(j) ;
      }
      sol->set_item( i, xx ) ;
   }
   bool ok = true ;
   for( size_t i=0 ; i<n ; i++ )
   {
      ok &= LA_Matrix_TEST::double_equality( sol->item(i), rhs->item(i) ) ;
   }
   notify_one_test_result( "solve_LU", ok ) ;
   
   rhs->destroy() ; rhs = 0 ;
}

//-------------------------------------------------------------------------
void
LA_SeqMatrix_TEST:: fill_sparse_matrix( LA_SeqMatrix* mat,
                                        doubleArray2D& values ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: fill_sparse_matrix" ) ;
   
   PEL_CHECK( mat != 0 ) ;

   int const n = mat->nb_rows() ;
   int const m = mat->nb_cols() ;

   // Nullify:
   values.re_initialize( n, m ) ;
   mat->nullify() ;

   // Set random values :
   for( int k=0 ; k<(n/4)*(m+n) ; ++k )
   {
      int const i1 = PEL::rand()%n ;
      int const j1 = PEL::rand()%m ;
      int const i2 = PEL::rand()%n ;
      int const j2 = PEL::rand()%m ;
      if( !SYM || ( i1 <= j1 && i2<= j2 ) )
      {
         double v = 1.*( PEL::rand()%100 ) ;
         mat->add_to_item( i1, j1, v ) ;
         values(i1,j1) += v ;
         double vinv = v ;
         if( v > 0. )
         {
            vinv = 1.0/v ;
         }
         mat->add_to_item( i2, j2, vinv ) ;
         values(i2,j2) += vinv ;
         if( SYM ) 
         {
            if( i1!=j1 )
               values(j1,i1) += v ;
            if( i2!=j2 )
               values(j2,i2) += vinv ;
         }
      }      
   }
   
   for( size_t i=0 ; i<(size_t) n ; ++i )
   {
      for( size_t j=0 ; j<(size_t) m ; ++j )
      {
         if( SYM )
            PEL_ASSERT( mat->item(i,j)==mat->item(j,i) ) ;
         PEL_ASSERT( values(i,j)==mat->item(i,j) ) ;
      }
   }
}

//-------------------------------------------------------------------------
LA_SeqMatrix*
LA_SeqMatrix_TEST:: create_matrix( PEL_Object* a_owner,
                                   size_t n, size_t m ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix_TEST:: create_matrix" ) ;
   
   LA_Matrix* mat = LA_Matrix_TEST:: create_matrix( a_owner, n, m ) ;
   LA_SeqMatrix* result = dynamic_cast< LA_SeqMatrix* >( mat ) ;
   if( result == 0 )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_SeqMatrix_TEST error:\n"
         "    matrix prototype is not LA_SeqMatrix" ) ;
   }
   PEL_ASSERT( result != 0 ) ;
   
   return( result ) ;
}

