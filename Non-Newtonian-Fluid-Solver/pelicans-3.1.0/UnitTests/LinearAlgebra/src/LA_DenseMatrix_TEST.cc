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

#include <LA_DenseMatrix_TEST.hh>

#include <LA_SeqMatrix.hh>
#include <LA_DenseMatrix.hh>
#include <LA_SymmetricMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>

#include <ios>
#include <iomanip>
#include <iostream>


using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ;

LA_DenseMatrix_TEST const* 
LA_DenseMatrix_TEST::PROTOTYPE = new LA_DenseMatrix_TEST()  ;

//-------------------------------------------------------------------------
LA_DenseMatrix_TEST:: LA_DenseMatrix_TEST( void )
//-------------------------------------------------------------------------
   : LA_SeqMatrix_TEST( "LA_DenseMatrix", "LA_DenseMatrix_TEST" )
{
}

//-------------------------------------------------------------------------
LA_DenseMatrix_TEST:: ~LA_DenseMatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: test( PEL_ModuleExplorer const* exp,
                            LA_Matrix const* matrix,
                            std::string const& test_matrix )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: test" ) ;
   LA_SeqMatrix_TEST::test(exp,matrix,test_matrix);
   
   doBlas3FullTests() ;
   doDeterminantFullTests(exp->string_data("det_test"),exp->string_data("eig_test")) ;
   doStatisticalTests() ;
}

//-------------------------------------------------------------------------
LA_DenseMatrix*
LA_DenseMatrix_TEST:: create_matrix(
                                  PEL_Object* a_owner, size_t n, size_t m ) 
//-------------------------------------------------------------------------
{
   return( LA_DenseMatrix::create( a_owner, n, m ) ) ;
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doBasicTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doBasicTests" ) ;
   size_t const n = 100 ;
   size_t m = 200 ;
   new int ;
   double petit = 1e-10 ;
   
   LA_DenseMatrix* mat = create_matrix( this, n, m ) ;

   notify_one_test_result( "nb_rows", mat->nb_rows()==n ) ;
   notify_one_test_result( "nb_cols", mat->nb_cols()==m ) ;
   
   for( size_t i = 0 ; i<n ; i++ )
   {
      for( size_t j = 0 ; j<m ; j++ )
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
         ok = ok && PEL::equal( v, v0 ) ;
      }
   }
   notify_one_test_result( "set_item,add_to_item,() ", ok ) ;
   new int ;
   
   mat->nullify() ;
   ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      for( size_t i = 0 ; i<n ; i++ )
      {
         double v = mat->item( i, j ) ;
         ok = ok && v==0.0 ;
      }
   }
   notify_one_test_result( "nullify ", ok ) ;
   new int [10] ;
   
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doBlockVectorTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doBlockVectorTests" ) ;
   const size_t n = 10 ;
   size_t m = 20 ;
   size_t iRow = n/2 ;
   size_t iCol = m/2 ;
   LA_SeqVector * row = LA_SeqVector::create( this, m ) ;
   LA_SeqVector * col = LA_SeqVector::create( this, n ) ;
   new int [10] ;

   LA_DenseMatrix * mat = create_matrix( this, n, m ) ;
   
   for( size_t k=0 ; k<m*n/2 ; k++ )
   {
      size_t i = PEL::rand()%n ;
      size_t j = PEL::rand()%m ;
      
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
   }
   mat->synchronize() ;
   LA_SeqVector* rowExtracted = LA_SeqVector::create( this, m ) ;
   mat->extract_row( iRow, rowExtracted ) ;
   bool ok = true ;
   for( size_t j=0 ; j<m ; j++ )
   {
      bool okt = PEL::equal( rowExtracted->item(j),
                             row->item(j) ) ;
      if( !okt )
      {
         out() << "Row extracted " << rowExtracted->item(j)
               << " differs from exact " << row->item(j) << std::endl ;
      }
      ok = ok && okt ;
   }
   notify_one_test_result( "extract_row", ok ) ;

   ok = true ;
   LA_SeqVector * colExtracted = LA_SeqVector::create( this, n ) ;
   mat->extract_col( iCol, colExtracted ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      bool okt = PEL::equal( colExtracted->item(i),
                             col->item(i) ) ;
      if(  !okt )
      {
         out() << "Col extracted " << colExtracted->item(i)
               << " differs from exact " << col->item(i) << std::endl ;
      }
      ok = ok && okt ;
   }
   notify_one_test_result( "extract_col", ok ) ; 
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doBlas2Tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doBlas2FullTests" ) ;
   const size_t n = 100 ;
   const size_t m = 200 ;
   double petit = 1e-10 ;
   double alpha = 0.5 ;
   double beta = 0.11 ;
   

   LA_DenseMatrix* mat = create_matrix( 0, n, m ) ;
   LA_SeqVector* y = LA_SeqVector::create( mat, n ) ;
   LA_SeqVector* x = LA_SeqVector::create( mat, m ) ;
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

      double v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1e-3/v ;
      }

      if( v > petit )
      {
         mat->add_to_item( i, j, v ) ;
         
         double r = (v*alpha*1.1*j)/m ;
         res[ i ] += r ;
         
         double r_t = (v*alpha*i)/n ;
         res_t[ j ] += r_t ;

      }
   }
   mat->synchronize() ;
   x->synchronize() ;
   y->synchronize() ;
   
   mat->multiply_vec_then_add( x, y, alpha, beta ) ;
   bool ok = true ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      bool okt = PEL::equal( res[ i ], y->item( i ) ) ;
      if(  !okt )
      {
         out() << res[ i ] << " ?=5=? " << y->item( i ) << " : " << okt << std::endl ;
      }
      ok = ok && okt ;
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
      bool okt = PEL::equal( res_t[ j ], x->item( j ) ) ;
      if(  !okt )
      {
         out().precision(10) ;
         out() << res_t[ j ] << " ?=6=? " << x->item( j ) << " : " << okt << std::endl ;
      }
      ok = ok && okt ;
   }
   notify_one_test_result( "tr_multiply_vec_then_add", ok ) ;

   LA_DenseMatrix* mat2 = create_matrix( mat, n, m ) ;
   double m2[n][m] ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         m2[i][j] = 0.0 ;
      }
   }
   for( size_t k=0 ; k<(n/2)*(m+n) ; k++ )
   {
      size_t i1 = PEL::rand()%n ;
      size_t j1 = PEL::rand()%m ;

      double v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1.0/v ;
      }
      mat2->add_to_item( i1, j1, v ) ;
      m2[i1][j1] += v ;
   }

   double alpha0 = 3.5 ;
   mat2->synchronize() ;
   
   mat2->add_Mat( mat, alpha0 ) ;
   ok = true ;
   for( size_t i=0 ; ok && i<n ; i++ )
   {
      for( size_t j=0 ; ok && j<m ; j++ )
      {
         ok &= PEL::equal( mat2->item(i,j), m2[i][j]+alpha0*mat->item(i,j)  ) ;
      }
   }
   notify_one_test_result( "add_Mat", ok ) ;

   LA_DenseMatrix* mat3 = create_matrix( mat, m, n ) ;
   double m3[m][n] ;
   for( size_t i=0 ; i<m ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         m3[i][j] = 0.0 ;
      }
   }
   for( size_t k=0 ; k<(n/2)*(m+n) ; k++ )
   {
      size_t i1 = PEL::rand()%m ;
      size_t j1 = PEL::rand()%n ;

      double v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1.0/v ;
      }
      mat3->add_to_item( i1, j1, v ) ;
      m3[i1][j1] += v ;
   }
   
   double alpha1 = -2. ;
   mat3->synchronize() ;
   mat3->add_tMat( mat, alpha1 ) ;
   ok = true ;
   for( size_t i=0 ; ok && i<m ; i++ )
   {
      for( size_t j=0 ; ok && j<n ; j++ )
      {
         ok &= PEL::equal( mat3->item(i,j), m3[i][j]+alpha1*mat->item(j,i)  ) ;
      }
   }
   notify_one_test_result( "add_tMat", ok ) ;
   
   mat->destroy() ; mat = 0 ; x = 0 ; y = 0 ; mat2 = 0 ; mat3 = 0 ;
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doBlas3FullTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doBlas3FullTests" ) ;
   
   const size_t n = 10 ;
   const size_t m = 20 ;
   double alpha = 0.5 ;   

   LA_DenseMatrix * matA = LA_DenseMatrix::create( this, n, m ) ;
   LA_DenseMatrix * matB= LA_DenseMatrix::create( this, m, n ) ;
   LA_DenseMatrix * matC= LA_DenseMatrix::create( this, m, n ) ;

   double mA[n][m] ;
   double mB[m][n] ;
   double mC[m][n] ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         mA[i][j] = 0.0 ;
         mB[j][i] = 0.0 ;
         mC[j][i] = 0.0 ;
      }
   }
   
   for( size_t k=0 ; k<(n/2)*(m+n) ; k++ )
   {
      size_t i1 = PEL::rand()%n ;
      size_t j1 = PEL::rand()%m ;

      double v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1.0/v ;
      }

      matA->add_to_item( i1, j1, v ) ;

      mA[i1][j1] += v ;
      
      size_t i2 = PEL::rand()%n ;
      size_t j2 = PEL::rand()%m ;

      v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1.0/v ;
      }

      matB->add_to_item( j2, i2, v ) ;

      mB[j2][i2] += v ;
      
      size_t i3 = PEL::rand()%n ;
      size_t j3 = PEL::rand()%m ;

      v = PEL::rand() ;
      if( v > 0 )
      {
         v = 1.0/v ;
      }

      matC->add_to_item( j3, i3, v ) ;

      mC[j3][i3] += v ;
   }
   
   LA_DenseMatrix * mat = LA_DenseMatrix::create( this, n, n ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         mat->set_item( i, j, (i+j)*1.0/(n) ) ;
      }
   }

   mat->synchronize() ;
   matA->synchronize() ;
   matB->synchronize() ;
   
   mat->add_Mat_Mat( matA, matB, alpha ) ;
   
   bool ok = true ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         double r = (i+j)*1.0/(n) ;
         for( size_t k=0 ; k<m ; k++ )
         {
            r += alpha * mA[i][k] * mB[k][j] ;
         }
         bool okt = PEL::equal( r, mat->item(i, j ) ) ;
         if(  !okt  )
         {
            out() << r << " ?=3=? " << mat->item(i, j ) << " : " << okt << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "add_Mat_Mat", ok ) ;
   
   mat->nullify() ;
   matC->synchronize() ;
   matB->synchronize() ;
   mat->add_tMat_Mat( matB, matC, alpha ) ;
   ok = true ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         double r = 0. ;
         for( size_t k=0 ; k<m ; k++ )
         {
            r += alpha * mB[k][i] * mC[k][j] ;
         }
         bool okt = PEL::equal( r, mat->item(i, j ) ) ;
         if( !okt )
         {
            out() << r << " ?=4=? " << mat->item(i, j ) << " : " << okt << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "add_tMat_Mat", ok ) ;
   
   mat->nullify() ;
   matA->synchronize() ;
   matB->synchronize() ;
   mat->add_tMat_tMat( matB, matA, alpha ) ;
   
   ok = true ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         double r = 0. ;
         for( size_t k=0 ; k<m ; k++ )
         {
            r += alpha * mB[k][i] * mA[j][k] ;
         }
         bool okt = PEL::equal( r, mat->item(i, j ) ) ;
         if( !okt )
         {
            out() << r << " ?=4=? " << mat->item(i, j ) << " : " << okt << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "add_tMat_tMat", ok ) ;
   
   mat->re_initialize( m, m ) ;
   for( size_t i=0 ; i<m ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         mat->set_item( i, j, (i+j)*1.0/(n) ) ;
      }
   }
   
   mat->synchronize() ;
   mat->add_Mat_tMat( matC, matB, alpha ) ;
   
   ok = true ;
   for( size_t i=0 ; i<m ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         double r = (i+j)*1.0/(n) ;
         for( size_t k=0 ; k<n ; k++ )
         {
            r += alpha * mC[i][k] * mB[j][k] ;
         }
         bool okt = PEL::equal( r, mat->item(i, j ) ) ;
         if( !okt )
         {
            out() << r << " ?=4=? " << mat->item(i, j ) << " : " << okt << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "add_Mat_tMat", ok ) ;
   
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doDeterminantFullTests( std::string const& fileNameDet,
                                              std::string const& fileNameEig )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doDeterminantFullTests" ) ;
   
   LA_DenseMatrix* det= LA_DenseMatrix::create( this, 0, 0 ) ;
   LA_DenseMatrix* detm1= LA_DenseMatrix::create( this, 0, 0 ) ;
   LA_DenseMatrix* eig= LA_DenseMatrix::create( this, 0, 0 ) ;
      
   det->readMM( fileNameDet ) ;
      
   eig->readMM( fileNameEig ) ;
      
   double d = det->determinant() ;
   double d0 = 1.0 ;
   for( size_t i=0 ; i<eig->nb_cols() ; i++ )
   {
      d0 *= eig->item( 0, i ) ;
      
   }
   bool ok = PEL::toler( PEL::relative( d, d0 ), 1e-10 ) ;
   if( !ok )
   {
      out().precision( 15 ) ;
      out() << "Calculated determinant : "
            << d
            << "Given determinant : "
            << d0 << std::endl ;
   }
   string resu = "getDeterminant " ;
   notify_one_test_result( resu, ok ) ;

   detm1->readMM( fileNameDet ) ;
   double d2 ;
   detm1->invert(d2) ;
   ok = PEL::toler( PEL::relative( d0, d2 ), 1e-10 ) ;
   notify_one_test_result( "determinant (invert) ", ok ) ;
      
   if( !ok )
   {
      out() << "Determinant : invert->" << d2
            << " eig->" << d0 << std::endl ;
   }
   LA_DenseMatrix * Id = LA_DenseMatrix::create( this,
                                                 detm1->nb_rows(), 
                                                 detm1->nb_rows() ) ;
   Id->synchronize() ;
   Id->add_Mat_Mat( det, detm1 ) ;
   ok = true ;
   for( size_t i= 0 ; i<detm1->nb_rows() ; i++ )
   {
      for( size_t j=0 ; j<detm1->nb_rows() ; j++ )
      {
         d = 0 ;
         if( i==j )
         {
            d=1 ;
         }
         bool okt = PEL::equal( Id->item( i, j ), d ) ;
         if( !okt )
         {
            out() << "Not identity M(" <<i<<","<<j<<") = " << Id->item( i, j ) << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "invert ", ok ) ;
   
}

//-------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: doStatisticalTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DenseMatrix_TEST:: doStatisticalTests" ) ;
   
   const size_t n = 10 ;
   const size_t m = 15 ;

   LA_DenseMatrix  * mat = LA_DenseMatrix::create( this, n, m ) ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      for( size_t j = 0 ; j<m ; j++ )
      {
         mat->set_item( i, j, i+j ) ;
      }
   }
   mat->synchronize() ;
   LA_SeqVector * means = mat->col_mean( this ) ;
   bool ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      bool okt = PEL::equal( means->item(j), j+(n-1)/2.0 ) ;
      ok = ok && okt ;
   }
   notify_one_test_result( "compute_mean", ok ) ;
   
   double v2 = 0 ;
   for( size_t k=0 ; k<n ; k++ )
   {
      v2 += ( k-(n-1)/2.0 )*( k-(n-1)/2.0 ) ;
   }
   v2=v2/n ;
   double v = PEL::sqrt( v2 ) ;
   
   LA_SeqVector * dev = mat->col_standard_deviation( this ) ;
   ok = true ;
   for( size_t j = 0 ; j<m ; j++ )
   {
      bool okt = PEL::equal( dev->item(j), v ) ;
      if(  !okt )
      {
         out() << "compute_standard_deviation comp. : " << dev->item(j)
               << "  theo. : " << v << std::endl ;
      }
      ok = ok && okt ;
   }
   notify_one_test_result( "getStandardDeviationOfColumns", ok ) ;
   
   LA_SymmetricMatrix * variance = mat->col_variance( this ) ;
   for( size_t i=0 ; i<m ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         bool okt = PEL::equal( variance->item( i, j ), v2 ) ;
         if(  !okt )
         {
            out() << "variance(" << i << "," << j << ") :  comp. : " << variance->item( i, j )
                  << "  theo. : " << v2 << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "getVarianceMatrix", ok ) ;
   
   LA_SymmetricMatrix * correlation = mat->col_correlation( this ) ;
   
   for( size_t i=0 ; i<m ; i++ )
   {
      for( size_t j=0 ; j<m ; j++ )
      {
         bool okt = PEL::equal( correlation->item( i, j ), 1.0 ) ;
         if(  !okt )
         {
            out() << "correlation(" << i << "," << j << ") :  comp. : " << correlation->item( i, j )
                  << "  theo. : " << 1.0 << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   notify_one_test_result( "getCorrelationMatrix", ok ) ;
   
}

//----------------------------------------------------------------------------
void
LA_DenseMatrix_TEST:: display_error( double xx_1, double xx_2 ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = out().flags() ;
   out().setf( ios_base::uppercase | ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;

   out() << std::setw( 20 ) << xx_1
         << std::setw( 20 ) << xx_2
         << endl ;

   out().flags( original_flags ) ;
}
