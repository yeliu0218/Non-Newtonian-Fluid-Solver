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

#include <LA_DenseMatrix.hh>

#include <LA_DenseMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_SymmetricMatrix.hh>
#include <LA_SymmetricMatrix_TEST.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>

#include <iostream>

//-------------------------------------------------------------------------
LA_SymmetricMatrix_TEST*
LA_SymmetricMatrix_TEST::unique_instance = new LA_SymmetricMatrix_TEST() ;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
LA_SymmetricMatrix_TEST:: LA_SymmetricMatrix_TEST( void )
//-------------------------------------------------------------------------
   : LA_SeqMatrix_TEST( "LA_SymmetricMatrix", "LA_SymmetricMatrix_TEST" )
{
}



//-------------------------------------------------------------------------
LA_SymmetricMatrix_TEST:: ~LA_SymmetricMatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_SymmetricMatrix_TEST:: test( PEL_ModuleExplorer const* exp,
                                LA_Matrix const* matrix,
                                std::string const& test_matrix )
//-------------------------------------------------------------------------
{
   LA_SeqMatrix_TEST::test(exp,matrix,test_matrix);
   
   doSymmetricTest( exp->string_data( "det_test" ),
                    exp->string_data( "eig_test" ) ) ;
}

//-------------------------------------------------------------------------
void
LA_SymmetricMatrix_TEST:: doSymmetricTest(
                        std::string const& mat, std::string const& eigmat )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SymmetricMatrix_TEST:: doSymmetricTest" ) ;
   
   LA_DenseMatrix * Af = LA_DenseMatrix::create( this, 0, 0 ) ;
   LA_DenseMatrix * eig = LA_DenseMatrix::create( this, 0, 0 ) ;
   
   Af->readMM( mat ) ;
   eig->readMM( eigmat ) ;
   
   size_t n = Af->nb_rows() ;
   
   LA_SymmetricMatrix * As = LA_SymmetricMatrix::create( this, n ) ;

   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=i ; j<n ; j++ )
      {
         As->set_item( i, j, Af->item( i, j ) ) ;

         PEL_ASSERT( PEL::equal( Af->item( i, j ), Af->item( j, i ) ) ) ;
      }
   }
   
   double d0 = 1.0 ;
   for( size_t i=0 ; i<n ; i++ )
   {
      d0 *= eig->item( 0, i ) ;
   }
   PEL_ASSERT( PEL::abs(d0) > 1e-10 ) ;

   As->synchronize() ;
   double d1 = PEL::bad_double() ;
   As->invert( d1 ) ;
   out() << " d1 " << d1 << " d0 " << d0 << std::endl ;
   
   notify_one_test_result( "determinant",
                           PEL::toler( PEL::relative( d0, d1 ), 1e-10 ) ) ;
   
   LA_DenseMatrix * Afm1 = LA_DenseMatrix::create( this, n, n ) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         Afm1->set_item( i, j, As->item( i, j ) ) ;
      }
   }
   LA_DenseMatrix * Id = LA_DenseMatrix::create( this, n, n ) ;
   Id->synchronize() ;
   Afm1->synchronize() ;
   Id->add_Mat_Mat( Af, Afm1 ) ;
   bool ok = true ;
   for( size_t i= 0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<n ; j++ )
      {
         double d = 0 ;
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
   notify_one_test_result( "invert", ok ) ;

   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=i ; j<n ; j++ )
      {
         As->set_item( i, j, Af->item( i, j ) ) ;
      }
   }
   As->synchronize() ;
   LA_SeqVector * eigVals = LA_SeqVector::create( this, n ) ;
   LA_DenseMatrix* eigVecs = LA_DenseMatrix::create( this, n, n ) ;
   
   As->eigen_reduce( n, eigVals, eigVecs, ok ) ;
   if( ok )
   {
      LA_SeqVector * tmp = LA_SeqVector::create( this, n ) ;
      LA_SeqVector * eigVec = LA_SeqVector::create( this, n ) ;
      for( size_t i=0 ; i<n ; i++ )
      {
         if( PEL::abs( eigVals->item(i) ) > 0.0 )
         {
            eigVecs->extract_col( i, eigVec ) ;
            tmp->set( eigVec ) ;
            As->multiply_vec_then_add( eigVec, tmp, -1.0/eigVals->item(i), 1.0 ) ;
            double r = tmp->two_norm() / eigVec->two_norm() ;
            bool okt = r < 1e-10 ;
            if( ! okt )
            {
               out() << "Residu : " << r << std::endl << " with vector " << eigVec
                     << " and value " << eigVals->item(i) << std::endl ;
            }
            ok = ok && okt ;
         }
      }
   }
   notify_one_test_result( "eigen_reduce", ok ) ;
   
            
}

//-------------------------------------------------------------------------
LA_SymmetricMatrix *
LA_SymmetricMatrix_TEST:: newMat( PEL_Object* a_owner, size_t n, size_t m ) 
//-------------------------------------------------------------------------
{
   PEL_ASSERT( n == m ) ;
   return(  LA_SymmetricMatrix::create( a_owner, n ) ) ;
}


