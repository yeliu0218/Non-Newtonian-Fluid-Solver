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

#include <GE_Matrix_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_List.hh> 
#include <PEL_Vector.hh> 
#include <doubleVector.hh> 

#include <GE_Matrix.hh>
#include <GE_Vector.hh>

#include <PEL_ModuleExplorer.hh>
 
#include <doubleArray2D.hh>
#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
GE_Matrix_TEST*
GE_Matrix_TEST::unique_instance = new GE_Matrix_TEST() ;
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
GE_Matrix_TEST:: GE_Matrix_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Matrix", "GE_Matrix_TEST" )
{
}


//----------------------------------------------------------------------------
GE_Matrix_TEST:: ~GE_Matrix_TEST( void )
//----------------------------------------------------------------------------
{
}


//----------------------------------------------------------------------------
void
GE_Matrix_TEST:: process_one_test( PEL_ModuleExplorer const* texp ) 
//----------------------------------------------------------------------------
{
   doubleArray2D const& matrix = texp->doubleArray2D_data( "matrix" ) ;
   doubleArray2D const& B = texp->doubleArray2D_data( "B" ) ;
   
   size_t n = matrix.index_bound(0) ;
   size_t m = matrix.index_bound(1) ;
   doubleArray2D X(B.index_bound(0),B.index_bound(1)) ;
   
   GE_Matrix* gemat = GE_Matrix::create( this, n, m ) ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<n ; j++ )
         gemat->set_item( i, j, matrix(i,j) ) ;
   
   notify_one_test_result( "nb_rows",
                        gemat->nb_rows()==n ) ;
   notify_one_test_result( "nb_cols",
                        gemat->nb_cols()==m ) ;
   notify_one_test_result( "determinant_is_computed1",
                        !gemat->determinant_is_computed() ) ;
   if( n==m )
   {
      gemat->compute_determinant() ;
      notify_one_test_result( "determinant_is_computed2",
                           gemat->determinant_is_computed() ) ;
      notify_one_test_result( "determinant",
                           PEL::equal( gemat->determinant(),
                                       texp->double_data( "determinant" ) ) ) ;
      gemat->invert( B, X ) ;
      bool ok = true ;
      for( size_t k=0 ; k<B.index_bound(0) ; k++ )
         for( size_t i=0 ; i<n ; i++ )
         {
            double r = -B(k,i) ;
            for( size_t j=0 ; j<m ; j++ )
            {
               r+=gemat->item(i,j)*X(k,j) ;
            }
            ok = ok && PEL::equal( r, 0.0 ) ;
         }
      notify_one_test_result( "invert",ok ) ;
      
      
   }
   gemat->nullify() ;
   bool ok = true ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
      {
         gemat->set_item( i, j, matrix(i,j) ) ;
         gemat->add_to_item( i, j, -matrix(i,j) ) ;
         ok = ok && PEL::equal( gemat->item(i,j), 0.0 ) ;
      }
   notify_one_test_result( "add_item",ok ) ;
  
   
}
