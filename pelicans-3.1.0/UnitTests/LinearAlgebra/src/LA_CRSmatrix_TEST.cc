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

#include <LA_CRSmatrix_TEST.hh>

#include <LA_CRSmatrix.hh>
#include <LA_PelMatrix.hh>
#include <LA_MatrixIterator.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>
#include <doubleArray2D.hh>

#include <iostream>

LA_CRSmatrix_TEST const* 
LA_CRSmatrix_TEST::PROTOTYPE = new LA_CRSmatrix_TEST()  ;

//-------------------------------------------------------------------------
LA_CRSmatrix_TEST:: LA_CRSmatrix_TEST( void )
//-------------------------------------------------------------------------
   : LA_SeqMatrix_TEST( "LA_CRSmatrix", "LA_CRSmatrix_TEST" )
{
}

//-------------------------------------------------------------------------
LA_CRSmatrix_TEST:: ~LA_CRSmatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_CRSmatrix_TEST:: doBlas3Tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix_TEST:: doBlas3Tests" ) ;

   LA_SeqMatrix_TEST::doBlas3Tests() ;
}

//-------------------------------------------------------------------------
void
LA_CRSmatrix_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix_TEST:: process_all_tests" ) ;
   
   LA_SeqMatrix_TEST::process_all_tests() ;
   
   do_distributed_tests() ;
}

//-------------------------------------------------------------------------
void
LA_CRSmatrix_TEST:: do_distributed_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CRSmatrix_TEST:: do_distributed_tests" ) ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   size_t const n = 10 ;
   size_t const m = 20 ;
   LA_CRSmatrix* mat = 0 ;
      
   if( com->rank()==0 )
   {
      LA_PelMatrix* pmat = LA_PelMatrix::create( this, n, m ) ;
      for( size_t i=0 ; i<n ; i++ )
         for( size_t j=0 ; j<m ; j++ )
            pmat->set_item( i, j, i+j ) ;
      pmat->synchronize() ;
      
      mat = LA_CRSmatrix::create( this, pmat ) ;
      
      for( size_t dest=1 ; dest<com->nb_ranks() ; dest++ )
         mat->send( com, dest ) ;
   }
   else
      mat = LA_CRSmatrix::receive( this, com, 0 ) ;
   
   bool ok = true ;
      
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<m ; j++ )
         ok &= mat->item(i,j) == (i+j) ;
   notify_one_test_result( "send/receive", ok ) ;
    
}
