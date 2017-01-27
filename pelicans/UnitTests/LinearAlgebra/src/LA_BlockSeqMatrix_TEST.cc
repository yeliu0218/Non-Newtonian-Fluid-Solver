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

#include <LA_BlockSeqMatrix_TEST.hh>

#include <LA_BlockSeqMatrix.hh>
#include <LA_PelMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_Vector.hh>

#include <PEL_Root.hh>
#include <PEL_System.hh>
#include <PEL_Timer.hh>

#include <size_t_vector.hh>

#include <iostream>

using std::cout ;
using std::endl ;

LA_BlockSeqMatrix_TEST const*
LA_BlockSeqMatrix_TEST::PROTOTYPE = new LA_BlockSeqMatrix_TEST() ;

//-------------------------------------------------------------------------
LA_BlockSeqMatrix_TEST:: LA_BlockSeqMatrix_TEST( void )
//-------------------------------------------------------------------------
   : LA_SeqMatrix_TEST( "LA_BlockSeqMatrix", "LA_BlockSeqMatrix_TEST" )
{
}

//-------------------------------------------------------------------------
LA_BlockSeqMatrix_TEST:: ~LA_BlockSeqMatrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
LA_BlockSeqMatrix_TEST:: test( PEL_ModuleExplorer const* exp,
                               LA_Matrix const* matrix,
                               std::string const& test_matrix )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix_TEST:: test" ) ;
   
   LA_SeqMatrix_TEST::test( exp, matrix, test_matrix ) ;

   doBlocTests() ;

   PEL_Timer * timer = PEL_Timer::create( 0 ) ;
   const int m = 3000 ;

   LA_SeqMatrix* sparse = create_matrix( timer, m, m ) ;
   LA_SeqMatrix* sparse2 = LA_PelMatrix::create( timer, m, m ) ;
   LA_SeqMatrix* sparsePEL = LA_PelMatrix::create( timer, m, m ) ;

   for( int i=0 ; i<100*m ; i++ )
   {
      int u = PEL::rand()%m ;
      int v = PEL::rand()%m ;
      double val = 1.0/( PEL::rand()+1 ) ;
      sparse->add_to_item( u, v, val ) ;
      sparse2->add_to_item( u, v, val ) ;
   }

   LA_SeqVector* vec = LA_SeqVector::create( timer, m ) ;
   LA_SeqVector* y = LA_SeqVector::create( timer, m ) ;
   for( int i=0 ; i<m ; i++ )
   {
      vec->set_item( i, PEL::rand() ) ;
   }
   sparse->synchronize() ;
   sparse2->synchronize() ;
   vec->synchronize() ;
   
   timer->start() ;
   sparsePEL->set( sparse2 ) ;
   timer->stop() ;
   out() << "Time to set SparsePEL matrix : " << timer->time() << endl ;
   timer->reset() ;
   sparsePEL->nullify() ;

   timer->start() ;
   sparsePEL->set( sparse ) ;
   timer->stop() ;
   out() << "Time to set SparseBLOC matrix : " << timer->time() << endl ;

   timer->reset() ;
   timer->start() ;
   sparse2->multiply_vec_then_add( vec, y ) ;
   timer->stop() ;
   out() << "Time to multiply_vec_then_add with PEL : " << timer->time() << endl ;

   timer->reset() ;
   timer->start() ;
   sparse->multiply_vec_then_add( vec, y ) ;
   timer->stop() ;
   out() << "Time to multiply_vec_then_add with Bloc : " << timer->time() << endl ;

   timer->destroy() ;
}

//-------------------------------------------------------------------------
LA_BlockSeqMatrix*
LA_BlockSeqMatrix_TEST:: create_matrix(
                            PEL_Object* a_owner, size_t n, size_t m ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix_TEST:: create_matrix" ) ;

   size_t_vector row( 2 ) ;
   row( 0 ) = n/2 ;
   row( 1 ) = n-n/2 ;

   size_t_vector col( 2 ) ;
   col( 0 ) = m/2 ;
   col( 1 ) = m-m/2 ;

   LA_PelMatrix* proto = LA_PelMatrix::create( 0, 0, 0 ) ;
   LA_BlockSeqMatrix* result =
                 LA_BlockSeqMatrix::create( a_owner, row, col, proto ) ;
   proto->destroy() ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
LA_BlockSeqMatrix_TEST:: doBlocTests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockSeqMatrix_TEST:: doBlocTests" ) ;

   size_t_vector rowPart( 3 ) ;
   size_t_vector colPart( 2 ) ;
   rowPart(0) = 2 ;
   rowPart(1) = 4 ;
   rowPart(2) = 8 ;
   colPart(0) = 5 ;
   colPart(1) = 4 ;

   LA_PelMatrix* proto = LA_PelMatrix::create( 0, 0, 0 ) ;
   LA_BlockSeqMatrix* mat = LA_BlockSeqMatrix::create(
                                             0, rowPart, colPart, proto ) ;
   proto->destroy() ; proto = 0 ;

   notify_one_test_result( "row_partition",
                           mat->row_partition()==rowPart ) ;
   notify_one_test_result( "col_partition",
                           mat->col_partition()==colPart ) ;

   mat->set_item( 7, 7, 1.0 ) ;

   mat->synchronize() ;
   mat->writeMM( "blocMatrix.mtx" ) ;
   size_t_vector rowPart2( 2 ) ;
   size_t_vector colPart2( 2 ) ;
   rowPart2(0) = 1 ;
   rowPart2(1) = 3 ;
   colPart2(0) = 3 ;
   colPart2(1) = 1 ;
   mat->set_partitioning_then_nullify( rowPart2, colPart2 ) ;

   notify_one_test_result( "reInit",
                               mat->row_partition()==rowPart2 &&
                               mat->col_partition()==colPart2 &&
                               mat->nb_stored_items()== 0 ) ;

   PEL_System::erase( "blocMatrix.mtx" ) ;

   mat->destroy() ; mat = 0 ;
}
