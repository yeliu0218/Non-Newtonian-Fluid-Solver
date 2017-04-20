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

#include <LA_MatrixSynchronization_TEST.hh>

#include <LA_DistImplementation.hh>
#include <LA_DenseMatrix.hh>
#include <LA_Matrix.hh>
#include <LA_PelMatrix.hh>
#include <LA_Scatter.hh>
#include <LA_SeqImplementation.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>

#include <string>
#include <iostream>
#include <iomanip>

LA_MatrixSynchronization_TEST const*
LA_MatrixSynchronization_TEST::PROTOTYPE =
                               new LA_MatrixSynchronization_TEST()  ;

//-------------------------------------------------------------------------
LA_MatrixSynchronization_TEST:: LA_MatrixSynchronization_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "LA_Matrix", "LA_MatrixSynchronization_TEST" )
   , MAT_PROTO( 0 )
   , FIRST_ROW( PEL::bad_index() )
   , LAST_ROW( PEL::bad_index() )
   , COM( 0 )
   , DIM( 3 )
   , DmatLHS( 0 )
   , matLHS( 0 )
   , vecSOL( 0 )
   , sumSOL( PEL::bad_double() )
   , matONES( 0 )
   , vecRHS( 0 )
   , sumRHS( PEL::bad_double() )
{
}

//-------------------------------------------------------------------------
LA_MatrixSynchronization_TEST:: ~LA_MatrixSynchronization_TEST( void )
//-------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: process_one_test(
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_MatrixSynchronization_TEST:: process_one_test" ) ;

   COM = PEL_Exec::communicator() ;

   DIM = (size_t) exp->int_data( "matrices_dimensions" ) ;

   EXP_PROTO = exp->create_subexplorer( this, "PROTOTYPE" ) ;
   MAT_PROTO = LA_Matrix::make( EXP_PROTO, EXP_PROTO ) ;
   MAT_PROTO->re_initialize( DIM, DIM ) ;
   NAME = MAT_PROTO->name() ;
   FIRST_ROW = MAT_PROTO->row_distribution()->first_local_index() ;
   LAST_ROW = MAT_PROTO->row_distribution()->local_index_limit() ;
   
   reference_matrices_initialization() ;

   //-----------------------------------------------------------
   // Tests of matrices creation

   //-- Does a matrix contain only zeros after creation using "make" ?
   test_zero_after_make() ;

   //-- Does a matrix contain only zeros after creation using "create_matrix" ?
   test_zero_after_create_matrix() ;

   //-- Does a matrix created with "make "contain only zeros
   //   after a call of "nullify" method ?
   test_zero_after_nullify_1() ;

   //-- Does a matrix created with "make "contain only zeros
   //   after a call of "nullify" method ?
   test_zero_after_nullify_2() ;

   //-- Does creation on only one process is allowed ?
   //   only for "!is_distributed()" matrix
   test_creation_on_only_one_process() ;

   //-----------------------------------------------------------
   /* The following tests verify the validity of synchronization
    * after call to different methods of LA_Matrix class.
    * Verification is perfomed using product matrix-vector
    * owing to the following equalities :
    *     " matLHS * vecSOL = vecRHS " and
    *     " [matONES *matLHS * vecSOL](i) = sumRHS, for all i=0..DIM-1 "
    */

   //-- set all matrix items on all process
   test_set_item_1() ;

   if( MAT_PROTO->implementation() == LA_DistImplementation::object() ||
       MAT_PROTO->implementation() == LA_SeqImplementation::object() )
   {
      //-- set only local items using `start/stop_local_modifs()'
      //   to avoid useless synchronization
      test_set_item_2() ;

      //-- mix add and set_item
      test_mix_add_and_set() ;
   }

   //-- use add_Mat without synchronization after
   test_add_Mat_1() ;

   //-- use add_Mat with an unsynchronized matrix
   test_add_Mat_2() ;

   //-- use add_tMat without synchronization after
   test_add_tMat_1() ;

   //-- use add_tMat with an unsynchronized matrix
   test_add_tMat_2() ;

   //-- use add_Mat_Mat without synchronization after
   test_add_Mat_Mat_1() ;

   //-- use add_Mat_Mat with an unsynchronized matrix
   test_add_Mat_Mat_2() ;
}

//----------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: reference_matrices_initialization( void )
//----------------------------------------------------------------------
{
   //-- ONES initialization
   /*
    *           | 1 1 1 |
    * matONES = | 1 1 1 |  matONES_{ij} = 1
    *           | 1 1 1 |
    */
   matONES =  MAT_PROTO->create_matrix( this ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         matONES->set_item( i, j, 1. ) ;
      }
   }
   matONES->synchronize() ;

   //-- Dense LHS initialization
   /*
    *           | 1 2 3 |
    * DmatLHS = | 4 5 6 |  DmatLHS_{ij} = i*DIM + j+1
    *           | 7 8 9 |
    */
   DmatLHS = LA_DenseMatrix::create( this, DIM, DIM ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         DmatLHS->set_item( i, j, (double)(i*DIM + j+1) ) ;
      }
   }
   DmatLHS->synchronize() ;

   //-- LHS initialization (same as DmatLHS )
   /*
    *                    | 1 2 3 |
    * matLHS = DmatLHS = | 4 5 6 |  matLHS_{ij} = i*DIM + j+1
    *                    | 7 8 9 |
    */
   matLHS = MAT_PROTO->create_matrix( this ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM; ++j )
      {
         matLHS->set_item( i, j, DmatLHS->item( i, j ) ) ;
      }
   }
   matLHS->synchronize() ;

   //-- SOL initialization
   /*
    *          | 1 |
    * vecSOL = | 2 |     vecSOL_{i} = i+1 ;
    *          | 3 |
    */
   vecSOL = MAT_PROTO->create_vector( this ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      vecSOL->set_item( i, (double)i + 1. )  ;
   }
   vecSOL->synchronize() ;

   //-- sum of vecSOL items
   /*
    * sumSOL = sum_{i=0}^{DIM-1} vecSOL(i)
    *        = DIM*(DIM+1)/2.
    *        = matONES * vecSOL
    */
   sumSOL =  (double)(DIM*(DIM+1))/2. ;

   //-- RHS initialization
   /*
    *          | 14 |
    * vecRHS = | 32 |  vecRHS_{i} = i*DIM*DIM*(DIM+1)/2
    *          | 50 |               + DIM*(DIM+1)*(2*DIM+1)/6
    */
   vecRHS = LA_SeqVector::create( this, DIM ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      double yi = (double)( i*DIM*DIM*(DIM+1))/2.
                + (double)(DIM*(DIM+1)*(2*DIM+1))/6. ;
      vecRHS->set_item( i , yi ) ;
   }
   vecRHS->synchronize() ;

   //-- sum of vecRHS items
   /*
    * sumRHS = sum_{i=0}^{DIM-1} vecRHS(i)
    *        = DIM*DIM*DIM*(DIM+1.)*(DIM-1.)/4.
    *        + DIM*DIM*(DIM+1.)*(2*DIM+1.)/6.
    *        = matONES * vecRHS
    */
   sumRHS = (double)(DIM*DIM*DIM*(DIM+1.)*(DIM-1.))/4.
          + (double)(DIM*DIM*(DIM+1.)*(2*DIM+1.))/6. ;

   //-- vector used product matrix-vector
   vecX = MAT_PROTO->create_vector( this ) ;
   scatterX = create_identity_scatter( vecX, vecX ) ;
   vecSX = LA_SeqVector::create( vecX, DIM ) ;
 }

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_zero_after_make( void )
//-------------------------------------------------------------------------
{
   matA = LA_Matrix::make( 0, EXP_PROTO ) ;
   matA->re_initialize( DIM, DIM ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "                    make",
                           double_equality( vecX->max_norm(), 0. ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_zero_after_create_matrix( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "            create_matrix",
                           double_equality( vecX->max_norm(), 0. ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_zero_after_nullify_1( void )
//-------------------------------------------------------------------------
{
   matA = LA_Matrix::make( 0, EXP_PROTO ) ;
   matA->re_initialize( DIM, DIM ) ;
   matA->add_Mat( matONES ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "            add_Mat-make",
                           double_equality( vecX->max_norm(), sumSOL ) ) ;
   matA->nullify() ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "            nullify-make",
                           double_equality( vecX->max_norm(), 0. ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_zero_after_nullify_2( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->add_Mat( matONES ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "     add_Mat-create_matrix",
                           double_equality( vecX->max_norm(), sumSOL ) ) ;
   matA->nullify() ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   notify_one_test_result( NAME + "      nullify-create_matrix",
                           double_equality( vecX->max_norm(), 0. ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_creation_on_only_one_process( void )
//-------------------------------------------------------------------------
{
   bool ok = true ;
   if( MAT_PROTO->distribution_strategy() == LA::NoDistribution && 
       COM->rank() == 0 )
   {
      matA = MAT_PROTO->create_matrix( 0 ) ;
      matA->add_Mat( matONES ) ;
      matA->synchronize() ;
      matA->multiply_vec_then_add( vecSOL, vecX ) ;
      ok = double_equality( vecX->max_norm(), sumSOL ) ;
      matA->destroy() ; matA = 0 ;
   }
   notify_one_test_result( NAME + "         creation only on proc 0", ok ) ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_set_item_1( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM; ++j )
      {
         matA->set_item( i, j, DmatLHS->item( i, j ) ) ;
      }
   }
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   notify_one_test_result( NAME + "               set (global)",
                           vector_equality( vecSX, vecRHS ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_set_item_2( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->start_local_modifs() ;
   for( size_t i=FIRST_ROW ; i<LAST_ROW ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         matA->set_item( i, j, DmatLHS->item( i, j ) ) ;
      }
   }
   matA->stop_local_modifs() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   notify_one_test_result( NAME + "    set (only_local_modifs)",
                           vector_equality( vecSX, vecRHS ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_mix_add_and_set( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->start_local_modifs() ;
   for( size_t i=FIRST_ROW ; i<LAST_ROW ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         matA->set_item( i, j, DmatLHS->item( i, j ) - 1. ) ;
         matA->add_to_item( i, j, 1. ) ;
      }
   }
   matA->stop_local_modifs() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   notify_one_test_result( NAME + "            mix add and set",
                           vector_equality( vecSX, vecRHS )  ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_Mat_1( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->add_Mat( matLHS ) ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   notify_one_test_result( NAME + "                 add_Mat(1)",
                           vector_equality( vecSX, vecRHS ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_Mat_2( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         matA->add_to_item( i, j, DmatLHS->item( i, j ) ) ;
      }
   }
   matA->add_Mat( matLHS, 1. ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   double alpha = ( matA->distribution_strategy() != LA::NoDistribution ) ? 
                     COM->nb_ranks() + 1. : 2. ;
   vecRHS->scale( alpha ) ;
   notify_one_test_result( NAME + "                 add_Mat(2)",
                           vector_equality( vecSX, vecRHS ) ) ;
   vecRHS->scale( 1. / alpha ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_tMat_1( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->add_tMat( matLHS ) ;
   matA->tr_multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   notify_one_test_result( NAME + "                add_tMat(1)",
                           vector_equality( vecSX, vecRHS ) ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_tMat_2( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      for( size_t j=0 ; j<DIM ; ++j )
      {
         matA->add_to_item( i, j, DmatLHS->item( j, i ) ) ;
      }
   }
   matA->add_tMat( matLHS, 1. ) ;
   matA->synchronize() ;
   matA->tr_multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   double alpha = ( matA->distribution_strategy() != LA::NoDistribution ) ? 
                       COM->nb_ranks() + 1. : 2. ;
   vecRHS->scale( alpha ) ;
   notify_one_test_result( NAME + "                add_tMat(2)",
                           vector_equality( vecSX, vecRHS ) ) ;
   vecRHS->scale( 1. / alpha ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_Mat_Mat_1( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   matA->synchronize() ;
   matA->add_Mat_Mat( matONES, matLHS ) ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   bool ok = true ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      ok = double_equality( vecSX->item( i ), sumRHS ) && ok;
   }
   notify_one_test_result( NAME + "             add_Mat_Mat(1)", ok ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
void
LA_MatrixSynchronization_TEST:: test_add_Mat_Mat_2( void )
//-------------------------------------------------------------------------
{
   matA = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t i=FIRST_ROW ; i<LAST_ROW ; ++i )
   {
      for( size_t j=0 ; j<DIM; ++j )
      {
         matA->add_to_item( i, j, DmatLHS->item( i, j ) ) ;
      }
   }
   matA->add_Mat_Mat( matONES, matLHS ) ;
   matA->synchronize() ;
   matA->multiply_vec_then_add( vecSOL, vecX ) ;
   scatterX->get( vecX, vecSX ) ;
   bool ok = true ;
   for( size_t i=0 ; i<DIM ; ++i )
   {
      ok = double_equality( vecSX->item( i ), sumRHS+vecRHS->item( i ) ) && ok ;
   }
   notify_one_test_result( NAME + "             add_Mat_Mat(2)", ok ) ;
   matA->destroy() ; matA = 0 ;
}

//-------------------------------------------------------------------------
bool
LA_MatrixSynchronization_TEST:: double_equality( double x,
                                                 double y ) const
//-------------------------------------------------------------------------
{
   bool result = PEL::double_equality( x, y, 1.E-8, 1.E-10 ) ;
   if( !result )
   {
      PEL::out() << std::setprecision( 10 )
                 << std::setiosflags( std::ios::scientific )
                 << "failed : " << x << " <-> " << y << std::endl ;
   }
   return( result ) ;
}


//-------------------------------------------------------------------------
bool
LA_MatrixSynchronization_TEST:: vector_equality(
                                           LA_SeqVector const* vec1,
                                           LA_SeqVector const* vec2 ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_MatrixSynchonization_TEST:: vector_equality" ) ;
   PEL_CHECK_PRE( vec1 != 0 ) ;
   PEL_CHECK_PRE( vec2 != 0 ) ;

   bool ok = vec1->nb_rows() == vec2->nb_rows() ;
   if( ok )
   {
      for( size_t i=0 ; i<vec1->nb_rows() ; i++ )
      {
         bool okt = double_equality( vec1->item(i), vec2->item(i) ) ;
         if( !ok )
         {
            out() << " vector_equality diff " << i << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   return ok ;
}

//--------------------------------------- ----------------------------------
LA_Scatter*
LA_MatrixSynchronization_TEST:: create_identity_scatter(
                        PEL_Object* a_owner, LA_Vector const* vec ) const
//-------------------------------------------------------------------------
{
   size_t n = vec->nb_rows() ;

   size_t_vector from( n ) ;
   size_t_vector to( n ) ;

   for( size_t i=0 ; i<n ; i++)
   {
      // Identity
      from( i ) = i ;
      to( i ) = i ;
   }
   LA_Scatter* result = vec->create_scatter( a_owner, from, to ) ;
   return( result ) ;
}
