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

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>
#include <LA_PelMatrix.hh>

#include <EXT_PETScMatrix.hh>
#include <EXT_PETScAPI.hh>
#include <EXT_PETScImplementation.hh>

#include <iomanip>
#include <iostream>
#include <fstream>

//                                      sequential  symmetric   block
EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_MPIAIJ =
   new EXT_PETScMatrix("PETSc_MPIAIJ",      false,     false,   false ) ;

EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_MPISBAIJ =
   new EXT_PETScMatrix("PETSc_MPISBAIJ",    false,      true,    true ) ;

EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_MPIBAIJ =
   new EXT_PETScMatrix("PETSc_MPIBAIJ",     false,     false,    true ) ;

EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_SeqSBAIJ =
   new EXT_PETScMatrix("PETSc_SeqSBAIJ",     true,      true,    true ) ;

EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_SeqAIJ =
   new EXT_PETScMatrix("PETSc_SeqAIJ",       true,     false,   false ) ;

EXT_PETScMatrix const*
EXT_PETScMatrix:: PROTOTYPE_AIJ =
   new EXT_PETScMatrix("PETSc_AIJ",          true,     false,   false ) ;

//----------------------------------------------------------------------
EXT_PETScMatrix:: EXT_PETScMatrix( std::string const& a_name,
                                   bool sequential,
                                   bool symmetric,
                                   bool block )
//----------------------------------------------------------------------
   : LA_Matrix( a_name )
   , EXP( 0 )
   , DESTROY_ON_EXIT( true )
   , SYMMETRIC( symmetric )
   , SEQ( sequential )
   , BLOCK( block )
   , HAS_OPT( false )
   , VERB( false )
   , MATRIX( 0 )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_DIST( 0 )
   , COL_DIST( 0 )
   , UPPER( false )
{
   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
EXT_PETScMatrix*
EXT_PETScMatrix:: create_replica( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;

   EXT_PETScMatrix* result = new EXT_PETScMatrix( a_owner, exp, this ) ;

   if( this==PROTOTYPE_AIJ )
   {
      if( PEL_Exec::communicator()->nb_ranks() > 1 )
      {
         result->SEQ = false ;
      }
   }

   if( SEQ )
   {
      PEL_ASSERT( result->distribution_strategy() == LA::NoDistribution ) ;
   }
   else
   {
      PEL_ASSERT( result->distribution_strategy() == LA::FromGlobalSize ) ;
      //la strategie LA::FromLocalSize est disponible dans PETSc
      //mais pas implemente ici
   }

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScMatrix:: EXT_PETScMatrix( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp,
                                   EXT_PETScMatrix const* other  )
//----------------------------------------------------------------------
   : LA_Matrix( a_owner, other->name(),
                other->SEQ ? LA::NoDistribution : LA::FromGlobalSize )
                // the stategy LA::FromLocalSize is available in PETSc
                // but not implemented here
   , EXP( exp->create_clone( this ) )
   , DESTROY_ON_EXIT(
      exp->has_entry( "destroy") ?  exp->bool_data( "destroy" ) : true )
   , SYMMETRIC( other->SYMMETRIC )
   , SEQ( other->SEQ )
   , BLOCK( other->BLOCK )
   , HAS_OPT( false )
   , VERB( exp->has_entry( "verbose" ) ? exp->bool_data( "verbose" ) : false )
   , MATRIX( 0 )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_DIST( PEL_DistributedPartition::create( this ) )
   , COL_DIST( PEL_DistributedPartition::create( this ) )
   , UPPER( false )
{
   PEL_LABEL( "EXT_PETScMatrix:: EXT_PETScMatrix" ) ;

   HAS_OPT = EXT_PETScAPI::parse_options( exp, VERB ) ;
   ROW_DIST->set_local_number( 0 ) ;
   COL_DIST->set_local_number( 0 ) ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   PEL_CHECK_POST( state() == LA::NotSync_undef ) ;
//   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
EXT_PETScMatrix*
EXT_PETScMatrix:: create_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: create_matrix" ) ;
   PEL_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;

   EXT_PETScMatrix* result = new EXT_PETScMatrix( a_owner, this ) ;

   PEL_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScMatrix:: EXT_PETScMatrix( PEL_Object* a_owner,
                                   EXT_PETScMatrix const* other )
//----------------------------------------------------------------------
   : LA_Matrix( a_owner, other->name(), other->distribution_strategy() )
   , EXP( other->EXP->create_clone( this ) )
   , DESTROY_ON_EXIT( other->DESTROY_ON_EXIT )
   , SYMMETRIC( other->SYMMETRIC )
   , SEQ( other->SEQ )
   , BLOCK( other->BLOCK )
   , HAS_OPT( other->HAS_OPT )
   , VERB( other->VERB )
   , MATRIX( 0 )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_DIST( PEL_DistributedPartition::create( this ) )
   , COL_DIST( PEL_DistributedPartition::create( this ) )
   , UPPER( false )
{
   PEL_LABEL( "EXT_PETScMatrix:: EXT_PETScMatrix" ) ;
   PEL_CHECK_PRE( !other->is_a_prototype() ) ;

   re_initialize( other->nb_rows(), other->nb_cols() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( state() == LA::NotSync_undef ) ;
   PEL_CHECK_POST( distribution_strategy() == other->distribution_strategy() ) ;
}

//----------------------------------------------------------------------
EXT_PETScMatrix:: ~EXT_PETScMatrix( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: ~EXT_PETScMatrix" ) ;
   destroy_matrix() ;
}

//----------------------------------------------------------------------
EXT_PETScVector*
EXT_PETScMatrix:: create_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: create_vector" ) ;
   PEL_CHECK_PRE( create_vector_PRE( a_owner ) ) ;

   EXT_PETScVector* result =
                     EXT_PETScVector::create( a_owner, SEQ, NB_ROWS ) ;

   PEL_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return result ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: re_initialize( size_t a_nb_rows, size_t a_nb_cols,
                                   size_t a_nb_local_rows,
                                   size_t a_nb_local_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: re_initialize" ) ;
   PEL_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_cols,
                                     a_nb_local_rows, a_nb_local_cols ) ) ;

   NB_ROWS = a_nb_rows ;
   NB_COLS = a_nb_cols ;
   build( intVector( 0 ) ) ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_cols,
                                       a_nb_local_rows, a_nb_local_cols ) ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
EXT_PETScMatrix:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: implementation" ) ;

   LA_Implementation const* result = EXT_PETScImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScMatrix:: is_symmetric( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: is_symmetric" ) ;

   bool result = SYMMETRIC ;

   PEL_CHECK_POST( is_symmetric_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: is_desynchronizable" ) ;

   bool result = true ;

   PEL_CHECK_POST( result == true ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: start_local_modifs" ) ;
   PEL_CHECK_PRE( start_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScMatrix:: start_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   PEL_Error::object()->raise_internal( mesg ) ;

   PEL_CHECK_POST( start_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: stop_local_modifs" ) ;
   PEL_CHECK_PRE( stop_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScMatrix:: stop_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   PEL_Error::object()->raise_internal( mesg ) ;

   PEL_CHECK_POST( stop_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: synchronize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: synchronize" ) ;
   PEL_CHECK_PRE( synchronize_PRE() ) ;

   PETSc_do( MatAssemblyBegin( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   PETSc_do( MatAssemblyEnd( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   LA_Matrix::synchronize() ;

   PEL_CHECK_POST( is_assembled() ) ;
   PEL_CHECK_POST( synchronize_POST() ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
EXT_PETScMatrix:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: row_distribution" ) ;
   PEL_CHECK_PRE( row_distribution_PRE() ) ;

   PEL_DistributedPartition const* result = ROW_DIST ;

   PEL_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
EXT_PETScMatrix:: col_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: col_distribution" ) ;
   PEL_CHECK_PRE( col_distribution_PRE() ) ;

   PEL_DistributedPartition const* result = COL_DIST ;

   PEL_CHECK_POST( col_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix*
EXT_PETScMatrix:: create_local_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: create_local_matrix" ) ;
   PEL_CHECK_PRE( create_local_matrix_PRE(a_owner) ) ;

   LA_PelMatrix* result = LA_PelMatrix::create( a_owner, NB_ROWS, NB_COLS ) ;

   PetscInt i0, i1 ;
   PETSc_do( MatGetOwnershipRange( matrix(), &i0, &i1 ) ) ;

   PetscInt nbcols = 0 ;
   const PetscInt* cols = 0 ;
   const PetscScalar* values = 0 ;
   if( UPPER ) PETSc_do( MatGetRowUpperTriangular( matrix() ) ) ;
   for( PetscInt i=i0 ; i<i1 ; ++i )
   {
      PETSc_do( MatGetRow( matrix(), i, &nbcols, &cols, &values ) ) ;
      for( PetscInt j=0 ; j<nbcols ; ++j )
      {
         result->add_to_item( (size_t) i, (size_t) cols[j],
                              (double) values[j] ) ;
      }
      PETSc_do( MatRestoreRow( matrix(), i, &nbcols, &cols, &values ) ) ;
   }
   if( UPPER ) PETSc_do( MatRestoreRowUpperTriangular( matrix() ) ) ;

   result->synchronize() ;

   PEL_CHECK_POST( create_local_matrix_POST(result, a_owner)  ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
EXT_PETScMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL("EXT_PETScMatrix:: nb_stored_items") ;
   PEL_CHECK_PRE( nb_stored_items_PRE() ) ;

   MatInfo info ;
   PETSc_do( MatGetInfo(matrix(),MAT_LOCAL,&info) ) ;

   return( (size_t)info.nz_allocated ) ;
}

//----------------------------------------------------------------------
size_t
EXT_PETScMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
EXT_PETScMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: extract_diag( LA_Vector* diag ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: extract_diag" ) ;
   PEL_CHECK_PRE( extract_diag_PRE( diag ) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector* >( diag ) != 0 ) ;
   EXT_PETScVector* pvec = static_cast<EXT_PETScVector* >( diag ) ;
   PETSc_do( MatGetDiagonal( matrix(), pvec->vector() ) ) ;
   diag->synchronize() ;

   PEL_CHECK_POST( extract_diag_POST( diag )  ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
                                         double alpha, double beta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const* >( x ) != 0 ) ;
   EXT_PETScVector const* px = static_cast<EXT_PETScVector const* >( x ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector* >( y ) != 0 ) ;
   EXT_PETScVector* py = static_cast<EXT_PETScVector* >( y ) ;

   if( beta==0.0 )
   {
      py->synchronize() ;//??
      PETSc_do( MatMult( matrix(), px->vector(), py->vector() ) ) ;
      py->scale( alpha ) ;
   }
   else
   {
      py->scale( beta/alpha ) ;
      PETSc_do( MatMultAdd( matrix(), px->vector(), py->vector(),
                            py->vector() ) ) ;
      py->scale( alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: tr_multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
                                            double alpha, double beta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: tr_multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const* >( x ) != 0 ) ;
   EXT_PETScVector const* px = static_cast<EXT_PETScVector const* >( x ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector* >( y ) != 0 ) ;
   EXT_PETScVector* py = static_cast<EXT_PETScVector* >( y ) ;

   if( beta==0.0 )
   {
      py->synchronize() ;//???
      PETSc_do( MatMultTranspose( matrix(), px->vector(), py->vector() ) ) ;
      py->scale( alpha ) ;
   }
   else
   {
      py->scale( beta/alpha ) ;
      PETSc_do( MatMultTransposeAdd( matrix(), px->vector(), py->vector(),
                                     py->vector() ) ) ;
      py->scale( alpha ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: scale_as_diag_mat_mat( LA_Vector const* lvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: scale_as_diag_mat_mat LA_Vector" ) ;
   PEL_CHECK_PRE( scale_as_diag_mat_mat_PRE( lvec ) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const* >( lvec ) != 0 ) ;
   EXT_PETScVector const* plvec =
                       static_cast<EXT_PETScVector const* >( lvec ) ;

   PETSc_do( MatDiagonalScale( matrix(), plvec->vector(), PETSC_NULL ) ) ;

   PEL_CHECK_POST( scale_as_diag_mat_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: scale_as_mat_diag_mat( LA_Vector const* rvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: scale_as_mat_diag_mat " ) ;
   PEL_CHECK_PRE( scale_as_mat_diag_mat_PRE( rvec ) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const* >( rvec ) != 0 ) ;
   EXT_PETScVector const* prvec =
                       static_cast<EXT_PETScVector const* >( rvec ) ;

   PETSc_do( MatDiagonalScale( matrix(), PETSC_NULL, prvec->vector() ) ) ;

   PEL_CHECK_POST( scale_as_mat_diag_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_to_diag( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_to_diag LA_Vector" ) ;
   PEL_CHECK_PRE( add_to_diag_PRE( vec ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const* >( vec ) != 0 ) ;
   EXT_PETScVector const* pvec = static_cast<EXT_PETScVector const* >( vec ) ;

   PETSc_do( MatDiagonalSet( matrix(), pvec->vector(), ADD_VALUES ) ) ;

   PEL_CHECK_POST( add_to_diag_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: set( LA_SeqMatrix const* A )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: set LA_SeqMatrix" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   NB_ROWS = (int) A->nb_rows() ;
   NB_COLS = (int) A->nb_cols() ;

   // Number of no-zero values par line :
   LA_MatrixIterator* it = A->create_stored_item_iterator(0) ;
   intVector nnz( nb_rows() ) ;
   for( size_t i=0 ; i<nb_rows() ; i++ )
   {
      int nb = 0 ;
      for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
      {
         if( !SYMMETRIC || it->col()>= (size_t) i ) nb++ ;
      }
      nnz( i ) = nb ;
   }

   build( nnz ) ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   // Copy A into the new PETSc matrix
   int first = (int) row_distribution()->first_local_index() ;
   int last = (int) row_distribution()->local_index_limit() ;
   for( int i=first ; i<last ; ++i)
   {
      int nb = nnz( i ) ;
      double* val = new double [nb] ;
      int* col = new int [nb] ;
      nb=0 ;
      for( it->start_row_items( i ) ; it->is_valid() ; it->go_next() )
      {
         if( !SYMMETRIC || it->col() >= (size_t) i )
         {
            val[nb] = it->item() ;
            col[nb] = it->col() ;
            nb++ ;
         }
      }
      PEL_CHECK( nb == nnz( i ) ) ;
      MatSetValues( matrix(), 1, &i, nb, col, val, INSERT_VALUES ) ;
      delete [] val ;
      delete [] col ;
   }
   PETSc_do( MatAssemblyBegin( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   PETSc_do( MatAssemblyEnd( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   it->destroy() ; it = 0 ;

   PEL_CHECK_POST( set_POST( A ) ) ;
   PEL_CHECK_POST( nb_rows()==A->nb_rows() ) ;
   PEL_CHECK_POST( nb_cols()==A->nb_cols() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScMatrix const* >( A ) != 0 ) ;
   EXT_PETScMatrix const* pA = static_cast<EXT_PETScMatrix const* >( A ) ;

   if( same_pattern && ( ! BLOCK ) )
   {
      PETSc_do( MatCopy( pA->matrix(), matrix(),
                         SAME_NONZERO_PATTERN ) ) ;
   }
   else
   {
      PETSc_do( MatCopy( pA->matrix(), matrix(),
                         DIFFERENT_NONZERO_PATTERN ) ) ;
   }
   PETSc_do( MatAssemblyBegin( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   PETSc_do( MatAssemblyEnd( matrix(), MAT_FINAL_ASSEMBLY ) ) ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( A ) ) ;

}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_Mat( LA_Matrix const* A,
                           double alpha,
                           bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScMatrix const* >( A ) != 0 ) ;
   EXT_PETScMatrix const* pA = static_cast<EXT_PETScMatrix const* >( A ) ;
   PEL_CHECK( pA->is_assembled() ) ;

   add_Mat_IMP( pA->matrix(), alpha, same_pattern ) ;

   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_tMat( LA_Matrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_tMat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScMatrix const* >( A ) != 0 ) ;
   EXT_PETScMatrix const* pA = static_cast<EXT_PETScMatrix const* >( A ) ;

   Mat tmp ;
   PETSc_do( MatTranspose( pA->matrix(), &tmp ) ) ;

   add_Mat_IMP( tmp, alpha, false ) ;

   PETSc_do( MatDestroy( tmp ) ) ;
   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                               double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_Mat_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScMatrix const* >( A ) != 0 ) ;
   EXT_PETScMatrix const* pA = static_cast<EXT_PETScMatrix const* >( A ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScMatrix const* >( B ) != 0 ) ;
   EXT_PETScMatrix const* pB = static_cast<EXT_PETScMatrix const* >( B ) ;

   // Initialement fill = 2.0 : pourquoi ?
   Mat tmp ;
   PETSc_do( MatMatMult( pA->matrix(), pB->matrix(),
                         MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp ) ) ;

   add_Mat_IMP( tmp, alpha, false ) ;
   PETSc_do( MatDestroy( tmp ) ) ;

   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: set_item" ) ;
   PEL_CHECK_PRE( set_item_PRE( i, j ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }

   MatSetValue( matrix(), i, j, x, INSERT_VALUES ) ;


   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_to_item" ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }
   MatSetValue( matrix(), i, j, x, ADD_VALUES ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: nullify( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: nullify" ) ;

   set_unsynchronized_state( LA::NotSync_undef ) ;
   PETSc_do( MatZeroEntries( matrix() ) ) ;

   PEL_CHECK_POST( nullify_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;

   PETSc_do( MatScale( matrix(), alpha ) ) ;

   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: readMM( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: readMM" ) ;
   PEL_CHECK_PRE( readMM_PRE( file ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   // Avoid unefficient call to insert function:
   LA_PelMatrix* m = LA_PelMatrix::create( 0, 0, 0 ) ;
   m->readMM( file ) ;
   set( m ) ;
   m->destroy() ; m = 0 ;
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( readMM_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: print" ) ;

   LA_Matrix::print( os, indent_width ) ;

   if( EXP != 0 )
   {
      std::string const& mat_type = name() ;
      std::string s( indent_width, ' ' ) ;
      if( SEQ )
      {
         if( EXP->has_entry( "nz" ) )
         {
            os << s << "   nz : " << EXP->int_data( "nz" ) << std::endl ;
         }
         if( mat_type == "PETSc_SeqSBAIJ" )
         {
            os << s << "   block_size : "
               << EXP->int_data( "block_size" ) << std::endl ;
         }
      }
      else
      {
         os << s << "   d_nz : " << EXP->int_data( "d_nz" ) << std::endl ;
         os << s << "   o_nz : " << EXP->int_data( "o_nz" ) << std::endl ;
         if( mat_type == "PETSc_MPISBAIJ" || mat_type == "PETSc_MPIBAIJ" )
         {
            os << s << "   block_size : "
               << EXP->int_data( "block_size" ) << std::endl ;
         }
      }
      if( EXP->has_entry( "subtype" ) )
      {
         os << s << "  subtype : \""
            << EXP->string_data( "subtype" ) << "\"" << std::endl ;
      }
   }
}

//----------------------------------------------------------------------
Mat const&
EXT_PETScMatrix:: matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: matrix" ) ;
   PEL_CHECK( MATRIX !=0 ) ;
   return( MATRIX ) ;
}


//----------------------------------------------------------------------
bool
EXT_PETScMatrix:: is_assembled( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: is_assembled" ) ;

   PetscTruth assembled ;
   PETSc_do( MatAssembled( matrix(), &assembled ) ) ;
   bool result = assembled==PETSC_TRUE ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: build( intVector const& nnz )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: build_sequential_matrix" ) ;

   destroy_matrix() ;

   UPPER = SYMMETRIC ;

   std::string const& mat_type = name()  ;

   if( SYMMETRIC && NB_ROWS!=NB_COLS )
   {
      PEL_Error::object()->raise_plain(
         "*** EXT_PETScMatrix error:\n"
         "    Unable to create a PETSc symmetric matrix with\n"
         "    a no square matrix" ) ;
   }


   if( SEQ )
   {
      int const*  nnzptr = ( nnz.size()==NB_ROWS ? nnz.data() : PETSC_NULL ) ;
      int nz = ( EXP->has_entry( "nz" ) && nnzptr == PETSC_NULL ?
                                     EXP->int_data( "nz" ) : PETSC_DEFAULT ) ;
      if( mat_type == "PETSc_SeqSBAIJ" )
      {
         int block_size = EXP->int_data( "block_size" ) ;
         PETSc_do(
            MatCreateSeqSBAIJ(
               PETSC_COMM_WORLD,
               block_size, NB_ROWS, NB_COLS,
               nz, nnzptr,
               &MATRIX ) );
      }
      else if( mat_type == "PETSc_SeqAIJ" || mat_type == "PETSc_AIJ" )
      {
         PETSc_do(
            MatCreateSeqAIJ(
               PETSC_COMM_WORLD,
               NB_ROWS, NB_COLS,
               nz, nnzptr,
               &MATRIX ) );
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value(
            EXP, "type",
            "  - \"PETSc_AIJ\"\n"
            "  - \"PETSc_SeqSBAIJ\"\n"
            "  - \"PETSc_SeqAIJ\"" ) ;
      }
   }
   else
   {
      int d_nz = EXP->int_data( "d_nz" ) ;
      int o_nz = EXP->int_data( "o_nz" ) ;
      int* d_nnz = PETSC_NULL ;
      int* o_nnz = PETSC_NULL ;

      if( mat_type == "PETSc_MPIAIJ" || mat_type == "PETSc_AIJ" )
      {
         PETSc_do(
            MatCreateMPIAIJ(
               PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
               NB_ROWS, NB_COLS,
               d_nz, d_nnz, o_nz, o_nnz,
               &MATRIX ) ) ;
      }
      else if( mat_type == "PETSc_MPISBAIJ" )
      {
         int block_size = EXP->int_data( "block_size" ) ;
         PETSc_do(
            MatCreateMPISBAIJ(
               PETSC_COMM_WORLD, block_size, PETSC_DECIDE, PETSC_DECIDE,
               NB_ROWS, NB_COLS,
               d_nz, d_nnz, o_nz, o_nnz,
               &MATRIX ) ) ;
      }
      else if( mat_type == "PETSc_MPIBAIJ" )
      {
         int block_size = EXP->int_data( "block_size" ) ;
         PETSc_do(
            MatCreateMPIBAIJ(
               PETSC_COMM_WORLD, block_size, PETSC_DECIDE, PETSC_DECIDE,
               NB_ROWS, NB_COLS,
               d_nz, d_nnz, o_nz, o_nnz,
               &MATRIX ) ) ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value(
            EXP, "type",
            "   - \"PETSc_MPIAIJ\"\n"
            "   - \"PETSc_MPISBAIJ\"\n"
            "   - \"PETSc_MPIBAIJ\"\n" ) ;
      }

   }

   if( EXP->has_entry( "subtype" ) )
   {
      std::string const& subtype = EXP->string_data( "subtype" ) ;
      if( subtype!=MATSUPERLU_DIST && subtype!=MATUMFPACK && subtype!=MATAIJSPOOLES )
      {
         PEL_Error::object()->raise_bad_data_value(
            EXP, "subtype",
            "   - \"" MATSUPERLU_DIST "\"\n"
            "   - \"" MATUMFPACK "\"\n"
            "   - \"" MATAIJSPOOLES "\"" ) ;
      }
      PETSc_do( MatSetType( matrix(), subtype.c_str() ) ) ;
   }
//    if( EXP->has_entry( "symmetric" ) && EXP->bool_data( "symmetric" ) )
//    {
//       PETSc_do( MatSetOption( matrix(), MAT_SYMMETRIC ) ) ;
//    }

   if( HAS_OPT ) PETSc_do( MatSetFromOptions( matrix() ) ) ;

   PetscInt m, n ;
   PETSc_do( MatGetLocalSize( matrix(), &m, &n ) ) ;
   ROW_DIST->set_local_number( m ) ;
   COL_DIST->set_local_number( n ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: destroy_matrix( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: destroy_matrix" ) ;

   if( MATRIX!=0 )
   {
      if( DESTROY_ON_EXIT ) PETSc_do( MatDestroy( MATRIX ) ) ;
      MATRIX = 0 ;
   }

   PEL_CHECK_POST( MATRIX == 0 ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScMatrix:: add_Mat_IMP( Mat A, double alpha, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScMatrix:: add_Mat_IMP" ) ;
   PEL_CHECK_PRE( state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( A != 0 ) ;

   MatInfo info ;
   PETSc_do( MatGetInfo( MATRIX, MAT_GLOBAL_SUM, &info ) ) ;
   double nb_MATRIX = info.nz_used ;

   PETSc_do( MatGetInfo( A, MAT_GLOBAL_SUM, &info ) ) ;
   double nb_source = info.nz_used ;

   PetscScalar ps = alpha ;

   if( nb_MATRIX == nb_source && same_pattern && ( !BLOCK ) )
   {
      PETSc_do( MatAXPY( MATRIX, ps, A, SAME_NONZERO_PATTERN ) ) ;
   }
   else if( nb_MATRIX >= nb_source )
   {
      PETSc_do( MatAXPY( MATRIX, ps, A, DIFFERENT_NONZERO_PATTERN ) ) ;
   }
   else
   {
      LA::SyncState old_state = state() ;
      if( !is_synchronized() )
      {
         synchronize() ;
      }
      Mat tmp ;
      PETSc_do( MatDuplicate( A, MAT_COPY_VALUES, &tmp ) ) ;
      PETSc_do( MatAYPX( tmp, ps, MATRIX, DIFFERENT_NONZERO_PATTERN ) ) ;
      PETSc_do( MatDestroy( MATRIX ) ) ;
      MATRIX = tmp ;
      if( old_state != LA::Sync ) set_unsynchronized_state( old_state ) ;
   }
}

//----------------------------------------------------------------------
bool
EXT_PETScMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::invariant() ) ;
   PEL_ASSERT( !only_local_modifs() ) ;
   PEL_ASSERT( IMPLIES( MATRIX==0,  ( NB_ROWS == 0 && NB_COLS == 0 ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScMatrix:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::implementation_POST( result ) ) ;
   PEL_ASSERT( result == EXT_PETScImplementation::object() ) ;
   return( true ) ;
}
