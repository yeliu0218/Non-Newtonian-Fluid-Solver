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

#include <LA_SeqMatrix.hh>

#include <LA_MatrixIterator.hh>
#include <LA_PelMatrix.hh>
#include <LA_SeqVector.hh>
#include <LA_SeqImplementation.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <fstream>
#include <sstream>
#include <string>

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_SeqMatrix:: make( PEL_Object* a_owner,
                     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   LA_SeqMatrix const* proto =
      static_cast<LA_SeqMatrix const*>(
                               sparse_mat_pluggin_map()->item( name ) ) ;
   PEL_ASSERT( proto!=0 ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   LA_SeqMatrix* result = proto->create_replica( a_owner, exp ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( IMPLIES( result->is_desynchronizable(),
                            result->state() == LA::NotSync_undef ) ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix:: LA_SeqMatrix( PEL_Object* a_owner,
                             std::string const& a_name )
//----------------------------------------------------------------------
   : LA_Matrix( a_owner, a_name, LA::NoDistribution )
   , ROW_DIST( 0 )
   , COL_DIST( 0 )
{
   PEL_LABEL( "LA_SeqMatrix:: LA_SeqMatrix" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::NoDistribution ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix:: LA_SeqMatrix( std::string const& a_name )
//----------------------------------------------------------------------
   : LA_Matrix( a_name )
   , ROW_DIST( 0 )
   , COL_DIST( 0 )
{
   PEL_LABEL( "LA_SeqMatrix:: LA_SeqMatrix( prototype )" ) ;
   sparse_mat_pluggin_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( !is_resizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix:: ~LA_SeqMatrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_SeqMatrix:: create_copy( PEL_Object* a_owner,
                            LA_SeqMatrix const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: create_copy" ) ;
   PEL_CHECK_PRE( create_copy_PRE( a_owner, other ) ) ;

   if( !is_resizable() )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_SeqMatrix:: create_copy error:\n"
         "    create copy of non resizable matrix is not implemented" ) ;
   }

   LA_SeqMatrix* result = create_matrix( a_owner ) ;
   result->re_initialize_with_global_sizes( other->nb_rows(),
                                            other->nb_cols() ) ;
   result->set( other ) ;

   PEL_CHECK_POST( create_copy_POST( result, a_owner, other ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
void
LA_SeqMatrix:: re_initialize( size_t a_nb_rows, size_t a_nb_cols,
                              size_t a_nb_local_rows, size_t a_nb_local_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: re_initialize" ) ;
   PEL_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_cols,
                                     a_nb_local_rows, a_nb_local_cols ) ) ;

   re_initialize_with_global_sizes( a_nb_rows, a_nb_cols ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_cols,
                                       a_nb_local_rows, a_nb_local_cols ) ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqMatrix:: create_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: create_vector" ) ;
   PEL_CHECK_PRE( create_vector_PRE( a_owner ) ) ;

   if( !is_resizable() )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_SeqMatrix:: create_vector error:\n"
         "    create vector of non resizable matrix is not implemented" ) ;
   }

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_rows() ) ;
   if( is_desynchronizable() ) result->make_desynchronizable() ;

   PEL_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return result ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_SeqMatrix:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: implementation" ) ;

   LA_Implementation const* result = LA_SeqImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
LA_SeqMatrix:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: row_distribution" ) ;
   PEL_CHECK_PRE( row_distribution_PRE() ) ;

   if( ROW_DIST == 0 )
   {
      ROW_DIST = PEL_DistributedPartition::create(
                                const_cast<LA_SeqMatrix*>( this ) ) ;
   }
   if( ROW_DIST->global_number() != nb_rows() )
   {
      ROW_DIST->set_global_number( nb_rows() ) ;
   }

   PEL_DistributedPartition const* result = ROW_DIST ;

   PEL_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
LA_SeqMatrix:: col_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: col_distribution" ) ;
   PEL_CHECK_PRE( col_distribution_PRE() ) ;

   if( COL_DIST == 0 )
   {
      COL_DIST = PEL_DistributedPartition::create(
                                const_cast<LA_SeqMatrix*>( this ) ) ;
   }
   if( COL_DIST->global_number() != nb_cols() )
   {
      COL_DIST->set_global_number( nb_cols() ) ;
   }

   PEL_DistributedPartition const* result = COL_DIST ;

   PEL_CHECK_POST( col_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_SeqMatrix:: create_local_matrix( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: create_local_matrix" ) ;
   PEL_CHECK_PRE( create_local_matrix_PRE( a_owner ) ) ;

   LA_SeqMatrix* result = create_matrix( a_owner ) ;
   result->set( this ) ;
   result->synchronize() ;

   PEL_CHECK_POST( create_local_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: extract_diag( LA_Vector* diag ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: extract_diag" ) ;
   PEL_CHECK_PRE( extract_diag_PRE( diag ) ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      diag->set_item( i, item( i, i ) ) ;
   }
   diag->synchronize() ;

   PEL_CHECK_POST( extract_diag_POST( diag ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: extract_row( size_t i, LA_Vector* row ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: extract_row" ) ;
   PEL_CHECK_PRE( extract_row_PRE( i, row ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t j=0 ; j<nb_cols() ; ++j )
   {
      row->set_item( j, item( i, j ) ) ;
   }
   row->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extract_row_POST( i, row ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: extract_col( size_t j, LA_Vector* col ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: extract_col" ) ;
   PEL_CHECK_PRE( extract_col_PRE( j, col ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      col->set_item( i, item( i, j ) ) ;
   }
   col->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extract_col_POST( j, col ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: extract_lump( LA_Vector* lump ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: extract_lump" ) ;
   PEL_CHECK_PRE( extract_lump_PRE( lump ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( size_t i=0 ; i<nb_rows() ; i++ )
   {
      double ss = 0.0 ;
      for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
      {
         ss += PEL::abs( it->item() ) ;
      }
      lump->set_item( i, ss ) ;
   }
   it->destroy() ; it = 0 ;
   lump->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( extract_lump_POST( lump ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: nullify( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: nullify" ) ;
   PEL_CHECK_PRE( nullify_PRE() ) ;

   set_stored_items( 0.0 ) ;

   PEL_CHECK_POST( nullify_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: nullify_row( size_t i )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: nullify_row" ) ;
   PEL_CHECK_PRE( nullify_row_PRE( i ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
   {
      it->set_item( 0. ) ;
   }
   it->destroy() ; it = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nullify_row_POST( i ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( alpha != 1. )
   {
      LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
      start_local_modifs() ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         it->set_item( it->item()*alpha ) ;
      }
      stop_local_modifs() ;
      it->destroy() ; it = 0 ;
   }

   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha, double beta ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;

   if( beta != 1. ) y->scale( beta ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      y->add_to_item( it->row(),
                      alpha * it->item() * bx->item( it->col() ) ) ;
   }
   it->destroy() ; it = 0 ;
   y->synchronize() ;

   PEL_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;

}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: tr_multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha, double beta )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: tr_multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;

   if( beta != 1. ) y->scale( beta ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      y->add_to_item( it->col(),
                      alpha * it->item() * bx->item( it->row() ) ) ;
   }
   it->destroy() ; it = 0 ;
   y->synchronize() ;

   PEL_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: scale_as_diag_mat_mat( LA_Vector const* lvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: scale_as_diag_mat_mat" ) ;
   PEL_CHECK_PRE( scale_as_diag_mat_mat_PRE( lvec ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( lvec ) != 0 ) ;
   LA_SeqVector const* dx = static_cast<LA_SeqVector const* >( lvec ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      it->set_item( it->item()*dx->item( it->row() ) ) ;
   }
   stop_local_modifs() ;
   it->destroy() ; it = 0 ;

   PEL_CHECK_POST( scale_as_diag_mat_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: scale_as_mat_diag_mat( LA_Vector const* rvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: scale_as_mat_diag_mat" ) ;
   PEL_CHECK_PRE( scale_as_mat_diag_mat_PRE( rvec ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( rvec ) != 0 ) ;
   LA_SeqVector const* dx = static_cast<LA_SeqVector const* >( rvec ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      it->set_item( it->item()*dx->item( it->col() ) ) ;
   }
   stop_local_modifs() ;
   it->destroy() ; it = 0 ;

   PEL_CHECK_POST( scale_as_mat_diag_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_to_diag( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_to_diag LA_Vector" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_to_diag_PRE( vec ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( vec ) != 0 ) ;
   LA_SeqVector const* dx = static_cast<LA_SeqVector const* >( vec ) ;

   start_local_modifs() ;
   for( size_t i=0 ; i<dx->nb_rows() ; ++i )
   {
      add_to_item( i, i, dx->item(i) ) ;
   }
   stop_local_modifs() ;

   PEL_CHECK_POST( add_to_diag_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   set_IMP( AS, same_pattern ) ;

   PEL_CHECK_POST( set_POST(A) ) ;
}


//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat( LA_Matrix const* A, double alpha, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   add_Mat_IMP( AS, alpha, same_pattern ) ;

   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat( LA_Matrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( A->is_symmetric() )
   {
      add_Mat( A, alpha ) ;
   }
   else
   {
      PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
      LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;
      add_tMat_IMP( AS, alpha ) ;
   }

   PEL_CHECK_POST( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                            double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat_Mat" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( B ) != 0 ) ;
   LA_SeqMatrix const* BS = static_cast<LA_SeqMatrix const*>(B) ;

   add_Mat_Mat_IMP( AS, BS, alpha ) ;

   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                             double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat_Mat" ) ;
   PEL_CHECK_PRE( add_tMat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( B ) != 0 ) ;
   LA_SeqMatrix const* BS = static_cast<LA_SeqMatrix const*>(B) ;

   add_tMat_Mat_IMP( AS, BS, alpha ) ;

   PEL_CHECK_POST( add_tMat_Mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                             double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat_tMat" ) ;
   PEL_CHECK_PRE( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( B ) != 0 ) ;
   LA_SeqMatrix const* BS = static_cast<LA_SeqMatrix const*>(B) ;

   add_Mat_tMat_IMP( AS, BS, alpha ) ;

   PEL_CHECK_POST( add_Mat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                              double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat_tMat" ) ;
   PEL_CHECK_PRE( add_tMat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
   LA_SeqMatrix const* AS = static_cast<LA_SeqMatrix const*>(A) ;

   PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( B ) != 0 ) ;
   LA_SeqMatrix const* BS = static_cast<LA_SeqMatrix const*>(B) ;

   add_tMat_tMat_IMP( AS, BS, alpha ) ;

   PEL_CHECK_POST( add_tMat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: writeMM( std::string const& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: writeMM" ) ;
   PEL_CHECK_PRE( !file.empty() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::ofstream out( file.c_str() ) ; //??????
   if( !out )
   {
      PEL_Error::object()->raise_file_handling( file, "open" ) ;
   }
   out << "%%MatrixMarket matrix coordinate real general" << std::endl ;
   out << nb_rows() << " " << nb_cols()
       << " " << nb_stored_items() << std::endl ;
   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   out.precision( 16 ) ;
   for( size_t i=0 ; i<nb_rows() ; i++ )
   {
      for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
      {
         out << it->row() + 1 << " "
             << it->col() + 1 << " "
             << it->item() << std::endl ;
      }
   }
   it->destroy() ; it=0 ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: save( PEL_Module* mod ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: save" ) ;
   PEL_CHECK_PRE( save_PRE( mod ) ) ;

   mod->add_entry( "sparse",    PEL_Bool::create( mod, false ) ) ;
   mod->add_entry( "symmetric", PEL_Bool::create( mod, is_symmetric() ) ) ;
   mod->add_entry( "nb_rows",   PEL_Int::create( mod, nb_rows() ) ) ;
   mod->add_entry( "nb_cols",   PEL_Int::create( mod, nb_cols() ) ) ;
   doubleVector V(0) ;
   if( !is_symmetric() )
   {
      V.re_initialize( nb_rows()*nb_cols() ) ;
      size_t k = 0 ;
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         for( size_t j=0 ; j<nb_cols() ; ++j)
         {
            V( k++ ) = item( i, j ) ;
         }
      }
   }
   else
   {
      V.re_initialize( nb_rows()*(nb_rows()+1)/2 ) ;
      size_t k = 0 ;
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         for( size_t j=i ; j<nb_cols() ; ++j)
         {
            V( k++ ) = item( i, j ) ;
         }
      }
   }
   mod->add_entry( "row_values",
                   PEL_DoubleVector::create( mod, V ) ) ;

   PEL_CHECK_POST( save_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: restore( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: restore" ) ;
   PEL_CHECK_PRE( restore_PRE( exp ) ) ;

   bool sym = exp->bool_data( "symmetric" ) ;
   bool sparse = exp->bool_data( "sparse" ) ;
   if( sparse )
   {
      PEL_Error::object()->raise_internal(
         "Sparse format not implemented" ) ;
   }
   if( sym != is_symmetric() )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_SeqMatrix:: restore error:\n"
         "    matrix internal structure is no compatible\n"
         "    with symmetry saving" ) ;
   }
   if( !is_resizable() )
   {
      PEL_Error::object()->raise_internal(
         "Reading of non resizable matrix is not implemented" ) ;
   }

   doubleVector const& values = exp->doubleVector_data( "row_values" ) ;
   size_t n = (size_t) exp->int_data( "nb_rows" ) ;
   size_t m = (size_t) exp->int_data( "nb_cols" ) ;

   if( sym )
   {
      PEL_ASSERT( n == m ) ;
      PEL_ASSERT( values.size() == n*(n+1)/2 ) ;
   }
   else
   {
      PEL_ASSERT( values.size()==n*m ) ;
   }

   re_initialize_with_global_sizes( n, m ) ;

   if( !sym )
   {
      size_t k=0 ;
      for( size_t i=0 ; i<n ; ++i )
      {
         for( size_t j=0 ; j<m ; ++j )
         {
            set_item( i, j, values( k++ ) );
         }
      }
   }
   else
   {
      size_t k=0 ;
      for( size_t i=0 ; i<n ; ++i )
      {
         for( size_t j=i ; j<m ; ++j )
         {
            set_item( i, j, values( k++ ) );
         }
      }
   }
   synchronize() ;

   PEL_CHECK_POST( restore_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: solve_LU( LA_SeqVector const* rhs, LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: solve_LU" ) ;
   PEL_CHECK_PRE( solve_LU_PRE( rhs, sol ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;

   size_t n = nb_rows() ;
   size_t i, j ;

   for( i=0 ; i<n ; i++ )
   {
      double r = rhs->item( i ) ;
      for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
      {
         j = it->col() ;
         if( j<i )
         {
            r -= it->item() * sol->item( j ) ;
         }
      }
      sol->set_item( i, r ) ;
   }

   for( i = n-1 ; i<n ; i-- )
   {
      double r = sol->item( i ) ;
      double d = 1. ;
      for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
      {
         j = it->col() ;
         if( j==i )
         {
            d = it->item() ;
         }
         else if( j>i )
         {
            r -= it->item() * sol->item( j ) ;
         }
      }
      sol->set_item( i, r/d ) ;
   }

   it->destroy() ;  it = 0 ;
   sol->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( solve_LU_POST( rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: relax( double omega,
                      LA_SeqMatrix::relaxation_mode mode,
                      LA_SeqVector const* omega_inv_diag,
                      LA_SeqVector const* rhs,
                      LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: relax" ) ;
   PEL_CHECK_PRE( relax_PRE( omega, mode, omega_inv_diag, rhs, sol ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;

   size_t const n = nb_rows() ;

   if( mode == forward || mode == symmetric )
   {
      for( size_t i=0 ; i<n ; ++i )
      {
         double xt =  rhs->item(i) ;
         for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
         {
            xt -= it->item()*sol->item( it->col() ) ;
         }
         sol->add_to_item( i, xt*omega_inv_diag->item(i) ) ;
      }
   }

   if( mode == backward || mode == symmetric )
   {
      for( size_t i=n-1 ; i<n ; --i )
      {
         double xt =  rhs->item(i) ;
         for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
         {
            xt -= it->item()*sol->item( it->col() ) ;
         }
         sol->add_to_item( i, xt*omega_inv_diag->item(i) ) ;
      }
   }

   it->destroy() ;  it = 0 ;
   sol->synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( relax_POST( omega, mode, omega_inv_diag, rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: factorize_MILU0( bool modified, double piv_min )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: factorize_MILU0" ) ;
   PEL_CHECK_PRE( factorize_MILU0_PRE( modified, piv_min ) ) ;

   size_t const n = nb_rows() ;

   LA_MatrixIterator* itIK = create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itIJ = create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itKJ = create_stored_item_iterator( 0 ) ;

   if( PEL::abs( item( 0, 0 ) )<piv_min )
   {
      nullify_row( 0 ) ;
      set_item( 0, 0, 1. ) ;
      raise_MILU0_zero_pivot( 0 ) ;
   }

   for( size_t i=1 ; i<n ; ++i )
   {
      double drop = 0.0 ;
      for( itIK->start_row_items(i) ; itIK->is_valid() ; itIK->go_next() )
      {
         size_t const k = itIK->col() ;
         if( k<i )
         {
            double const akk = item( k, k ) ;
            double const aik = itIK->item()/akk ;
            itIK->set_item( aik ) ;
            for( itKJ->start_row_items(k) ; itKJ->is_valid() ; itKJ->go_next() )
            {
               size_t const j = itKJ->col() ;
               if( j>k )
               {
                  double const akj = itKJ->item() ;
                  bool found = false ;
                  for( itIJ->start_row_items(i) ; itIJ->is_valid() ; itIJ->go_next() )
                  {
                     if( itIJ->col() == j )
                     {
                        found = true ;
                        itIJ->add_to_item( -aik*akj ) ;
                        break ;
                     }
                  }
                  if( !found ) drop -= aik*akj ;
               }
            }
         }
      }
      double aii = item( i, i ) ;
      if( modified )
      {
         aii -= drop ;
         add_to_item( i, i, -drop ) ;
      }
      if( PEL::abs( aii )<piv_min )
      {
         nullify_row( i ) ;
         set_item( i, i, 1. ) ;
         raise_MILU0_zero_pivot( i ) ;
      }
   }
   itIK->destroy() ; itIK = 0 ;
   itIJ->destroy() ; itIJ = 0 ;
   itKJ->destroy() ; itKJ = 0 ;
   synchronize() ;

   PEL_CHECK_POST( factorize_MILU0_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: raise_MILU0_zero_pivot( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: raise_MILU0_zero_pivot" ) ;
   PEL_CHECK_PRE( i<nb_rows() ) ;

   std::ostringstream mesg ;
   mesg << "*** ILU(0) factorization :" << std::endl ;
   mesg << "***   zero pivot at step " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::invariant() ) ;
   PEL_ASSERT( distribution_strategy() == LA::NoDistribution ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: create_copy_PRE( PEL_Object* a_owner,
                                LA_SeqMatrix const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( other != 0 ) ;
   PEL_ASSERT( other->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: create_copy_POST( LA_SeqMatrix* result,
                                 PEL_Object* a_owner,
                                 LA_SeqMatrix const* other ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->is_synchronized() ) ;
   PEL_ASSERT( result->name() == name() ) ;
   PEL_ASSERT( result->nb_rows() == other->nb_rows() ) ;
   PEL_ASSERT( result->nb_cols() == other->nb_cols() ) ;
   PEL_ASSERT( result->nb_stored_items() == other->nb_stored_items() ) ;
   PEL_ASSERT( result->is_resizable() ) ;
   PEL_ASSERT( result->distribution_strategy() == LA::NoDistribution ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: item_PRE( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( j < nb_cols() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_row_PRE( size_t i, LA_Vector* row ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( row != 0 ) ;
   PEL_ASSERT( row->nb_rows() == nb_cols() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_row_POST( size_t i, LA_Vector* row ) const
//----------------------------------------------------------------------
{

   PEL_ASSERT( row->is_synchronized() ) ;
   PEL_ASSERT(
      FORMAL( FORALL( ( size_t j=0 ; j<nb_cols() ; ++j ),
                      row->item( j ) == item( i, j ) ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_col_PRE( size_t j, LA_Vector* col ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( j < nb_cols() ) ;
   PEL_ASSERT( col != 0 ) ;
   PEL_ASSERT( col->nb_rows() == nb_rows() ) ;
   return( true ) ;
}


//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_col_POST( size_t j, LA_Vector* col ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( col->is_synchronized() ) ;
   PEL_ASSERT(
      FORMAL( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                      col->item( i ) == item( i, j ) ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_lump_PRE( LA_Vector* lump ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( lump != 0 ) ;
   PEL_ASSERT( lump->nb_rows() == nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: extract_lump_POST( LA_Vector* lump ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( lump->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: nullify_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::nullify_POST() ) ;

   LA_MatrixIterator* it_mat = create_stored_item_iterator( 0 ) ;
   PEL_ASSERT(
      FORALL( ( it_mat->start_all_items() ;
                it_mat->is_valid() ;
                it_mat->go_next() ),
              it_mat->item() == 0. ) ) ;
   it_mat->destroy() ; it_mat = 0 ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: nullify_row_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: nullify_row_POST( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( FORALL( ( size_t j=0 ; j<nb_cols() ; j++ ),
                       item( i, j ) == 0.0 ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        state() == LA::NotSync_undef ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        !is_synchronized() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_tMat_Mat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->nb_cols() == nb_rows() ) ;
   PEL_ASSERT( B != 0 ) ;
   PEL_ASSERT( B != this ) ;
   PEL_ASSERT( B->implementation() == implementation() ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( B->nb_cols() == nb_cols() ) ;
   PEL_ASSERT( A->nb_rows() == B->nb_rows() ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_tMat_Mat_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_Mat_tMat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( B != 0 ) ;
   PEL_ASSERT( B != this ) ;
   PEL_ASSERT( B->implementation() == implementation() ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( B->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == B->nb_cols() ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_Mat_tMat_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_tMat_tMat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                  double alpha ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->nb_cols() == nb_rows() ) ;
   PEL_ASSERT( B != 0 ) ;
   PEL_ASSERT( B != this ) ;
   PEL_ASSERT( B->implementation() == implementation() ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( B->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( A->nb_rows() == B->nb_cols() ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_tMat_tMat_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: set_stored_items_POST( double val ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_desynchronizable(), state() == LA::NotSync_undef ) ) ;
   LA_MatrixIterator* it_mat = create_stored_item_iterator( 0 ) ;
   PEL_ASSERT(
      FORALL( ( it_mat->start_all_items() ;
                it_mat->is_valid() ;
                it_mat->go_next() ),
              it_mat->item() == val ) ) ;
   it_mat->destroy() ; it_mat = 0 ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: set_PRE( LA_Matrix const* A, bool same_pattern ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::set_PRE( A, same_pattern )  ) ;
   PEL_ASSERT( !same_pattern || has_same_pattern( A ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: add_Mat_PRE( LA_Matrix const* A,
                            double alpha,
                            bool same_pattern ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_ASSERT( !same_pattern || has_same_pattern( A ) ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: factorize_MILU0_PRE(
                                   bool modified, double piv_min ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( nb_rows() !=0 ) ;
   PEL_ASSERT( nb_cols() == nb_rows() ) ;
   PEL_ASSERT( piv_min >= 0. ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: factorize_MILU0_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: solve_LU_PRE(
                LA_SeqVector const* rhs, LA_SeqVector const* sol ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( nb_cols() == nb_rows() ) ;
   PEL_ASSERT( nb_rows() !=0 ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
                       item( i, i ) != 0.0 ) ) ;
   PEL_ASSERT( rhs != 0 ) ;
   PEL_ASSERT( rhs->is_synchronized() ) ;
   PEL_ASSERT( rhs->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( sol != 0 ) ;
   PEL_ASSERT( sol->nb_rows() == nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: solve_LU_POST(
                LA_SeqVector const* rhs, LA_SeqVector const* sol ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( sol->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: relax_PRE(
                       double omega,
                       LA_SeqMatrix::relaxation_mode mode,
                       LA_SeqVector const* omega_inv_diag,
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( nb_cols() == nb_rows() ) ;
   PEL_ASSERT( nb_rows() !=0 ) ;
   PEL_ASSERT( rhs != 0 ) ;
   PEL_ASSERT( rhs->is_synchronized() ) ;
   PEL_ASSERT( rhs->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( sol != 0 ) ;
   PEL_ASSERT( sol->is_synchronized() ) ;
   PEL_ASSERT( sol->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( omega>0. && omega<2. ) ;
   PEL_ASSERT( omega_inv_diag!= 0 ) ;
   PEL_ASSERT( omega_inv_diag->is_synchronized() ) ;
   PEL_ASSERT( omega_inv_diag->nb_rows() == nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: relax_POST(
                       double omega,
                       LA_SeqMatrix::relaxation_mode mode,
                       LA_SeqVector const* omega_inv_diag,
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( sol->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: set_IMP( LA_SeqMatrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: set_IMP" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( this != A )
   {
      bool not_sym = !is_symmetric() ;
      nullify() ;
      LA_MatrixIterator* it = A->create_stored_item_iterator( 0 ) ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         size_t i = it->row() ;
         size_t j = it->col() ;
         if( not_sym || i<=j ) set_item( i, j, it->item() ) ;
      }
      it->destroy() ; it = 0 ;
   }
   synchronize() ;

   PEL_CHECK_POST( set_POST(A) ) ;
}


//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat_IMP( LA_SeqMatrix const* A, double alpha,
                            bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* itA = A->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itA->start_all_items() ; itA->is_valid() ; itA->go_next() )
   {
      add_to_item( itA->row(), itA->col(), alpha*itA->item() ) ;
   }
   stop_local_modifs() ;
   itA->destroy() ; itA = 0 ;

   PEL_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat_IMP( LA_SeqMatrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* itA = A->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itA->start_all_items() ; itA->is_valid() ; itA->go_next() )
   {
      add_to_item( itA->col(), itA->row(), alpha*itA->item() ) ;
   }
   stop_local_modifs() ;
   itA->destroy() ; itA = 0 ;

   PEL_CHECK_POST( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat_Mat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                                double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat_Mat_IMP" ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* itA = A->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itB = B->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itA->start_all_items() ; itA->is_valid() ; itA->go_next() )
   {
      size_t const i = itA->row() ;
      size_t const k = itA->col() ;
      for( itB->start_row_items( k ) ; itB->is_valid() ; itB->go_next() )
      {
         size_t const j = itB->col() ;
         add_to_item( i, j, alpha * itA->item() * itB->item() ) ;
      }
   }
   stop_local_modifs() ;
   itA->destroy() ; itA = 0 ;
   itB->destroy() ; itB = 0 ;

   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat_Mat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                                 double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat_Mat_IMP" ) ;
   PEL_CHECK_PRE( add_tMat_Mat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* itA = A->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itB = B->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itA->start_all_items() ; itA->is_valid() ; itA->go_next() )
   {
      size_t const i = itA->col() ;
      size_t const k = itA->row() ;
      for( itB->start_row_items( k ) ; itB->is_valid() ; itB->go_next() )
      {
         size_t const j = itB->col() ;
         add_to_item( i, j, alpha*itA->item()*itB->item() ) ;
      }
   }
   stop_local_modifs() ;
   itA->destroy() ; itA = 0 ;
   itB->destroy() ; itB = 0 ;

   PEL_CHECK_POST( add_tMat_Mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_Mat_tMat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                                 double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_Mat_tMat_IMP" ) ;
   PEL_CHECK_PRE( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_PelMatrix* tB = LA_PelMatrix::create( 0, B->nb_cols(), A->nb_rows() ) ;
   tB->add_tMat( B ) ;
   tB->synchronize() ;
   add_Mat_Mat_IMP( A, tB, alpha ) ;
   tB->destroy() ; tB = 0 ;

   PEL_CHECK_POST( add_Mat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqMatrix:: add_tMat_tMat_IMP( LA_SeqMatrix const* A,
                                  LA_SeqMatrix const* B,
                                  double alpha  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: add_tMat_tMat_IMP" ) ;
   PEL_CHECK_PRE( add_tMat_tMat_PRE( A, B, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_MatrixIterator* itA = A->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itB = B->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itB->start_all_items() ; itB->is_valid() ; itB->go_next() )
   {
      size_t const j = itB->row() ;
      size_t const k = itB->col() ;
      for( itA->start_row_items( k ) ; itA->is_valid() ; itA->go_next() )
      {
         size_t const i = itA->col() ;
         add_to_item( i, j, alpha*itA->item()*itB->item() ) ;
      }
   }
   stop_local_modifs() ;
   itA->destroy() ; itA = 0 ;
   itB->destroy() ; itB = 0 ;

   PEL_CHECK_POST( add_tMat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: has_same_pattern( LA_Matrix const* A ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqMatrix:: has_same_pattern" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( A->is_synchronized() ) ;
   PEL_CHECK_PRE( A->implementation() == implementation() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( A == this ) ;

   if( !result )
   {
      PEL_CHECK( dynamic_cast<LA_SeqMatrix const* >( A ) != 0 ) ;
      LA_SeqMatrix const* sA = static_cast<LA_SeqMatrix const* >( A ) ;

      result = ( A->nb_rows()==nb_rows() && A->nb_cols() == nb_cols() ) ;

      LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
      LA_MatrixIterator* itA = sA->create_stored_item_iterator( 0 ) ;
      for( it->start_all_items(), itA->start_all_items()  ;
           result && it->is_valid() && itA->is_valid() ;
           it->go_next(), itA->go_next() )
      {
         result &= it->row()==itA->row() && it->col()==itA->col() ;
      }
      result &= ( ! it->is_valid() ) &&  ( ! itA->is_valid() ) ;
      it->destroy()  ; it = 0 ;
      itA->destroy() ; itA = 0 ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Matrix::implementation_POST( result ) ) ;
   PEL_ASSERT( result == LA_SeqImplementation::object() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: re_initialize_with_global_sizes_PRE( size_t a_nb_rows,
                                                    size_t a_nb_cols ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_resizable() ) ;
   PEL_ASSERT( a_nb_rows != PEL::bad_index() ) ;
   PEL_ASSERT( a_nb_cols != PEL::bad_index() ) ;
   PEL_ASSERT( IMPLIES( is_symmetric(), a_nb_rows == a_nb_cols ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: re_initialize_with_global_sizes_POST( size_t a_nb_rows,
                                                     size_t a_nb_cols ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( nb_rows() == a_nb_rows ) ;
   PEL_ASSERT( nb_cols() == a_nb_cols ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        state() == LA::NotSync_undef ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        !is_synchronized() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: save_PRE( PEL_Module const* mod ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( mod != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: save_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: restore_PRE( PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqMatrix:: restore_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
LA_SeqMatrix:: sparse_mat_pluggin_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "LA_SeqMatrix descendant" ) ;
   return( result ) ;
}
