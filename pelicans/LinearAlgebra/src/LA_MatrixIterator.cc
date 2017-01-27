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

#include <LA_MatrixIterator.hh>

#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>

#include <LA.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_MatrixIterator:: LA_MatrixIterator( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//----------------------------------------------------------------------
LA_MatrixIterator:: ~LA_MatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_MatrixIterator:: unsynchronized_matrix_state_for_set( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_MatrixIterator:: unsynchronized_matrix_state_for_set" ) ;
   PEL_CHECK_PRE( matrix()->is_desynchronizable() ) ;
   PEL_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   PEL_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   PEL_CHECK_PRE( matrix()->state() != LA::NotSync_add ) ;

   LA_Matrix* m = const_cast<LA_Matrix*>( matrix() ) ;
   m->set_unsynchronized_state( LA::NotSync_set ) ;

   PEL_CHECK_POST( matrix()->state() == LA::NotSync_set ) ;
}

//----------------------------------------------------------------------
void
LA_MatrixIterator:: unsynchronized_matrix_state_for_add( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_MatrixIterator:: unsynchronized_matrix_state_for_add" ) ;
   PEL_CHECK_PRE( matrix()->is_desynchronizable() ) ;
   PEL_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   PEL_CHECK_PRE( matrix()->state() != LA::NotSync_set ) ;

   LA_Matrix* m = const_cast<LA_Matrix*>( matrix() ) ;
   m->set_unsynchronized_state( LA::NotSync_add ) ;

   PEL_CHECK_POST( matrix()->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: start_all_items_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( matrix()->state() ==  LA::Sync ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: start_row_items_PRE( size_t i_row ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( matrix()->state() == LA::Sync ) ;
   PEL_ASSERT( i_row < nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: go_next_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: matrix_POST( LA_Matrix const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: row_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: row_POST( size_t result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result < matrix()->nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: col_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: col_POST( size_t result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result < matrix()->nb_cols() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
double
LA_MatrixIterator:: item_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: set_item_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT(
       IMPLIES( matrix()->only_local_modifs(),
                row() >= matrix()->row_distribution()->first_local_index() &&
                row() < matrix()->row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( matrix()->state() != LA::NotSync_add ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: set_item_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
   IMPLIES( matrix()->is_desynchronizable() && !matrix()->only_local_modifs(),
            matrix()->state() == LA::NotSync_set ) ) ;
   PEL_ASSERT(
   IMPLIES( !matrix()->is_desynchronizable() || matrix()->only_local_modifs(),
            matrix()->state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: add_to_item_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT(
       IMPLIES( matrix()->only_local_modifs(),
                row() >= matrix()->row_distribution()->first_local_index() &&
                row() < matrix()->row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( matrix()->state() != LA::NotSync_set ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: add_to_item_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      IMPLIES( matrix()->is_desynchronizable() && !matrix()->only_local_modifs(),
               matrix()->state() == LA::NotSync_add ) ) ;
   PEL_ASSERT(
   IMPLIES( !matrix()->is_desynchronizable() || matrix()->only_local_modifs(),
            matrix()->state() == old_state ) ) ;
   return( true ) ;
}



