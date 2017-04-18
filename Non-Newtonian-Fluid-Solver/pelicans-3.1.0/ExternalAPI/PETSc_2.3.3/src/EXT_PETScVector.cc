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

#include <EXT_PETScVector.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>

#include <EXT_PETScImplementation.hh>

#include <iomanip>
#include <fstream>

//----------------------------------------------------------------------
EXT_PETScVector*
EXT_PETScVector:: create( PEL_Object* a_owner,
                          bool sequential,
                          size_t a_nb_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: create" ) ;

   EXT_PETScVector* result = new EXT_PETScVector( a_owner, sequential,
                                                  a_nb_rows ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;
   PEL_CHECK_POST( result->distribution_strategy() !=
                   LA::InvalidDistribution ) ;
   PEL_CHECK_POST( 
       EQUIVALENT( result->distribution_strategy() == LA::NoDistribution, 
                   sequential ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScVector:: EXT_PETScVector( PEL_Object* a_owner,
                                   bool sequential,
                                   size_t a_nb_rows )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, 0 )
   , VECTOR( 0 )
   , DIST( PEL_DistributedPartition::create( this ) )
   , SEQ( sequential )
{
   PEL_LABEL( "EXT_PETScVector:: EXT_PETScVector" ) ;

   if( sequential )
   {
      set_distribution_strategy( LA::NoDistribution ) ;
   }
   else
   {
      set_distribution_strategy( LA::FromGlobalSize) ;
      //la strategie LA::FromLocalSize est disponible dans PETSc
      //mais pas implemente ici
   }
   re_initialize( a_nb_rows ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( 
       EQUIVALENT( distribution_strategy() == LA::NoDistribution, 
                   sequential ) ) ;
   PEL_CHECK_POST( nb_rows() == a_nb_rows ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: re_initialize" ) ;
   PEL_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;

   int nloc = PETSC_DECIDE ;
   int nglob = a_nb_rows ;

   if( VECTOR ==0 || a_nb_rows != nb_rows() )
   {
      if( VECTOR!=0 ) PETSc_do( VecDestroy( VECTOR ) ) ;

      if( SEQ )
      {
         PETSc_do( VecCreateSeq( PETSC_COMM_SELF, a_nb_rows, &VECTOR ) ) ;
      }
      else
      {
         PETSc_do( VecCreateMPI( PETSC_COMM_WORLD, nloc, nglob, &VECTOR ) ) ;
      }
      int s ;
      PETSc_do( VecGetLocalSize( VECTOR, &s ) ) ;
      DIST->set_local_number( s ) ;
      int i0, i1 ;
      PETSc_do( VecGetOwnershipRange( VECTOR, &i0, &i1 ) ) ;
      FIRST = i0 ;
      LAST = i1 ;
      set_rows_number( a_nb_rows ) ;
   }
   set( 0.0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=row_distribution()->first_local_index() ;
                             i<row_distribution()->local_index_limit() ; ++i ),
                             item( i ) == 0.0 ) ) ;
}

//----------------------------------------------------------------------
EXT_PETScVector*
EXT_PETScVector:: create_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: create_vector" ) ;
   PEL_CHECK_PRE( create_vector_PRE( a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   EXT_PETScVector* result =
      new EXT_PETScVector( a_owner, SEQ, nb_rows() ) ;

   PEL_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScScatter*
EXT_PETScVector:: create_scatter( PEL_Object* a_owner,
                                  size_t_vector const& from,
                                  size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: create_scatter" ) ;
   PEL_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;

   PEL_CHECK_INV( invariant() ) ;

   EXT_PETScScatter* result =
      EXT_PETScScatter::create( a_owner, this, from, to ) ;

   PEL_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
EXT_PETScVector:: create_local_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: create_local_vector" ) ;

   PEL_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_rows() ) ;
   PetscScalar* values = 0 ;
   PETSc_do( VecGetArray( VECTOR, &values ) ) ;

   for( size_t i=FIRST ; i<LAST; i++ )
   {
      result->set_item(i, (double) values[i-FIRST]  ) ;
   }
   PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   result->synchronize() ;

   PEL_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
EXT_PETScVector:: item( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: item" ) ;
   PEL_CHECK_PRE( item_PRE( i ) ) ;

   PetscScalar values = PEL::bad_index() ;
   PetscInt ni = 1 ;
   PetscInt ix = i ;
   PETSc_do( VecGetValues( VECTOR, ni, &ix, &values ) ) ;

   return( values ) ;
}

//----------------------------------------------------------------------------
LA_Implementation const*
EXT_PETScVector:: implementation( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: implementation" ) ;

   LA_Implementation const* result = EXT_PETScImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScVector:: ~EXT_PETScVector( void )
//----------------------------------------------------------------------
{
   if( VECTOR!=0 ) PETSc_do( VecDestroy( VECTOR ) ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: is_desynchronizable" ) ;

   bool result = true ;

   PEL_CHECK_POST( result == true ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
EXT_PETScVector:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: row_distribution" ) ;
   PEL_CHECK_PRE( row_distribution_PRE() ) ;

   PEL_DistributedPartition const* result = DIST ;

   PEL_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: start_local_modifs" ) ;
   PEL_CHECK_PRE( start_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScVector:: start_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   PEL_Error::object()->raise_internal( mesg ) ;

   PEL_CHECK_POST( start_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: stop_local_modifs" ) ;
   PEL_CHECK_PRE( stop_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScVector:: stop_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   PEL_Error::object()->raise_internal( mesg ) ;

   PEL_CHECK_POST( stop_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: set( double value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: set double" ) ;
   PEL_CHECK_PRE( set_PRE( value ) ) ;

   PETSc_do( VecSet( vector(), value ) ) ;
   synchronize() ;

   PEL_CHECK_POST( set_POST( value ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: set LA_Vector" ) ;
   PEL_CHECK_PRE( set_PRE( a ) ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   PEL_CHECK_INV( invariant() ) ;

   PETSc_do( VecCopy( pa->vector(), vector() ) ) ;
   synchronize() ;

   PEL_CHECK_POST( set_POST( a ) ) ;
}


//----------------------------------------------------------------------
void
EXT_PETScVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PETSc_do( VecScale( VECTOR, alpha ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: sum LA_Vector" ) ;
   PEL_CHECK_PRE( sum_PRE( a, alpha ) ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;

   PETSc_do( VecAXPY( VECTOR, alpha, pa->vector() ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( sum_POST( a, alpha ) ) ;

}

//----------------------------------------------------------------------
double
EXT_PETScVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: dot LA_Vector" ) ;
   PEL_CHECK_PRE( dot_PRE( a ) ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   double result ;
   PETSc_do( VecDot( VECTOR, pa->vector(), &result ) ) ;

   PEL_CHECK_POST( dot_POST( result ) ) ;
   return result ;
}

//----------------------------------------------------------------------
double
EXT_PETScVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: two_norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( two_norm_PRE() ) ;

   double result ;
   PETSc_do( VecNorm( VECTOR, NORM_2, &result ) ) ;

   PEL_CHECK_POST( two_norm_POST( result ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
double
EXT_PETScVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   PEL_LABEL( "EXT_PETScVector:: max_norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( max_norm_PRE() ) ;

   double result ;
   PETSc_do( VecNorm( VECTOR, NORM_INFINITY, &result ) ) ;

   PEL_CHECK_POST( max_norm_POST( result ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: synchronize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: synchronize" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PETSc_do( VecAssemblyBegin( VECTOR ) ) ;
   PETSc_do( VecAssemblyEnd( VECTOR ) ) ;

   LA_Vector::synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
}


//----------------------------------------------------------------------
Vec&
EXT_PETScVector:: vector( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: vector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( VECTOR!=0 ) ;

   return const_cast<EXT_PETScVector*>(this)->VECTOR ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: print_items" ) ;
   PEL_CHECK_PRE( print_items_PRE( os, indent_width ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;

   std::string const space( indent_width, ' ' ) ;
   os << space << "nb_rows:" << nb_rows() << std::endl ;
   if( nb_rows()>0 )
   {
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      std::streamsize p = os.precision() ;
      os << std::setprecision( 7 ) ;
      double * values ;
      PETSc_do( VecGetArray( VECTOR, &values ) ) ;
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = values[iRow-FIRST] ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
      PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   }
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: set_as_reciprocal( LA_Vector const* a,
                                     double smallest_inverted_item,
                                     double default_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: set_as_reciprocal LA_Vector" ) ;
   PEL_CHECK_PRE( set_as_reciprocal_PRE(a,smallest_inverted_item,default_value) ) ;

   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;

   PETSc_do( VecCopy( pa->vector(), VECTOR ) ) ;

   if( smallest_inverted_item!=0.0 )
   {
      double * values ;
      PETSc_do( VecGetArray( pa->vector(), &values ) ) ;
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         if( PEL::abs( values[iRow-FIRST] )<smallest_inverted_item )
         {
            VecSetValue( VECTOR, iRow, 1.0/default_value, INSERT_VALUES ) ;
         }
      }
      PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   }

   PETSc_do( VecReciprocal( VECTOR ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_as_reciprocal_POST(a) ) ;

}

//----------------------------------------------------------------------
void
EXT_PETScVector:: set_as_v_product( LA_Vector const* a,
                                         LA_Vector const* b )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: set_as_v_product LA_Vector" ) ;
   PEL_CHECK_PRE( set_as_v_product_PRE(a,b) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( b ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   EXT_PETScVector const* pb = static_cast<EXT_PETScVector const*>( b ) ;

   PETSc_do( VecPointwiseMult( VECTOR, pa->vector(), pb->vector() ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: set_item( size_t i, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: set_item" ) ;
   PEL_CHECK_PRE( set_item_PRE( i ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }
   VecSetValue( VECTOR, i, x, INSERT_VALUES ) ;

   PEL_CHECK_POST( set_item_POST( i, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: add_to_item( size_t i, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: add_to_item" ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }
   if( x != 0.0 )
   {
      VecSetValue( VECTOR, i, x, ADD_VALUES ) ;
   }
   PEL_CHECK_POST( add_to_item_POST( i, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::invariant() ) ;
   PEL_ASSERT( distribution_strategy() != LA::InvalidDistribution ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   PEL_ASSERT( result == EXT_PETScImplementation::object() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScVector:: write" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;

   size_t rank = PEL_Exec::communicator()->rank() ;
   size_t size = PEL_Exec::communicator()->nb_ranks() ;
   int dummy ;

   if( rank>0 ) PEL_Exec::communicator()->receive( rank-1, dummy ) ;

   std::ofstream out( filename.c_str(), (rank==0 ? std::ios::out : std::ios::out|std::ios::app) ) ;
   if( rank==0 )
      out << nb_rows() ;

   double * values ;
   PETSc_do( VecGetArray( VECTOR, &values ) ) ;
   for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
   {
      out << std::endl << values[iRow-FIRST]  ;
   }
   PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   out.close() ;

   if( rank!=size-1 )
   {
      PEL_Exec::communicator()->send( rank+1, dummy ) ;
      PEL_Exec::communicator()->receive( size-1, dummy ) ;
   }
   else
   {
      for(size_t i=0 ; i<size-1 ; i++ )
         PEL_Exec::communicator()->send( i, dummy ) ;
   }
}
