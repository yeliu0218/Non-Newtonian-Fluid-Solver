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

#include <LA_DistVector.hh>

#include <LA_DistImplementation.hh>

#include <iomanip>
#include <fstream>
#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Int.hh>
#include <PEL_Vector.hh>

#ifdef OUTLINE
   #define inline
   #include <LA_DistVector.icc>
   #undef inline
#endif

//----------------------------------------------------------------------
LA_DistVector*
LA_DistVector:: create( PEL_Object* a_owner,
                        size_t a_nb_rows, size_t a_nb_local_rows,
                        LA::DistributionStrategy dist_strat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: create" ) ;
   PEL_CHECK_PRE( dist_strat == LA::FromLocalSize ||
                  dist_strat == LA::FromGlobalSize ) ;

   LA_DistVector* result = new LA_DistVector( a_owner,
                                              a_nb_rows, a_nb_local_rows,
                                              dist_strat ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST(
       FORALL( ( size_t i=result->row_distribution()->first_local_index() ;
                 i<result->row_distribution()->local_index_limit() ; ++i ),
                 result->item(i) == 0.0 ) ) ;
   PEL_CHECK_POST( result->distribution_strategy() == dist_strat ) ;
   PEL_CHECK_POST( IMPLIES( result->distribution_strategy() == LA::FromGlobalSize,
                            result->nb_rows() == a_nb_rows  ) ) ;
   PEL_CHECK_POST( IMPLIES( result->distribution_strategy() == LA::FromLocalSize,
                            result->row_distribution()->local_index_limit()
                          - result->row_distribution()->first_local_index()
                          == a_nb_local_rows ) ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistVector:: LA_DistVector( PEL_Object* a_owner,
                               size_t a_nb_rows,
                               size_t a_nb_local_rows,
                               LA::DistributionStrategy dist_strat )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, 0 )
   , DIST( PEL_DistributedPartition::create( this ) )
   , GLOBAL_INDEXES( PEL_Exec::communicator()->nb_ranks(), intVector( 0 ) )
   , GLOBAL_VALUES( PEL_Exec::communicator()->nb_ranks(), doubleVector( 0 ) )
   , LOCAL_VECTOR( LA_SeqVector::create( this, 0 ) )
   , VERBOSE( false )
   , INITIALIZED( false )
   , FIRST( 0 )
   , LAST( 0 )
   , NB_RANKS( PEL_Exec::communicator()->nb_ranks() )
   , HAS_BEEN_NULLIFIED( true )
   , NB_EXTRA(0)
{
   PEL_LABEL( "LA_DistVector:: LA_DistVector" ) ;
   PEL_CHECK_PRE( dist_strat == LA::FromLocalSize ||
                  dist_strat == LA::FromGlobalSize ) ;

   set_distribution_strategy( dist_strat ) ;

   re_initialize( a_nb_rows, a_nb_local_rows ) ;
   INITIALIZED = true ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( distribution_strategy() == dist_strat ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: re_initialize" ) ;
   PEL_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( distribution_strategy() == LA::FromLocalSize ||
       DIST->global_number() != a_nb_rows )
   {
      if( distribution_strategy() == LA::FromLocalSize )
      {
         DIST->set_local_number( a_nb_local_rows ) ;
      }
      else if( distribution_strategy() == LA::FromGlobalSize )
      {
         DIST->distribute_global_number( a_nb_rows ) ;
      }
      LOCAL_VECTOR->re_initialize( DIST->local_number() ) ;
      FIRST = DIST->first_local_index() ;
      LAST  = DIST->local_index_limit() ;
      set_rows_number( DIST->global_number() ) ;
   }
   set( 0.0 ) ;

   PEL_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
}

//----------------------------------------------------------------------
LA_DistVector*
LA_DistVector:: create_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: create_vector" ) ;
   PEL_CHECK_PRE( create_vector_PRE( a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_DistVector* result =
      new LA_DistVector( a_owner,
                         DIST->global_number(),
                         DIST->local_number(),
                         distribution_strategy() ) ;

   PEL_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistScatter*
LA_DistVector:: create_scatter( PEL_Object* a_owner,
                                size_t_vector const& from,
                                size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: create_scatter" ) ;
   PEL_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;

   PEL_CHECK_INV( invariant() ) ;

   LA_DistScatter* result =
      LA_DistScatter::create( a_owner, row_distribution() ) ;
   result->set_unsorted( from, to ) ;

   PEL_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DistVector:: create_local_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: create_local_vector" ) ;

   PEL_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_rows() ) ;

   for( size_t i=FIRST ; i<LAST; i++ )
   {
      result->set_item( i, LOCAL_VECTOR->item( i-FIRST ) ) ;
   }
   result->synchronize() ;

   PEL_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistVector:: ~LA_DistVector( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_DistVector:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: implementation" ) ;

   LA_Implementation const* result = LA_DistImplementation::object() ;


   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_synchronized( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set_synchronized" ) ;
   PEL_CHECK_INV( invariant() ) ;
   NB_EXTRA=0 ;
   for( size_t i=0 ; i<NB_RANKS ; i++ )
   {
      GLOBAL_VALUES[i].re_initialize( 0 ) ;
      GLOBAL_INDEXES[i].re_initialize( 0 ) ;
   }
   LA_Vector::synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
}

//----------------------------------------------------------------------
bool
LA_DistVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: is_desynchronizable" ) ;

   bool result = true ;

   PEL_CHECK_POST( result == true ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( double value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set LA_Vector" ) ;
   PEL_CHECK_PRE( set_PRE( value ) ) ;

   LOCAL_VECTOR->set( value ) ;
   set_synchronized() ;

   PEL_CHECK_POST( set_POST( value ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set LA_Vector" ) ;
   PEL_CHECK_PRE( set_PRE( a ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->set( da->LOCAL_VECTOR ) ;
   set_synchronized() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( a ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( LA_SeqVector const* a,
                     size_t_vector const& local_to_global )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set" ) ;
   PEL_CHECK_PRE( a->nb_rows() == local_to_global.size() ) ;
   PEL_CHECK_PRE( state() != LA::NotSync_add ) ;
   PEL_CHECK_INV( invariant() ) ;

   set_unsynchronized_state( LA::NotSync_set ) ;
   for( size_t i=0 ; i<local_to_global.size() ; i++ )
   {
      set_item( local_to_global(i), a->item(i) ) ;
   }

   PEL_CHECK_POST( state() == LA::NotSync_set ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: add( LA_SeqVector const* a,
                     size_t_vector const& local_to_global )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: add" ) ;
   PEL_CHECK_PRE( a->nb_rows() == local_to_global.size() ) ;
   PEL_CHECK_PRE( state()!=LA::NotSync_set ) ;
   PEL_CHECK_INV( invariant() ) ;

   set_unsynchronized_state( LA::NotSync_add ) ;
   for( size_t i=0 ; i<local_to_global.size() ; i++ )
      if( a->item(i)!=0.0 )
         add_to_item( local_to_global(i), a->item(i) ) ;

   PEL_CHECK_POST( state() == LA::NotSync_add ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LOCAL_VECTOR->scale( alpha ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: sum LA_Vector" ) ;
   PEL_CHECK_PRE( sum_PRE( a, alpha ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->sum( da->LOCAL_VECTOR, alpha ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( sum_POST( a, alpha ) ) ;

}

//----------------------------------------------------------------------
double
LA_DistVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: dot LA_Vector" ) ;
   PEL_CHECK_PRE( dot_PRE( a ) ) ;

   PEL_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;

   double result =
      ( DIST->local_number() > 0 ?
                      LOCAL_VECTOR->dot( da->LOCAL_VECTOR ) : 0.0 ) ;
   result = PEL_Exec::communicator()->sum( result ) ;

   PEL_CHECK_POST( dot_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_DistVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: two_norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( two_norm_PRE() ) ;

   double result =
      ( DIST->local_number() > 0 ? LOCAL_VECTOR->two_norm() : 0.0 ) ;
   result = PEL::sqrt( PEL_Exec::communicator()->sum( PEL::sqr( result ) ) ) ;

   PEL_CHECK_POST( two_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_DistVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   PEL_LABEL( "LA_DistVector:: max_norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( max_norm_PRE() ) ;

   double result =
      ( DIST->local_number() > 0 ?
                    LOCAL_VECTOR->max_norm() : -PEL::max_double() ) ;
   result = PEL_Exec::communicator()->max( result ) ;

   PEL_CHECK_POST( max_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: recover_global_vector( LA_SeqVector* vec ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: recover_global_vector" ) ;
   PEL_CHECK_PRE( vec!=0 ) ;
   PEL_CHECK_PRE( vec->nb_rows()==nb_rows() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t dec = FIRST ;
   size_t first = FIRST ;
   size_t last = LAST ;

   for( size_t i=first ; i<last; i++ )
   {
      vec->set_item( i, LOCAL_VECTOR->item( i - dec ) ) ;
   }

   PEL_Exec::communicator()->all_gather_v( LOCAL_VECTOR->data(),
                                           DIST->local_number(),
                                           const_cast< double* >( vec->data() ),
                                           DIST->partitioning(),
                                           DIST->start_of_partition() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: synchronize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: synchronize" ) ;
   PEL_CHECK_PRE( synchronize_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK(
      PEL_Exec::communicator()->boolean_and( state()!=LA::NotSync_add ) ||
      PEL_Exec::communicator()->boolean_and( state()!=LA::NotSync_set ) ) ;

   LA::SyncState effective_mode =
      ( PEL_Exec::communicator()->boolean_and( state() != LA::NotSync_set ) ?
        LA::NotSync_add : LA::NotSync_set ) ;

   if( VERBOSE )
   {
      PEL::out() << "LA_DistVector:: synchronize" << std::endl ;
   }

   size_t const nb_ranks = PEL_Exec::communicator()->nb_ranks() ;
   size_t const rank = PEL_Exec::communicator()->rank() ;

   // C'est la que les choses serieuses vont commencer ...

   // Tout d'abords, on va commencer par s'envoyer le nombre d'elements
   //  a echanger entre process
   intVector exported_items( nb_ranks ) ;
   for( size_t i=0 ; i<nb_ranks ; ++i )
   {
      if( i != rank )
      {
         exported_items( i ) = GLOBAL_VALUES[i].size() ;
      }
   }

   intVector exchanged_items( nb_ranks ) ;
   PEL_Exec::communicator()->all_to_all( exported_items, exchanged_items ) ;

   // Maintenant, on va les envoyer et recevoir
   for( size_t i=0 ; i<rank ; ++i )
   {
      size_t N = (size_t) exchanged_items( i ) ;
      if( N > 0 )
      {
         receive_subvector( N, i, effective_mode ) ;
      }
   }

   for( size_t i=0 ; i<nb_ranks ; ++i )
   {
      int N = exported_items( i ) ;
      if( N > 0 )
      {
         send_subvector( N, i ) ;
      }

   }

   for( size_t i=rank+1 ; i<nb_ranks ; ++i )
   {
      size_t N = (size_t) exchanged_items( i ) ;
      if( N > 0 )
      {
         receive_subvector( N, i, effective_mode ) ;
      }
   }

   set_synchronized() ;

   HAS_BEEN_NULLIFIED = false ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
   PEL_CHECK_POST( local_vector()->is_synchronized() ) ;
}


//----------------------------------------------------------------------
void
LA_DistVector:: send_subvector( size_t N, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: send_subvector" ) ;
   PEL_CHECK( N > 0 ) ;
   PEL_CHECK( i < PEL_Exec::communicator()->nb_ranks() ) ;
   PEL_CHECK( i != PEL_Exec::communicator()->rank() ) ;

   if( VERBOSE )
   {
      PEL::out() << "Process " << PEL_Exec::communicator()->rank()
                 << " send " << N << " values to proccess " << i ;
   }

   PEL_Exec::communicator()->send( i, GLOBAL_INDEXES[i].data(), N ) ;
   PEL_Exec::communicator()->send( i, GLOBAL_VALUES[i].data(), N ) ;

   if( VERBOSE )
   {
      PEL::out() << " ... OK " << std::endl ;
      PEL::out().flush() ;
   }
}

//----------------------------------------------------------------------
void
LA_DistVector:: receive_subvector( size_t N, size_t i,
                                   LA::SyncState effective_mode )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: receive_subvector" ) ;
   PEL_CHECK( N > 0 ) ;
   PEL_CHECK( i < PEL_Exec::communicator()->nb_ranks() ) ;
   PEL_CHECK( i != PEL_Exec::communicator()->rank() ) ;

   doubleVector values( N ) ;
   double* ptr_values = const_cast<double*>( values.data() ) ;

   intVector idx( N ) ;
   int* ptr_idx = const_cast<int*>( idx.data() ) ;

   double* ptr_vec = LOCAL_VECTOR->data() ;

   if( VERBOSE )
   {
      PEL::out() << "Process " << PEL_Exec::communicator()->rank()
                 << " receive " << N << " values from " << i ;
   }

   PEL_Exec::communicator()->receive( i, ptr_idx, N ) ;
   PEL_Exec::communicator()->receive( i, ptr_values, N ) ;

   if( effective_mode == LA::NotSync_add )
   {
      for( size_t j=0 ; j<N ; j++ )
      {
         PEL_CHECK( idx(j)>=(int)FIRST && idx(j)<(int)LAST ) ;
         ptr_vec[ptr_idx[j]-FIRST] += ptr_values[j] ;
      }

   }
   else
   {
      for( size_t j=0 ; j<N ; j++ )
      {
         PEL_CHECK( idx(j)>=(int)FIRST && idx(j)<(int)LAST ) ;
         PEL_CHECK(
            IMPLIES( HAS_BEEN_NULLIFIED ,
                     ( LOCAL_VECTOR->item(idx(j)-FIRST)==0.0
                       || LOCAL_VECTOR->item(idx(j)-FIRST)==values(j) ) ) ) ;
         ptr_vec[ptr_idx[j]-FIRST] = ptr_values[j] ;
      }
   }
   if( VERBOSE )
   {
      PEL::out() << " ... OK " << std::endl ;
      PEL::out().flush() ;
   }
}

//----------------------------------------------------------------------
LA_SeqVector const*
LA_DistVector:: local_vector( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: local_vector const" ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = LOCAL_VECTOR ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->nb_rows()==row_distribution()->local_number() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DistVector:: local_vector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: local_vector" ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = LOCAL_VECTOR ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->nb_rows()==row_distribution()->local_number() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: print_items" ) ;
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
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = item(iRow) ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
   }
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_as_reciprocal( LA_Vector const* a,
                                   double smallest_inverted_item,
                                   double default_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set_as_reciprocal LA_Vector" ) ;
   PEL_CHECK_PRE( set_as_reciprocal_PRE(
                         a, smallest_inverted_item, default_value ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->set_as_reciprocal(
           da->LOCAL_VECTOR, smallest_inverted_item, default_value ) ;
   set_synchronized() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_as_reciprocal_POST(a) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_as_v_product( LA_Vector const* a, LA_Vector const* b )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set_as_v_product LA_Vector" ) ;
   PEL_CHECK_PRE( set_as_v_product_PRE(a,b) ) ;

   PEL_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   PEL_CHECK( dynamic_cast<LA_DistVector const*>( b ) !=0 ) ;
   LA_DistVector const* db = static_cast<LA_DistVector const*>( b ) ;

   LOCAL_VECTOR->set_as_v_product( da->LOCAL_VECTOR, db->LOCAL_VECTOR ) ;
   set_synchronized() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_item( size_t i, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: set_item" ) ;
   PEL_CHECK_PRE( set_item_PRE( i ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }

   if( i<FIRST || i>=LAST )
   {
      size_t r = DIST->rank_of( i ) ;
      PEL_CHECK( r != PEL_Exec::communicator()->rank() ) ;
      intVector& idx = GLOBAL_INDEXES[r] ;
      doubleVector& val = GLOBAL_VALUES[r] ;
      bool found = false ;
      size_t n = idx.size() ;
      for( size_t j=0 ; j<n ; j++ )
      {
         if( idx(j)==(int)i )
         {
            val(j) = x ;
            found = true ;
            break ;
         }
      }
      if( !found )
      {
         idx.append( i ) ;
         val.append( x ) ;
         NB_EXTRA++ ;
      }
   }
   else
   {
      LOCAL_VECTOR->data()[i-FIRST] = x ;
   }
   PEL_CHECK_POST( set_item_POST( i, OLD( state ) ) ) ;

}

//----------------------------------------------------------------------
void
LA_DistVector:: add_to_item( size_t i, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: add_to_item" ) ;
   PEL_CHECK_PRE( add_to_item_PRE( i ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }

   if( x != 0.0 )
   {
      if( i<FIRST || i>=LAST )
      {
         size_t r = DIST->rank_of( i ) ;
         PEL_CHECK( r != PEL_Exec::communicator()->rank() ) ;
         intVector& idx = GLOBAL_INDEXES[r] ;
         doubleVector& val = GLOBAL_VALUES[r] ;
         bool found = false ;
         size_t n = idx.size() ;
         for( size_t j=0 ; j<n ; j++ )
         {
            if( idx(j)==(int)i )
            {
               val(j) += x ;
               found = true ;
               break ;
            }
         }
         if( !found )
         {
            idx.append( i ) ;
            val.append( x ) ;
            NB_EXTRA++ ;
         }
      }
      else
      {
         LOCAL_VECTOR->data()[i-FIRST] += x ;
      }
   }

   PEL_CHECK_POST( add_to_item_POST( i, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_DistVector:: write" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;

   size_t rank = PEL_Exec::communicator()->rank() ;
   size_t size = PEL_Exec::communicator()->nb_ranks() ;
   int dummy ;

   if( rank>0 ) PEL_Exec::communicator()->receive( rank-1, dummy ) ;

   std::ofstream out( filename.c_str(),
                      ( rank==0 ? std::ios::out
                                : std::ios::out|std::ios::app ) ) ;
   if( rank==0 ) out << nb_rows() ;

   for( size_t i=0 ; i<LOCAL_VECTOR->nb_rows() ; i++ )
   {
      out << std::endl << LOCAL_VECTOR->item(i) ;
   }

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

//----------------------------------------------------------------------
bool
LA_DistVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::invariant() ) ;
   if( INITIALIZED )
   {
      PEL_ASSERT( distribution_strategy() == LA::FromLocalSize ||
                  distribution_strategy() == LA::FromGlobalSize ) ;
      PEL_ASSERT( PEL_Exec::communicator()!=0 ) ;
      PEL_ASSERT( DIST!=0 ) ;
      PEL_ASSERT( LOCAL_VECTOR!=0 ) ;
      PEL_ASSERT( LOCAL_VECTOR->distribution_strategy() == 
                  LA::NoDistribution ) ;
      PEL_ASSERT( GLOBAL_VALUES.size() ==
                             PEL_Exec::communicator()->nb_ranks() ) ;
      PEL_ASSERT( GLOBAL_INDEXES.size() ==
                             PEL_Exec::communicator()->nb_ranks() ) ;
      PEL_ASSERT( DIST->local_number()==LOCAL_VECTOR->nb_rows() ) ;
      PEL_ASSERT( ( state() != LA::Sync ) ||
                  FORALL( (size_t i=0 ; i<GLOBAL_VALUES.size() ; i++),
                          GLOBAL_VALUES[i].size()==0 ) ) ;
      PEL_ASSERT( FORALL( (size_t i=0 ; i<GLOBAL_VALUES.size() ; i++),
                          GLOBAL_VALUES[i].size()==GLOBAL_INDEXES[i].size() ) ) ;
      PEL_ASSERT( FIRST == DIST->first_local_index() ) ;
      PEL_ASSERT( LAST == DIST->local_index_limit() ) ;
   }
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_DistVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   PEL_ASSERT( result == LA_DistImplementation::object() ) ;
   return( true ) ;
}
