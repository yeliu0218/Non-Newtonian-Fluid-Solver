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

#include <LA_Vector.hh>

#include <LA_SeqVector.hh>
#include <LA_Scatter.hh>

#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <size_t_vector.hh>

#ifdef OUTLINE
#define inline
   #include <LA_Vector.icc>
#undef inline
#endif

#include <iostream>
#include <fstream>

//----------------------------------------------------------------------
LA_Vector:: LA_Vector( PEL_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_ROWS( a_nb_rows )
   , DIST_STATUS( LA::Sync )
   , RESIZABLE( true )
   , ONLY_LOCAL_MODIFS( false )
   , DIST_STRAT( LA::InvalidDistribution )
{
   PEL_LABEL( "LA_Vector:: LA_Vector" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( nb_rows() == a_nb_rows ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
LA_Vector:: ~LA_Vector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: ~LA_Vector" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: is_resizable( void ) const
//----------------------------------------------------------------------
{
   return( RESIZABLE ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: make_non_resizable( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: make_non_resizable" ) ;

   RESIZABLE = false ;

   PEL_CHECK_POST( !is_resizable() ) ;
}

//----------------------------------------------------------------------
size_t
LA_Vector:: nb_local_rows( void ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: nb_local_rows" ) ;
   
   size_t result = row_distribution()->local_number() ;
   
   PEL_CHECK_POST( result == row_distribution()->local_number() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: only_local_modifs( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: only_local_modifs" ) ;

   bool result = ONLY_LOCAL_MODIFS ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: set_only_local_modifs_state( bool only_local )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: set_only_local_modifs_state" ) ;

   ONLY_LOCAL_MODIFS = only_local ;

   PEL_CHECK_POST( only_local_modifs() == only_local ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: start_local_modifs" ) ;
   PEL_CHECK_PRE( start_local_modifs_PRE() ) ;

   set_only_local_modifs_state( true ) ;

   PEL_CHECK_POST( start_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: stop_local_modifs" ) ;
   PEL_CHECK_PRE( stop_local_modifs_PRE() ) ;

   set_only_local_modifs_state( false ) ;

   PEL_CHECK_POST( stop_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: set_distribution_strategy(
                              LA::DistributionStrategy dist_strat ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: set_distribution_strategy" ) ;
   PEL_CHECK_PRE( distribution_strategy() == LA::InvalidDistribution ) ;
   PEL_CHECK_PRE( dist_strat != LA::InvalidDistribution ) ;

   DIST_STRAT = dist_strat ;

   PEL_CHECK_POST( distribution_strategy() == dist_strat ) ;
}

//----------------------------------------------------------------------
LA::DistributionStrategy
LA_Vector:: distribution_strategy( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: distribution_strategy" ) ;

   return( DIST_STRAT ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: is_synchronized( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: is_synchronized" ) ;
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;

   bool result = ( is_desynchronizable() ?
                      PEL_Exec::communicator()->boolean_and(
                                                        state() == LA::Sync ) :
                      true ) ;

   PEL_CHECK_POST(
      IMPLIES( !is_desynchronizable(), result == true ) ) ;
   PEL_CHECK_POST(
      IMPLIES( is_desynchronizable(),
               result == PEL_Exec::communicator()->boolean_and(
                                               state() == LA::Sync ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: nullify( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: nullify" ) ;
   PEL_CHECK_PRE( nullify_PRE() ) ;

   set( 0.0 ) ;

   PEL_CHECK_POST( nullify_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "vector: \"" << type_name() << "\"" << std::endl ;
}
//----------------------------------------------------------------------
void
LA_Vector:: read( std::string const& filename )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: read" ) ;
   PEL_CHECK_PRE( read_PRE( filename ) ) ;

   if( !is_resizable() )
   {
      PEL_Error::object()->raise_internal(
         "Reading of non resizable vector is not implemented" ) ;
   }

   size_t rank = PEL_Exec::communicator()->rank() ;

   std::ifstream in( filename.c_str() ) ;
   if( !in )
   {
      PEL_Error::object()->raise_file_handling( filename, "open" ) ;
   }
   size_t nbrows ;

   in >> nbrows ;
   re_initialize( nbrows ) ;
   size_t i = 0 ;

   if( rank==0 || !is_desynchronizable() )
   {
      while( !in.eof() )
      {
         if( i>=nbrows )
         {
            PEL_Error::object()->raise_plain( "Invalid vector (rows>nb_rows) in "
                                              +filename ) ;
         }
         double v ;
         in >> v ;
         set_item( i, v ) ;
         i++ ;
      }
      if( i!=nbrows )
      {
         PEL_Error::object()->raise_plain( "Invalid vector (rows<nb_rows) in "
                                           +filename ) ;
      }
   }
   in.close() ;
   synchronize() ;

   PEL_CHECK_POST( read_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Vector:: set_rows_number( size_t a_nb_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Vector:: set_rows_number" ) ;
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;

   NB_ROWS = a_nb_rows ;

   PEL_CHECK_POST( nb_rows() == a_nb_rows ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: re_initialize_PRE( size_t a_nb_rows,
                               size_t a_nb_local_rows ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_resizable() ) ;

   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        a_nb_rows != PEL::bad_index() ) ) ;
   PEL_ASSERT(
     IMPLIES( distribution_strategy() != LA::FromLocalSize && is_desynchronizable(),
       PEL_Exec::communicator()->same_value_everywhere( (int) a_nb_rows ) ) ) ;

   PEL_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                           a_nb_local_rows != PEL::bad_index() ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: re_initialize_POST( size_t a_nb_rows,
                                size_t a_nb_local_rows ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        nb_rows() == a_nb_rows ) ) ;
   PEL_ASSERT(
       IMPLIES( distribution_strategy() == LA::FromLocalSize,
                row_distribution()->local_number() ==  a_nb_local_rows ) ) ;

   PEL_ASSERT( state() == LA::Sync ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_vector_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_vector_POST(
                          LA_Vector* result, PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->implementation() == implementation() ) ;
   PEL_ASSERT( result->is_resizable() == is_resizable() ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( result->is_desynchronizable() == is_desynchronizable() ) ;
   PEL_ASSERT( result->state() == LA::Sync ) ;
   PEL_ASSERT( result->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( result->is_synchronized() ) ;
   PEL_ASSERT( result->distribution_strategy() == distribution_strategy() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: start_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: start_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: stop_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: stop_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: synchronize_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: synchronize_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: row_distribution_PRE( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: row_distribution_POST(
                          PEL_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( result->partitioning().size()==
                               PEL_Exec::communicator()->nb_ranks() ) ;
   PEL_ASSERT( result->global_number()==nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_scatter_PRE( PEL_Object* a_owner,
                                size_t_vector const& from,
                                size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( from.size() == to.size() ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<from.size() ; ++i ),
                       from(i) < nb_rows() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_scatter_POST( LA_Scatter const* result,
                                 PEL_Object const* a_owner,
                                 size_t_vector const& from,
                                 size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->implementation() == implementation() ) ;
   PEL_ASSERT( result->size() == from.size() ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<result->size() ; ++i ),
              result->repatriated_items().has( from(i) ) ) ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<result->size() ; ++i ),
              result->local_indices()(i) == to(
                 from.index_of( result->repatriated_items()(i) ) ) ) ) ;
   PEL_ASSERT( result->distribution()->is_compatible(
                                                row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_local_vector_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: create_local_vector_POST(
                       LA_SeqVector* result, PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( result->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < row_distribution()->local_index_limit()  &&
               i >= row_distribution()->first_local_index() ) ;
   PEL_ASSERT( state() == LA::Sync ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( state() != LA::NotSync_add ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_item_POST( size_t i, LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
   IMPLIES( is_desynchronizable() && !only_local_modifs(),
            state() == LA::NotSync_set ) ) ;
   PEL_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: add_to_item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: add_to_item_POST( size_t i, LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      IMPLIES( is_desynchronizable() && !only_local_modifs(),
               state() == LA::NotSync_add ) ) ;
   PEL_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: dot_PRE( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( a != 0 ) ;
   PEL_ASSERT( a->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( a->implementation() == implementation() ) ;
   PEL_ASSERT( a->is_synchronized() ) ;
   PEL_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: dot_POST( double result ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: two_norm_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( nb_rows() > 0 ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: two_norm_POST( double result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result >= 0.0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: max_norm_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( nb_rows() > 0 ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: max_norm_POST( double result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result >= 0.0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: nullify_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: nullify_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: scale_PRE( double value ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( value ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: scale_POST( double value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_PRE( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( a != 0 ) ;
   PEL_ASSERT( a->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( a->implementation() == implementation() ) ;
   PEL_ASSERT( a->is_synchronized() ) ;
   PEL_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_POST( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_PRE( double value ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( value ) ) )  ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_POST( double value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_as_v_product_PRE(
                          LA_Vector const* a, LA_Vector const* b ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( a != 0 ) ;
   PEL_ASSERT( a->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( a->implementation() == implementation() ) ;
   PEL_ASSERT( a->is_synchronized() ) ;
   PEL_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( b != 0 ) ;
   PEL_ASSERT( b->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( b->implementation() == implementation() ) ;
   PEL_ASSERT( b->is_synchronized() ) ;
   PEL_ASSERT( b->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_as_v_product_POST(
                          LA_Vector const* a, LA_Vector const* b ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( FORMAL( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                               item(i) == a->item(i)*b->item(i) ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_as_reciprocal_PRE( LA_Vector const* a,
                                   double smallest_inverted_item,
                                   double default_value ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( a != 0 ) ;
   PEL_ASSERT( a->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( a->implementation() == implementation() ) ;
   PEL_ASSERT( a->is_synchronized() ) ;
   PEL_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( smallest_inverted_item >= 0. ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( smallest_inverted_item ) ) ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( default_value ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: set_as_reciprocal_POST( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: sum_PRE( LA_Vector const* a, double alpha ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( a != 0 ) ;
   PEL_ASSERT( a->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( a->implementation() == implementation() ) ;
   PEL_ASSERT( a->is_synchronized() ) ;
   PEL_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: sum_POST(  LA_Vector const* a, double alpha ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: print_items_PRE( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: write_PRE( std::string const& filename ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( !filename.empty() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: write_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: read_PRE( std::string const& filename ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( !filename.empty() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Vector:: read_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}
