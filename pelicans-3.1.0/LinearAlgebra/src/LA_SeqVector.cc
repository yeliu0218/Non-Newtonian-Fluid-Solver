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

#include <LA_SeqVector.hh>

#include <LA_SeqImplementation.hh>
#include <LA_SeqScatter.hh>

#include <PEL.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Int.hh>
#include <PEL_Vector.hh>

#include <fstream>
#include <iomanip>

#ifdef OUTLINE
#define inline
#include <LA_SeqVector.icc>
#undef inline
#endif

//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create( PEL_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: create" ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, a_nb_rows ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_rows() == a_nb_rows ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                           result->item(i) == 0.0 ) ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( ! result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create( PEL_Object* a_owner, doubleVector const& dvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: create" ) ;
   PEL_CHECK_PRE( dvec.size()>0 ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, dvec.size() ) ;
   for( size_t i=0 ; i<dvec.size() ; ++i )
   {
      result->set_item( i, dvec(i) ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_rows() == dvec.size() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                           result->item(i) == dvec(i) ) ) ;
   PEL_CHECK_POST( result->state() == LA::Sync ) ;
   PEL_CHECK_POST( result->is_resizable() ) ;
   PEL_CHECK_POST( ! result->is_desynchronizable() ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;
   PEL_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector:: LA_SeqVector( PEL_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, a_nb_rows )
   , UNSYNCHRO( false )
   , ROW_DIST( 0 )
   , OWNS_DATA( true )
   , DATA( 0 )
{
   PEL_LABEL( "LA_SeqVector:: LA_SeqVector" ) ;

   if( a_nb_rows>0)
   {
      DATA = new double [ a_nb_rows ] ;
      for( size_t i=0 ; i<a_nb_rows ; ++i ) DATA[i] = 0.0 ;
   }

   set_distribution_strategy( LA::NoDistribution ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::NoDistribution ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: create_vector" ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, nb_rows() ) ;
   if( UNSYNCHRO ) result->make_desynchronizable() ;

   PEL_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           result->item(i) == 0. ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector:: ~LA_SeqVector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: ~LA_SeqVector" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( OWNS_DATA && DATA != 0 )
   {
      delete [] DATA ;
   }
   DATA = 0 ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_SeqVector:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: implementation" ) ;

   LA_Implementation const* result = LA_SeqImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: is_desynchronizable" ) ;

   bool result = UNSYNCHRO ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: make_desynchronizable( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: make_desynchronizable" ) ;
   PEL_CHECK_PRE( ! is_desynchronizable() ) ;

   UNSYNCHRO = true ;

   PEL_CHECK_POST( is_desynchronizable() ) ;
   PEL_CHECK_POST( state() == LA::Sync ) ;
   PEL_CHECK_POST( is_synchronized() ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
LA_SeqVector:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: row_distribution" ) ;
   PEL_CHECK_PRE( row_distribution_PRE() ) ;

   if( ROW_DIST == 0 )
   {
      ROW_DIST = PEL_DistributedPartition::create(
                                 const_cast<LA_SeqVector*>( this ) ) ;
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
void
LA_SeqVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: re_initialize" ) ;
   PEL_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;

   re_initialize_with_global_sizes( a_nb_rows ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           item(i) == 0.0 ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: re_initialize_with_global_sizes( size_t a_nb_rows )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: re_initialize_with_global_sizes(size_t)" ) ;
   PEL_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( OWNS_DATA ) ;

   if( a_nb_rows != nb_rows() )
   {
      if( DATA != 0 )
      {
         delete [] DATA ;
         DATA = 0 ;
      }
      if( a_nb_rows>0)
      {
         DATA = new double [ a_nb_rows ] ;
         for( size_t i=0 ; i<a_nb_rows ; ++i ) DATA[i] = 0.0 ;
      }
      set_rows_number( a_nb_rows ) ;
      if( UNSYNCHRO )
      {
         synchronize() ;
      }
   }
   else
   {
      nullify() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           item(i) == 0.0 ) ) ;
}

//----------------------------------------------------------------------
LA_SeqScatter*
LA_SeqVector:: create_scatter( PEL_Object* a_owner,
                               size_t_vector const& from,
                               size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: create_scatter" ) ;
   PEL_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   LA_SeqScatter* result =
                LA_SeqScatter::create( a_owner, nb_rows(), from, to ) ;

   PEL_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create_local_vector( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: create_local_vector" ) ;
   PEL_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = create_vector( a_owner ) ;
   result->set( this ) ;

   PEL_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: set" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( set_PRE( a ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   if( this != a )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = ba->DATA[i] ;
   }
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( a ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: set( double value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: set double" ) ;
   PEL_CHECK_PRE( set_PRE( value ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = value ;
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_POST( value ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: nullify( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: nullify" ) ;
   PEL_CHECK_PRE( nullify_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = 0.0 ;
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nullify_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: scale" ) ;
   PEL_CHECK_PRE( scale_PRE( alpha ) ) ;

   if( alpha == 0. )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         DATA[i] = 0. ;
      }
   }
   else if( alpha != 1. )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         DATA[i] *= alpha ;
      }
   }

   PEL_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: sum" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( sum_PRE( a, alpha ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( a ) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>( a ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      DATA[i] += alpha*ba->DATA[i] ;
   }

   PEL_CHECK_POST( sum_POST( a, alpha ) ) ;
}

//----------------------------------------------------------------------
double
LA_SeqVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: dot" ) ;
   PEL_CHECK_PRE( dot_PRE(a) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   double result = 0. ;
   for( size_t i=0 ; i<a->nb_rows() ; ++i )
   {
      result += DATA[i] * ba->DATA[i] ;
   }

   PEL_CHECK_POST( dot_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_SeqVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: two_norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( two_norm_PRE() ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      result += DATA[i] * DATA[i] ;
   }
   result = PEL::sqrt( result ) ;

   PEL_CHECK_POST( two_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_SeqVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   PEL_LABEL( "LA_SeqVector:: max_norm" ) ;
   PEL_CHECK_PRE( max_norm_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      if( DATA[i]>result )
      {
         result = DATA[i] ;
      }
      else if( -DATA[i]>result )
      {
         result = -DATA[i] ;
      }
   }
   PEL_CHECK_POST( max_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_SeqVector:: mean( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: mean" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( nb_rows() > 0 ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      result += DATA[i] ;
   }
   return( result/double( nb_rows() ) ) ;
}

//----------------------------------------------------------------------
double
LA_SeqVector:: standard_deviation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: standard_deviation" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( nb_rows() > 1 ) ;
   double result = 0.0 ;

   if( nb_rows() > 1 )
   {
      double variance = 0.0 ;
      double moyenne =  mean() ;
      double elem ;
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         elem = DATA[i] - moyenne ;
         variance += elem*elem ;
      }
      variance /= double( nb_rows() )  ;
      result = PEL::sqrt( variance ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: print_items" ) ;
   PEL_CHECK_PRE( print_items_PRE( os, indent_width ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const space( indent_width, ' ' ) ;
   os << space << "nb_rows:" << nb_rows() << std::endl ;
   if( nb_rows()>0 )
   {
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      std::streamsize p = os.precision() ;
      os << std::setprecision( 7 ) ;
      for( size_t iRow = 0 ; iRow<nb_rows() ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = DATA[iRow] ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
   }
}

//----------------------------------------------------------------------
void
LA_SeqVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: write" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::ofstream file( filename.c_str() ) ;
   PEL_ASSERT( file ) ;
   file.precision( 15 ) ;

   file << nb_rows() ;
   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      file << std::endl << DATA[iRow] ;
   }
}

//----------------------------------------------------------------------
void
LA_SeqVector:: save( PEL_Module* mod ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: save" ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( mod != 0 ) ;

   mod->add_entry( "nb_rows",
                   PEL_Int::create( mod, (int)nb_rows() ) ) ;
   doubleVector dvec( nb_rows() );
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      dvec(i)=DATA[i] ;
   }

   mod->add_entry( "data", PEL_DoubleVector::create( mod, dvec ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: restore( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: restore" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   re_initialize_with_global_sizes( exp->int_data("nb_rows") ) ;
   doubleVector const& v = exp->doubleVector_data("data");
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      DATA[i]=v(i) ;
   }
   synchronize() ;

   PEL_CHECK_POST( is_synchronized() ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: set_as_v_product( LA_Vector const* a, LA_Vector const* b )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: set_as_v_product" ) ;
   PEL_CHECK_PRE( set_as_v_product_PRE( a, b ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector const*>(b) != 0 ) ;
   LA_SeqVector const* bb = static_cast<LA_SeqVector const*>(b) ;

   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      DATA[iRow] = ba->DATA[iRow] * bb->DATA[iRow] ;
   }
   synchronize() ;

   PEL_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqVector:: set_as_reciprocal( LA_Vector const* a,
                                  double smallest_inverted_item,
                                  double default_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqVector:: set_as_reciprocal" ) ;
   PEL_CHECK_PRE( set_as_reciprocal_PRE( a, smallest_inverted_item, default_value ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      DATA[iRow] = ( PEL::abs( ba->DATA[iRow] ) < smallest_inverted_item ?
                     default_value :
                     1.0 / ba->DATA[iRow] ) ;
   }
   synchronize() ;

   PEL_CHECK_POST( set_as_reciprocal_POST( a ) ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqVector:: re_initialize_with_global_sizes_PRE( size_t a_nb_rows ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_resizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqVector:: re_initialize_with_global_sizes_POST( size_t a_nb_rows ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( nb_rows() == a_nb_rows ) ;
   PEL_ASSERT( state() == LA::Sync ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::invariant() ) ;
   PEL_ASSERT( distribution_strategy() == LA::NoDistribution ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_SeqVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   PEL_ASSERT( result == LA_SeqImplementation::object() ) ;
   return( true ) ;
}
