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

#include <PEL_Communicator.hh>

#include <PEL.hh>
#include <PEL_DoubleComparator.hh>
#include <PEL_Error.hh>
#include <PEL_NumberedDoubleVectors.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <boolArray2D.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <stringVector.hh>

#include <string>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

struct PEL_Communicator_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( void ) ;
   static void n3( void ) ;
} ;

bool PEL_Communicator:: TRACE = false ;

//----------------------------------------------------------------------
PEL_Communicator:: PEL_Communicator( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , MY_NAME( a_name )
{
   PEL_LABEL( "PEL_Communicator:: PEL_Communicator" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_under_ownership_of( PEL_Root::object() ) ) ;
}

//----------------------------------------------------------------------
PEL_Communicator:: ~PEL_Communicator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PEL_Communicator const*
PEL_Communicator:: object( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: object" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_Communicator const* result =
      static_cast<PEL_Communicator const*>(
                               plugins_map()->item( a_name ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Communicator:: name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: name" ) ;
   return MY_NAME ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, size_t value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send size_t" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int val = (int) value ;
   send( dest, &val, 1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, size_t& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive size_t" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int val = PEL::bad_int() ;
   receive( src, &val, 1 ) ;
   if( val != (int) PEL::bad_index() && val < 0 )
      PEL_Communicator_ERROR:: n2() ;
   value = (size_t) val ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, size_t_vector const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send size_t_vector" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=0 ; i<value.size() ; ++i ),
         value(i)==PEL::bad_index() || (int) value(i)<PEL::max_int() ) ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, size_t_vector& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive size_t_vector" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, size_t_array2D const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send size_t_array2D" ) ;
   PEL_CHECK_PRE( dest<nb_ranks() ) ;
   PEL_CHECK_PRE( dest!=rank() ) ;
   PEL_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<value.index_bound(0) ; ++i ),
         FORALL(
            ( size_t j=0 ; j<value.index_bound(1) ; ++j ),
            value(i,j)==PEL::bad_index() || (int) value(i,j)<PEL::max_int() ) ) ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, size_t_array2D& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive size_t_array2D" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, int value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send int" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   send( dest, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, int& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive int" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   receive( src, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, intVector const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send intVector" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      send( dest, value.data(), dim ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, intVector& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive intVector" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = PEL::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) PEL_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) dim ) ;
   if( dim>0 )
   {
      receive( src, const_cast<int*>( value.data() ), dim ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, intArray2D const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send intArray2D" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] ;
   D[0] = (int) value.index_bound( 0 ) ;
   D[1] = (int) value.index_bound( 1 ) ;
   send( dest, D, 2 ) ;
   if( D[0]*D[1]>0 )
   {
      send( dest, value.data(), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, intArray2D& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive intArray2D" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] = { PEL::bad_int(), PEL::bad_int() } ;
   receive( src, D, 2 ) ;
   if( D[0]<0 || D[1]<0 ) PEL_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) D[0], (size_t) D[1] ) ;
   if( D[0]*D[1]>0 )
   {
      receive( src, const_cast<int*>( value.data() ) , D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send double" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   send( dest, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, double& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive double" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   receive( src, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, doubleVector const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send doubleVector" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      send( dest, value.data(), dim ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, doubleVector& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive doubleVector" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = PEL::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) PEL_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) dim ) ;
   if( dim>0 )
   {
      receive( src, const_cast<double*>( value.data() ), dim ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, doubleArray2D const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send doubleArray2D" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] ;
   D[0] = (int) value.index_bound( 0 ) ;
   D[1] = (int) value.index_bound( 1 ) ;
   send( dest, D, 2 ) ;
   if( D[0]*D[1]>0 )
   {
      send( dest, value.data(), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, doubleArray2D& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive doubleArray2D" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] = { PEL::bad_int(), PEL::bad_int() } ;
   receive( src, D, 2 ) ;
   if( D[0]<0 || D[1]<0 ) PEL_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) D[0], (size_t) D[1] ) ;
   if( D[0]*D[1]>0 )
   {
      receive( src, const_cast<double*>( value.data() ), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, bool value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send bool" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int i = ( value ? 1 : 0 ) ;
   send( dest, &i, 1 ) ;

}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, bool& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive bool" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int i ;
   receive( src, &i, 1 ) ;
   if( i != 0 && i != 1 ) PEL_Communicator_ERROR:: n3() ;
   value = ( i==1 ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, boolVector const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send boolVector" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, boolVector& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive boolVector" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, boolArray2D const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send boolArray2D" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, boolArray2D& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive boolArray2D" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, std::string const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send string" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) (value.size()+1) ;
   send( dest, &dim, 1 ) ;
   send( dest, value.c_str(), dim ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, std::string& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive string" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = PEL::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<1 ) PEL_Communicator_ERROR::n0() ;
   char* val = new char[dim] ;
   receive( src, val, dim ) ;
   if( val[dim-1] != '\0' ) PEL_Communicator_ERROR::n1() ;
   value = std::string( val ) ;
   delete [] val ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: send( size_t dest, stringVector const& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: send stringVector" ) ;
   PEL_CHECK_PRE( dest < nb_ranks() ) ;
   PEL_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      for( size_t i=0 ; i<value.size() ; ++i )
      {
         send( dest, value(i) ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: receive( size_t src, stringVector& value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: receive stringVector" ) ;
   PEL_CHECK_PRE( src < nb_ranks() ) ;
   PEL_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = PEL::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) PEL_Communicator_ERROR::n0() ;
   value.re_initialize( dim ) ;
   if( dim>0 )
   {
      for( size_t i=0 ; i<value.size() ; ++i )
      {
         receive( src, value(i) ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: wait( void* request ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( size_t& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(size_t)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      int dummy = ( rank() != root ? PEL::bad_int() : (int) value ) ;
      broadcast( &dummy, (size_t)1, root ) ;
      if( dummy != (int) PEL::bad_index() && dummy < 0 )
         PEL_Communicator_ERROR:: n2() ;
      value = (size_t) dummy ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( size_t_vector& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(size_t_vector)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intVector dummy(0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( size_t_array2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(size_t_array2D)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intArray2D dummy(0,0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( int& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(int)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   broadcast( &value, (size_t)1, root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( intVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(intVector)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim = ( rank() != root ? PEL::bad_int() : (int) values.size() ) ;
   broadcast( &dim , (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim ) ;
   broadcast( const_cast<int*>( values.data() ), (size_t) dim, root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( intArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(intArray2D)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim0 = ( rank() != root ? PEL::bad_int() : (int) values.index_bound(0) ) ;
   int dim1 = ( rank() != root ? PEL::bad_int() : (int) values.index_bound(1) ) ;
   broadcast( &dim0, (size_t) 1, root ) ;
   broadcast( &dim1, (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim0, dim1 ) ;
   broadcast( const_cast<int*>( values.data() ), (size_t) dim0*dim1, root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( double& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(double)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   broadcast( &value, (size_t)1, root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( doubleVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(doubleVector)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim = ( rank() != root ? PEL::bad_int() : (int) values.size() ) ;
   broadcast( &dim , (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim ) ;
   broadcast( const_cast<double*>( values.data() ), (size_t) dim, root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( doubleArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(doubleArray2D)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim0 = ( rank() != root ? PEL::bad_int() : (int) values.index_bound(0) ) ;
   int dim1 = ( rank() != root ? PEL::bad_int() : (int) values.index_bound(1) ) ;
   broadcast( &dim0, (size_t) 1, root ) ;
   broadcast( &dim1, (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim0, dim1 ) ;
   broadcast( const_cast<double*>( values.data() ), (size_t) dim0*dim1, root ) ;
}


//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( bool& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(bool)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      int dummy = ( rank() != root ? PEL::bad_int() : ( value ? 1 : 0 ) ) ;
      broadcast( &dummy, (size_t)1, root ) ;
      if( dummy != 0 && dummy != 1 )
         PEL_Communicator_ERROR:: n2() ;
      value = ( dummy == 1 ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( boolVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(boolVector)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intVector dummy(0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( boolArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(boolArray2D)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intArray2D dummy(0,0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( std::string& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(string)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim = PEL::bad_int() ;
   if( rank() == root ) dim = int(value.size()+1) ;
   broadcast( dim, root ) ;

   char* tab = 0 ;
   if( rank() == root )
   {
      tab = const_cast<char*>( value.c_str() ) ;
   }
   else
   {
      tab = new char[dim] ;
   }
   broadcast( tab, dim, root ) ;
   if( rank() != root )
   {
      value = std::string( tab ) ;
      delete [] tab ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: broadcast( stringVector& value, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: broadcast(stringVector)" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;

   int dim = PEL::bad_int() ;
   if( rank() == root ) dim = (int) value.size() ;
   broadcast( dim, root ) ;
   if( rank() != root ) value.re_initialize( (size_t) dim ) ;
   for( size_t i=0 ; i<value.size() ; ++i )
   {
      broadcast( value(i), root ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: gather( double value, doubleVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: gather(double)" ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;
   PEL_CHECK_PRE( root!=rank() || result.size()==nb_ranks() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   gather( &value, 1, const_cast<double*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: gather( int value, intVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: gather(double)" ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;
   PEL_CHECK_PRE( root!=rank() || result.size()==nb_ranks() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   gather( &value, 1, const_cast<int*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: gather( doubleVector const& values,
                           doubleVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: gather(intVector)" ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;
   PEL_CHECK_PRE( root!=rank() || result.size()==nb_ranks()*values.size() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   gather( values.data(), values.size(),
           const_cast<double*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: gather( intVector const& values,
                           intVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: gather(intVector)" ) ;
   PEL_CHECK_PRE( root < nb_ranks() ) ;
   PEL_CHECK_PRE( root!=rank() || result.size()==nb_ranks()*values.size() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   gather( values.data(), values.size(),
           const_cast<int*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: boolean_and( bool value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: boolean_and" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   bool result = value ;
   if( rank()==0 )
   {
      bool recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result &= recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: boolean_or( bool value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: boolean_or" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   bool result = value ;
   if( rank()==0 )
   {
      bool recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result |= recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
double
PEL_Communicator:: sum( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: sum" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank()==0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result += recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: sum_vector( doubleVector& vec ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: sum_vector" ) ;
   PEL_CHECK_PRE( same_value_everywhere( (int) vec.size() ) ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_SAVEOLD( size_t, vec_size, vec.size() ) ;

   if( trace() ) print_method() ;

   int const nb = (int) vec.size() ;
   if( nb_ranks()>1 && nb>0 )
   {
      sum_vector( const_cast<double*>( vec.data() ), nb ) ;
   }
   
   PEL_CHECK_POST( vec.size() == OLD( vec_size ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: sum_array( doubleArray2D& array ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: sum_array" ) ;
   PEL_CHECK_PRE( same_value_everywhere( (int) array.index_bound(0) ) ) ;
   PEL_CHECK_PRE( same_value_everywhere( (int) array.index_bound(1) ) ) ;
   PEL_SAVEOLD( size_t, array_index_bound0, array.index_bound(0) ) ;
   PEL_SAVEOLD( size_t, array_index_bound1, array.index_bound(1) ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   int const nb = (int) array.index_bound(0)*array.index_bound(1) ;
   if( nb_ranks()>1 && nb>0 )
   {
      sum_vector( const_cast<double*>( array.data() ), nb ) ;
   }

   PEL_CHECK_POST( array.index_bound(0) == OLD( array_index_bound0 ) ) ;
   PEL_CHECK_POST( array.index_bound(1) == OLD( array_index_bound1 ) ) ;
}


//----------------------------------------------------------------------
void
PEL_Communicator:: merge( PEL_DoubleComparator const* dbl_comp,
                          doubleArray2D& coord,
                          size_t_vector& idx ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: merge" ) ;
   PEL_CHECK_PRE( dbl_comp != 0 ) ;
   PEL_SAVEOLD( doubleArray2D, coord, coord ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   size_t dim = coord.index_bound( 0 ) ;

   size_t_array2D idx_all( 0, nb_ranks() ) ;
   size_t_vector nbvert( nb_ranks() ) ;
   doubleArray2D coords( 0, 0 ) ;

   if( rank()>0 )
   {
      send( 0, coord ) ;
      receive( 0, idx ) ;
   }
   else
   {
      PEL_NumberedDoubleVectors* numvec =
                  PEL_NumberedDoubleVectors::create( 0, dbl_comp, dim ) ;
      doubleVector vec( dim ) ;
      for( size_t i=0 ; i<nb_ranks() ; i++ )
      {
         if( i>0 )
         {
            receive( i, coords ) ;
         }
         else
         {
            coords = coord ;
         }
         PEL_CHECK( coords.index_bound( 0 ) == dim ) ;
         size_t n = nbvert( i ) = coords.index_bound( 1 ) ;
         if( idx_all.index_bound( 0 )<nbvert( i ) )
            idx_all.raise_first_index_bound( nbvert( i ) ) ;
         for( size_t j=0 ; j<n; j++ )
         {
            for( size_t d=0 ; d<dim; d++ )
               vec( d ) = coords( d, j ) ;
            numvec->extend( vec ) ;
            idx_all( j, i ) = numvec->index( vec ) ;
         }
      }

      // Re-ordering nodes to be independent from splitting
      size_t_vector const& perm = numvec->order() ;

      for( size_t i=nb_ranks()-1 ; i<nb_ranks() ; i-- )
      {
         idx.re_initialize( nbvert( i ) ) ;
         for( size_t j=0 ; j<idx.size() ; j++ )
         {
            idx( j ) = perm( idx_all( j, i ) ) ;
         }
         if( i>0 )
         {
            send( i, idx ) ;
         }
      }
      coord = numvec->ordered_items() ;
      numvec->destroy() ;
   }

   PEL_CHECK_POST( IMPLIES( rank()!=0, coord==OLD(coord) ) ) ;
   PEL_CHECK_POST( rank() != 0 ||
      FORALL( ( size_t i=0 ; i<coord.index_bound(1)-1 ; ++i ),
         dbl_comp->three_way_comparison( coord(0,i), coord(0,i+1) ) <= 0 ) ) ;
   PEL_CHECK_POST( coord.index_bound( 0 ) == OLD(coord).index_bound( 0 ) ) ;
   PEL_CHECK_POST( idx.size() == OLD(coord).index_bound( 1 ) ) ;
}

//----------------------------------------------------------------------
double
PEL_Communicator:: min( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: min" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank() == 0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = PEL::min( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   PEL_CHECK_POST( min_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
size_t
PEL_Communicator:: min( size_t value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: min" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   size_t result = value ;
   if( rank() == 0 )
   {
      size_t recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = PEL::min( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   PEL_CHECK_POST( min_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
double
PEL_Communicator:: max( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: max" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank() == 0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = PEL::max( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   PEL_CHECK_POST( max_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
size_t
PEL_Communicator:: max( size_t value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: max" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   size_t result = value ;
   if( rank() == 0 )
   {
      size_t recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = PEL::max( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   PEL_CHECK_POST( max_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
void
PEL_Communicator:: all_gather( int value, intVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: all_gather" ) ;
   PEL_CHECK_PRE( result.size()==nb_ranks() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   all_gather( &value, 1, const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: all_gather( intVector const& values, intVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: all_gather" ) ;
   PEL_CHECK_PRE( values.size()!=0 ) ;
   PEL_CHECK_PRE( result.size()==nb_ranks()*values.size() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   all_gather( values.data(), values.size(), const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: all_to_all( intVector const& values,
                               intVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: all_to_all" ) ;
   PEL_CHECK_PRE( values.size() % nb_ranks() == 0 ) ;
   PEL_CHECK_PRE( result.size()==values.size() ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   int N = values.size() / nb_ranks() ;

   all_to_all( values.data(), N, const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_Communicator:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
           PEL_ObjectRegister::create( PEL_Root::object(),
                                       "PEL_Communicator descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: trace( void )
//----------------------------------------------------------------------
{
   return( TRACE ) ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: print_method( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: print_method" ) ;
   if( PEL_Marker::nb_labels()>1 )
   {
      PEL::out() << "Methods : " ;
      for( size_t i=0 ; i<PEL_Marker::nb_labels()-1 ; ++i )
      {
         PEL::out() << PEL_Marker::label(i) << ";" ;
      }
      PEL::out() << std::endl ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: do_trace( void )
//----------------------------------------------------------------------
{
   TRACE = true ;
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                       size_t_vector const& src, intVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      dest(i) = (int) src(i) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                       intVector const& src, size_t_vector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      if( src(i)!= (int) PEL::bad_index() && src(i) < 0 )
         PEL_Communicator_ERROR:: n2() ;
      dest(i) = (size_t) src(i) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                     size_t_array2D const& src, intArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         dest(i,j) = (int) src(i,j) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                     intArray2D const& src, size_t_array2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         if( src(i,j) != (int) PEL::bad_index() && src(i,j) < 0 )
            PEL_Communicator_ERROR:: n2() ;
         dest(i,j) = (size_t) src(i,j) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                          boolVector const& src, intVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      dest(i) = ( src(i) ? 1 : 0 ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                          intVector const& src, boolVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      if( src(i) != 0 && src(i) != 1 )
         PEL_Communicator_ERROR:: n3() ;
      dest(i) = ( src(i) == 1 ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                        boolArray2D const& src, intArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         dest(i,j) = ( src(i,j) ? 1 : 0 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: convert(
                        intArray2D const& src, boolArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         if( src(i,j) != 0 && src(i,j) != 1 )
            PEL_Communicator_ERROR:: n3() ;
         dest(i,j) = ( src(i,j) == 1 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_Communicator:: sum_vector( double* values, int nb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator:: sum_vector(double*)" ) ;
   PEL_CHECK_PRE( sum_vector_PRE( values, nb ) ) ;
   
   if( rank()==0 )
   {
      int s = nb ;
      double* recep = new double[nb] ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep, s ) ;
         PEL_ASSERT( s == nb ) ;
         for( int j=0 ; j<nb ; ++j )
            values[j] += recep[j] ;
      }
      delete [] recep ;
   }
   else
   {
      send( 0, values, nb ) ;
   }
   broadcast( values, nb, 0 ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: nb_ranks_POST( size_t result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result > 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: rank_POST( size_t result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result < nb_ranks() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: send_PRE( size_t dest, int const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( dest < nb_ranks() ) ;
   PEL_ASSERT( dest != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: send_PRE( size_t dest, double const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( dest < nb_ranks() ) ;
   PEL_ASSERT( dest != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: send_PRE( size_t dest, char const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( dest < nb_ranks() ) ;
   PEL_ASSERT( dest != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: receive_PRE( size_t src, int const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( src < nb_ranks() ) ;
   PEL_ASSERT( src != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: receive_PRE( size_t src, double const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( src < nb_ranks() ) ;
   PEL_ASSERT( src != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: receive_PRE( size_t src, char const* value, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( src < nb_ranks() ) ;
   PEL_ASSERT( src != rank() ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: min_POST( double result, double value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result <= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: min_POST( size_t result, size_t value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result <= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: sum_vector_PRE( double const* values, int nb ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( same_value_everywhere( nb ) ) ;
   PEL_ASSERT( nb > 0 ) ;
   PEL_ASSERT( nb_ranks() > 1 ) ;
   PEL_ASSERT( values != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: max_POST( double result, double value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result >= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_Communicator:: max_POST( size_t result, size_t value ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result >= value ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void
PEL_Communicator_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** PEL_Communicator:" << std::endl
        << "    negative table bound received" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PEL_Communicator_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** PEL_Communicator:" << std::endl
        << "    pointer to a no null-terminated array of characters received." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PEL_Communicator_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** PEL_Communicator:" << std::endl
        << "    negative size_t received." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PEL_Communicator_ERROR:: n3( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** PEL_Communicator:" << std::endl
        << "    bad boolean value received." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
