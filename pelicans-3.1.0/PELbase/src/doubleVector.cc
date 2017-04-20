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

#include <doubleVector.hh>

#ifdef OUTLINE
#define inline
#include <doubleVector.icc>
#undef inline
#endif

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_System.hh>

#include <iostream>

struct doubleVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
doubleVector:: doubleVector( size_t dim, double val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "doubleVector:: doubleVector( size_t )" ) ;
   
   allocate( LENGTH ) ;
   set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
doubleVector:: doubleVector( int dim, double val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( static_cast<size_t>( dim ) )
   , CAPACITY( 0 )
{
   PEL_LABEL( "doubleVector:: doubleVector( int )" ) ;
   PEL_CHECK_PRE( dim >= 0 ) ;

   allocate( LENGTH ) ;
   set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == static_cast<size_t>(dim) ) ;
   PEL_CHECK_POST( capacity() == size() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
doubleVector:: ~doubleVector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: ~doubleVector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
doubleVector:: doubleVector( doubleVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   PEL_LABEL( "doubleVector:: doubleVector( doubleVector const& )" ) ;

   allocate( other.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = other.VECTOR[i] ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == other.size() ) ;
   PEL_CHECK_POST( capacity() == other.capacity() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == other(i) ) ) ;
}

//----------------------------------------------------------------------
doubleVector const&
doubleVector:: operator=( doubleVector const& other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: operator=" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   if( this != &other )
   {
      LENGTH = other.LENGTH ;
      if( LENGTH > CAPACITY )
      {
         deallocate() ;
         allocate( LENGTH ) ;
      }
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = other.VECTOR[i] ;
      }
   }
   doubleVector const& result = *this ;
  
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result.size() == other.size() ) ;
   PEL_CHECK_POST( IMPLIES( other.size() <= OLD(capacity),
                            result.capacity() == OLD(capacity) ) ) ;
   PEL_CHECK_POST( IMPLIES( other.size() > OLD(capacity),
                            result.capacity() == other.size() ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   result(i) == other(i) ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
doubleVector:: re_initialize( size_t dim, double val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: re_initialize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   LENGTH = dim ;
   if( CAPACITY != dim )
   {
      deallocate() ;
      allocate( LENGTH ) ;
   }
   set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
            operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
bool
doubleVector:: operator==( doubleVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = ( LENGTH==other.LENGTH ) ;
      for( size_t i=0 ; result && i<size() ; ++i )
      {
         result = ( operator()(i)==other(i) ) ;
      }
   }

   PEL_CHECK_POST( IMPLIES( result, size()==other.size() ) ) ;
   PEL_CHECK_POST( !result || FORALL( ( size_t i=0 ; i<size() ; ++i ),
                                       operator()(i) == other(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
doubleVector:: operator!=( doubleVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   PEL_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
doubleVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   PEL_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
doubleVector:: index_of( double val, double a_dbl_eps, double a_dbl_min ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: index_of" ) ;
   
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      if( PEL::double_equality( VECTOR[i], val, a_dbl_eps, a_dbl_min ) )
      {
         result = i ;
         break ;
      }
   }
   
   PEL_CHECK_POST( 
      result == PEL::bad_index() ||
      PEL::double_equality( operator()(result) , val, 
                            a_dbl_eps, a_dbl_min ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
doubleVector:: has( double val, double a_dbl_eps, double a_dbl_min ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: has" ) ;
   
   size_t ifound = index_of( val, a_dbl_eps, a_dbl_min ) ;
   bool result = ( ifound != PEL::bad_index() ) ;

   PEL_CHECK_POST(
      IMPLIES( 
         result, 
         PEL::double_equality( operator()(
                                  index_of( val, a_dbl_eps, a_dbl_min ) ), 
                               val, a_dbl_eps, a_dbl_min ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
doubleVector:: resize( size_t dim, double val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: resize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = PEL_System::new_block_size( CAPACITY, dim ) ;
      double const* oldVector = VECTOR ;
      VECTOR = 0 ;
      allocate( newCapacity ) ;
      if( oldVector!=0 )
      {
         for( size_t i=0 ; i<LENGTH ; ++i )
         {
            VECTOR[i] = oldVector[i] ;
         }
         delete[] oldVector ;
      }
   }
   if( dim > LENGTH )
   {
      for( size_t i=LENGTH ; i<dim ; ++i )
      {
         VECTOR[i] = val ;
      }
   }
   LENGTH = dim ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( dim <= OLD(capacity) ,
                            capacity() == OLD(capacity) ) ) ;
   PEL_CHECK_POST( 
      IMPLIES( dim > OLD(capacity) ,
               capacity() == PEL_System::new_block_size( OLD(capacity),dim ))) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=OLD(size) ; i<size() ; ++i ),
                           operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
void
doubleVector:: set( double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: set" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = val ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
void
doubleVector:: append( double val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: append" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t n = LENGTH ;
   resize( n+1 ) ;
   operator()( n ) = val ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( 
      IMPLIES( OLD(size)+1 <= OLD(capacity),
               capacity() == OLD( capacity ) ) ) ;
   PEL_CHECK_POST( 
      IMPLIES( OLD(size)+1 > OLD(capacity),
               capacity() == PEL_System::new_block_size( OLD(capacity),
                                                         OLD(size)+1 ) ) ) ;
   PEL_CHECK_POST( size() == OLD(size)+1 ) ;
   PEL_CHECK_POST( operator()(size()-1) == val ) ;
}

//----------------------------------------------------------------------
void
doubleVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleVector:: remove_at" ) ;
   PEL_CHECK_PRE( idx<size() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   for( size_t j=idx ; j<size()-1 ; j++ )
      VECTOR[j]=VECTOR[j+1] ;
   LENGTH-- ;

   PEL_CHECK_POST( size()==OLD(size)-1 ) ;
   PEL_CHECK_POST( capacity() == OLD(capacity) ) ;
}

//----------------------------------------------------------------------
std::ostream& operator<<( std::ostream& out, 
                          doubleVector const& vec ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, doubleVector const& )" ) ;

   size_t j=0 ;
   out << "< " ;
   for( size_t i=0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      PEL::print_double( out, vec(i) ) ;
      out << " " ;
   }
   out << ">" ;
   return out ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, doubleVector& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator>>( std::ostream&, doubleVector const& )" ) ;

   doubleVector aux( 0 ) ;
   
   std::string c ;
   double val = 0. ;
   
   if( !( in >> c ) || c != "<" )
   {
      doubleVector_ERROR::n0() ;
   }
   while( in >> val )
   {
      aux.append( val ) ;
   }
   in.clear() ;
   if( !( in >> c ) || c != ">" )
   {
      doubleVector_ERROR::n0() ;
   }

   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
doubleVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   PEL_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new double [a_size] ;
      PEL_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
doubleVector:: deallocate( void )
//----------------------------------------------------------------------
{
   if( VECTOR != 0 )
   {
      delete [] VECTOR ;
      VECTOR = 0 ;
   }
   CAPACITY = 0 ;
}

//-----------------------------------------------------------------------
bool
doubleVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   PEL_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
doubleVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"double\":\n"
      "    < val1 val2 ... > is expected." ) ;
}

