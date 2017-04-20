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

#include <stringVector.hh>

#ifdef OUTLINE
#define inline
#include <stringVector.icc>
#undef inline
#endif

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_System.hh>
#include <iostream>

struct stringVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
stringVector:: stringVector( size_t dim )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "stringVector:: stringVector( size_t )" ) ;

   allocate( LENGTH ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; i++ ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( int dim )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "stringVector:: stringVector( int )" ) ;
   PEL_CHECK_PRE( dim >= 0 ) ;

   allocate( LENGTH ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == static_cast<size_t>(dim) ) ;
   PEL_CHECK_POST( capacity() == size() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; i++ ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( std::string const& elements,
                             char const separator )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( 0 )
   , CAPACITY( 0 )
{
   PEL_LABEL( "stringVector:: stringVector( string, char )" ) ;

   build_from_string( elements, separator ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( char const* elements,
                             char const separator )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( 0 )
   , CAPACITY( 0 )
{
   PEL_LABEL( "stringVector:: stringVector( const char*, char )" ) ;

   build_from_string( elements, separator ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( stringVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   PEL_LABEL( "stringVector:: stringVector( stringVector const& )" ) ;

   allocate( other.capacity() ) ;
   if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = other(i) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == other.size() ) ;
   PEL_CHECK_POST( capacity() == other.capacity() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
                           operator()(i) == other(i) ) ) ;
}

//----------------------------------------------------------------------
stringVector:: ~stringVector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: ~stringVector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
stringVector const&
stringVector:: operator=( stringVector const& other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: operator=" ) ;
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
      if( LENGTH != 0 )
      {
         for( size_t i=0 ; i<LENGTH ; ++i )
         {
            VECTOR[i] = other(i) ;
         }
      }
   }
   stringVector const& result = *this ;

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
stringVector:: re_initialize( size_t dim )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: re_initialize" ) ;
   PEL_CHECK_INV( invariant() ) ;

   LENGTH = dim ;
   if( CAPACITY != dim )
   {
      deallocate() ;
      allocate( LENGTH ) ;
   }
   else if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = "" ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator==( stringVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = (LENGTH == other.LENGTH) ;
      for( size_t i=0 ; result && i<size() ; i++ )
      {
         result = operator()(i)==other(i) ;
      }
   }

   PEL_CHECK_POST( IMPLIES( result, size()==other.size() ) ) ;
   PEL_CHECK_POST( !result || FORALL( ( size_t i=0 ; i<size() ; ++i ),
                                       operator()(i) == other(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator!=( stringVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   PEL_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
int
stringVector:: three_way_comparison( stringVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: three_way_comparison" ) ;

   size_t min_size = PEL::min( size(), other.size() ) ;
   int result = 0 ;

   size_t i=0 ;
   for( ; i < min_size ; ++i )
   {
      int comp = operator()(i).compare( other.operator()(i) ) ;
      if( comp != 0 )
      {
         result = comp ;
         break ;
      }
   }

   if( result == 0 )
   {
      if( i == size() && i< other.size() )
      {
         result = -1 ;
      }
      else if( i < size() && i == other.size() )
      {
         result = 1 ;
      }
   }

  return( result ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator>( stringVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: operator>" ) ;

  return( three_way_comparison( other ) > 0 ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator<( stringVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: operator<" ) ;

   return( three_way_comparison( other ) < 0 ) ;
}

//----------------------------------------------------------------------
void stringVector:: build_from_string( std::string const& elements,
                                       char const separator )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: build_from_string" ) ;

   if( !elements.empty() )
   {
      size_t start = 0 ;
      size_t idx = 0 ;
      size_t end = elements.length()  ;
      
      while( ( idx = elements.find( separator, start ) ) < end )
      {
         std::string a_choice = elements.substr( start, idx-start ) ;
         PEL::remove_enclosing_characters(a_choice,'\"') ;
         append( a_choice ) ;
         start = idx+1 ;
      }
      std::string a_choice = elements.substr( start, elements.length()-start ) ;
      PEL::remove_enclosing_characters(a_choice,'\"') ;
      append( a_choice ) ;
   }
}

//----------------------------------------------------------------------
size_t
stringVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: capacity" ) ;

   size_t result = CAPACITY ;

   PEL_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
stringVector:: index_of( std::string const& val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: index_of" ) ;

   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; i++ )
   {
      if ( VECTOR[i] == val )
      {
         result = i ;
         break ;
      }
   }

   PEL_CHECK_POST( result == PEL::bad_index() || operator()(result) == val ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: has( std::string const& val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: has" ) ;

   size_t ifound = index_of( val ) ;
   bool result = ( ifound != PEL::bad_index() ) ;

   PEL_CHECK_POST( IMPLIES( result, operator()(index_of(val))==val ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringVector:: resize( size_t dim, std::string const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: resize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY )
   {
      size_t const newCapacity = PEL_System::new_block_size( CAPACITY, dim ) ;
      std::string const* oldVector = VECTOR ;
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
      for( size_t i=LENGTH ; i<dim ; i++ )
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
stringVector:: set( std::string const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: set" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<LENGTH ; i++ )
   {
      VECTOR[i] = val ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
void
stringVector:: append( std::string const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: append" ) ;
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
stringVector:: extend( std::string const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: extend" ) ;
   PEL_SAVEOLD( bool, has, has(val) ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t i = index_of( val ) ;
   if( i >= LENGTH )
   {
      append( val ) ;
   }

   PEL_CHECK_POST( has( val ) ) ;
   PEL_CHECK_POST( IMPLIES( OLD(has), size()==OLD(size) ) ) ;
   PEL_CHECK_POST( IMPLIES( OLD(has), capacity() == OLD(capacity) ) ) ;
   PEL_CHECK_POST(
         IMPLIES( !OLD(has),
                   size()==OLD(size)+1 && operator()(size()-1)==val ) ) ;
   PEL_CHECK_POST(
      IMPLIES( !OLD(has) && ( OLD(size)+1 <= OLD(capacity) ),
               capacity() == OLD( capacity ) ) ) ;
   PEL_CHECK_POST(
      IMPLIES( !OLD(has) && ( OLD(size)+1 > OLD(capacity) ),
               capacity() == PEL_System::new_block_size( OLD(capacity),
                                                         OLD(size)+1 ) ) ) ;
}

//----------------------------------------------------------------------
void
stringVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: remove_at" ) ;
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
void
stringVector:: sort( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: sort" ) ;

   // Tri "à bulles"
   for( size_t i=0 ; i<size() ; ++i )
   {
      for( size_t j=i+1 ; j<size() ; ++j )
      {
         if( operator()(i).compare( operator()(j) )>0 )
         {
            std::string dummy = operator()(i) ;
            operator()(i) = operator()(j) ;
            operator()(j) = dummy ;
         }
      }
   }
}

//----------------------------------------------------------------------
stringVector stringVector:: sorted_copy( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringVector:: sorted_copy" ) ;

   stringVector result( *this ) ;
   result.sort() ;

   return( result ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, stringVector const& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, stringVector const& )" ) ;

   size_t j=0 ;
   out << "< " ;
   for( size_t i = 0 ; i<vec.size() ; i++ )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << "\"" << vec( i ) << "\" " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, stringVector& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator>>( std::ostream&, stringVector const& )" ) ;

   stringVector aux( 0 ) ;

   std::string c ;

   if( !( in >> c ) || c != "<" )
   {
      stringVector_ERROR::n0() ;
   }
   while( in >> c && c != ">" )
   {
      if( c[0] != '\"' || c[c.size()-1] != '\"' )
      {
         stringVector_ERROR::n0() ;
      }
      std::string c2( c, 1, c.size()-2 ) ;
      aux.append( c2 ) ;
   }
   if( c != ">" )
   {
      stringVector_ERROR::n0() ;
   }
   in.clear() ;

   vec = aux ;

   return( in ) ;
}

//----------------------------------------------------------------------
void
stringVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   PEL_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new std::string [a_size] ;
      PEL_CHECK( VECTOR != 0 ) ;
      for( size_t i=0 ; i<a_size ; ++i )
      {
         VECTOR[i] = "" ;
      }
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
stringVector:: deallocate( void )
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
stringVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   PEL_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void
stringVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"string\":\n"
      "    < val1 val2 ... > is expected." ) ;
}

