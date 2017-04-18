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

#include <boolVector.hh>

#ifdef OUTLINE
#define inline
#include <boolVector.icc>
#undef inline
#endif

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_System.hh>

#include <iostream>

struct boolVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
boolVector:: boolVector( size_t dim, bool val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "boolVector:: boolVector( size_t )" ) ;
   
   allocate( LENGTH ) ;
   set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
boolVector:: boolVector( int dim, bool val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "boolVector:: boolVector( int )" ) ;
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
boolVector:: boolVector( boolVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   PEL_LABEL( "boolVector:: boolVector( boolVector const& )" ) ;
   
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
boolVector:: ~boolVector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: ~boolVector" ) ;
   PEL_CHECK_INV( invariant() ) ;   
   deallocate() ;
}

//----------------------------------------------------------------------
boolVector const&
boolVector:: operator=( boolVector const& other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: operator=" ) ;
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
   boolVector const& result = *this ;
  
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
boolVector:: re_initialize( size_t dim, bool val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: re_initialize" ) ;
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
boolVector:: operator==( boolVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = ( LENGTH == other.LENGTH ) ;
      for( size_t i=0 ; result && i<size() ; ++i )
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
boolVector:: operator!=( boolVector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   PEL_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
boolVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   PEL_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
boolVector:: has( bool val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: has" ) ;
   bool result = false ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      if ( VECTOR[i] == val )
      {
         result = true ;
         break ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
void
boolVector:: resize( size_t dim, bool val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: resize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = PEL_System::new_block_size( CAPACITY, dim ) ;
      bool const* oldVector = VECTOR ;
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
boolVector:: set( bool val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: set" ) ;
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
boolVector:: append( bool val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: append" ) ;
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
boolVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "boolVector:: remove_at" ) ;
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
std::ostream&
operator<<( std::ostream& out, boolVector const& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, boolVector const& )" ) ;
   
   size_t j=0 ;
   out << "< " ;
   for( size_t i = 0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << ( vec( i ) ? "true" : "false" ) << " " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, boolVector& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator>>( std::ostream&, boolVector const& )" ) ;

   boolVector aux( 0 ) ;
   
   std::string c ;
   
   if( !( in >> c ) || c != "<" )
   {
      boolVector_ERROR::n0() ;
   }
   while( in >> c && c != ">" )
   {
      if( c == "true" )
      {
         aux.append( true ) ;
      }
      else if( c == "false" )
      {
         aux.append( false ) ;
      }
      else
      {
         boolVector_ERROR::n0() ;
      }
   }
   if( c != ">" )
   {
      boolVector_ERROR::n0() ;
   }
   in.clear() ;

   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
boolVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   PEL_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new bool [a_size] ;
      PEL_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}


//----------------------------------------------------------------------
void
boolVector:: deallocate( void )
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
boolVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   PEL_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
boolVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"bool\":\n"
      "    < val1 val2 ... > is expected." ) ;
}
