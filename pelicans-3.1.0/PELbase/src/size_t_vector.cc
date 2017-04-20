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

#include <size_t_vector.hh>

#ifdef OUTLINE
#define inline
#include <size_t_vector.icc>
#undef inline
#endif

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_System.hh>
#include <intVector.hh>

#include <iostream>

struct size_t_vector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( size_t dim, size_t val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   PEL_LABEL( "size_t_vector:: size_t_vector( size_t )" ) ;

   allocate( LENGTH ) ;
   set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == dim ) ;
   PEL_CHECK_POST( capacity() == dim ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( int dim, size_t val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( static_cast<size_t>( dim ) )
   , CAPACITY( 0 )
{
   PEL_LABEL( "size_t_vector:: size_t_vector( int )" ) ;
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
size_t_vector:: size_t_vector( intVector const& ivec )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( ivec.size() )
   , CAPACITY( 0 )
{
   PEL_LABEL( "size_t_vector:: size_t_vector( intVector )" ) ;
   PEL_CHECK_PRE( FORALL( (size_t i=0 ;i<ivec.size();++i ),
                          ivec(i)>=0 ) ) ;

   allocate( ivec.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = (size_t)ivec(i) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == ivec.size() ) ;
   PEL_CHECK_POST( capacity() == ivec.capacity() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == (size_t)ivec(i) ) ) ;
}

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( size_t_vector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   PEL_LABEL( "size_t_vector:: size_t_vector( size_t_vector const& )" ) ;

   allocate( other.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = other.VECTOR[i] ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
size_t_vector:: ~size_t_vector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: ~size_t_vector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
size_t_vector const&
size_t_vector:: operator=( size_t_vector const& other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: operator=" ) ;
   PEL_CHECK_INV( invariant() ) ;

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
   size_t_vector const& result = *this ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result == other ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
size_t_vector:: re_initialize( size_t dim, size_t val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: re_initialize" ) ;
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
size_t_vector:: operator==( size_t_vector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = (LENGTH == other.LENGTH) ;
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
size_t_vector:: operator!=( size_t_vector const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   PEL_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: capacity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   PEL_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: index_of( size_t val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: index_of" ) ;
   
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; ++i )
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
size_t_vector:: has( size_t val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: has" ) ;
   
   size_t ifound = index_of( val ) ;
   bool result = ( ifound != PEL::bad_index() ) ;

   PEL_CHECK_POST( IMPLIES( result, operator()(index_of(val))==val ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: sum( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: sum" ) ;
   size_t resu = 0 ;
   for( size_t i = 0 ; i<LENGTH ; ++i )
   {
      resu += operator()( i ) ;
   }
   return resu ;
}

//----------------------------------------------------------------------
void
size_t_vector:: resize( size_t dim, size_t val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: resize" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;
    
   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = PEL_System::new_block_size( CAPACITY, dim ) ;
      size_t const* oldVector = VECTOR ;
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
size_t_vector:: set( size_t val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: set" ) ;
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
size_t_vector:: append( size_t val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: append" ) ;
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
size_t_vector:: extend( size_t val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: extend" ) ;
   PEL_SAVEOLD( bool, has, has(val) ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   PEL_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t i = index_of( val ) ;
   if( i == PEL::bad_index() )
   {
      append( val ) ;
   }

   PEL_CHECK_POST( has( val ) ) ;
   PEL_CHECK_POST( IMPLIES( OLD(has), size() == OLD(size) ) ) ;
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
size_t_vector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: remove_at" ) ;
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
size_t_vector:: sort_increasingly( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "size_t_vector:: sort_increasingly" ) ;

   // BUBBLE SORT
   bool modif = true ;
   while(modif)
   {
      modif = false ;
      for( size_t i=0 ; i<LENGTH-1 ; ++i )
      {
         if(VECTOR[i]>VECTOR[i+1])
         {
            size_t tmp = VECTOR[i] ;
            VECTOR[i] = VECTOR[i+1] ;
            VECTOR[i+1] = tmp ;
            modif = true ;
         }
      }
   }

   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<size()-1 ; ++i ),
            operator()(i) <= operator()(i+1) ) ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, size_t_vector const& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, size_t_vector const& )" ) ;
   
   size_t j=0 ;
   out << "< " ;
   for( size_t i = 0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << vec( i ) << " " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, size_t_vector& vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator>>( std::istream&, size_t_vector& )" ) ;

   size_t_vector aux( 0 ) ;
   
   std::string c ;
   int val = 0 ;
   
   if( !( in >> c ) || c != "<" )
   {
      size_t_vector_ERROR::n0() ;
   }
   while( in >> val )
   {
      if( val<0 ) size_t_vector_ERROR::n0() ;
      aux.append( (size_t) val ) ;
   }
   in.clear() ;
   if( !( in >> c ) || c != ">" )
   {
      size_t_vector_ERROR::n0() ;
   }
   
   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
size_t_vector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   PEL_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new size_t [a_size] ;
      PEL_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
size_t_vector:: deallocate( void )
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
size_t_vector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   PEL_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
size_t_vector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"size_t\":\n"
      "    < val1 val2 ... > is expected." ) ;
}


