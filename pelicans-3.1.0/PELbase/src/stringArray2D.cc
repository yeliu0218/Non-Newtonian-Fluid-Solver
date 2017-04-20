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

#include <stringArray2D.hh>

#ifdef OUTLINE
#define stringArray2D_HH
#include <stringArray2D.icc>
#undef inline
#endif

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

//----------------------------------------------------------------------
stringArray2D:: stringArray2D( size_t dim0, size_t dim1, string val )
//----------------------------------------------------------------------
   : vector( dim0*dim1 )
   , d0( dim0 )
   , d1( dim1 )
{
   PEL_LABEL( "stringArray2D:: stringArray2D" ) ;
   vector.set( val ) ;
   PEL_CHECK_POST( index_bound(0) == dim0 ) ;
   PEL_CHECK_POST( index_bound(1) == dim1 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
stringArray2D:: ~stringArray2D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
stringArray2D:: stringArray2D( stringArray2D const& other )
//----------------------------------------------------------------------
   : vector( other.vector )
   , d0( other.d0 )
   , d1( other.d1 )
{
   PEL_LABEL( "stringArray2D:: stringArray2D( stringArray2D const& )" ) ;
   PEL_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
stringArray2D:: operator=( stringArray2D const& other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: operator=" ) ;
   
   if( this != &other )
   {
      vector = other.vector ;
      d0 = other.d0 ;
      d1 = other.d1 ;
   }

   stringArray2D const& result = *this ;

   PEL_CHECK_POST( result == other ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: re_initialize( size_t dim0, size_t dim1, string val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;
   vector.re_initialize( d0*d1 ) ;
   vector.set( val ) ;

   PEL_CHECK_POST( index_bound(0) == dim0 ) ;
   PEL_CHECK_POST( index_bound(1) == dim1 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
bool
stringArray2D:: operator==( stringArray2D const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: operator==" ) ;

   bool result = ( d0==other.d0 &&
                   d1==other.d1 &&
                   vector==other.vector ) ;

   PEL_CHECK_POST(
      IMPLIES( result,
               index_bound(0)==other.index_bound(0) &&
               index_bound(1)==other.index_bound(1) ) ) ;
   PEL_CHECK_POST(
      !result ||
      FORALL( ( size_t i=0 ; i<index_bound(0) ; ++i ),
              FORALL( ( size_t j=0 ; j<index_bound(1) ; ++j ),
                      operator()(i,j) == other(i,j) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
stringArray2D:: operator!=( stringArray2D const& other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: operator!=" ) ;

   bool result = !operator==( other ) ;

   PEL_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: raise_first_index_bound( size_t dim0, string val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: raise_first_index_bound" ) ;
   PEL_CHECK_PRE( dim0 > index_bound(0) ) ;
   PEL_CHECK_PRE( index_bound(1) > 0 ) ;
   PEL_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;
   vector.resize( d0*d1, val ) ;

   PEL_CHECK_POST( index_bound(0) == dim0 ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
         FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
            operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: set( string val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: set" ) ;
   vector.set( val ) ;

   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: set_section( size_t an_index, size_t index_value, 
                            stringVector const& x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: set_section" ) ;
   PEL_CHECK_PRE( an_index < 2 ) ;
   PEL_CHECK_PRE( index_value < index_bound(an_index) ) ;
   PEL_CHECK_PRE( IMPLIES( an_index==0, x.size() == index_bound(1) ) ) ;
   PEL_CHECK_PRE( IMPLIES( an_index==1, x.size() == index_bound(0) ) ) ;
         
   if( an_index == 0 )
   {
      for( size_t k=0 ; k<d1 ; k++)
      {
         operator()(index_value,k) = x(k) ;
      }
   }
   else if( an_index == 1 )
   {
      for( size_t k=0 ; k<d0 ; k++ )
      {
         operator()(k,index_value) = x(k) ;
      }
   }

   PEL_CHECK_POST( an_index==1 ||
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                      operator()(index_value,i1) == x(i1) ) ) ;
   PEL_CHECK_POST( an_index==0 ||
                   FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                      operator()(i0,index_value) == x(i0) ) ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: extract_section( size_t an_index, size_t index_value,
                                 stringVector& x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "stringArray2D:: extract_section" ) ;
   PEL_CHECK_PRE( an_index < 2 ) ;
   PEL_CHECK_PRE( index_value < index_bound(an_index) ) ;
   PEL_CHECK_PRE( IMPLIES( an_index==0, x.size() == index_bound(1) ) ) ;
   PEL_CHECK_PRE( IMPLIES( an_index==1, x.size() == index_bound(0) ) ) ;
   
   if( an_index == 0 )
   {
      for( size_t k=0 ; k<d1 ; k++)
      {
         x(k) = operator()(index_value,k) ;
      }
   }
   else if( an_index == 1 )
   {
      for( size_t k=0 ; k<d0 ; k++ )
      {
         x(k) = operator()(k,index_value) ;
      }
   }

   PEL_CHECK_POST( an_index==1 || 
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                           x(i1) == operator()(index_value,i1) ) ) ;
   PEL_CHECK_POST( an_index==0 || 
                   FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                           x(i0) == operator()(i0,index_value) ) ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, stringArray2D const& a )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, stringArray2D const& )" ) ;

   ios_base::fmtflags original_flags = out.flags() ;
   out.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = out.precision() ;
   out << setprecision( 6 )  ;

   out << setw( 5 ) << "i" << setw( 5 ) << "j" 
       << setw( 15 ) << "item(i,j)" << endl ;
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         out << setw( 5 ) << i << setw( 5 ) << j ;
         out << setw( 15 ) << a( i, j ) << endl ;
      }
   }

   out << std::setprecision(p) ;
   out.flags( original_flags ) ;
   return( out ) ;
}
