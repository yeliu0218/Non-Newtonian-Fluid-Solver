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

#include <doubleArray5D.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <doubleArray5D.icc>
#undef inline
#endif

//----------------------------------------------------------------------
doubleArray5D:: doubleArray5D( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, size_t dim4, double val )
//----------------------------------------------------------------------
   : vector( dim0*dim1*dim2*dim3*dim4, val )
   , d0( dim0 )
   , d1( dim1 )
   , d2( dim2 )
   , d3( dim3 )
   , d4( dim4 )
{
   PEL_LABEL( "doubleArray5D:: doubleArray5D" ) ;
   PEL_CHECK_POST( index_bound(0) == dim0 ) ;
   PEL_CHECK_POST( index_bound(1) == dim1 ) ;
   PEL_CHECK_POST( index_bound(2) == dim2 ) ;
   PEL_CHECK_POST( index_bound(3) == dim3 ) ;
   PEL_CHECK_POST( index_bound(4) == dim4 ) ;
#ifndef _WIN32
   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                      operator()(i0,i1,i2,i3,i4) == val )))))) ;
#endif
}

//----------------------------------------------------------------------
doubleArray5D:: ~doubleArray5D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
doubleArray5D:: re_initialize( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, size_t dim4, double val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleArray5D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;
   d2 = dim2 ;
   d3 = dim3 ;
   d4 = dim4 ;
   vector.re_initialize( d0*d1*d2*d3*d4, val ) ;

   PEL_CHECK_POST( index_bound(0) == dim0 ) ;
   PEL_CHECK_POST( index_bound(1) == dim1 ) ;
   PEL_CHECK_POST( index_bound(2) == dim2 ) ;
   PEL_CHECK_POST( index_bound(3) == dim3 ) ;
   PEL_CHECK_POST( index_bound(4) == dim4 ) ;
#ifndef _WIN32
   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                      operator()(i0,i1,i2,i3,i4) == val )))))) ;
#endif
}

//----------------------------------------------------------------------
void
doubleArray5D:: raise_first_index_bound( size_t dim0, double val ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleArray5D:: raise_first_index_bound" ) ;
   PEL_CHECK_PRE( dim0 > index_bound(0) ) ;
   PEL_CHECK_PRE( index_bound(1) > 0 &&
                  index_bound(2) > 0 &&
                  index_bound(3) > 0 &&
                  index_bound(4) > 0 ) ;
   PEL_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;
   vector.resize( d0*d1*d2*d3*d4, val ) ;

   PEL_CHECK_POST( index_bound(0)==dim0 ) ;
#ifndef _WIN32
   PEL_CHECK_POST( 
         FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                      operator()(i0,i1,i2,i3,i4) == val )))))) ;
#endif
}

//----------------------------------------------------------------------
size_t
doubleArray5D:: index_bound( size_t an_index ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleArray5D:: index_bound" ) ;
   PEL_CHECK_PRE( an_index < 5 ) ;
   size_t result ;
   
   switch( an_index )
   {
      case 0 : 
         result = ( d0 ) ;
         break ;
      case 1 : 
         result = ( d1 ) ;
         break ;
      case 2 : 
         result = ( d2 ) ;
         break ;
      case 3 : 
         result = ( d3 ) ;
         break ;
      default : 
         result = ( d4 ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
doubleArray5D:: set( double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "doubleArray5D:: set" ) ;

   vector.set( val ) ;

#ifndef _WIN32
   PEL_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                      operator()(i0,i1,i2,i3,i4) == val )))))) ;
#endif
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, doubleArray5D const& a )
//----------------------------------------------------------------------
{
   PEL_LABEL( "operator<<( std::ostream&, doubleArray5D const& )" ) ;
   out << a.d0 << " " << a.d1 << " " << a.d2
       << " " << a.d3 << " " << a.d4 << " : " << std::endl ;
   
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         for( size_t k = 0 ; k<a.d2 ; k++ )
         {
            for( size_t l = 0 ; l<a.d3 ; l++ )
            {
               for( size_t m = 0 ; m<a.d4 ; m++ )
               {
                  out << i << " " << j << " " << k
                      << " " << l << " " << m
                      << " " << a( i, j, k, l, m ) << std::endl ;
               }
            }
         }
      }
   }
   return( out ) ;
}
