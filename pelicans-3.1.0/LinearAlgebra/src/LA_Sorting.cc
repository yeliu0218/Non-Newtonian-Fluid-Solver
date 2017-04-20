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

#include <LA_Sorting.hh>

#include <LA_DenseMatrix.hh>
#include <LA_SeqVector.hh>
#include <size_t_array2D.hh> 
#include <size_t_vector.hh>
#include <intArray2D.hh>
#include <iostream>

//----------------------------------------------------------------------
void
LA_Sorting:: exchange( LA_SeqVector* sortedValues,
                       size_t a,
                       size_t b )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Sorting:: exchange" ) ;
   
   double tmp = sortedValues->item(a) ;
   sortedValues->set_item(a, sortedValues->item(b) ) ;
   sortedValues->set_item(b, tmp ) ;
}

   
//----------------------------------------------------------------------
void
LA_Sorting:: quick_sort( LA_SeqVector* sortedValues,
                         size_t start,
                         size_t end )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Sorting:: quick_sort" ) ;
   PEL_CHECK_PRE( sortedValues!=0 ) ;
   PEL_CHECK_PRE( start<=end ) ;
   
   if( end==start+1 )
   {
      if( sortedValues->item(start)>sortedValues->item(end ))
      {
         exchange( sortedValues, start, end ) ;
      }
   }
   else if( end > start )
   {
      size_t a = start ;
      size_t b = end ;
      size_t c = (b+a)/2 ;
      double piv = sortedValues->item( c ) ;
      
      while( a<b )
      {
         double va = sortedValues->item(a) ;
         double vb = sortedValues->item(b) ;
         
         if( va<piv || ( va==vb && va==piv ) )
         {
            a++ ;
         }
         else if( vb>piv )
         {
            b--;
         }
         else
         {
            exchange( sortedValues, a, b ) ;
         }
      }
      PEL_CHECK( a==b ) ;
      PEL_CHECK( sortedValues->item(a)==piv ) ;
      PEL_CHECK( FORALL(
         ( size_t i=start;i<a;i++), 
         sortedValues->item(i)<=piv)) ;
      PEL_CHECK( FORALL(
         ( size_t i=a+1;i<=end;i++), 
         sortedValues->item(i)>=piv)) ;
      if( a>start+1 ) quick_sort( sortedValues, start, a-1 ) ;
      if( a<end-1 ) quick_sort( sortedValues, a+1, end ) ;
   }
   PEL_CHECK_POST( FORALL(
      ( size_t i=start;i<end-1;i++),
      sortedValues->item(i)<=sortedValues->item(i+1))) ;
}



//----------------------------------------------------------------------
void
LA_Sorting:: sort( LA_SeqVector const* values,
                   LA_SeqVector* sortedValues,
                   size_t_vector& rank )
//----------------------------------------------------------------------
// insertSort
{
   PEL_LABEL( "LA_Sorting:: sort" ) ;
   
   PEL_CHECK_PRE ( values->nb_rows() > 0 ) ;
   PEL_CHECK_PRE( sortedValues->nb_rows() == values->nb_rows() ) ;
   PEL_CHECK ( rank.size() == values->nb_rows()  ) ;
   const size_t nSize= values->nb_rows() ;

   size_t_vector inda( nSize ) ;
   size_t_vector indb( nSize ) ;
   
   sortedValues->set( values ) ;

   for( size_t i =0 ; i< nSize ; i++ )
   {
      inda( i ) = i ;
   }

   for( size_t i=1 ; i<nSize ; i++ )
   {
      double elem=sortedValues->item( i ) ;
      indb = inda ;
      int j ;
      for( j=i-1 ; j >= 0 ; j-- )
      {
         if( sortedValues->item( j ) > elem )
         {
            sortedValues->set_item( j+1, sortedValues->item( j ) ) ;
            inda( j+1 ) =  inda( j ) ;
         }
         else
         { 
            break ;
         }
      }
      sortedValues->set_item( j+1, elem ) ;
      inda( j+1 ) = indb( i ) ;
   }


   for( size_t i=0 ; i<nSize ; i++ )
   {
      rank( inda( i ) )= i ;
   } 
   
   PEL_CHECK_POST( FORALL(
      ( size_t i=0;i<values->nb_rows();i++),
      sortedValues->item(rank(i))==values->item(i))) ;
   PEL_CHECK_POST( FORALL(
      ( size_t i=0;i<values->nb_rows()-1;i++),
      sortedValues->item(i)<=sortedValues->item(i+1))) ;
   
}



//----------------------------------------------------------------------
void
LA_Sorting:: sort( intVector const& values,
                   intVector& sortedValues,
                   size_t_vector& rank )
//----------------------------------------------------------------------
// insertSort
{
   PEL_LABEL( "LA_Sorting:: sort" ) ;
   
   PEL_CHECK_PRE ( values.size() > 0 ) ;
   PEL_CHECK_PRE( sortedValues.size() == values.size() ) ;
   PEL_CHECK ( rank.size() == values.size()  ) ;

   const size_t nSize= values.size() ;

   size_t_vector inda( nSize ) ;
   size_t_vector indb( nSize ) ;
   

   sortedValues = values ;

   for( size_t i =0 ; i<nSize ; i++ )
   {
      inda( i ) = i ;
   }

   for( size_t i=1 ; i<nSize ; i++ )
   {
      int elem=sortedValues( i ) ;
      indb = inda ;
      int j ;
      for( j=i-1 ; j>=0 ; j-- )
      {
         if( sortedValues( j ) > elem )
         {
            sortedValues( j+1 ) = sortedValues( j ) ;
            inda( j+1 ) = inda( j ) ;
         }
         else
         {
            break ;
         }
      }
      sortedValues( j+1 ) = elem ;
      inda( j+1 ) = indb(i) ;
   }

   for( size_t i=0 ; i<nSize ; i++ )
   {
      rank( inda( i ) ) = i ;
   }

   PEL_CHECK_POST( FORALL(
      ( size_t i=0;i<values.size();i++),
      sortedValues(rank(i))==values(i))) ;
   PEL_CHECK_POST( FORALL(
      ( size_t i=0;i<values.size()-1;i++),
      sortedValues(i)<=sortedValues(i+1))) ;
}



//----------------------------------------------------------------------
void
LA_Sorting:: column_rank( LA_SeqMatrix const* mat,
                          size_t_array2D& rank )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Sorting:: column_rank" ) ;
   PEL_CHECK_PRE( mat!=0 ) ;
   PEL_CHECK_PRE( mat->nb_rows()>0 ) ;
   PEL_CHECK_PRE( mat->nb_cols()>0 ) ;

   rank.re_initialize( mat->nb_rows(), mat->nb_cols() ) ;
   
   size_t n = mat->nb_rows() ;
   LA_SeqVector* vecr = LA_SeqVector::create( 0, n ) ;
   LA_SeqVector* vSort = LA_SeqVector::create( 0, n ) ;
   size_t_vector iSort( n ) ;
   
   for( size_t j=0 ; j<mat->nb_cols() ; j++ )
   {
      mat->extract_col( j, vecr ) ;
      LA_Sorting::sort( vecr, vSort, iSort ) ;
      rank.set_section( 1, j, iSort ) ;
   }
   vecr->destroy() ;
   vSort->destroy() ;
   
   PEL_CHECK_POST( rank.index_bound(0)==mat->nb_rows() ) ;
   PEL_CHECK_POST( rank.index_bound(1)==mat->nb_cols() ) ;
   PEL_CHECK_POST( FORALL( (size_t j=0;j<mat->nb_cols();j++),
                           FORALL( (size_t i=0;i<mat->nb_rows()-1;i++),
                                   IMPLIES(
                                      rank(i,j)<=rank(i+1,j),
                                       mat->item(i,j)<=mat->item(i+1,j))))) ;
   
}

