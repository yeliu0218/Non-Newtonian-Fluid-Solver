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

#include <GE_Matrix.hh>

#include <memory.h>
#include <iostream>
#include <sstream>
#include <string>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <doubleArray2D.hh>

#ifdef OUTLINE
#define inline
#include <GE_Matrix.icc>
#undef inline
#endif

//----------------------------------------------------------------------
GE_Matrix*
GE_Matrix:: create( PEL_Object* a_owner,
                    size_t row_nb,
                    size_t col_nb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: create" ) ;
   PEL_CHECK_PRE( row_nb<=3 ) ;
   PEL_CHECK_PRE( col_nb<=3 ) ;

   GE_Matrix* result = new GE_Matrix( a_owner, row_nb, col_nb ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_rows() == row_nb ) ;
   PEL_CHECK_POST( result->nb_cols() ==  col_nb) ;
   PEL_CHECK_POST( !result->determinant_is_computed() ) ;

   return( result ) ;
}
   
//----------------------------------------------------------------------
GE_Matrix:: GE_Matrix( PEL_Object* a_owner,
                       size_t row_nb,
                       size_t col_nb )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     rows( row_nb ),
     cols( col_nb ),
     det( 0.0 ),
     det_comp( false )
{
   PEL_LABEL( "GE_Matrix:: GE_Matrix" ) ;
   nullify() ;
}

//----------------------------------------------------------------------
GE_Matrix:: ~GE_Matrix( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: ~GE_Matrix" ) ;
}

//----------------------------------------------------------------------
void 
GE_Matrix:: nullify( void  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: nullify" ) ;
   size_t size = 3*nb_cols() ;
   for( size_t i=0 ; i<size ; i++ )
      mat[i]=0.0 ;
   det = 0.0 ;
   det_comp = false ;
}

//----------------------------------------------------------------------
void  
GE_Matrix:: compute_determinant( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: compute_determinant" ) ;
//   PEL_CHECK_PRE( !determinant_is_computed() ) ;
   PEL_CHECK_PRE( nb_rows()==nb_cols() ) ;
   
   if( rows==1 )
   {
      det = item(0,0) ;
   }
   else if( rows==2 )
   {
      det = item(0,0)*item(1,1) - item(1,0)*item(0,1) ;
   }
   else
   {
      det = item(0,0)* ( item(1,1)*item(2,2) - item(1,2)*item(2,1) ) 
         - item(0,1)* ( item(1,0)*item(2,2) - item(1,2)*item(2,0) )
         + item(0,2)* ( item(1,0)*item(2,1) - item(1,1)*item(2,0) ) ;
      
   }
   det_comp = true ;
}

//----------------------------------------------------------------------
void  
GE_Matrix:: invert( doubleArray2D const& B,
                    doubleArray2D& X ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: gauss_multi" ) ;
   PEL_CHECK_PRE( B.index_bound(0)==X.index_bound(0) ) ;
   PEL_CHECK_PRE( B.index_bound(1)==X.index_bound(1) ) ;
   PEL_CHECK_PRE( B.index_bound(1)==nb_rows() ) ;
   PEL_CHECK_PRE( determinant_is_computed() ) ;
   PEL_CHECK_PRE( determinant()!=0 ) ;
   
   switch( rows ) 
   {
      case 1 :
         gauss1x1_multi(B,X) ;
         break ;
      case 2 :
         gauss2x2_multi(B,X) ;
         break ;
      case 3 :
         gauss3x3_multi(B,X) ;
         break ;
   }   
}

//----------------------------------------------------------------------
bool
GE_Matrix:: determinant_is_computed( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: determinant_is_computed" ) ;
   return det_comp ;
}

//----------------------------------------------------------------------
double
GE_Matrix:: determinant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: determinant" ) ;
   PEL_CHECK_PRE( determinant_is_computed() ) ;
   return det ;
}

//----------------------------------------------------------------------
void  
GE_Matrix:: gauss1x1_multi( doubleArray2D const& B,
                            doubleArray2D& X )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: gauss1x1_multi" ) ;
   size_t const nbSystems = B.index_bound(0) ;
   
   for( size_t k=0 ; k<nbSystems ; ++k )
   {
      X( k, 0 ) = B( k, 0 )/det ;
   }
}

//----------------------------------------------------------------------
void  
GE_Matrix:: gauss2x2_multi( doubleArray2D const& B,
                            doubleArray2D& X )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: gauss2x2_multi" ) ;

   size_t const nbSystems = B.index_bound(0) ;
   
   double const A00 = item( 0, 0 ) ;
   double const A10 = item( 1, 0 ) ;
   double const A11 = item( 1, 1 ) ;
   double const A01 = item( 0, 1 ) ;
   
   for( size_t k=0 ; k<nbSystems ; ++k )
   {
      double const x0 =
         B( k, 0 )*A11-B( k, 1 )*A01 ;
      double const x1 =
         A00*B( k, 1 )-A10*B( k, 0 ) ;
      X( k, 0 ) = x0/det ;
      X( k, 1 ) = x1/det ;
   }
}

//----------------------------------------------------------------------
void  
GE_Matrix:: gauss3x3_multi( doubleArray2D const& B,
                            doubleArray2D& X ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: gauss3x3_multi" ) ;

   size_t const nbSystems = B.index_bound(0) ;
   double const A00 = item( 0, 0 ) ;
   double const A10 = item( 1, 0 ) ;
   double const A20 = item( 2, 0 ) ;
   double const A01 = item( 0, 1 ) ;
   double const A11 = item( 1, 1 ) ;
   double const A21 = item( 2, 1 ) ;
   double const A02 = item( 0, 2 ) ;
   double const A12 = item( 1, 2 ) ;
   double const A22 = item( 2, 2 ) ;
   
   for( size_t k=0 ; k<nbSystems ; ++k )
   {
      double const x0 =
         B( k, 0 )*( A11*A22 - A12*A21 )
         - A01*( B( k, 1 )*A22
                 -A12*B( k, 2 ) )
         + A02*( B( k, 1 )*A21
                 -A11*B( k, 2 ) ) ;
      double const x1 =
         A00*( B( k, 1 )*A22
               -A12*B( k, 2 ) )
         - B( k, 0 )*( A10*A22 - A12*A20 )
         + A02*( A10*B( k, 2 )
                 -B( k, 1 )*A20 ) ;
      double const x2 =
         A00*( A11*B( k, 2 )
               -B( k, 1 )*A21 )
         - A01*( A10*B( k, 2 )
                 -B( k, 1 )*A20 )
         + B( k, 0 )*( A10*A21 - A11*A20 ) ;
   
      X( k, 0 ) = x0/det ;
      X( k, 1 ) = x1/det ;
      X( k, 2 ) = x2/det ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Matrix:: print( std::ostream& os, size_t indent_width ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Matrix:: print" ) ;

   std::string const space( indent_width, ' ' ) ;
   os << space << " GE_Matrix : ( " ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      os << "( " ;
      for( size_t j=0 ; j<nb_cols()-1 ; ++j )
      {
         os << item( i, j ) << " , " ;
      }
      os << item( i, nb_cols()-1 ) << " ) " ;
      if( i<nb_rows()-1 )
         os << " , " ;
      else 
         os << " ) " << std::endl ;
   }
}
