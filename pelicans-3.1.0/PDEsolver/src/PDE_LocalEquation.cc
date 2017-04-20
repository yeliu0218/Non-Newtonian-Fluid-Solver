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

#include <PDE_LocalEquation.hh>

#ifdef OUTLINE
   #define inline
   #include <PDE_LocalEquation.icc>
   #undef inline
#endif

#include <iostream>

//----------------------------------------------------------------------
PDE_LocalEquation*
PDE_LocalEquation:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: create" ) ;
   PDE_LocalEquation* result = new PDE_LocalEquation( a_owner ) ;

   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_rows() == 0 ) ;
   PEL_CHECK_POST( result->nb_columns() == 0 ) ;
   PEL_CHECK_POST( result->nb_row_sub_indices() == 0 ) ;
   PEL_CHECK_POST( result->nb_column_sub_indices() == 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LocalEquation:: PDE_LocalEquation( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , nbRows( 0 )
   , nbCols( 0 )
   , nbComps1( 0 )
   , nbComps2( 0 )
   , A_IS_SET( 0 )
   , A( 0 )
   , b( 0 )
   , SIZE_A( 0 )
   , SIZE_b( 0 )
   , loc2node1( 0 )
   , loc2node2( 0 )
{
}

//----------------------------------------------------------------------
PDE_LocalEquation:: ~PDE_LocalEquation( void )
//----------------------------------------------------------------------
{
   if( A != 0 )
   {
      delete [] A ;
      A = 0 ;
   }
   if( A_IS_SET != 0 )
   {
      delete [] A_IS_SET ;
      A_IS_SET = 0 ;
   }
   if( b != 0 )
   {
      delete [] b ;
      b = 0 ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: initialize( size_t_vector const& row_node_connectivity,
                                size_t a_nb_row_sub_indices,
                                size_t_vector const& column_node_connectivity,
                                size_t a_nb_col_sub_indices )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: initialize" ) ;

   nbComps1 = a_nb_row_sub_indices ;
   nbComps2 = a_nb_col_sub_indices ;

   nbRows = row_node_connectivity.size() ;
   nbCols = column_node_connectivity.size() ;

   loc2node1 = row_node_connectivity ;
   loc2node2 = column_node_connectivity ;
   
   size_t old_size_A = SIZE_A ;
   SIZE_A = nbRows * nbCols * nbComps1 * nbComps2 ;
   if( SIZE_A > old_size_A )
   {
      if( A != 0 ) delete [] A ;
      A = new double[ SIZE_A ] ;
      if( A_IS_SET != 0 ) delete [] A_IS_SET ;
      A_IS_SET = new bool[ SIZE_A ] ;
   }
   
   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      A[ i ] = 0.0 ;
      A_IS_SET[ i ] = false ;
   }
   
   size_t old_size_b = SIZE_b ;
   SIZE_b = nbRows * nbComps1 ;
   if( SIZE_b > old_size_b )
   {
      if( b != 0 ) delete [] b ;
      b = new double[ SIZE_b ] ;
   }
   
   for( size_t i=0 ; i<SIZE_b ; ++i )
   {
      b[ i ] = 0.0 ;
   }

   PEL_CHECK_POST( nb_rows() == row_node_connectivity.size() ) ;
   PEL_CHECK_POST( nb_columns() == column_node_connectivity.size() ) ;
   PEL_CHECK_POST( nb_row_sub_indices() == a_nb_row_sub_indices ) ;
   PEL_CHECK_POST( nb_column_sub_indices() == a_nb_col_sub_indices ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
                      row_node(i) == row_node_connectivity(i) ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_columns() ; i++ ),
                      column_node(i) == column_node_connectivity(i) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
            vector_item(i,ic) == 0.0 ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                  !matrix_item_is_set(i,j,ic,jc) ) ) ) ) ) ;
#ifndef _WIN32
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                  matrix_item(i,j,ic,jc) == 0.0 ) ) ) ) ) ;
#endif
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: set_as_transpose( PDE_LocalEquation const* other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: set_as_transpose" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   
   nbComps1 = other->nbComps2 ;
   nbComps2 = other->nbComps1 ;
   
   nbRows = other->nbCols ;
   nbCols = other->nbRows ;
   
   loc2node1 = other->loc2node2 ;
   loc2node2 = other->loc2node1 ;
   
   if( other->SIZE_A > SIZE_A )
   {
      if( A != 0 ) delete [] A ;
      A = new double[ other->SIZE_A ] ;
      if( A_IS_SET != 0 ) delete [] A_IS_SET ;
      A_IS_SET = new bool[ other->SIZE_A ] ;
   }
   SIZE_A = other->SIZE_A ;
   
   for( size_t i=0 ; i<nbRows ; ++i )
   {
      for( size_t j=0 ; j<nbCols ; ++j )
      {
         for( size_t ic=0 ; ic<nbComps1 ; ++ic )
         {
            for( size_t jc=0 ; jc<nbComps2 ; ++jc )
            {
               size_t idx1 = mat_idx( i, j, ic, jc ) ;
               size_t idx2 = other->mat_idx( j, i, jc, ic ) ;
               if( other->A_IS_SET[ idx2 ] )
               {
                  A[ idx1 ] = other->A[ idx2 ] ;
                  A_IS_SET[ idx1 ] = true ;
               }
               else
               {
                  A_IS_SET[ idx1 ] = false ;
               }
            }
         }
      }
   }
   
   size_t old_size_b = SIZE_b ;
   SIZE_b = nbRows * nbComps1 ;
   if( SIZE_b > old_size_b )
   {
      if( b != 0 ) delete [] b ;
      b = new double[ SIZE_b ] ;
   }
   
   for( size_t i=0 ; i<SIZE_b ; ++i )
   {
      b[ i ] = 0.0 ;
   }
   
   PEL_CHECK_POST( nb_rows() == other->nb_columns() ) ;
   PEL_CHECK_POST( nb_columns() == other->nb_rows() ) ;
   PEL_CHECK_POST( nb_row_sub_indices() ==  other->nb_column_sub_indices() ) ;
   PEL_CHECK_POST( nb_column_sub_indices() == other->nb_row_sub_indices() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
                      row_node(i) == other->column_node(i) ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_columns() ; i++ ),
                           column_node(i) == other->row_node(i) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
            vector_item(i,ic) == 0.0 ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                  matrix_item(i,j,ic,jc) == 
                  other->matrix_item(j,i,jc,ic) ) ) ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   for ( size_t i=0 ; i<nbRows ; i++ )
   {
      os << space ;
      for ( size_t j=0 ; j<nbCols ; j++ )
      {  
         os << "(" << i << ";" << j << ")" 
            << "   (" << loc2node1(i) << ";" << loc2node2(j) << ")" 
            << std::endl << space;
         for ( size_t iComp1=0 ; iComp1<nbComps1 ; iComp1++ )         {
            for ( size_t iComp2=0 ; iComp2<nbComps2 ; iComp2++ )
            {
               os << "   " ;
               if( matrix_item_is_set(i,j,iComp1,iComp2) )
               {
                  os << matrix_item(i,j,iComp1,iComp2) ;
               }
               else
               {
                  os << "not_set" ;
               }
            }
            os << std::endl << space ;
         }
      }  
      os << "Second membre : " << std::endl ;
      for ( size_t iComp1=0 ; iComp1<nbComps1 ; iComp1++ )
      {  
        os << space << "   " <<  vector_item(i,iComp1) << std::endl ;
      }
      os << std::endl <<  std::endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: scale_matrix( double coeff )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: scale_matrix" ) ;

   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      if( A_IS_SET[ i ] ) A[ i ] *= coeff ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: scale_vector( double coeff )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: scale_vector" ) ;
   
   for( size_t i=0 ; i<SIZE_b ; ++i )
   {
      b[ i ] *= coeff ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: nullify_then_unset_small_matrix_items( double relative_tolerance )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: nullify_then_unset_small_matrix_items" ) ;
   PEL_CHECK_PRE( relative_tolerance > 0.0 ) ;

   double val_max = 0. ;
   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      if( A_IS_SET[ i ] )
      {
         val_max = PEL::max( val_max, PEL::abs( A[ i ] ) ) ;
      }
   }

   double const epsilon = relative_tolerance*val_max ;
   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      if( A_IS_SET[ i ] && PEL::abs( A[ i ] )<epsilon )
      {
         A_IS_SET[ i ] = false ;
         A[ i ] = 0.0 ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: nullify_then_unset_all_matrix_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: nullify_then_unset_all_matrix_items" ) ;

   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      A_IS_SET[ i ] = false ;
      A[ i ] = 0.0 ;
   }
   
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                       !matrix_item_is_set(i,j,ic,jc) ) ) ) ) ) ;
#ifndef _WIN32
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                  matrix_item(i,j,ic,jc) == 0.0 ) ) ) ) ) ;
#endif
}

//----------------------------------------------------------------------
void
PDE_LocalEquation:: mark_all_matrix_items_as_set( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalEquation:: mark_all_matrix_items_as_set" ) ;

   for( size_t i=0 ; i<SIZE_A ; ++i )
   {
      A_IS_SET[ i ] = true ;
   }
   
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_rows() ; i++ ),
         FORALL( ( size_t j=0 ; j<nb_columns() ; j++ ),
            FORALL( ( size_t ic=0 ; ic<nb_row_sub_indices() ; ic++ ),
               FORALL( ( size_t jc=0 ; jc<nb_column_sub_indices() ; jc++ ),
                       matrix_item_is_set(i,j,ic,jc) ) ) ) ) ) ;
}
