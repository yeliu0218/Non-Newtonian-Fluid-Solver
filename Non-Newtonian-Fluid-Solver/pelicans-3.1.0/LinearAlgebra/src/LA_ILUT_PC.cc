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

#include <LA_ILUT_PC.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqImplementation.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <sstream>

struct LA_ILUT_PC_ERROR
{
   static void n0( int i ) ;
} ;

LA_ILUT_PC const* LA_ILUT_PC::PROTOTYPE = new LA_ILUT_PC() ;

//----------------------------------------------------------------------
LA_ILUT_PC:: LA_ILUT_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_ILUT_PC" )
   , RULE( -PEL::max_double() )
   , LFIL( PEL::bad_index() )
   , MODIFIED( false )
   , PIV_MIN( -PEL::max_double() )
   , LU( 0 )
   , COL_P_LARGEST( 0 )
   , VAL_P_LARGEST( 0 )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILUT_PC*
LA_ILUT_PC:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double tolerance = exp->double_data( "tolerance" ) ;
   if( tolerance<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "tolerance", "greater than zero" ) ;
   }
   int order = exp->int_data( "order" ) ;
   if( order<0 )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "order", "greater than zero" ) ;
   }
   bool diagonal_compensation = exp->bool_data( "diagonal_compensation" ) ;
   double smallest_nonzero_pivot = exp->double_data( "smallest_nonzero_pivot" ) ;
   if( smallest_nonzero_pivot<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "smallest_nonzero_pivot", "greater than zero" ) ;
   }
   PEL_ModuleExplorer const* se =
                       exp->create_subexplorer( 0, "LA_SeqMatrix" ) ;
   LA_Matrix const* matrix_prototype = LA_Matrix::make( 0, se ) ;
   if( matrix_prototype->implementation() != LA_SeqImplementation::object() )
   {
      PEL_Error::object()->raise_plain( "LA_ILUT_PC:: create_replica::invalid" ) ;
   }
   PEL_CHECK( dynamic_cast<LA_SeqMatrix const*>( matrix_prototype ) != 0 ) ;
   LA_SeqMatrix const* smatrix_prototype = static_cast<LA_SeqMatrix const*>( matrix_prototype ) ;
   
   LA_ILUT_PC* result = new LA_ILUT_PC( a_owner,
                                        (size_t) order, tolerance,
                                        diagonal_compensation,
                                        smallest_nonzero_pivot,
                                        smatrix_prototype ) ;

   matrix_prototype->destroy() ; matrix_prototype = 0 ;
   se->destroy() ; se = 0 ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ILUT_PC* 
LA_ILUT_PC:: create( PEL_Object* a_owner,
                     size_t order,
                     double tolerance,
                     bool diagonal_compensation,
                     double smallest_nonzero_pivot,
                     LA_SeqMatrix const* matrix_prototype )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: create" ) ;
   PEL_CHECK_PRE( tolerance > 0 ) ;
   PEL_CHECK_PRE( smallest_nonzero_pivot > 0 ) ;
   PEL_CHECK_PRE( matrix_prototype != 0 ) ;

   LA_ILUT_PC* result = new LA_ILUT_PC( a_owner, 
                                        order, 
                                        tolerance, 
                                        diagonal_compensation,
                                        smallest_nonzero_pivot,
                                        matrix_prototype ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
LA_ILUT_PC:: LA_ILUT_PC( PEL_Object* a_owner,
                         size_t order,
                         double tolerance,
                         bool diagonal_compensation,
                         double smallest_nonzero_pivot,
                         LA_SeqMatrix const* matrix_prototype )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , RULE( tolerance )
   , LFIL( order )
   , MODIFIED( diagonal_compensation )
   , PIV_MIN( smallest_nonzero_pivot )
   , LU( matrix_prototype->create_matrix( this ) )
   , COL_P_LARGEST(0)
   , VAL_P_LARGEST(0)
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILUT_PC:: ~LA_ILUT_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_ILUT_PC*
LA_ILUT_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: create_clone" ) ;

   LA_ILUT_PC* result =
          new LA_ILUT_PC( a_owner, LFIL, RULE, MODIFIED, PIV_MIN, LU ) ;

   PEL_CHECK_POST( create_clone_POST( result,  a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_ILUT_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: is_valid" ) ;
   return( BUILD_OK ) ;
}

//----------------------------------------------------------------------
size_t
LA_ILUT_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( LU->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUT_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;
   
   if( dynamic_cast<LA_SeqMatrix const*>(mat) == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_ILUT_PC error:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected" ) ;
   }
   
   size_t const length = mat->nb_rows() ;
   
   LU->re_initialize( length, length ) ;
   LU->set( mat ) ;
   COL_P_LARGEST.re_initialize( 2*length+1 ) ;
   VAL_P_LARGEST.re_initialize( 2*length+1 ) ;
   
   LA_MatrixIterator* itI = LU->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itK = LU->create_stored_item_iterator( 0 ) ;
   
   LA_SeqVector* diag = LA_SeqVector::create( 0, length ) ;
   double d00 = LU->item(0,0) ;
   if( PEL::abs( d00 )<PIV_MIN )
   {
      d00 = 1. ;
      LU->nullify_row( 0 ) ;
      LU->set_item( 0, 0, d00 ) ;
      LA_ILUT_PC_ERROR::n0( 0 ) ;
   }
   diag->set_item( 0, d00 ) ;
   
   for( size_t i=1 ; i<length ; ++i ) 
   {
      //Compute the 1-norm of the i-th line
      double norm1 = 0. ;
      size_t nb = 0 ;
      size_t lowEl = 0 ;
      size_t upEl = 0 ;
      for( itI->start_row_items(i) ; itI->is_valid() ; itI->go_next() )
      {
         norm1 += PEL::abs( itI->item() ) ;
         size_t j = itI->col() ;
         if( i>j )
         {
            lowEl++ ;
         }
         else if( j>i )
         {
            upEl++ ;
         }
         nb++ ;
      }
      norm1 /= nb ;
      double relativeEps = norm1*RULE ;
      double drop = 0.0 ;
      for( itI->start_row_items(i) ; itI->is_valid() ; itI->go_next() )
      {
         size_t k = itI->col() ;
         if( k>=i ) break ;
         PEL_CHECK( k<i ) ;
         // aik <- aik/akk ;
         double aik =  itI->item()/diag->item( k ) ;
         itI->set_item( aik ) ;
         // rule on aik :
         if( PEL::abs(aik)>relativeEps )
         {
            for( itK->start_row_items(k) ; itK->is_valid() ; itK->go_next() )
            {
               size_t j = itK->col() ;
               if( j>k )
               {
                  double akj = itK->item() ;
                  LU->add_to_item( i, j, -aik*akj ) ;
               }
            }
         }
      }
      LU->synchronize() ;
      drop += keepLUplargest( i, relativeEps, lowEl, upEl ) ;
      double d = LU->item( i, i ) ;
      if( MODIFIED )
      {
         LU->add_to_item( i, i, -drop ) ;
         d = LU->item( i, i ) ;
      }
      if( PEL::abs( d )<PIV_MIN ) 
      {
         d = 1.0 ;
         LU->nullify_row( i ) ;
         LU->set_item( i, i, d ) ;
         LA_ILUT_PC_ERROR::n0( i ) ;
      }
      diag->set_item( i, d ) ;
   }

   diag->destroy() ; diag = 0 ;
   itI->destroy() ; itI = 0 ;
   itK->destroy() ; itK = 0 ;

   LU->synchronize() ;
   
   BUILD_OK = true ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUT_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_ILUT_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   BUILD_OK = false ;
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUT_PC:: solve( LA_Vector const* rhs,
                    LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( rhs ) != 0 ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( sol ) != 0 ) ;

   LA_SeqVector const* brhs = static_cast<LA_SeqVector const*>( rhs ) ;
   LA_SeqVector* bsol = static_cast<LA_SeqVector*>( sol ) ;
   
   LU->solve_LU( brhs, bsol ) ;
   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;   
}

//----------------------------------------------------------------------
bool
LA_ILUT_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
double
LA_ILUT_PC:: keepLUplargest( int i, double norm, int lowEl, int upEl )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: keepLUplargest" ) ;
   PEL_CHECK_PRE( LU->is_synchronized() ) ;
   
   int nbL = lowEl+LFIL ;
   int nbU = upEl+LFIL ;
   int med = nbL ;
   int nbTot = nbL+nbU+1 ;
   double drop = 0.0 ;
   
   for( int j=0 ; j<nbTot ; j++ )
   {
      COL_P_LARGEST( j ) = -1 ;
   }
   
   
   // COL_P_LARGEST.set( -1 ) ;
   int compteurLow = 0 ;
   int compteurUp = 0 ;
   LA_MatrixIterator* itI = LU->create_stored_item_iterator( 0 ) ;
   for( itI->start_row_items(i) ; itI->is_valid() ; itI->go_next() )
   {
      double x = itI->item() ;
      int k = itI->col() ;
      if( k==i )
      {
         COL_P_LARGEST( med ) = i ;
         VAL_P_LARGEST( med ) = x ;
      }
      else
      {
         if( PEL::abs(x) > norm )
         {
            if( k<i )
            {
               if( compteurLow<nbL )
               {
                  COL_P_LARGEST( compteurLow ) =  k ;
                  VAL_P_LARGEST( compteurLow ) =  x ;
                  compteurLow++ ;
               }
               else
               {
                  // seach the inf of valPlargest
                  int ip = 0 ;
                  int idx = 0 ;
                  double minLow = 1.E30 ;
                  while( ip<nbL && COL_P_LARGEST(ip)!=-1 )
                  {
                     if( PEL::abs(VAL_P_LARGEST(ip))<minLow )
                     {
                        minLow = PEL::abs(VAL_P_LARGEST(ip)) ;
                        idx = ip ;
                     }
                     ip++ ;
                  }
                  if( PEL::abs(x)>minLow )
                  {
                     drop += VAL_P_LARGEST( idx ) ;
                     COL_P_LARGEST(idx ) = k ;
                     VAL_P_LARGEST(idx ) = x ;
                  }
               }
            }
            else if( k>i )
            {
               if( compteurUp<nbU )
               {
                  COL_P_LARGEST( med+1+compteurUp ) = k ;
                  VAL_P_LARGEST( med+1+compteurUp ) = x ;
                  compteurUp++ ; 
               }
               else
               {
                  int ip = med+1 ;
                  int idx = med+1 ;
                  double minUp = 1.E30 ;
                  while( ip<nbTot && COL_P_LARGEST(ip)!=-1 )
                  {
                     if( PEL::abs(VAL_P_LARGEST(ip))<minUp )
                     {
                        minUp = PEL::abs(VAL_P_LARGEST(ip)) ;
                        idx = ip ;
                     }
                     ip++ ;
                  }
                  if( PEL::abs(x)>minUp )
                  {
                     drop += VAL_P_LARGEST( idx ) ;
                     COL_P_LARGEST( idx ) = k ;
                     VAL_P_LARGEST( idx ) = x ;
                  }
               }
            }
         }
         else
         {
            drop += x ;
         }
      }
   }
   
   LU->nullify_row( i ) ;
   
   for( int ip=0 ; ip<nbTot ; ip++ )
   {
      if( COL_P_LARGEST(ip)!=-1 )
      {
         LU->add_to_item( i, 
                         COL_P_LARGEST(ip), 
                         VAL_P_LARGEST(ip) ) ;
      }
   }

   itI->destroy() ;
   itI = 0 ;
   return( drop ) ;
}

//----------------------------------------------------------------------
void
LA_ILUT_PC:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUT_PC:: print_more" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string s( indent_width, ' ') ;

   os << s << "smallest_nonzero_pivot : " << PIV_MIN << std::endl ;
   os << s << "diagonal_compensation  : "
      << ( MODIFIED ? "true" : "false" ) << std::endl ;
   os << s << "order : " << LFIL << std::endl ;
   os << s << "tolerance : " << RULE << std::endl ;
}

//internal--------------------------------------------------------------
void 
LA_ILUT_PC_ERROR:: n0( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** ILUT factorization :" << std::endl ;
   mesg << "***   zero pivot at step " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
