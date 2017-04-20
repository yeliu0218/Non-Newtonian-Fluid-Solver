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

#include <LA_ILUp_PC.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqImplementation.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <sstream>

struct LA_ILUp_PC_ERROR
{
   static void n0( int i ) ;
} ;

LA_ILUp_PC const* LA_ILUp_PC::PROTOTYPE = new LA_ILUp_PC() ;

//----------------------------------------------------------------------
LA_ILUp_PC:: LA_ILUp_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_ILUp_PC" )
   , LEVEL_ORDER( PEL::bad_index() )
   , MODIFIED( false )
   , PIV_MIN( -PEL::max_double() )
   , LU( 0 )
   , LEVELS( 0 )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILUp_PC*
LA_ILUp_PC:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_ModuleExplorer const* se =
                       exp->create_subexplorer( 0, "LA_SeqMatrix" ) ;
   LA_Matrix const* matrix_prototype = LA_Matrix::make( 0, se ) ;
   if( matrix_prototype->implementation() != LA_SeqImplementation::object() )
   {
      PEL_Error::object()->raise_plain( "LA_ILU0_PC:: create_replica::invalid" ) ;
   }
   PEL_CHECK( dynamic_cast<LA_SeqMatrix const*>( matrix_prototype ) != 0 ) ;
   LA_SeqMatrix const* smatrix_prototype = static_cast<LA_SeqMatrix const*>( matrix_prototype ) ;
   
   bool diagonal_compensation = exp->bool_data( "diagonal_compensation" ) ;
   int order = exp->int_data( "order" ) ;
   if( order<0 )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "order", "greater than zero" ) ;
   }
   double smallest_nonzero_pivot = exp->double_data( "smallest_nonzero_pivot" ) ;
   if( smallest_nonzero_pivot<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "smallest_nonzero_pivot", "greater than zero" ) ;
   }

   LA_ILUp_PC* result = new LA_ILUp_PC( a_owner,
                                        (size_t) order,
                                        diagonal_compensation,
                                        smallest_nonzero_pivot,
                                        smatrix_prototype ) ;

   matrix_prototype->destroy() ; matrix_prototype = 0 ;
   se->destroy() ; se = 0 ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ILUp_PC*
LA_ILUp_PC:: create( PEL_Object* a_owner,
                     size_t order,
                     bool diagonal_compensation,
                     double smallest_nonzero_pivot,
                     LA_SeqMatrix const* matrix_prototype )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: create" ) ;
   PEL_CHECK_PRE( smallest_nonzero_pivot > 0. ) ;
   PEL_CHECK_PRE( matrix_prototype != 0 ) ;
   PEL_CHECK_PRE( matrix_prototype->nb_rows()==0 &&
                  matrix_prototype->nb_cols()==0 ) ;

   LA_ILUp_PC* result = new LA_ILUp_PC( a_owner, 
                                        order, 
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
LA_ILUp_PC:: LA_ILUp_PC( PEL_Object* a_owner,
                         size_t order,
                         bool diagonal_compensation,
                         double smallest_nonzero_pivot,
                         LA_SeqMatrix const* matrix_prototype )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , LEVEL_ORDER( order )
   , MODIFIED( diagonal_compensation )
   , PIV_MIN( smallest_nonzero_pivot )
   , LU( matrix_prototype->create_matrix( this ) )
   , LEVELS( matrix_prototype->create_matrix( this ) )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILUp_PC:: ~LA_ILUp_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_ILUp_PC*
LA_ILUp_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: create_clone" ) ;

   LA_ILUp_PC* result = new LA_ILUp_PC(
                         a_owner, LEVEL_ORDER, MODIFIED, PIV_MIN, LU ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_ILUp_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: is_valid" ) ;
   return( BUILD_OK ) ;
}

//----------------------------------------------------------------------
size_t
LA_ILUp_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( LU->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;
   
   LA_SeqMatrix const* smat = dynamic_cast<LA_SeqMatrix const*>(mat) ;
   if( smat == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_ILUp_PC error:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected" ) ;
   }
   
   LU->re_initialize( mat->nb_rows(), mat->nb_rows() ) ;
   LU->set( mat ) ; 

   LEVELS->re_initialize( mat->nb_rows(), mat->nb_rows() ) ;
   
   computeLevel( LU ) ;
   factorizeLU() ;

   // Free memory :
   LEVELS->nullify() ;

   BUILD_OK = true ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_ILUp_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   BUILD_OK=false ;
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: solve( LA_Vector const* rhs,
                    LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: solve" ) ;
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
LA_ILUp_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: factorizeLU( void )
//----------------------------------------------------------------------
{   
   size_t n = LU->nb_rows() ;
   
   LA_MatrixIterator* itI = LU->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itK = LU->create_stored_item_iterator( 0 ) ;

   if( PEL::abs( LU->item( 0, 0 ) )<PIV_MIN )
   {
      LU->nullify_row( 0 ) ;
      LU->set_item( 0, 0, 1. ) ;
      LA_ILUp_PC_ERROR::n0( 0 ) ;
   }
   
   for( size_t i=1 ; i<n ; i++ ) 
   {
      double drop = 0.0 ;
      for( itI->start_row_items(i) ; itI->is_valid() ; itI->go_next() )
      {
         size_t k = itI->col() ;
         if( k<i )
         {
            double akk = LU->item( k, k ) ;
            double aik = itI->item()/akk ;
            itI->set_item( aik ) ;
            for( itK->start_row_items(k) ; itK->is_valid() ; itK->go_next() )
            {
               size_t j = itK->col() ;
               if( j>k )
               {
                  int l = static_cast<int>( LEVELS->item( i, j ) ) ;
                  double akj =  itK->item() ;
                  if( j==i || l <= (int) LEVEL_ORDER )
                  {  
                     LU->add_to_item( i, j, - aik*akj ) ;
                  }
                  else if( MODIFIED )
                  {
                     drop -= aik*akj ;
                  }
               }
            }
         }
      }
      double aii = LU->item( i, i ) ;
      if( MODIFIED )
      {
         LU->add_to_item( i, i, -drop ) ;
         aii = LU->item( i, i ) ;
      }
      if( PEL::abs( aii )<PIV_MIN ) 
      {
         LU->nullify_row( i ) ;
         LU->set_item( i, i, 1. ) ;
         LA_ILUp_PC_ERROR::n0( i ) ;
      }
   }
   itI->destroy() ;
   itI = 0 ;
   itK->destroy() ;
   itK = 0 ;

   LU->synchronize() ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: computeLevel( LA_SeqMatrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: computeLevel" ) ;
   PEL_CHECK_PRE( mat != 0 && mat->is_synchronized() ) ;
   
   LEVELS->nullify() ;
   int n = mat->nb_rows() ;
   double const zero = 1e-5 ;
   
   for( int i=0 ; i<n ; i++ )
   {
      LEVELS->set_item( i, i, zero ) ;
   }
   LA_MatrixIterator* it = mat->create_stored_item_iterator( 0 ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      int i = it->row() ;
      int j = it->col() ;
      if( i!=j )
      {
         LEVELS->set_item( i, j, zero ) ;
      }
   }
   it->destroy() ;
   it = 0 ;

   LEVELS->synchronize() ;
   LA_MatrixIterator* itI = LEVELS->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itJ = LEVELS->create_stored_item_iterator( 0 ) ;
   for( int i=1 ; i<n ; i++ ) 
   {
      for( itI->start_row_items(i) ; itI->is_valid() ; itI->go_next() )
      {
         int k = itI->col() ;
         if( k<i )
         {
            
            int aik = (int) itI->item() ;
            for( itJ->start_row_items(k) ; itJ->is_valid() ; itJ->go_next() )
            {
               int j = itJ->col() ;
               if( j>k )
               {
                  int akj =  (int) itJ->item() ;
                  if( akj <= (int) LEVEL_ORDER )
                  {
                     double aij_d = LEVELS->item( i, j ) ;
                     int aij = (int) aij_d ;
                     if( aij_d == 0.0 )
                     {
                        aij = 1000+LEVEL_ORDER ;
                     }
                     if( aij != 0 )
                     {
                        int lev = (int) PEL::min( aij, aik+akj+1 ) ;
                        LEVELS->set_item( i, j, lev ) ;
                     }
                  }
               }
            }
         }
      }
   }
   itI->destroy() ; itI = 0 ;
   itJ->destroy() ; itJ = 0 ;
   
   LEVELS->synchronize() ;
}

//----------------------------------------------------------------------
void
LA_ILUp_PC:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILUp_PC:: print_more" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string s( indent_width, ' ') ;

   os << s << "smallest_nonzero_pivot : " << PIV_MIN << std::endl ;
   os << s << "diagonal_compensation  : "
      << ( MODIFIED ? "true" : "false" ) << std::endl ;
   os << s << "order : " << LEVEL_ORDER << std::endl ;
   
}

//internal--------------------------------------------------------------
void 
LA_ILUp_PC_ERROR:: n0( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** ILUp factorization :" << std::endl ;
   mesg << "***   zero pivot at step " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
