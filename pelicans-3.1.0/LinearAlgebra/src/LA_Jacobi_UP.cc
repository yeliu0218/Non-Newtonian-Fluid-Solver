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

#include <LA_Jacobi_UP.hh>

#include <LA_MatrixIterator.hh>
#include <LA_Matrix.hh>
#include <LA_Scatter.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <boolVector.hh>
#include <size_t_vector.hh>

#include <sstream>

struct LA_Jacobi_UP_ERROR
{
   static void n0( int i ) ;
} ;

//-----------------------------------------------------------------------------
LA_Jacobi_UP*
LA_Jacobi_UP:: create( PEL_Object* a_owner,
                       double smallest_inverted_item )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: create" ) ;
   PEL_CHECK_PRE( smallest_inverted_item > 0. ) ;
   
   LA_Jacobi_UP* result = new LA_Jacobi_UP( a_owner,
                                            smallest_inverted_item ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
LA_Jacobi_UP:: LA_Jacobi_UP( PEL_Object* a_owner,
                             double smallest_inverted_item )
//-----------------------------------------------------------------------------
   : LA_UzawaPreconditioner( a_owner )
   , MIN_DIAG( smallest_inverted_item )
   , DIAG_A( LA_SeqVector::create( this, 0 ) )
   , INV_PREC( 0 )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//-----------------------------------------------------------------------------
LA_Jacobi_UP:: ~LA_Jacobi_UP( void )
//-----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Jacobi_UP*
LA_Jacobi_UP:: create( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: create" ) ;
   PEL_CHECK( exp != 0 ) ;

   double a_smallest_item = exp->double_data( "smallest_inverted_item" ) ;
   if( a_smallest_item<=0. )
   {
      PEL_Error::object()->raise_bad_data_value( exp, 
                                                 "smallest_inverted_item", 
                                                 "greater than zero" ) ;
   }
   
   LA_Jacobi_UP* result = LA_Jacobi_UP::create( a_owner, a_smallest_item ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
LA_Jacobi_UP*
LA_Jacobi_UP:: create_clone( PEL_Object* a_owner ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: create_clone" ) ;

   LA_Jacobi_UP* result = new LA_Jacobi_UP( a_owner, MIN_DIAG ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Jacobi_UP:: is_valid( void ) const
//-----------------------------------------------------------------------------
{
   return( BUILD_OK ) ;
}

//-----------------------------------------------------------------------------
size_t
LA_Jacobi_UP:: dimension( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( INV_PREC->nb_rows() ) ;
}

//-----------------------------------------------------------------------------
LA_Implementation const*
LA_Jacobi_UP:: implementation( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: implementation" ) ;
   PEL_CHECK_PRE( implementation_PRE() ) ;

   return( INV_PREC->implementation() ) ;
}

//-----------------------------------------------------------------------------
void
LA_Jacobi_UP:: build( LA_Matrix const* A,
                      LA_Matrix const* B,
                      LA_Matrix const* C )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: build" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, C ) ) ;

   if( INV_PREC != 0 &&
       INV_PREC->implementation() != A->implementation() )
   {
      destroy_possession( INV_PREC ) ;
      INV_PREC = 0 ;
   }

   // Get A diagonal:
   get_diagA( A, B, DIAG_A ) ;
   
   if( INV_PREC == 0 ) INV_PREC = A->create_vector( this ) ;
   build_prec( A, B, C, DIAG_A, INV_PREC ) ;

   BUILD_OK = true ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( A, B, C ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//-----------------------------------------------------------------------------
void
LA_Jacobi_UP:: solve( LA_Vector const* rhs,
                      LA_Vector* sol )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   sol->set_as_v_product( rhs, INV_PREC ) ;
   SOLVE_OK = true ;
   
   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_Jacobi_UP:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------------
void
LA_Jacobi_UP::  get_diagA( LA_Matrix const* A,
                           LA_Matrix const* B,
                           LA_SeqVector* DiagA ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: get_diagA" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, 0 ) ) ;
   PEL_CHECK_PRE( DiagA != 0 ) ;

   DiagA->re_initialize( A->nb_rows() ) ;

   LA_SeqMatrix const* seq_mat = dynamic_cast<LA_SeqMatrix const*>( A ) ;
   if( seq_mat != 0 )
   {
      for( size_t i=0 ; i<A->nb_rows() ; ++i )
      {
         DiagA->set_item( i, seq_mat->item(i,i) ) ;
      }
   }
   else
   {
      // Global diagonal vector:
      LA_Vector* D = A->create_vector( 0 ) ;
      D->re_initialize( A->nb_rows() ) ;
      A->extract_diag( D ) ;

      boolVector needed( A->nb_rows() ) ;
      LA_MatrixIterator* B_it = B->create_stored_item_iterator( 0 ) ;
      for( B_it->start_all_items() ; B_it->is_valid() ; B_it->go_next() )
      {
         needed( B_it->col() ) = true ;
      }
      B_it->destroy() ; B_it = 0 ;
      size_t_vector elm_idx(0) ;
      for( size_t i=0 ; i<A->nb_rows() ; ++i )
      {
         if( needed(i) ) elm_idx.append(i) ;
      }

      LA_Scatter* scatter = D->create_scatter( 0, elm_idx, elm_idx ) ;
      scatter->get( D, DiagA ) ;
      scatter->destroy() ; scatter = 0 ;
      D->destroy() ; D = 0 ;
   }
}

//----------------------------------------------------------------------------
void
LA_Jacobi_UP::  build_prec( LA_Matrix const* A,
                            LA_Matrix const* B,
                            LA_Matrix const* C,
                            LA_SeqVector const* DiagA,
                            LA_Vector* mres ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: build_prec" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, C ) ) ;
   PEL_CHECK_PRE( DiagA != 0 ) ;
   PEL_CHECK_PRE( DiagA->nb_rows() == A->nb_rows() ) ;
   PEL_CHECK_PRE( mres != 0 ) ;
   PEL_CHECK_PRE( mres->implementation() == A->implementation() ) ;
   
   // Compute mres :
   mres->re_initialize( B->nb_rows() ) ;
   LA_MatrixIterator* B_it = B->create_stored_item_iterator( 0 ) ;
   for( B_it->start_all_items() ; B_it->is_valid() ; B_it->go_next() )
   {
      size_t const i = B_it->row() ;
      size_t const k = B_it->col() ;
      double akk = DiagA->item(k) ;
      if( PEL::abs( akk )<MIN_DIAG )
      {
         akk = 1. ;
         LA_Jacobi_UP_ERROR:: n0( k ) ;
      }
      double d = 1./akk ;
      double const v = B_it->item() ;
      mres->add_to_item( i, v*d*v ) ;
   }
   B_it->destroy() ; B_it = 0 ;
   
   if( C!=0 )
   {
      LA_MatrixIterator* C_it = C->create_stored_item_iterator( 0 ) ;
      for( size_t i=0 ; i<C->nb_rows() ; ++i )
      {
         double cii = 0. ;
         for( C_it->start_row_items(i) ; C_it->is_valid() ; C_it->go_next() )
         {
            if( C_it->col() == i )
            {
               cii = C_it->item() ;
               break ;
            }
         }
         mres->add_to_item( i, cii ) ;
      }
      C_it->destroy() ; C_it = 0 ;
   }
   mres->synchronize() ;
   
   mres->synchronize() ;
   mres->set_as_reciprocal( mres, MIN_DIAG, 1. ) ;
   mres->synchronize() ;
}

//----------------------------------------------------------------------
void
LA_Jacobi_UP:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_UP:: print" ) ;

   LA_UzawaPreconditioner:: print( os, indent_width ) ;
   std::string s( indent_width+3, ' ' ) ;
   os << s << "smallest_inverted_item:" << MIN_DIAG << std::endl ;
}

//internal--------------------------------------------------------------
void 
LA_Jacobi_UP_ERROR:: n0( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** Uzawa Jacobi preconditioner:" << std::endl ;
   mesg << "    vanishing diagonal term at line " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
