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

#include <LA_Cholesky_DS.hh>

#include <LA_SymmetricMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

LA_Cholesky_DS const* LA_Cholesky_DS::PROTOTYPE = new LA_Cholesky_DS() ;

struct LA_Cholesky_DS_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
LA_Cholesky_DS:: LA_Cholesky_DS( void )
//----------------------------------------------------------------------
   : LA_Solver( "LA_Cholesky_DS" )
   , PIV_MIN_VAL( PEL::bad_double() )
   , MAT( 0 )
{
}


//----------------------------------------------------------------------
LA_Cholesky_DS*
LA_Cholesky_DS:: create( PEL_Object* a_owner,
                        double const a_pivot_minimal_value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: create" ) ;
   PEL_CHECK_PRE( a_pivot_minimal_value>0. ) ;
   
   LA_Cholesky_DS* result = new LA_Cholesky_DS( a_owner,
                                              a_pivot_minimal_value ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Cholesky_DS*
LA_Cholesky_DS:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double const a_piv_min_val = exp->double_data( "pivot_minimal_value" ) ;
   if( a_piv_min_val<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "pivot_minimal_value", "a value greater than 0 is expected" ) ;
   }
   
   LA_Cholesky_DS* result = new LA_Cholesky_DS( a_owner, a_piv_min_val ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Cholesky_DS:: LA_Cholesky_DS( PEL_Object* a_owner,
                                 double const a_pivot_minimal_value )
//----------------------------------------------------------------------
   : LA_Solver( a_owner )
   , PIV_MIN_VAL( a_pivot_minimal_value )
   , MAT( LA_SymmetricMatrix::create( this, 0 ) )
{
}

//----------------------------------------------------------------------
LA_Cholesky_DS*
LA_Cholesky_DS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: create_clone" ) ;
   
   LA_Cholesky_DS* result = new LA_Cholesky_DS( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
LA_Cholesky_DS:: LA_Cholesky_DS( PEL_Object* a_owner,
                               LA_Cholesky_DS const* other )
//----------------------------------------------------------------------
   : LA_Solver( a_owner )
   , PIV_MIN_VAL( other->PIV_MIN_VAL )
   , MAT( other->MAT->create_matrix( this ) )
{
}

//----------------------------------------------------------------------
LA_Cholesky_DS:: ~LA_Cholesky_DS( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
LA_Cholesky_DS:: set_matrix_self( LA_Matrix const* mat,
                                  bool &ok, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   LA_SeqMatrix const* smat =
      dynamic_cast<LA_SeqMatrix const*>( mat ) ;
   
   ok = true ;
   if( smat!=0 && smat->is_symmetric() )
   {
      factorize_LDLt( smat, ok ) ;
   }
   else
   {
      PEL_Error::object()->raise_plain(
         "*** LA_Cholesky_DS :\n"
         "    only a sequential symmetric matrix is allowed" ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------
void
LA_Cholesky_DS:: unset_matrix_self( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: unset_matrix_self" ) ;
   MAT->re_initialize( 0, 0 ) ;
}

//----------------------------------------------------------------------
void
LA_Cholesky_DS:: solve_self(  LA_Vector const* b, 
                              LA_Vector* x,
                              size_t &nb_iter,
                              bool &ok )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   LA_SeqVector const* bb = static_cast<LA_SeqVector const*>( b ) ;
   LA_SeqVector* bx = static_cast<LA_SeqVector*>( x ) ;
   
   const size_t n = b->nb_rows() ;
   
   for( size_t i=0 ; i<n ; i++ )
   {
      double c = 0.0 ;
      for( size_t j=0 ; j<i ; j++ )
      {
         c += MAT->item( i, j )*bx->item(j) ;
      }
      bx->set_item( i, ( bb->item( i ) - c )/MAT->item( i, i ) ) ;
   }
   for( int i=n-1 ; i>=0 ; i-- )
   {
      double c = 0.0 ;
      for( size_t j=i+1 ; j<n ; j++ )
      {
         c += MAT->item( i, j )*bx->item(j) ;
      }
      bx->set_item( i, ( bx->item(i) - c )/MAT->item( i, i ) ) ;
   }   
   ok = true ;
   nb_iter = 1 ;

   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------
LA_SymmetricMatrix const*
LA_Cholesky_DS:: factorized_LDLt_matrix( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: factorized_LDLt_matrix" ) ;
   PEL_CHECK_PRE( matrix_is_set() ) ;
   
   return( MAT ) ;
}

//----------------------------------------------------------------------
void
LA_Cholesky_DS:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "Direct solver : \"LA_Cholesky_DS\"" << std::endl ;
   os << s << "   pivot_minimal_value = " << PIV_MIN_VAL << std::endl ;
}

//----------------------------------------------------------------------
void
LA_Cholesky_DS:: factorize_LDLt( LA_SeqMatrix const* smat, bool &ok ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Cholesky_DS:: factorize_LDLt" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t  const n = smat->nb_rows() ;
   MAT->re_initialize( n, n ) ;
   double aux ,p=0.0 ;

   for( size_t i=0; i<n ; ++i )
   {   
      for( size_t j=i; j<n ; ++j )
      {
         aux = smat->item( i, j );
         for( size_t k=0; k<i ; ++k )
         {
            aux -= MAT->item( k, j )*MAT->item( k, i ) ;
         }
         if( i==j )
         {
            if( aux<PIV_MIN_VAL )
            {
               LA_Cholesky_DS_ERROR::n0() ;
               ok = false ;
               return ;
            }
            p = 1.0/PEL::sqrt( aux );
         }   
         MAT->set_item( i, j ,aux*p );
      }
   }

}

//----------------------------------------------------------------------
bool
LA_Cholesky_DS:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( EQUIVALENT( is_a_prototype(), MAT==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_a_prototype(), !matrix_is_set() ) ) ;
   PEL_ASSERT( IMPLIES( MAT!=0, MAT->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( MAT!=0, MAT->nb_rows() == MAT->nb_cols() ) ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
LA_Cholesky_DS_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** LA_Cholesky_DS : " << std::endl ;
   mesg << "    null pivot found for matrix : " << std::endl ;
   mesg << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
