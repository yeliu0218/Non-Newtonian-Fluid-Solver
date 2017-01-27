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

#include <EXT_UMFsolver_DS.hh>

#include <LA_Matrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

extern "C" {
#include <umfpack.h>
}

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

EXT_UMFsolver_DS const* EXT_UMFsolver_DS::PROTOTYPE = new EXT_UMFsolver_DS() ;

//----------------------------------------------------------------------
EXT_UMFsolver_DS:: EXT_UMFsolver_DS( void )
//----------------------------------------------------------------------
  : LA_Solver( "EXT_UMFsolver_DS" )
  , SIZE( PEL::bad_index() )
  , NB_NO_ZERO( PEL::bad_index() )
  , UMF_Ap( 0 )
  , UMF_Ai( 0 )
  , UMF_Ax( 0 )
  , UMF_b( 0 )
  , UMF_x( 0 )
  , UMF_Numeric( 0 )
  , UMF_Control( 0 )
  , UMF_Info( 0 )
  , VERBOSE( false )
{
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_UMFPACK" ), val ) ;
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS*
EXT_UMFsolver_DS:: create( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   EXT_UMFsolver_DS* result = new EXT_UMFsolver_DS( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS*
EXT_UMFsolver_DS:: create_replica( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   EXT_UMFsolver_DS* result = new EXT_UMFsolver_DS( a_owner, exp ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS:: EXT_UMFsolver_DS( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
  : LA_Solver( a_owner )
  , EXP( exp->create_clone( this ) )
  , SIZE( PEL::bad_index() )
  , NB_NO_ZERO( PEL::bad_index() )
  , UMF_Ap( 0 )
  , UMF_Ai( 0 )
  , UMF_Ax( 0 )
  , UMF_b( 0 )
  , UMF_x( 0 )
  , UMF_Numeric( 0 )
  , UMF_Control( new double[UMFPACK_CONTROL] )
  , UMF_Info( new double[UMFPACK_INFO] )
  , VERBOSE( exp->has_entry( "verbose" ) ? exp->bool_data( "verbose" )
                                         : false )
{
   PEL_LABEL( "EXT_UMFsolver_DS:: EXT_UMFsolver_DS" ) ;
   initialize() ;
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS*
EXT_UMFsolver_DS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: create_clone" ) ;
   
   EXT_UMFsolver_DS* result = new EXT_UMFsolver_DS( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS:: EXT_UMFsolver_DS( PEL_Object* a_owner,
                                     EXT_UMFsolver_DS const* other )
//----------------------------------------------------------------------
  : LA_Solver( a_owner )
  , EXP( other->EXP->create_clone( this ) )
  , SIZE( PEL::bad_index() )
  , NB_NO_ZERO( PEL::bad_index() )
  , UMF_Ap( 0 )
  , UMF_Ai( 0 )
  , UMF_Ax( 0 )
  , UMF_b( 0 )
  , UMF_x( 0 )
  , UMF_Numeric( 0 )
  , UMF_Control( new double[UMFPACK_CONTROL] )
  , UMF_Info( new double[UMFPACK_INFO] )
  , VERBOSE( other->VERBOSE )
{
   PEL_LABEL( "EXT_UMFsolver_DS:: EXT_UMFsolver_DS" ) ;
   initialize() ;

}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: initialize( void )
//----------------------------------------------------------------------
{
   umfpack_di_defaults( UMF_Control ) ;
   
   if( VERBOSE )
   {
      UMF_Control[UMFPACK_PRL] = (double) 2 ;
   }

   if( EXP->has_entry( "UMFPACK_STRATEGY" ) )
   {
      std::string const& strategy = EXP->string_data( "UMFPACK_STRATEGY" ) ;
      if( strategy == "UMFPACK_STRATEGY_AUTO" )
      {
         UMF_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_AUTO ;
      }
      else if( strategy == "UMFPACK_STRATEGY_UNSYMMETRIC" )
      {
         UMF_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC ;
      }
      else if( strategy == "UMFPACK_STRATEGY_SYMMETRIC" )
      {
         UMF_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC ;
      }
      else if( strategy == "UMFPACK_STRATEGY_2BY2" )
      {
         UMF_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_2BY2 ;
      }
      else
      {
         std::string  allowed_values =
            "   - \"UMFPACK_STRATEGY_AUTO\"\n"
            "   - \"UMFPACK_STRATEGY_UNSYMMETRIC\"\n"
            "   - \"UMFPACK_STRATEGY_SYMMETRIC\"\n"
            "   - \"UMFPACK_STRATEGY_2BY2\"" ;
         PEL_Error::object()->raise_bad_data_value( EXP, "UMFPACK_STRATEGY",
                                                    allowed_values ) ;
      }
   }
   if( EXP->has_entry( "UMFPACK_DENSE_COL" ) )
   {
      UMF_Control[UMFPACK_DENSE_COL] = EXP->double_data( "UMFPACK_DENSE_COL" ) ;
   }
   if( EXP->has_entry( "UMFPACK_DENSE_ROW" ) )
   {
      UMF_Control[UMFPACK_DENSE_ROW] = EXP->double_data( "UMFPACK_DENSE_ROW" ) ;
   }
   if( EXP->has_entry( "UMFPACK_AMD_DENSE" ) )
   {
      UMF_Control[UMFPACK_AMD_DENSE] = EXP->double_data( "UMFPACK_AMD_DENSE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_BLOCK_SIZE" ) )
   {
      UMF_Control[UMFPACK_BLOCK_SIZE] = EXP->double_data( "UMFPACK_BLOCK_SIZE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_2BY2_TOLERANCE" ) )
   {
      UMF_Control[UMFPACK_2BY2_TOLERANCE] = EXP->double_data( "UMFPACK_2BY2_TOLERANCE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_SCALE" ) )
   {
      std::string const& scale = EXP->string_data( "UMFPACK_SCALE" ) ;
      if( scale == "UMFPACK_SCALE_NONE" )
      {
         UMF_Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE ;
      }
      else if( scale == "UMFPACK_SCALE_SUM" )
      {
         UMF_Control[UMFPACK_SCALE] = UMFPACK_SCALE_SUM ;
      }
      else if( scale == "UMFPACK_SCALE_MAX" )
      {
         UMF_Control[UMFPACK_SCALE] = UMFPACK_SCALE_MAX ;
      }
      else
      {
         std::string  allowed_values =
            "   - \"UMFPACK_SCALE_NONE\"\n"
            "   - \"UMFPACK_SCALE_SUM\"\n"
            "   - \"UMFPACK_SCALE_MAX\"" ;
         PEL_Error::object()->raise_bad_data_value( EXP, "UMFPACK_SCALE",
                                                    allowed_values ) ;
      }
   }
   if( EXP->has_entry( "UMFPACK_FIXQ" ) )
   {
      UMF_Control[UMFPACK_FIXQ] = EXP->double_data( "UMFPACK_FIXQ" ) ;
   }
   if( EXP->has_entry( "UMFPACK_AGGRESSIVE" ) )
   {
      UMF_Control[UMFPACK_AGGRESSIVE] = EXP->double_data( "UMFPACK_AGGRESSIVE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_PIVOT_TOLERANCE" ) )
   {
      UMF_Control[UMFPACK_PIVOT_TOLERANCE] = EXP->double_data( "UMFPACK_PIVOT_TOLERANCE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_SYM_PIVOT_TOLERANCE" ) )
   {
      UMF_Control[UMFPACK_SYM_PIVOT_TOLERANCE] = EXP->double_data( "UMFPACK_SYM_PIVOT_TOLERANCE" ) ;
   }
   if( EXP->has_entry( "UMFPACK_DROPTOL" ) )
   {
      UMF_Control[UMFPACK_DROPTOL] = EXP->double_data( "UMFPACK_DROPTOL" ) ;
   }
   if( EXP->has_entry( "UMFPACK_IRSTEP" ) )
   {
      UMF_Control[UMFPACK_IRSTEP] = EXP->double_data( "UMFPACK_IRSTEP" ) ;
   }
   if( EXP->has_entry( "UMFPACK_PRL" ) )
   {
      int const prl = EXP->int_data( "UMFPACK_PRL" ) ; ;
      UMF_Control[UMFPACK_PRL] = (double) prl ;
      VERBOSE = ( prl>=2 ) ;
   }
   set_iterative( false );
}

//----------------------------------------------------------------------
EXT_UMFsolver_DS:: ~EXT_UMFsolver_DS( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: ~EXT_UMFsolver_DS" ) ;

   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   else
   {
      delete [] UMF_Control ; UMF_Control = 0 ;
      delete [] UMF_Info ; UMF_Info = 0 ;
   }
   if( matrix_is_set() )
   {
      unset_matrix() ;
   }
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: set_matrix_self( LA_Matrix const* mat,
                                    bool &ok, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   ok = true ;
   
   LA_SeqMatrix const* smat = dynamic_cast<LA_SeqMatrix const*>( mat ) ;
   if( smat == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** EXT_UMFsolver_DS:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected" ) ;
   }
   
   // Printing :
   if( VERBOSE )
   {
      PEL::out() << "Matrix factorization : " << std::endl ;
   }

   SIZE = smat->nb_rows() ;
   NB_NO_ZERO = smat->nb_stored_items() ;

   if( SIZE==0 || NB_NO_ZERO==0 )
   {
      PEL_Error::object()->raise_plain(
         "*** EXT_UMFsolver_DS :\n"
         "    a no nil matrix is allowed" ) ;
   }

   UMF_Ap = new int[SIZE+1] ;
   for( size_t i=0 ; i<SIZE+1 ; i++ )
      UMF_Ap[i] = -1 ;
   UMF_Ai = new int[NB_NO_ZERO] ;
   UMF_Ax = new double[NB_NO_ZERO] ;
   UMF_x = new double[SIZE] ;
   UMF_b = new double[SIZE] ;

   // Build matrix :
   {
      int* Ti = new int[NB_NO_ZERO] ;
      int* Tj = new int[NB_NO_ZERO] ;
      double* Tx = new double[NB_NO_ZERO] ;
      int* Map = 0 ;
      
      LA_MatrixIterator* mat_it = smat->create_stored_item_iterator( 0 ) ;
      mat_it->start_all_items() ;
      size_t i = 0 ;
      for( ; mat_it->is_valid() ; ++i, mat_it->go_next() )
      {
         Ti[i] = mat_it->row() ;
         Tj[i] = mat_it->col() ;
         Tx[i] = mat_it->item() ;
      }
      mat_it->destroy() ; mat_it = 0 ;

      int s = umfpack_di_triplet_to_col( SIZE, SIZE, NB_NO_ZERO, Ti, Tj, Tx,
                                         UMF_Ap, UMF_Ai, UMF_Ax, Map ) ;

      delete [] Ti ; Ti = 0 ;
      delete [] Tj ; Tj = 0 ;
      delete [] Tx ; Tx = 0 ;

      if( s != UMFPACK_OK )
      {
         raise_fatal_error( "umfpack_di_triplet_to_col", s ) ;
         ok = false ;
      }
   }
   if( VERBOSE )
   {
      umfpack_di_report_matrix( SIZE, SIZE, UMF_Ap, UMF_Ai, UMF_Ax, 0,
                                UMF_Control ) ;
   }
   
   if(ok) inverse(ok) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: inverse( bool& ok )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: inverse" ) ;
   
   // Decompose matrix :
   {
      void *UMF_Symbolic = 0 ;
      int const s1 = umfpack_di_symbolic( SIZE, SIZE,
                                          UMF_Ap, UMF_Ai, UMF_Ax,
                                          &UMF_Symbolic,
                                          UMF_Control, UMF_Info ) ;
      if( s1!=UMFPACK_OK )
      {
         raise_fatal_error( "umfpack_di_symbolic", s1 ) ;
      }
      int const s2 = umfpack_di_numeric( UMF_Ap, UMF_Ai, UMF_Ax,
                                         UMF_Symbolic,
                                         &UMF_Numeric,
                                         UMF_Control, UMF_Info ) ;
      if( s2==UMFPACK_OK )
      {
         // Success
      }
      else
      {
         raise_fatal_error( "umfpack_di_numeric", s2 ) ;
         ok = false ;
      }
      umfpack_di_free_symbolic( &UMF_Symbolic ) ;
      UMF_Symbolic = 0 ;
   }

   // Printing :
   if( VERBOSE )
   {
      umfpack_di_report_info( UMF_Control, UMF_Info ) ;
   }
      
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: unset_matrix_self( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: unset_matrix_self" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   SIZE = PEL::bad_index() ;
   delete [] UMF_Ap ; UMF_Ap = 0 ;
   delete [] UMF_Ai ; UMF_Ai = 0 ;
   delete [] UMF_Ax ; UMF_Ax = 0 ;
   delete [] UMF_b  ; UMF_b  = 0 ;
   delete [] UMF_x  ; UMF_x  = 0 ;
   if( UMF_Numeric!=0 )
   {
      umfpack_di_free_numeric( &UMF_Numeric ) ;
      UMF_Numeric = 0 ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   nb_iter=1 ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( b ) != 0 ) ;
   LA_SeqVector const* bb = static_cast<LA_SeqVector const* >( b ) ;
   
   // Printing :
   if( VERBOSE )
   {
      PEL::out() << "System resolution : " << std::endl ;
   }
   
   for( size_t i=0 ; i<SIZE ; ++i )
   {
      UMF_b[i] = bb->item(i) ;
   }

   int s = umfpack_di_solve( UMFPACK_A,
                             UMF_Ap, UMF_Ai, UMF_Ax, UMF_x, UMF_b,
                             UMF_Numeric,
                             UMF_Control, UMF_Info ) ;
   if( s!=UMFPACK_OK )
   {
      raise_fatal_error( "umfpack_di_solve", s ) ;
      ok = false ;
   }
   else
   {
      ok = true ;
      for( size_t i=0 ; i<SIZE ; ++i )
      {
         x->set_item( i, UMF_x[i] ) ;
      }
      
      // Printing :
      if( VERBOSE )
      {
         umfpack_di_report_info( UMF_Control, UMF_Info ) ;
         PEL::out() << "System resolution : successful" << std::endl ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_UMFsolver_DS:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "Direct solver : \"EXT_UMFsolver_DS\"" << std::endl ;
   os << s << "   UMFPACK options : " << std::endl ;

   double aux = UMF_Control[UMFPACK_PRL] ;
   UMF_Control[UMFPACK_PRL] = (double) 2 ;
   umfpack_di_report_control( UMF_Control ) ;
   UMF_Control[UMFPACK_PRL] = aux ;
}

//----------------------------------------------------------------------
void
EXT_UMFsolver_DS:: raise_fatal_error( std::string const& routine,
                                      int status ) const
//----------------------------------------------------------------------
{
   UMF_Control[UMFPACK_PRL] = 2 ;      
   umfpack_di_report_status( UMF_Control, status ) ;

   std::ostringstream mesg ;
   mesg << "*** EXT_UMFsolver_DS:" << endl << endl ;
   mesg << "    the status value: " << status << endl ;
   mesg << "    has been returned by: \"" << routine << "\"" << endl ;
   mesg << "    which means that something was not successful" << endl << endl ;
   mesg << "    see above for the associated output of" << endl ;
   mesg << "      \"umfpack_di_report_status\"" ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
