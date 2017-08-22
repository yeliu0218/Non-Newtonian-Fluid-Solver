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

#include <LA_CahouetChabard_UP.hh>

#include <LA_Solver.hh>
#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;

struct LA_CahouetChabard_UP_ERROR
{
   static void n0( void ) ;
   static void n1( int i ) ;
   static void n2( void ) ;
} ;

//-----------------------------------------------------------------------------
LA_CahouetChabard_UP:: LA_CahouetChabard_UP(
                            PEL_Object* a_owner,
                            LA_Matrix const* a_matrix_prototype,
                            LA_Solver* a_solver_of_M )
//-----------------------------------------------------------------------------
   : LA_UzawaPreconditioner( a_owner )
   , MIN_DIAG( 1.E-30 )
   , MAT_M( a_matrix_prototype->create_matrix( this ) )
   , BT( a_matrix_prototype->create_matrix( this ) )
   , SOLVER( a_solver_of_M )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
   a_solver_of_M->set_owner( this ) ;
}

//-----------------------------------------------------------------------------
LA_CahouetChabard_UP*
LA_CahouetChabard_UP:: create( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PEL_ModuleExplorer const* se1 =
                       exp->create_subexplorer( 0, "LA_Matrix" ) ;
   LA_Matrix const* matrix_prototype = LA_Matrix::make( 0, se1 ) ;
   se1->destroy() ; se1 = 0 ;
   
   PEL_ModuleExplorer* se2 =
                     exp->create_subexplorer( 0, "LA_Solver" );
   LA_Solver* prec_solver =  LA_Solver::make( 0, se2 ) ;
   se2->destroy() ; se2 = 0 ;

   LA_CahouetChabard_UP* result =
          new LA_CahouetChabard_UP( a_owner,
                                    matrix_prototype,
                                    prec_solver ) ;

   matrix_prototype->destroy() ; matrix_prototype = 0 ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
LA_CahouetChabard_UP:: ~LA_CahouetChabard_UP( void )
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
LA_CahouetChabard_UP*
LA_CahouetChabard_UP:: create_clone( PEL_Object* a_owner ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: create_clone" ) ;

   LA_CahouetChabard_UP* result =
      new LA_CahouetChabard_UP( a_owner,
                                MAT_M,
                                SOLVER->create_clone( 0 ) ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
LA_CahouetChabard_UP:: is_valid( void ) const
//-----------------------------------------------------------------------------
{
   return( BUILD_OK ) ;
}

//-----------------------------------------------------------------------------
size_t
LA_CahouetChabard_UP:: dimension( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( MAT_M->nb_rows() ) ;
}

//-----------------------------------------------------------------------------
LA_Implementation const*
LA_CahouetChabard_UP:: implementation( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: implementation" ) ;
   PEL_CHECK_PRE( implementation_PRE() ) ;

   return( MAT_M->implementation() ) ;
}

//-----------------------------------------------------------------------------
void
LA_CahouetChabard_UP:: build( LA_Matrix const* A,
                              LA_Matrix const* B,
                              LA_Matrix const* C )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: build" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, C ) ) ;

   if( A->implementation() != MAT_M->implementation() )
   {
      LA_CahouetChabard_UP_ERROR::n2() ;
   }
   build_prec( A, B, C, MAT_M ) ;
   if( SOLVER->matrix_is_set() ) SOLVER->unset_matrix() ;
   
   SOLVER->set_matrix( MAT_M ) ;
   BUILD_OK = SOLVER->matrix_is_set() ;
   SOLVE_OK = false ;
   
   PEL_CHECK_POST( build_POST( A, B, C ) ) ;
}

//-----------------------------------------------------------------------------
void
LA_CahouetChabard_UP:: solve( LA_Vector const* rhs,
                              LA_Vector* sol )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   SOLVER->set_initial_guess_nonzero( false ) ;
   SOLVER->solve( rhs, sol ) ;
   SOLVE_OK = SOLVER->solution_is_achieved() ;
   if( !SOLVE_OK )
   {
      LA_CahouetChabard_UP_ERROR::n0() ;
   }

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
}

//----------------------------------------------------------------------
bool
LA_CahouetChabard_UP:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------------
void
LA_CahouetChabard_UP:: build_prec( LA_Matrix const* A,
                                   LA_Matrix const* B,
                                   LA_Matrix const* C,
                                   LA_Matrix* mres ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: build_prec" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, C ) ) ;
   PEL_CHECK_PRE( mres != 0 ) ;
   PEL_CHECK_PRE( mres->implementation() == A->implementation() ) ;
   
   // Compute BT :
   BT->re_initialize( B->nb_cols(), B->nb_rows() ) ;
   BT->add_tMat( B ) ;
   BT->synchronize() ;
   
   // Compute mres :
   LA_MatrixIterator* BT_it1 = BT->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* BT_it2 = BT->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* A_it = A->create_stored_item_iterator( 0 ) ;
   mres->re_initialize( B->nb_rows(), B->nb_rows() ) ;
   for( BT_it1->start_all_items() ; BT_it1->is_valid() ; BT_it1->go_next() )
   {
      size_t const k = BT_it1->row() ;
      size_t const i = BT_it1->col() ;
      double const v1 = BT_it1->item() ;
      double akk = 0. ;
      // Akk is owned (stored) by current process:
      for( A_it->start_row_items(k) ; A_it->is_valid() ; A_it->go_next() )
      {
         if( A_it->col() == k )
         {
            akk = A_it->item() ;
            break ;
         }
      }
      if( PEL::abs( akk )<MIN_DIAG )
      {
         akk = 1. ;
         LA_CahouetChabard_UP_ERROR:: n1( k ) ;
      }
      double const d = 1./akk ;
      for( BT_it2->start_row_items(k) ; BT_it2->is_valid() ; BT_it2->go_next() )
      {
         size_t const j = BT_it2->col() ;
         double const v2 = BT_it2->item() ;
         mres->add_to_item( i, j, v1*d*v2 ) ;
      }
   }
   BT_it1->destroy() ; BT_it1 = 0 ;
   BT_it2->destroy() ; BT_it2 = 0 ;
   A_it->destroy() ; A_it = 0 ;

   // Free memory :
   BT->nullify() ;
   
   // Add C :
   if( C!=0 )
   {
      mres->add_Mat( C, -1.0 ) ;
   }
   mres->synchronize() ;
}

//----------------------------------------------------------------------
void
LA_CahouetChabard_UP:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_CahouetChabard_UP:: print" ) ;

   LA_UzawaPreconditioner:: print( os, indent_width ) ;
   std::string s( indent_width+3, ' ' ) ;
   os << s << "smallest_inverted_item:" << MIN_DIAG << std::endl ;
   SOLVER->print( os, indent_width+3 ) ;
   MAT_M->print( os, indent_width+3 ) ;
}

//internal--------------------------------------------------------------
void 
LA_CahouetChabard_UP_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** Cahouet and Chabard preconditioner:" << std::endl ;
   mesg << "    convergence failure in the solution of the linear" << endl ;
   mesg << "    system whose matrix is the preconditioner." << endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
LA_CahouetChabard_UP_ERROR:: n1( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** Cahouet and Chabard preconditioner:" << std::endl ;
   mesg << "    vanishing diagonal term at line " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
LA_CahouetChabard_UP_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** Cahouet and Chabard preconditioner:" << std::endl ;
   mesg << "    incompatible matrix types" << std::endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
