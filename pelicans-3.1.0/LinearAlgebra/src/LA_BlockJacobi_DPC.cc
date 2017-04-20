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

#include <LA_BlockJacobi_DPC.hh>

#include <LA_SeqVector.hh>
#include <LA_DistMatrix.hh>
#include <LA_DistVector.hh>
#include <LA_Solver.hh>
#include <LA_SeqMatrix.hh>

#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

#include <iostream>

LA_BlockJacobi_DPC const* 
LA_BlockJacobi_DPC:: PROTOTYPE = new LA_BlockJacobi_DPC() ;

//----------------------------------------------------------------------
LA_BlockJacobi_DPC:: LA_BlockJacobi_DPC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_BlockJacobi_DPC" )
   , SIZE( 0 )
   , SOLVE_OK( false )
   , SOLVER( 0 )
{
}

//----------------------------------------------------------------------
LA_BlockJacobi_DPC*
LA_BlockJacobi_DPC:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: create_replica" ) ;

   LA_BlockJacobi_DPC* result = new LA_BlockJacobi_DPC( a_owner, exp ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockJacobi_DPC:: LA_BlockJacobi_DPC( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , SIZE( 0 )
   , SOLVE_OK( false )
   , SOLVER( 0 )
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: LA_BlockJacobi_DPC" ) ;
   PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, sexp ) ;
   sexp->destroy() ;
}

//----------------------------------------------------------------------
LA_BlockJacobi_DPC:: ~LA_BlockJacobi_DPC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_BlockJacobi_DPC* 
LA_BlockJacobi_DPC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: create_clone" ) ;

   LA_BlockJacobi_DPC* result = new LA_BlockJacobi_DPC( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockJacobi_DPC:: LA_BlockJacobi_DPC( PEL_Object* a_owner, 
                                         LA_BlockJacobi_DPC const* other ) 
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , SIZE( 0 )
   , SOLVE_OK( false )
   , SOLVER( other->SOLVER->create_clone( this ) )
{
}

//----------------------------------------------------------------------
bool
LA_BlockJacobi_DPC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( SIZE != 0 ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockJacobi_DPC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( SIZE ) ;
}

//----------------------------------------------------------------------
void
LA_BlockJacobi_DPC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_BlockJacobi_DPC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   LA_DistMatrix const* dmat = dynamic_cast<LA_DistMatrix const*>( mat ) ;
   if( dmat == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_BlockJacobi_DPC:\n"
         "    a matrix of type \"LA_DistMatrix\" is expected" ) ;
   }
   
   if( dmat->row_distribution()->local_number() > 0 )
   {
      LA_SeqMatrix const* block = dmat->diagonal_block_matrix() ;
      if( SOLVER->matrix_is_set() ) SOLVER->unset_matrix() ;
      SOLVER->set_matrix( block ) ;
   }

   // since a preconditioner solves a system whose unknown is in fact the
   // residual of the original system, a zero initial guess may be pertinent
   SOLVER->set_initial_guess_nonzero( false ) ;
   
   SIZE = mat->nb_rows() ;
   SOLVE_OK = false ;
   
   PEL_CHECK_POST( build_POST( mat ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockJacobi_DPC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_BlockJacobi_DPC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   
   SOLVER->unset_matrix() ;
   SIZE=0 ;
   SOLVE_OK = false ;
   
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_BlockJacobi_DPC:: solve( LA_Vector const* rhs, LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;
   
   PEL_CHECK( dynamic_cast<LA_DistVector const*>( rhs ) != 0 ) ;
   LA_DistVector const* drhs = static_cast<LA_DistVector const*>( rhs ) ;
   PEL_CHECK( dynamic_cast<LA_DistVector*>( sol ) != 0 ) ;
   LA_DistVector* dsol = static_cast<LA_DistVector*>( sol ) ;
   
   if( drhs->row_distribution()->local_number() > 0 )
   {
      SOLVER->solve( drhs->local_vector(), dsol->local_vector() ) ;
      SOLVE_OK = SOLVER->solution_is_achieved() ;
   }
   else
   {
      SOLVE_OK = true ;
   }
   sol->synchronize() ;
   
   SOLVE_OK = PEL_Exec::communicator()->boolean_and( SOLVE_OK ) ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockJacobi_DPC:: print_more(
                           std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_BlockJacobi_DPC:: print_more" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   SOLVER->print( os, indent_width ) ;   
}

//----------------------------------------------------------------------
bool
LA_BlockJacobi_DPC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   return( SOLVE_OK ) ;
}
