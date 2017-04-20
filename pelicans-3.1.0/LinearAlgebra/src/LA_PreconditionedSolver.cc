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

#include <LA_PreconditionedSolver.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>

#include <LA_Preconditioner.hh>
#include <LA_SeqMatrix.hh>
#include <LA_IterativeSolver.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;

LA_PreconditionedSolver const*
LA_PreconditionedSolver::PROTOTYPE = new LA_PreconditionedSolver() ;

//----------------------------------------------------------------------------
LA_PreconditionedSolver:: LA_PreconditionedSolver( void ) 
//----------------------------------------------------------------------------
   : LA_Solver( "LA_PreconditionedSolver" )
   , SOLVER( 0 )
   , PRECOND( 0 )
   , MY_A( 0 )
{   
}

//----------------------------------------------------------------------------
LA_PreconditionedSolver:: ~LA_PreconditionedSolver( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: ~LA_PreconditionedSolver" ) ;
   
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   if( matrix_is_set() )
   {
      unset_matrix() ;
   }   
}

//----------------------------------------------------------------------------
LA_PreconditionedSolver*
LA_PreconditionedSolver:: create( PEL_Object* a_owner, 
                                  PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: create" ) ;

   LA_PreconditionedSolver* result = new LA_PreconditionedSolver( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->is_iterative() ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
LA_PreconditionedSolver*
LA_PreconditionedSolver:: create_replica( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   LA_PreconditionedSolver* result =
                                 new LA_PreconditionedSolver( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
LA_PreconditionedSolver:: LA_PreconditionedSolver(
                         PEL_Object* a_owner,  PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : LA_Solver( a_owner )
   , SOLVER( 0 )
   , PRECOND( 0 )
   , MY_A( 0 )
{
   PEL_LABEL( "LA_PreconditionedSolver:: LA_PreconditionedSolver" ) ;

   PEL_ModuleExplorer const* ee =
                exp->create_subexplorer( 0, "LA_IterativeSolver" ) ;
   SOLVER = LA_IterativeSolver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   ee = exp->create_subexplorer( 0, "LA_Preconditioner" ) ;
   PRECOND = LA_Preconditioner::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   set_iterative( true ) ;
}

//----------------------------------------------------------------------
LA_PreconditionedSolver*
LA_PreconditionedSolver:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: create_clone" ) ;
   
   LA_PreconditionedSolver* result =
                          new LA_PreconditionedSolver( a_owner, this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PreconditionedSolver:: LA_PreconditionedSolver(
             PEL_Object* a_owner, LA_PreconditionedSolver const* other )
//----------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , SOLVER( other->SOLVER->create_clone( this ) )
   , PRECOND( other->PRECOND->create_clone( this ) )
   , MY_A( 0 )
{
}

//----------------------------------------------------------------------------
void
LA_PreconditionedSolver:: set_matrix_self( LA_Matrix const* mat,
                                           bool &ok, bool same_pattern )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;
   
   MY_A = mat ;   
   PRECOND->build( MY_A ) ;
   ok = PRECOND->is_valid() ;

   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}


//----------------------------------------------------------------------------
void
LA_PreconditionedSolver:: unset_matrix_self( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: unset_matrix_self" ) ;
   PEL_CHECK( unset_matrix_self_PRE() ) ;
   
   PRECOND->unbuild() ;
   MY_A = 0 ;   
}

//----------------------------------------------------------------------------
void
LA_PreconditionedSolver:: solve_self( LA_Vector const* b, LA_Vector* x,
                                      size_t& nb_iter, bool& ok  ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   
   PEL_CHECK( MY_A!=0 ) ;

   SOLVER->solve( MY_A, b, PRECOND, zero_initial_guess(), x ) ;
   
   ok = SOLVER->convergence_achieved() ;
   nb_iter = SOLVER->nb_iterations_achieved() ;

   if( !ok )
   {
      LA_ConvergenceTest::ConvergedReason reason = SOLVER->converged_reason() ;
      PEL_ASSERT( reason<=0 ) ;
      std::ostringstream mesg ;
      mesg << "*** LA_PreconditionedSolver:" << endl ;
      mesg << "    convergence failure " << reason << endl ;
      PEL_Error::object()->display_info( mesg.str() ) ;
   }
      
   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------------
void
LA_PreconditionedSolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_PreconditionedSolver:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   SOLVER->print( os, indent_width ) ;
   PRECOND->print( os, indent_width ) ;
}
