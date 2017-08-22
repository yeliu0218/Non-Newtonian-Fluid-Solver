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

#include <LA_IterativeSolver.hh>

#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <stringVector.hh>

#include <LA_ConvergenceTest.hh>
#include <LA_DefaultConvergenceMonitor.hh>
#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_SkipConvergenceTest.hh>

#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::string ;

//----------------------------------------------------------------------
LA_IterativeSolver*
LA_IterativeSolver:: make( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   LA_IterativeSolver const* proto =
      static_cast<LA_IterativeSolver const*>( plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   LA_IterativeSolver* result = proto->create_replica( a_owner, exp ) ;
   result->SOLVER_NAME = name ;
   if( exp->has_entry( "name" ) )
   {
      result->SOLVER_NAME += "(" + exp->string_data( "name" ) +")" ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->converged_reason() ==
                   LA_ConvergenceTest::Undetermined ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_IterativeSolver:: LA_IterativeSolver( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , MAXITS( 100000 )
   , CVG_TEST( 0 )
   , MONITOR( 0 )
   , RNORM( PEL::bad_double() )
   , LAST_IT( PEL::bad_index() )
   , VERBOSE( false )
{
   PEL_LABEL( "LA_IterativeSolver:: LA_IterativeSolver" ) ;

   if( exp->has_entry( "verbose" ) )
   {
      VERBOSE = exp->bool_data( "verbose" ) ;
   }

   if( exp->has_entry( "nb_iterations_max" ) )
   {
      MAXITS = exp->int_data( "nb_iterations_max" ) ;
      exp->test_data( "nb_iterations_max", "nb_iterations_max>0" ) ;
   }

   if( exp->has_module( "LA_ConvergenceTest" ) )
   {
      PEL_ModuleExplorer const* ee =
                         exp->create_subexplorer( 0, "LA_ConvergenceTest" ) ;
      CVG_TEST = LA_ConvergenceTest::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   else
   {
      size_t a_nb_iterations = ( exp->has_entry( "nb_iterations_max" ) ?
                                    MAXITS : 1 ) ;
      CVG_TEST = LA_SkipConvergenceTest::create( this, a_nb_iterations ) ;
   }

   if( exp->has_module( "LA_ConvergenceMonitor" ) )
   {
      PEL_ModuleExplorer const* ee =
                exp->create_subexplorer( 0, "LA_ConvergenceMonitor" ) ;
      MONITOR = LA_ConvergenceMonitor:: make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   else
   {
      MONITOR = LA_DefaultConvergenceMonitor:: create( this ) ;
   }

   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
LA_IterativeSolver:: LA_IterativeSolver( PEL_Object* a_owner,
                                         LA_IterativeSolver const* other )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , SOLVER_NAME( other->SOLVER_NAME )  //??? is it a decent behavior ???
   , MAXITS( other->MAXITS )
   , CVG_TEST( other->CVG_TEST->create_clone( this ) )
   , MONITOR( other->MONITOR->create_clone( this ) )
   , RNORM( other->RNORM )
   , LAST_IT( other->LAST_IT )
   , VERBOSE( other->VERBOSE )
{
   PEL_LABEL( "LA_IterativeSolver:: LA_IterativeSolver" ) ;
   PEL_CHECK_PRE( !other->is_a_prototype() ) ;
}

//----------------------------------------------------------------------------
LA_IterativeSolver:: LA_IterativeSolver( std::string const& a_name )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "LA_IterativeSolver:: LA_IterativeSolver" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
LA_IterativeSolver:: ~LA_IterativeSolver( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: ~LA_IterativeSolver" ) ;
}

//----------------------------------------------------------------------
std::string const&
LA_IterativeSolver:: name( void ) const
//----------------------------------------------------------------------
{
   return( SOLVER_NAME ) ;
}

//----------------------------------------------------------------------
LA_ConvergenceTest::NormType
LA_IterativeSolver:: norm_type( void ) const
//----------------------------------------------------------------------
{
   return( CVG_TEST->norm_type() ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: solve( LA_Matrix const* A,
                            LA_Vector const* b,
                            LA_Preconditioner* prec,
                            bool zero_initial_guess,
                            LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: do_solve" ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( b != 0 ) ;
   PEL_CHECK_PRE( prec != 0 ) ;
   PEL_CHECK_PRE( x != 0 ) ;
   PEL_CHECK_PRE( A->nb_rows() == A->nb_cols() ) ;
   PEL_CHECK_PRE( A->nb_rows() == b->nb_rows() ) ;
   PEL_CHECK_PRE(
       A->row_distribution()->is_compatible( b->row_distribution() ) ) ;
   PEL_CHECK_PRE( A->nb_cols() == x->nb_rows() ) ;
   PEL_CHECK_PRE(
       A->col_distribution()->is_compatible( x->row_distribution() ) ) ;
   PEL_CHECK_PRE( prec->dimension() == A->nb_rows() ) ;
   PEL_CHECK_INV( invariant() ) ;

   LAST_IT = PEL::bad_index() ;
   CVG_TEST->set_converged_reason( LA_ConvergenceTest::ConvergedIterating ) ;

   if( VERBOSE )
   {
      PEL::out() << name() << std::endl ;
      MONITOR->display_at_entry( A, b, prec, zero_initial_guess, x,
                                 CVG_TEST ) ;
   }

   do_solve( A, b, prec, zero_initial_guess, x ) ;

   if( ( CVG_TEST->converged_reason() ==
         LA_ConvergenceTest::ConvergedIterating ) )
   {
      PEL_ASSERT( LAST_IT >= MAXITS ) ;
      CVG_TEST->set_converged_reason( LA_ConvergenceTest::DivergedIts ) ;
   }

   if( VERBOSE )
   {
      MONITOR->display_at_exit( A, b, prec, x, CVG_TEST ) ;
   }
}

//----------------------------------------------------------------------
bool
LA_IterativeSolver:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
LA_ConvergenceTest::ConvergedReason
LA_IterativeSolver:: converged_reason( void ) const
//----------------------------------------------------------------------
{
   return( CVG_TEST->converged_reason() ) ;
}

//----------------------------------------------------------------------
bool
LA_IterativeSolver:: convergence_achieved( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: convergence_achieved" ) ;

   bool result = ( CVG_TEST->converged_reason() > 0 ) ;

   PEL_CHECK_POST( EQUIVALENT( result, ( converged_reason() > 0 ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_IterativeSolver:: nb_iterations_achieved( void ) const
//----------------------------------------------------------------------
{
   return( LAST_IT ) ;
}

//----------------------------------------------------------------------
double
LA_IterativeSolver:: residual_norm_achieved( void ) const
//----------------------------------------------------------------------
{
   return( RNORM ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: print" ) ;
   std::string s( indent_width, ' ') ;

   os << s << "LA_IterativeSolver: \"" << name() << "\"" << endl ;
   os << s << "   maximum iteration: " << MAXITS << endl ;
   print_more( os, indent_width+3 ) ;
   CVG_TEST->print( os, indent_width+3 ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: set_norm_type( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: set_norm_type" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   if( exp->has_entry( "norm_type" ) )
   {
      std::string const& nn = exp->string_data( "norm_type" ) ;
      exp->test_data_in( "norm_type", "preconditioned,unpreconditioned" ) ;
      if( nn == "preconditioned" )
      {
         CVG_TEST->set_norm_type( LA_ConvergenceTest::Preconditioned ) ;
      }
      else if( nn == "unpreconditioned" )
      {
         CVG_TEST->set_norm_type( LA_ConvergenceTest::Unpreconditioned ) ;
      }
   }

   PEL_CHECK_POST( IMPLIES( !exp->has_entry( "norm_type" ),
                      norm_type() == LA_ConvergenceTest::Preconditioned ) ) ;
}

//----------------------------------------------------------------------
size_t
LA_IterativeSolver:: max_nb_iterations( void ) const
//----------------------------------------------------------------------
{
   return( MAXITS ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: apply_prec( size_t iter,
                                 LA_Preconditioner* prec,
                                 LA_Vector const* rhs,
                                 LA_Vector* sol,
                                 bool& success )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: apply_prec" ) ;
   PEL_CHECK_PRE( prec != 0 ) ;
   PEL_CHECK_PRE( rhs != 0 ) ;
   PEL_CHECK_PRE( sol != 0 ) ;
   PEL_CHECK_PRE( prec->is_valid() ) ;
   PEL_CHECK_PRE( rhs->nb_rows() == prec->dimension() ) ;
   PEL_CHECK_PRE( sol->nb_rows() == prec->dimension() ) ;

   prec->solve( rhs, sol ) ;
   success = prec->successful_solve() ;
   if( !success )
   {
      set_converged_reason( iter, LA_ConvergenceTest::PreconditionerFailure ) ;
   }

   PEL_CHECK_POST( EQUIVALENT( !success,
       ( converged_reason() == LA_ConvergenceTest::PreconditionerFailure ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: set_converged_reason(
                                   size_t iter,
                                   LA_ConvergenceTest::ConvergedReason reason )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: set_converged_reason" ) ;

   LAST_IT = iter ;
   CVG_TEST->set_converged_reason( reason ) ;

   PEL_CHECK_POST( converged_reason() == reason ) ;
}

//----------------------------------------------------------------------
void
LA_IterativeSolver:: test_convergence( size_t iter,
                                       double r_norm,
                                       LA_Matrix const* A,
                                       LA_Vector const* b,
                                       LA_Preconditioner* prec,
                                       bool zero_initial_guess )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_IterativeSolver:: test_convergence" ) ;
   PEL_CHECK_PRE( r_norm >= 0.0 ) ;
   PEL_CHECK_PRE( A != 0 ) ;
   PEL_CHECK_PRE( b != 0 ) ;
   PEL_CHECK_PRE( prec != 0 ) ;
   PEL_CHECK_PRE( A->nb_rows() == A->nb_cols() ) ;
   PEL_CHECK_PRE( A->nb_cols() == b->nb_rows() ) ;
   PEL_CHECK_PRE( prec->dimension() == A->nb_rows() ) ;
   PEL_CHECK_PRE(
    IMPLIES( ( iter == 0 ), ( nb_iterations_achieved() == PEL::bad_index() ) ) ) ;
   PEL_CHECK_PRE(
    IMPLIES( ( iter != 0 ), ( ( iter == nb_iterations_achieved() ) ||
                              ( iter == nb_iterations_achieved()+1 ) ) ) ) ;

   RNORM = r_norm ;
   LAST_IT = iter ;

   if( VERBOSE )
   {
      MONITOR->monitor( iter, r_norm ) ;
   }
   CVG_TEST->test_convergence( iter, r_norm,
                               A, b, prec, zero_initial_guess ) ;
}

//----------------------------------------------------------------------
bool
LA_IterativeSolver:: do_solve_PRE( LA_Matrix const* A,
                                   LA_Vector const* b,
                                   LA_Preconditioner const* prec,
                                   bool zero_initial_guess,
                                   LA_Vector* x ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( converged_reason() ==
               LA_ConvergenceTest::ConvergedIterating ) ;
   PEL_ASSERT( !convergence_achieved() ) ;
   PEL_ASSERT( nb_iterations_achieved() == PEL::bad_index() ) ;
   PEL_ASSERT( A!=0 && b!=0 && prec!=0 && x!=0 ) ;
   PEL_ASSERT( A->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == b->nb_rows() ) ;
   PEL_ASSERT( A->nb_cols() == x->nb_rows() ) ;
   PEL_ASSERT( prec->dimension() == A->nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_IterativeSolver:: create_replica_PRE( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_IterativeSolver:: create_replica_POST(
                                    LA_IterativeSolver const* result,
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
LA_IterativeSolver:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "LA_IterativeSolver descendant" ) ;
   return( result ) ;
}
