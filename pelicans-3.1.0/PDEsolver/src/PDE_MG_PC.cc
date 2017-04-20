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

#include <PDE_MG_PC.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_LinkDOF2Unknown.hh>

#include <LA_Solver.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>

#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::set ;
using std::string ; using std::ostringstream ;

PDE_MG_PC const*
PDE_MG_PC::PROTOTYPE = new PDE_MG_PC() ;

//----------------------------------------------------------------------
PDE_MG_PC:: PDE_MG_PC( void )
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( "PDE_MG_PC" )
{
}

//----------------------------------------------------------------------
PDE_MG_PC*
PDE_MG_PC:: create_replica( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MG_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PDE_MG_PC* result = new PDE_MG_PC( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_MG_PC:: PDE_MG_PC( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( a_owner, exp )
   , NB_CYCLES( exp->int_data( "nb_cycles" ) )
   , NB_PRESMOO( exp->int_data( "nb_presmoothing_steps" ) )
   , NB_POSTSMOO( exp->int_data( "nb_postsmoothing_steps" ) )
   , MC( exp->int_data( "cycling_strategy" ) )
   , I_CYCLE( 0 )
   , SOLVER( 0 )
   , SOLVE_OK( false )
   , VERBOSE( 0 )
{
   PEL_LABEL( "PDE_MG_PC:: PDE_MG_PC" ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "coarse_solver" ) ;
   SOLVER = LA_Solver::make( this, se ) ;
   se->destroy() ; se=0 ;

   if( exp->has_entry( "verbose_level" ) )
      VERBOSE = exp->int_data( "verbose_level" ) ;
}

//----------------------------------------------------------------------
PDE_MG_PC:: ~PDE_MG_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PDE_MG_PC*
PDE_MG_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MG_PC:: create_clone" ) ;
   return( new PDE_MG_PC( a_owner, this ) ) ;
}

//----------------------------------------------------------------------
PDE_MG_PC:: PDE_MG_PC( PEL_Object* a_owner,
                                   PDE_MG_PC const* other )
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( a_owner, other )
{
   PEL_ASSERT( false ) ; //???????????????????????????
}

//----------------------------------------------------------------------
void
PDE_MG_PC:: solve( LA_Vector const* rhs,  LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MG_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   double norm_res_0 = PEL::bad_double() ;
   LA_Vector* resid = res_of_level( nb_levels()-1 ) ;
   if( VERBOSE > 1 )
   {
      compute_residual( finest_mat(), rhs, sol, resid ) ;
      norm_res_0 = resid->two_norm() ;
   }

   sol->nullify() ;

   I_CYCLE = 0 ;
   do
   {
      ++I_CYCLE ;

      do_one_multigrid_iteration( nb_levels()-1, finest_mat(), rhs, sol ) ;

      if( VERBOSE > 1 )
      {
         compute_residual( finest_mat(), rhs, sol, resid ) ;
         double norm_res = resid->two_norm() ;
         print_residuals( INDENT, I_CYCLE, norm_res, norm_res_0 ) ;
      }
   } while( I_CYCLE != NB_CYCLES ) ;

   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
void
PDE_MG_PC:: do_one_multigrid_iteration( size_t level,
                                        LA_Matrix const* mat,
                                        LA_Vector const* rhs,
                                        LA_Vector* sol ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MG_PC:: do_one_multigrid_iteration" ) ;

   std::string mi( "      " ) ;
   if( level == 0 )
   {
      if( VERBOSE > 2 )
         PEL::out() << INDENT << mi << "level " << level
                    << " : solution on the coarsest grid..." << endl ;
      SOLVER->set_matrix( mat ) ;
      SOLVER->solve( rhs, sol ) ;
      SOLVER->unset_matrix() ;
   }
   else
   {
      LA_Matrix const* pr = coarse_to_fine( level-1 ) ;
      LA_Vector const* to_be_smoothed = unknowns_of_level( level-1 ) ;

      // *** presmoothing
      if( VERBOSE > 2 )
         PEL::out() << INDENT << mi << "level " << level
                                   << " : presmoothing..." << endl ;

      smooth_GaussSeidel( NB_PRESMOO, to_be_smoothed, mat, rhs, sol ) ;

      // *** residual computation
      LA_Vector* rr = res_of_level( level ) ;
      compute_residual( mat, rhs, sol, rr ) ;

      LA_Vector* bb = rhs_of_level( level-1 ) ;
      LA_Vector* xx = sol_of_level( level-1 ) ;

      // *** residual restriction
      // cB = tr(Pr) * RES
      pr->tr_multiply_vec_then_add( rr, bb ) ;  // bb = t(pr) * rr

      // *** coarse grid solution (coarse grid initial guess is zero)
      xx->nullify() ;
      size_t nb_c = ( level==1 ? 1 : MC ) ;
      for( size_t i=0 ; i<nb_c ; ++i )
      {
         do_one_multigrid_iteration( level-1, mat_of_level( level-1 ), 
                                     bb, xx ) ;
      }

      // *** interpolation an error correction
      pr->multiply_vec_then_add( xx, sol, 1.0, 1.0 ) ;  // sol = sol + pr*xx

      // *** postsmoothing
      if( VERBOSE > 2 )
         PEL::out() << INDENT << mi << "level " << level
                                   << " : postsmoothing..." << endl ;

      smooth_GaussSeidel( NB_POSTSMOO, to_be_smoothed, mat, rhs, sol ) ;
   }
}

//----------------------------------------------------------------------
bool
PDE_MG_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MG_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MG_PC:: nb_cycles_performed( void ) const
//----------------------------------------------------------------------
{
   return( I_CYCLE ) ;
}
