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

#include <LA_Solver_TEST.hh>

#include <LA_PreconditionedSolver.hh>
#include <LA_PelMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>

#include <intVector.hh>
#include <stringVector.hh>

#include <LA_Matrix_TEST.hh>

#include <iostream>
#include <sstream>

LA_Solver_TEST const*
LA_Solver_TEST::PROTOTYPE = new LA_Solver_TEST() ;

//----------------------------------------------------------------------
LA_Solver_TEST:: LA_Solver_TEST( void )
//----------------------------------------------------------------------
  : PEL_ObjectTest( "LA_Solver", "LA_Solver_TEST" )
{
}

//----------------------------------------------------------------------
LA_Solver_TEST:: ~LA_Solver_TEST( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_Solver_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Solver_TEST:: process_one_test" ) ;
   PEL_CHECK( texp!=0 ) ;
   
   std::string const& test_name = texp->name() ;
   PATH = "." ;
   PEL_ModuleExplorer const* dexp = data_deck_explorer() ;
   
   if( dexp->has_entry( "test_matrices_path" ) )
   {
      PATH = dexp->string_data( "test_matrices_path" ) ;
   }

   out() << "| ... " << test_name << " : " << std::endl ;

   // Build the solver :
   PEL_ModuleExplorer* solver_exp =
                texp->create_subexplorer( 0, "LA_Solver" ) ;
   LA_Solver* solver = LA_Solver::make( 0, solver_exp ) ;
   bool verbose = false ;
   if( solver_exp->has_entry( "verbose" ) )
   {
      verbose = solver_exp->bool_data( "verbose" ) ;
   }
   if( verbose )
   {
      solver->print( out(), 3 ) ;
   }
   solver_exp->destroy() ; solver_exp = 0 ;
   
   bool iterative = texp->bool_data( "is_iterative" ) ;
   
   notify_one_test_result( test_name+"/is_iterative", solver->is_iterative()==iterative ) ;   

   PEL_ModuleExplorer* mat_exp =
                texp->create_subexplorer( 0, "LA_Matrix" ) ;
   LA_Matrix* mat = LA_Matrix::make( 0, mat_exp ) ;
   mat_exp->destroy() ; mat_exp=0 ;
   
   // Tests :
   stringVector const& matrices = texp->stringVector_data( "tested_matrix" ) ;
   double error = texp->double_data( "error_bound" ) ;
   if( error<=0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         texp, "inversion_error_bound", "a positive value is expected" ) ;
   }
   
   for( size_t i=0 ; i<matrices.size() ; ++i )
   {
      run_one_test( mat, solver, PATH+PEL_System::path_name_separator()+matrices(i), error ) ;
   }
   solver->destroy() ; solver = 0 ;
   mat->destroy() ; mat = 0 ;
}

//-------------------------------------------------------------------------
void
LA_Solver_TEST:: run_one_test( LA_Matrix* mat,
                               LA_Solver* solver,
                               std::string const& matrix,
                               double error_bound )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Solver_TEST:: run_one_test" ) ;

   out() << "| ... " <<  matrix << std::endl ;
   mat->readMM( matrix ) ;
   LA_Vector* x = mat->create_vector( solver ) ;
   LA_Vector* b = x->create_vector( solver ) ;
   LA_Vector* x0 = x->create_vector( solver ) ;
   LA_Vector* residu = x->create_vector( solver )  ;
   if( mat->nb_rows() != mat->nb_cols() )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_Solver_TEST :\n"
         "    the matrix is not a square matrix" ) ;
   }
   for( size_t i=0 ; i<mat->nb_rows() ; i++ )
   {
      x->set_item( i, PEL::sin(6.28/(i+1)) ) ;
   }
   
   bool ok = true ;

   x->synchronize() ;
   b->synchronize() ;
   mat->multiply_vec_then_add( x, b ) ;
   
   notify_one_test_result( "b is synchonized", b->is_synchronized() ) ;

   solver->set_matrix( mat ) ;
   notify_one_test_result( "matrix is set", solver->matrix_is_set() ) ;

   solver->set_initial_guess_nonzero( false ) ;
   notify_one_test_result( "nonzero_initial_guess", 
              !solver->is_iterative() || solver->zero_initial_guess() ) ;
   solver->solve( b, x0 ) ;
   ok = solver->solution_is_achieved() ;
   size_t nit =  solver->nb_iterations_achieved() ;
   
   notify_one_test_result( "solution_is_achieved", ok ) ;
   if( !solver->is_iterative() )
      notify_one_test_result( "!iterative in 1 it", nit==1 ) ;

   solver->set_initial_guess_nonzero( true ) ;
   solver->solve( b, x0 ) ;
   ok = solver->nb_iterations_achieved()<=1 ;
   if( !ok ) out() << "   nb_iterations_achieved " << solver->nb_iterations_achieved() << std::endl ;
   
   notify_one_test_result( "re-solve", ok ) ;
   if( ok )
   {
      LA_Solver* cloned_solver = solver->create_clone(solver) ;
      ok = cloned_solver->is_iterative()==solver->is_iterative() ;
      notify_one_test_result( "clone iterative", ok ) ;
      
      cloned_solver->set_matrix( mat ) ;
      cloned_solver->set_initial_guess_nonzero( false ) ;
      cloned_solver->solve( b, x0 ) ;
      ok = cloned_solver->solution_is_achieved() &&
         cloned_solver->nb_iterations_achieved()==nit;
      cloned_solver->unset_matrix() ;
      if(!ok)
      {
         bool boo = cloned_solver->solution_is_achieved() ;
         out()<<" solution_is_achieved "<<boo<<std::endl ;
         if(boo)  out()<<" nb_iterations_achieved cloned "<<cloned_solver->nb_iterations_achieved()
                     << " <> "<<nit<<std::endl ;
      }
      
      notify_one_test_result( "create_clone", ok ) ;

      if( !ok  )
      {
         std::cout << " The solver failed !!!" << std::endl ;
      }
      else
      {
         residu->set( b ) ;
         mat->multiply_vec_then_add( x0, residu , -1.0, 1.0 ) ;
         ok = PEL::toler( residu->max_norm() / b->max_norm(), error_bound ) ;
         notify_one_test_result( "inversion verification", ok ) ;
         if( !ok )
         {
            out() << "|| Ax-b || : " << residu->max_norm() << std::endl ;
            out() << "|| Ax-b ||/||b|| : " << residu->max_norm() / b->max_norm()
                  << " ( error_bound : " << error_bound << " )" << std::endl ;
            out() << " The solver failed !!!" << std::endl ;
         }
      }
   }
   solver->unset_matrix() ;
   notify_one_test_result( "!matrix is set", !solver->matrix_is_set() ) ;
}
