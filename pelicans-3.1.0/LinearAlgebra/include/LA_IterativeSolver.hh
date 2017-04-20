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

#ifndef LA_ITERATIVE_SOLVER_HH
#define LA_ITERATIVE_SOLVER_HH

#include <PEL_Object.hh>

#include <LA_ConvergenceTest.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class LA_ConvergenceMonitor ;
class LA_Preconditioner ;
class LA_Matrix ;
class LA_Vector ;

/*
Preconditioned iterative solvers of linear systems.

FRAMEWORK INSTANTIATION
*/

class PEL_EXPORT LA_IterativeSolver : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_IterativeSolver* make( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;
      
      virtual LA_IterativeSolver* create_clone( 
                                         PEL_Object* a_owner ) const = 0 ;

   //-- Characteristics
         
      std::string const& name( void ) const ;
      
      // type of the norm that is passed to 
      // `LA_ConvergenceTest::test_convergence'
      LA_ConvergenceTest::NormType norm_type( void ) const ;
      
   //-- Iteration
      
      // Solve the system defined by `A' as the LHS matrix and by `b' as the
      // RHS, and copy the solution into `x' ; use preconditioner `prec'. If
      // `zero_initial_guess' is false, use incoming `x' as the initial 
      // guess.
      void solve( LA_Matrix const* A,
                  LA_Vector const* b,
                  LA_Preconditioner* prec,
                  bool zero_initial_guess,
                  LA_Vector* x ) ;

      LA_ConvergenceTest::ConvergedReason converged_reason( void ) const ;
      
      // Is convergence achieved at `::solve' completion ?
      bool convergence_achieved( void ) const ;

      // number of iterations performed at `::solve' completion
      size_t nb_iterations_achieved( void ) const ;
      
      // residual norm, as defined by `::norm_type()', 
      // achieved at `::solve' completion
      double residual_norm_achieved( void ) const ;
     
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      virtual void print_more( std::ostream& os, 
                               size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   //-- Plug in
      
      virtual ~LA_IterativeSolver( void ) ;

      LA_IterativeSolver( std::string const& a_name ) ;

      LA_IterativeSolver( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp ) ;

      LA_IterativeSolver( PEL_Object* a_owner, 
                          LA_IterativeSolver const* other ) ;

      virtual LA_IterativeSolver* create_replica( 
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;
      
   //-- Internals
      
      void set_norm_type( PEL_ModuleExplorer const* exp ) ;
      
      size_t max_nb_iterations( void ) const ;
      
      void apply_prec( size_t iter,
                       LA_Preconditioner* prec,
                       LA_Vector const* rhs,
                       LA_Vector* sol,
                       bool& success ) ;
      
      virtual void do_solve( LA_Matrix const* A,
                             LA_Vector const* b,
                             LA_Preconditioner* prec,
                             bool zero_initial_guess,
                             LA_Vector* x ) = 0 ;
      
      void set_converged_reason(
                       size_t iter,
                       LA_ConvergenceTest::ConvergedReason reason ) ;
                  
      void test_convergence( size_t iter, 
                             double r_norm,
                             LA_Matrix const* A, 
                             LA_Vector const* b,
                             LA_Preconditioner* prec,
                             bool zero_initial_guess ) ;
            
   //-- Preconditions, Postconditions, Invariant    

      bool do_solve_PRE( LA_Matrix const* A, 
                         LA_Vector const* b,
                         LA_Preconditioner const* prec,
                         bool zero_initial_guess,
                         LA_Vector* x ) const ;

      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( LA_IterativeSolver const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;
      
   private: //----------------------------------------------------------

      LA_IterativeSolver( void ) ;
      LA_IterativeSolver( LA_IterativeSolver const& other ) ; 
      LA_IterativeSolver& operator=( LA_IterativeSolver const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      bool const IS_PROTO ;
      std::string SOLVER_NAME ;
      size_t MAXITS ;
      LA_ConvergenceTest* CVG_TEST ;
      LA_ConvergenceMonitor* MONITOR ;
      double RNORM ;
      size_t LAST_IT ;
      bool VERBOSE ;
} ;

#endif
