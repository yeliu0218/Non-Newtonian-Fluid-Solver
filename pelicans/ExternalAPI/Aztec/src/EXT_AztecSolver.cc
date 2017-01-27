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

#include <EXT_AztecSolver.hh>

#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>
#include <iostream>

#include <stdlib.h> // for free operator


using std::cout ;
using std::endl ;
using std::string ;

EXT_AztecSolver const* EXT_AztecSolver::PROTOTYPE = new EXT_AztecSolver() ;

//----------------------------------------------------------------------------
EXT_AztecSolver:: EXT_AztecSolver( void ) 
//----------------------------------------------------------------------------
   : LA_Solver( "EXT_AztecSolver" )
   , EXP( 0 )
{   
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_Aztec" ), val ) ;
}

//----------------------------------------------------------------------------
EXT_AztecSolver:: ~EXT_AztecSolver( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: ~EXT_AztecSolver" ) ;
   
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
EXT_AztecSolver*
EXT_AztecSolver:: create( PEL_Object* a_owner, PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: create" ) ;

   EXT_AztecSolver* result = new EXT_AztecSolver( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_AztecSolver*
EXT_AztecSolver:: create_replica( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE(  a_owner, exp ) ) ;
   
   EXT_AztecSolver* result = new EXT_AztecSolver( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_AztecSolver:: EXT_AztecSolver( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner )
   , EXP( exp->create_clone( this ) )
   , update( 0 )
   , my_Ai( 0 )
   , my_Av( 0 )
{
   PEL_LABEL( "EXT_AztecSolver:: EXT_AztecSolver" ) ;
   raise_fatal_error_if_not_sequential() ;
   AZ_set_proc_config(proc_config, AZ_NOT_MPI);
   initialize() ;
}

//----------------------------------------------------------------------------
EXT_AztecSolver*
EXT_AztecSolver:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: create_clone" ) ;
   
   EXT_AztecSolver* result = new EXT_AztecSolver( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_AztecSolver:: EXT_AztecSolver( PEL_Object* a_owner,
                                   EXT_AztecSolver const* other ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , EXP( other->EXP->create_clone( this ) )
   , update( 0 )
   , my_Ai( 0 )
   , my_Av( 0 )
{
   PEL_LABEL( "EXT_AztecSolver:: EXT_AztecSolver" ) ;
   raise_fatal_error_if_not_sequential() ;
   initialize() ;
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: initialize( void ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: initialize" ) ;
   
   AZ_set_proc_config(proc_config, AZ_NOT_MPI);

   PEL_ModuleExplorer* sexp =
                  EXP->create_subexplorer( 0, "Aztec_IterativeSolver" ) ;
   build_ksp( sexp ) ;
   sexp->destroy() ; sexp = 0 ;
   
   sexp = EXP->create_subexplorer( 0, "Aztec_Preconditioner" ) ;
   build_pc( sexp ) ;
   sexp->destroy() ; sexp = 0 ;
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string s( indent_width, ' ') ;
   
   // Solver :
   std::string solver ;
   if( options[AZ_solver]==AZ_gmres )
   {
      solver = "AZ_GMRES" ;
   }
   else if( options[AZ_solver]==AZ_cg )
   {
      solver = "AZ_CG" ;
   }
   else if( options[AZ_solver]==AZ_cgs )
   {
      solver = "AZ_CGS" ;
   }
   else if( options[AZ_solver]==AZ_tfqmr )
   {
      solver = "AZ_TFQMR" ;
   }
   else if( options[AZ_solver]==AZ_bicgstab )
   {
      solver = "AZ_BICGSTAB" ;
   }
   else if( options[AZ_solver]==AZ_lu )
   {
      solver = "AZ_LU" ;
   }
   os << s << "linear solver : \"" << solver << "\"" << std::endl ;
   os << s << "linear solver relative error bound : "
      << params[AZ_tol] << std::endl ;
   os << s << "linear solver maximal number of iterations : "
      << options[AZ_max_iter] << std::endl ;
   if( options[AZ_solver]==AZ_gmres )
   {
      os << s << "iteration of restart : "
         << options[AZ_kspace] << std::endl ;
   }
   
   // Preconditioner :
   std::string prec = "none" ;
   if( options[AZ_precond]==AZ_Jacobi )
   {
      prec = "AZ_Jacobi" ;
   }
   else if( options[AZ_precond]==AZ_Neumann )
   {
      prec = "AZ_Neumann" ;
   }
   else if( options[AZ_precond]==AZ_ls )
   {
      prec = "AZ_ls" ;
   }
   else if( options[AZ_precond]==AZ_sym_GS )
   {
      prec = "AZ_symGS" ;
   }
   else if( options[AZ_precond]==AZ_dom_decomp &&
            options[AZ_subdomain_solve]==AZ_ilu )
   {
      prec = "AZ_LU" ;
   }
   else if( options[AZ_precond]==AZ_dom_decomp &&
            options[AZ_subdomain_solve]==AZ_icc )
   {
      prec = "AZ_ICC" ;
   }
   else if( options[AZ_precond]==AZ_dom_decomp &&
            options[AZ_subdomain_solve]==AZ_ilut )
   {
      prec = "AZ_ILUT" ;
   }
   else if( options[AZ_precond]==AZ_dom_decomp &&
            options[AZ_subdomain_solve]==AZ_rilu )
   {
      prec = "AZ_RILU" ;
   }
   os << s << "preconditioner : \"" << prec << "\"" << std::endl ;
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: build_ksp( PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: build_ksp" ) ;
   PEL_CHECK( exp != 0 ) ;
   
   string algo = exp->string_data( "type" ) ;
   double toler = exp->double_data( "relative_tolerance" ) ;
   int maxit = exp->int_data( "nb_iterations_max" ) ;
   AZ_defaults(options, params);
   set_iterative( true );
   if( algo=="GMRES" )
   {
      options[AZ_solver] = AZ_gmres ;
      options[AZ_kspace] = exp->int_data( "restart" ) ;
   }
   else if( algo=="CG" )
   {
      options[AZ_solver] = AZ_cg ;
   }
   else if( algo=="CGS" )
   {
      options[AZ_solver] = AZ_cgs ;
   }
   else if( algo=="TFQMR" )
   {
      options[AZ_solver] = AZ_tfqmr ;
   }
   else if( algo=="BICGSTAB" )
   {
      options[AZ_solver] = AZ_bicgstab ;
   }
   else if( algo=="LU" )
   {
      set_iterative( false );
      options[AZ_solver] = AZ_lu ;
   }
   else
   {
      PEL_Error::object()->raise_plain(
         "Unknown Aztec solver (not in CG,GMRES,CGS,TFQMR,BICGSTAB,LU) : "+algo ) ;
   }
   options[AZ_scaling] = AZ_none ;
   options[AZ_conv] = AZ_rhs ;
   if( exp->has_entry( "residual_norm" ) )
   {
      std::string const& conv = exp->string_data( "residual_norm" ) ;
      if( conv == "r0" )
      {
         options[AZ_conv] = AZ_r0 ;
      }
      else if( conv == "rhs" )
      {
         options[AZ_conv] = AZ_rhs ;
      }
      else if( conv == "Anorm" )
      {
         options[AZ_conv] = AZ_Anorm ;
      }
      else if( conv == "noscaled" )
      {
         options[AZ_conv] = AZ_noscaled ;
      }
      else if( conv == "sol" )
      {
         options[AZ_conv] = AZ_sol ;
      }
      else if( conv == "weighted" )
      {
         options[AZ_conv] = AZ_weighted ;
      }
      else
      {
         PEL_Error::object()->raise_plain(
            "residual_norm be in : r0, rhs, Anorm, noscaled, sol, weighted " ) ;
      }
   }
   if( exp->has_entry( "verbose" ) && exp->bool_data( "verbose" ) )
   {
      options[AZ_output] = 10 ;
   }
   else
   {
      options[AZ_output] = AZ_none ;
   }
   options[AZ_max_iter] = maxit ;
   params[AZ_tol] = toler ;
   params[AZ_orthog] = AZ_modified ;

   if( exp->has_entry( "scaling" ) )
   {
      std::string const& scal = exp->string_data( "scaling" ) ;
      if( scal == "none" )
      {
         options[AZ_scaling] = AZ_none ;
      }
      else if( scal == "jacobi" )
      {
         options[AZ_scaling] = AZ_Jacobi ;
      }
      else if( scal == "row_sum" )
      {
         options[AZ_scaling] = AZ_row_sum ;
      }
      else if( scal == "sym_diag" )
      {
         options[AZ_scaling] = AZ_sym_diag ;
      }
      else if( scal == "sym_row_sum" )
      {
         options[AZ_scaling] = AZ_sym_row_sum ;
      }
      else
      {
         PEL_Error::object()->raise_plain(
            "scaling must be in : none, jacobi, row_sum, sym_diag, sym_row_sum " ) ;
      }
      
   }
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: build_pc( PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: build_pc" ) ;
   PEL_CHECK( exp != 0 ) ;
   
   string precond = exp->string_data( "type" ) ;
   if( precond=="Jacobi" )
   {
      options[AZ_precond] = AZ_Jacobi ;
      options[AZ_poly_ord] = exp->int_data( "steps" ) ;
   }
   else if( precond=="Neumann" )
   {
      options[AZ_precond] = AZ_Neumann ;
      options[AZ_poly_ord] = exp->int_data( "order" ) ;
   }
   else if( precond=="ls" )
   {
      options[AZ_precond] = AZ_ls ;
      options[AZ_poly_ord] = exp->int_data( "order" ) ;
   }
   else if( precond=="symGS" )
   {
      options[AZ_precond] = AZ_sym_GS ;
      options[AZ_poly_ord] = exp->int_data( "steps" ) ;
   }
   else if( precond=="LU" )
   {
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilu;
      params[AZ_drop] = exp->double_data( "drop_tolerance" ) ; ;
   }
   else if( precond=="ILU" )
   {
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilu;
      options[AZ_graph_fill] = exp->int_data( "order" ) ; ;
      
   }
   else if( precond=="ICC" )
   {
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_icc;
      options[AZ_graph_fill] = exp->int_data( "order" ) ; ;
      
   }
   else if( precond=="ILUT" )
   {
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilut ;
      params[AZ_drop] = exp->double_data( "drop_tolerance" ) ; ;
      params[AZ_ilut_fill] = exp->double_data( "fillin" ) ; ;
   }
   else if( precond=="RILU" )
   {
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_rilu ;
      params[AZ_omega] = exp->double_data( "omega" ) ; ;
      params[AZ_ilut_fill] = exp->double_data( "fillin" ) ; ;
   }
   else if( precond!="none" )
   {
      PEL_Error::object()->raise_plain(
         "Unknown Aztec preconditioner (not in Jacobi,Neumann,ls,symGS,LU,ILU,ILUT,RILU,ICC) : "
         +precond ) ;
   }

}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: unset_matrix_self( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: unset_matrix_self" ) ;
   PEL_CHECK( unset_matrix_self_PRE() ) ;
   
   if( update!=0 )
   {
      free((void *) update);   free((void *) update_index);
      free((void *) external); free((void *) extern_index);
      free((void *) data_org);
      update = 0 ;
   }
   if( my_Ai!=0 )
   {
      delete [] my_Ai ; my_Ai=0 ;
      delete [] my_Av ; my_Av=0 ;
   }
      
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: set_matrix_self( LA_Matrix const* mat,
                                   bool &ok, bool same_pattern )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   if( dynamic_cast<LA_SeqMatrix const*>( mat ) == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** EXT_AztecSolver:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected\n" ) ;
   }
   
   LA_SeqMatrix const* smat = static_cast<LA_SeqMatrix const*>( mat ) ;
   
   nrow = smat->nb_rows() ;
   
   size_t nnz = smat->nb_stored_items() ;
   
   LA_MatrixIterator* it = smat->create_stored_item_iterator( 0 ) ;
  
   AZ_read_update(&N_update, &update, proc_config, nrow, 1, AZ_linear);

   my_Ai = new int [ nnz+1 ] ;  
   my_Av   = new double [ nnz+1 ] ;

   my_Ai[0] = N_update+1;
   
   for( int i=0 ; i< N_update; i++ )
   {
      int global_row = update[i] ;
      int k = my_Ai[i] ;
      for( it->start_row_items(i) ;
           it->is_valid() ;
           it->go_next() )
      {
         size_t j = it->col() ;
         double Aij = it->item() ;
         if( global_row==j )
         {
            my_Av[i] = Aij ;
         }
         else
         {
            my_Ai[k]=j ;
            my_Av[k++] = Aij ;
         }
      }
      my_Ai[i+1] = k ;
      PEL_ASSERT( k<=nnz+1 ) ;
   }
   it->destroy() ; it = 0 ;
   
   AZ_transform(proc_config, &external, my_Ai, my_Av, update, &update_index,
                &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
                AZ_MSR_MATRIX);
   ok = true ;

   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------------
void
EXT_AztecSolver:: solve_self( LA_Vector const* b, LA_Vector* x,
                              size_t &nb_iter,  bool &ok ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_AztecSolver:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   
   PEL_CHECK( my_Ai!=0 ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( b ) != 0 ) ;
   LA_SeqVector const* brhs = static_cast<LA_SeqVector const* >( b ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;
   
   double* my_rhs = new double [N_update] ;
   for( int i=0 ; i<N_update ; i++ )
   {
      my_rhs[i] = brhs->item(update[i]) ;
   }
   
   double * my_x = new double [ N_update ] ;
   if( zero_initial_guess() )
   {
      for( size_t i=0 ; i<N_update ; i++ )
         my_x[i] = 0.0 ;
   }
   else
   {
      for( size_t i=0 ; i<N_update ; i++ )
         my_x[i] = bx->item(i) ;
   }
   
   AZ_solve(my_x, my_rhs, options, params, NULL, my_Ai, NULL, NULL, NULL,
            my_Av, data_org, status, proc_config);
   for( size_t i=0 ; i<N_update ; i++ )
   {
      x->set_item( update[i], my_x[i] ) ;
   }
   delete [] my_x ;
   delete [] my_rhs ;
   
   nb_iter = (size_t) status[AZ_its];
   ok = status[AZ_why]==AZ_normal ;

   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}
