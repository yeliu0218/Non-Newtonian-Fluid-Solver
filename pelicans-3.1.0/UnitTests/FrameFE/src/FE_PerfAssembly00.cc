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

#include <FE_PerfAssembly00.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GeometricMultilevel_PC.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_LocalBCsBuilder.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TauStab.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <vector>

using std::endl ;
using std::string ;
using std::vector ;

FE_PerfAssembly00 const* 
FE_PerfAssembly00:: PROTOTYPE = new FE_PerfAssembly00() ;

//---------------------------------------------------------------------------
FE_PerfAssembly00:: FE_PerfAssembly00( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_PerfAssembly00" )
{
}

//---------------------------------------------------------------------------
FE_PerfAssembly00*
FE_PerfAssembly00:: create_replica( PEL_Object* a_owner,
                                    PDE_DomainAndFields const* dom,
                                    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_PerfAssembly00* result = 
                      new FE_PerfAssembly00( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_PerfAssembly00:: FE_PerfAssembly00( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , NB_DIMS( dom->nb_space_dimensions() )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data( "field" ) ) )
   , L_UU( exp->int_data( "level_of_field" ) )
   , GAMMA( 0 )
   , KAPPA( prms->item( exp->string_data( "param_diffusion" ) ) )
   , PI( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , ID_ALGO( exp->int_data( "loop_on_cells" ) )
   , NMB( 0 )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , PREC( 0 )
   , TIMER_SET_MATRIX( PEL_Timer::create( this ) )
   , TIMER_ELM( PEL_Timer::create( this ) )
   , TIMER_GLOB( PEL_Timer::create( this ) )
{
   PEL_LABEL( "FE_PerfAssembly00:: FE_PerfAssembly00" ) ;
   
   exp->test_data_in( "loop_on_cells", "1,2" ) ;
   
   check_field_storage_depth( UU, L_UU ) ;
   
   size_t nb_c = UU->nb_components() ;

   cFE->require_field_calculation( UU, PDE_LocalFE::N )  ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;

   check_param_nb_components( KAPPA, "param_diffusion", 1 ) ;
   KAPPA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   
   if( exp->has_entry( "param_source" ) )
   {
      PI = prms->item( exp->string_data( "param_source" ) ) ;
      check_param_nb_components( PI, "param_source", nb_c ) ;
      PI->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   }
   
   if( exp->has_entry( "param_reaction" ) )
   {
      GAMMA = prms->item( exp->string_data( "param_reaction" ) ) ;
      check_param_nb_components( GAMMA, "param_reaction", 1 ) ;
      GAMMA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   }

   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
   NMB = PDE_SystemNumbering::create( this, uu_link ) ;
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   std::string nn = "multilevel_preconditioner_name" ;
   configure_multilevel_preconditioner( exp, nn, dom, NMB ) ;

   // only for saving the number of iterations
   if( exp->has_entry( "multilevel_preconditioner_name" ) )
   {
      PREC = PDE_GeometricMultilevel_PC::object(
                 exp->string_data( "multilevel_preconditioner_name" ) ) ;
   }
}

//---------------------------------------------------------------------------
FE_PerfAssembly00:: ~FE_PerfAssembly00( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
FE_PerfAssembly00:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_PerfAssembly00:: do_one_inner_iteration" ) ;
   // --------------

   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   start_assembling_timer() ;
   // -------------------
   
   if( ID_ALGO == 1 )
   {
      loop_on_cells_1( t_it ) ;
   }
   else if( ID_ALGO == 2 )
   {
      loop_on_cells_2( t_it ) ;
   }
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ----------------
   
   SOLVER->set_initial_guess_nonzero( false ) ;
   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
   
   stop_solving_timer() ;
   // ---------------
   
   if( ! SOLVER->solution_is_achieved() )
   {
      PEL_Error::object()->display_info(
         "*** FE_PerfAssembly00 error\n"
         "    No convergence of the linear solver" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      LA_Scatter const* sca = NMB->scatter() ;
      sca->get( X, X_LOC ) ;
      UU->update_free_DOFs_value( L_UU, X_LOC, NMB->link() ) ;
      UU->enforce_constraints_for_DOFs( L_UU ) ;
   }
   
   stop_total_timer() ;
   // -------------
}

//----------------------------------------------------------------------
void
FE_PerfAssembly00:: print_additional_times( std::ostream& os,
                                            size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: print_additional_times" ) ;
   
   std::string space( indent_width, ' ' ) ;
   PEL::out() << space << "FE_RefinementAdvectionDiffusion" << std::endl ;
  
   PEL::out() << space << "   build local systems : " ;
   TIMER_ELM->print( PEL::out(), 0 ) ;
   PEL::out() << std::endl ;
   PEL::out() << space << "   assemble in global  : " ;
   TIMER_GLOB->print( PEL::out(), 0 ) ;
   PEL::out() << std::endl ;
   PEL::out() << space << "   set_matrix          : " ;
   TIMER_SET_MATRIX->print( PEL::out(), 0 ) ;
   PEL::out() << std::endl ;
}

//---------------------------------------------------------------------------
void
FE_PerfAssembly00:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: print" ) ;

   FE_OneStepIteration:: print( os, indent_width ) ;

   std::string const s( indent_width+3, ' ' ) ;

   if( ID_ALGO == 1 ) os << s << "assembly using FE::" << endl ;
   if( ID_ALGO == 2 ) os << s << "assembly without FE::" << endl ;
}

//------------------------------------------------------------------------
void
FE_PerfAssembly00:: loop_on_cells_1( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: loop_on_cells_1" ) ;
   
   double m_xx, kappa ;
   size_t nb_c = UU->nb_components() ;
   doubleVector rhs( nb_c ) ;
   doubleVector aa( NB_DIMS ) ;
   
   PEL_ASSERT( UU->nb_components() == 1 ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      TIMER_ELM->start() ;
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_c,
                              cFE->col_field_node_connectivity(), nb_c ) ;
      
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         compute_coefs_at_IP( t_it, m_xx, kappa, rhs );
         
         FE::add_row_col_S( ELEMENT_EQ, cFE, m_xx ) ;
         FE::add_row( ELEMENT_EQ, cFE, rhs ) ;
         FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, kappa ) ;         
      }
      TIMER_ELM->stop() ;
      
      TIMER_GLOB->start() ;
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
      TIMER_GLOB->stop() ;
   }
}

//------------------------------------------------------------------------
void
FE_PerfAssembly00:: loop_on_cells_2( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: loop_on_cells_2" ) ;
   
   double m_xx, kappa ;
   size_t nb_c = UU->nb_components() ;
   doubleVector rhs( nb_c ) ;
   doubleVector aa( NB_DIMS ) ;
   
   PEL_ASSERT( UU->nb_components() == 1 ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      TIMER_ELM->start() ;
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_c,
                              cFE->col_field_node_connectivity(), nb_c ) ;
      
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         compute_coefs_at_IP( t_it, m_xx, kappa, rhs );
         
         doubleVector const& N_row = cFE->Ns_at_IP( PDE_LocalFE::row ) ;
         doubleVector const& N_col = cFE->Ns_at_IP( PDE_LocalFE::col ) ;
         doubleArray2D const& dN_row = cFE->dNs_at_IP( PDE_LocalFE::row ) ;
         doubleArray2D const& dN_col = cFE->dNs_at_IP( PDE_LocalFE::col ) ;
         size_t nb_nodes = cFE->nb_basis_functions( PDE_LocalFE::row ) ;
         
         double ww = cFE->weight_of_IP() ;
         double m_xx_w  = ww * m_xx ;
         double kappa_w = ww * kappa ;
         
         for( size_t i=0 ; i<nb_nodes ; ++i )
         {
            for( size_t j=i ; j<nb_nodes ; ++j )
            {
               double xx = m_xx_w * N_row( i ) * N_col( j ) ;
               for( size_t d=0 ; d<NB_DIMS ; ++d )
               {
                  xx += kappa_w * dN_row( i, d ) * dN_col( j, d ) ;
               }

               for( size_t ic=0 ; ic<nb_c ; ++ic )
               {
                  ELEMENT_EQ->add_to_matrix( xx, i, j, ic, ic ) ;
                  if( i!=j ) ELEMENT_EQ->add_to_matrix( xx, j, i, ic, ic ) ;
               }
            }
            double yy = ww * N_row( i ) ;
            for( size_t ic=0 ; ic<nb_c ; ++ic )
            {
               ELEMENT_EQ->add_to_vector( yy*rhs( ic ), i, ic ) ;
            }
         }
      }
      TIMER_ELM->stop() ;
      
      TIMER_GLOB->start() ;
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
      TIMER_GLOB->stop() ;
   }
}

// coefs for the equation
//    m_xx u - div( kappa grad u ) = rhs
//---------------------------------------------------------------------------
void
FE_PerfAssembly00:: compute_coefs_at_IP( FE_TimeIterator const* t_it, 
                                         double& m_xx,
                                         double& kappa, 
                                         doubleVector& rhs ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PerfAssembly00:: compute_coefs_at_IP" ) ;
      
   size_t nb_c = UU->nb_components() ;
   
   m_xx = GAMMA->cell_value_at_IP( t_it, cFE ) ;
   kappa  = KAPPA->cell_value_at_IP( t_it, cFE ) ;
   rhs.set( 0.0 ) ;
   
   if( PI != 0 )
   {
      for( size_t ic=0 ; ic<nb_c ; ++ic )
      {
         rhs( ic ) += PI->cell_value_at_IP( t_it, cFE, ic ) ;
      }
   }
}
