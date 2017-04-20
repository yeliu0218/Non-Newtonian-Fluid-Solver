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

#include <AP_AdvectionDiffusion1G.hh>

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

using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_AdvectionDiffusion1G const* 
AP_AdvectionDiffusion1G::PROTOTYPE = new AP_AdvectionDiffusion1G() ;

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1G:: AP_AdvectionDiffusion1G( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_AdvectionDiffusion1G" )
{
}

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1G*
AP_AdvectionDiffusion1G:: create_replica( PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_AdvectionDiffusion1G* result = 
               new AP_AdvectionDiffusion1G( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1G:: AP_AdvectionDiffusion1G( 
                                               PEL_Object* a_owner,
                                               PDE_DomainAndFields const* dom,
                                               FE_SetOfParameters const* prms,
                                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , NB_DIMS( dom->nb_space_dimensions() )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data( "field" ) ) )
   , L_UU( exp->int_data( "level_of_field" ) )
   , UU_EXP( 0 )
   , L_UU_EXP( PEL::bad_index() )
   , UU_EXP_EXP( 0 )
   , L_UU_EXP_EXP( PEL::bad_index() )
   , TDISC( NoTime )
   , AA( 0 )
   , ALPHA( 0 )
   , GAMMA( 0 )
   , KAPPA( prms->item( exp->string_data( "param_diffusion" ) ) )
   , PI( 0 )
   , BCs( dom->set_of_boundary_conditions() )
   , LOCAL_BC( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , TAU( 0 )
   , STAB( NoStab )
   , SGNTAU( 0.0 )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , NMB( 0 )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , PREC( 0 )
   , TIMER_SET_MATRIX( PEL_Timer::create( this ) )
   , TIMER_ELM( PEL_Timer::create( this ) )
   , TIMER_GLOB( PEL_Timer::create( this ) )
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: AP_AdvectionDiffusion1G" ) ;
   
   check_field_storage_depth( UU, L_UU ) ;
   
   size_t nb_c = UU->nb_components() ;

   cFE->require_field_calculation( UU, PDE_LocalFE::N )  ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::N )  ;

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

   PDE_SetOfDiscreteFields const* sdf = dom->set_of_discrete_fields() ;
   if( exp->has_module( "time_discretization") )
   {
      PEL_ModuleExplorer const* se = 
                         exp->create_subexplorer( 0, "time_discretization") ;
      UU_EXP = sdf->item( se->string_data( "field_explicit" ) ) ;
      L_UU_EXP = se->int_data( "level_of_field_explicit" ) ;      
      check_field_nb_components( UU_EXP, nb_c ) ;
      check_field_storage_depth( UU_EXP, L_UU_EXP ) ;
      
      cFE->require_field_calculation( UU_EXP, PDE_LocalFE::N )  ;
      
      ALPHA = prms->item( se->string_data( "param_unsteady" ) ) ;
      check_param_nb_components( ALPHA, "param_unsteady", 1 ) ;
      ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      
      string const& tt = se->string_data( "type" ) ;
      if( tt == "Euler" )
      {
         TDISC = Euler ;
      }
      else if( tt == "BDF2" )
      {
         TDISC = BDF2 ;
         UU_EXP_EXP = sdf->item( 
                           se->string_data( "field_explicit_explicit" ) ) ;
         check_field_nb_components( UU_EXP_EXP, nb_c ) ;
         L_UU_EXP_EXP = se->int_data( "level_of_field_explicit_explicit" ) ;      
         check_field_storage_depth( UU_EXP_EXP, L_UU_EXP_EXP ) ;
         cFE->require_field_calculation( UU_EXP_EXP, PDE_LocalFE::N )  ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value( se, "type",
                                                    "\"Euler\" or \"BDF2\"" ) ;
      }
      se->destroy() ; se = 0 ;
   }
   
   if( exp->has_module( "advection" ) )
   {
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "advection" ) ;
      AA = prms->item( e->string_data( "param_advective_velocity" ) ) ;
      check_param_nb_components( AA, "param_advective_velocity", NB_DIMS ) ;
      AA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      
      if( e->has_module( "stabilization" ) )
      {
         PEL_ModuleExplorer* se = e->create_subexplorer( 0, "stabilization" ) ;
         string const& type = se->string_data( "type" ) ;
         if( type == "SUPG" ) 
         {
            STAB = SUPG ;
         }
         else if( type == "GLS" )
         {
            STAB = GLS ;
            SGNTAU = 1.0 ;      
         }
         else if( type == "USFEM" )
         {
            STAB = USFEM ;
            SGNTAU = -1.0 ; 
         }
         else
         {
            PEL_Error::object()->raise_bad_data_value( se, "type",
                                 "\"SUPG\" or \"GLS\" or \"USFEM\"" ) ;
         }

         PEL_ModuleExplorer const* sse = 
                            se->create_subexplorer( se, "MI_TauStab" ) ;
         TAU = FE_TauStab::create( this, sse ) ;
         cFE->require_field_calculation( UU, PDE_LocalFE::d2N ) ;
         se->destroy() ;
      }
      e->destroy() ;
   }

   if( exp->has_entry( "boundary_conditions_types" ) )
   {
      stringVector const& bt = 
                     exp->stringVector_data( "boundary_conditions_types" ) ;
      string const& fname = exp->string_data( "velocity" ) ;
      LOCAL_BC = FE_LocalBCsBuilder::create( this, dom, fname, bt, prms ) ;
      LOCAL_BC->transfer_calculation_requirements( bFE ) ;
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
AP_AdvectionDiffusion1G:: ~AP_AdvectionDiffusion1G( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_AdvectionDiffusion1G:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_AdvectionDiffusion1G:: do_one_inner_iteration" ) ;
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
   
   loop_on_cells( t_it ) ;
   loop_on_bounds( t_it ) ;
   
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
         "*** AP_AdvectionDiffusion1G error\n"
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
AP_AdvectionDiffusion1G:: save_other_than_time_and_fields(
                                            FE_TimeIterator const* t_it,
                                            PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   if( PREC != 0 )
   {
      start_total_timer(
        "FE_RefinementAdvectionDiffusion:: save_other_than_time_and_fields" ) ;

      rs->save_variable( (double)PREC->nb_cycles_performed(), "XNBC" ) ;

      stop_total_timer() ;
   }
}

//----------------------------------------------------------------------
void
AP_AdvectionDiffusion1G:: print_additional_times( std::ostream& os,
                                                  size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: print_additional_times" ) ;
   
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

//------------------------------------------------------------------------
void
AP_AdvectionDiffusion1G:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: loop_on_cells" ) ;
   
   double m_xx, norm_aa, kappa ;
   size_t nb_c = UU->nb_components() ;
   doubleVector rhs( nb_c ) ;
   doubleVector aa( NB_DIMS ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      TIMER_ELM->start() ;
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_c,
                              cFE->col_field_node_connectivity(), nb_c ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         compute_coefs_at_IP( t_it, m_xx, aa, norm_aa, kappa, rhs );
         
         if( AA != 0 )
         {
            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;
         }
         
         // *** assembly optimization
         // the following calls to FE member functions is replaced by
         // a fully reimplemented version which is more efficient
         //
         // FE::add_row_col_S( ELEMENT_EQ, cFE, m_xx ) ;
         // FE::add_row( ELEMENT_EQ, cFE, rhs ) ;
         // FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, kappa ) ;
         
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
         
         // *** end of assembly optimization
         
         double taustab = PEL::bad_double() ;
         if( TAU != 0 ) 
         {
            double h = cFE->polyhedron()->inter_vertices_maximum_distance() ;
            taustab = TAU->tau( h, m_xx, norm_aa, kappa ) ;
         }
         if( STAB == SUPG || STAB == GLS || STAB == USFEM )
         {
            FE::add_vvgrad_row_col( ELEMENT_EQ, cFE, aa, taustab*m_xx ) ;

            FE::add_vvgrad_row_vvgrad_col( ELEMENT_EQ, cFE, aa, taustab ) ;

            if( FE::geometry() != FE::axisymmetrical )
               FE::add_vvgrad_row_lapl_col( ELEMENT_EQ, cFE, aa,
                                            - taustab*kappa ) ;

            FE::add_vvgrad_row( ELEMENT_EQ, cFE, aa, rhs, taustab ) ;
         }

         if( STAB == GLS || STAB == USFEM  )
         {
            FE::add_row_col_S( ELEMENT_EQ, cFE,
                                SGNTAU * taustab * m_xx * m_xx ) ;

            FE::add_row( ELEMENT_EQ, cFE, rhs, SGNTAU * taustab * m_xx ) ;

            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa,
                                    SGNTAU * taustab * m_xx  ) ;

            FE::add_lapl_row_vvgrad_col( ELEMENT_EQ, cFE, aa,
                                         SGNTAU * taustab * kappa* m_xx  ) ; 

            FE::add_row_lapl_col( ELEMENT_EQ, cFE, 
                                  - SGNTAU * taustab * kappa* m_xx   ) ;

            FE::add_lapl_row_col( ELEMENT_EQ, cFE, 
                                  - SGNTAU * taustab * kappa* m_xx  ) ;

            FE::add_lapl_row_lapl_col( ELEMENT_EQ, cFE, 
                                       SGNTAU * taustab * kappa* kappa  ) ;

            FE::add_lapl_row( ELEMENT_EQ, cFE, rhs, -SGNTAU *taustab * m_xx ) ;
         }
      }
      TIMER_ELM->stop() ;
      
      TIMER_GLOB->start() ;
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
      TIMER_GLOB->stop() ;
   }
}

//---------------------------------------------------------------------------
void
AP_AdvectionDiffusion1G:: loop_on_bounds( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: loop_on_bounds" ) ;
   
   if( LOCAL_BC != 0 )
   {
      for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
      {
         GE_Color const* color = bFE->color() ;
         if( BCs->has_BC( color, UU ) )
         {
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;
            LOCAL_BC->set_current_BC_type( ee->string_data( "type" ) ) ;
            if( LOCAL_BC->current_BC_type_is_ok() )
            {
               bFE->set_row_and_col_fields( UU, UU ) ;
               ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 1,
                                       bFE->col_field_node_connectivity(), 1 );

               LOCAL_BC->build_current_BC( ELEMENT_EQ, bFE, t_it, QRP ) ;
 
               PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
            }
         }
      }
   }
}

// coefs for the equation
//    m_xx u + aa.grad u - div( kappa grad u ) = rhs
//---------------------------------------------------------------------------
void
AP_AdvectionDiffusion1G:: compute_coefs_at_IP( FE_TimeIterator const* t_it, 
                                               double& m_xx,
                                               doubleVector& aa,
                                               double& norm_aa,
                                               double& kappa, 
                                               doubleVector& rhs ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1G:: compute_coefs_at_IP" ) ;
      
   size_t nb_c = UU->nb_components() ;
   
   m_xx = 0.0 ;
   kappa  = KAPPA->cell_value_at_IP( t_it, cFE ) ;
   rhs.set( 0.0 ) ;
   
   if( PI != 0 )
   {
      for( size_t ic=0 ; ic<nb_c ; ++ic )
      {
         rhs( ic ) += PI->cell_value_at_IP( t_it, cFE, ic ) ;
      }
   }
   
   if( ALPHA != 0 )
   {
      double dt = t_it->time_step() ;   
      doubleVector uu_exp( nb_c, PEL::bad_double() ) ;
      if( TDISC == Euler || ( TDISC == BDF2 && t_it->iteration_number() == 1 ) )
      {
         for( size_t ic=0 ; ic<nb_c ; ++ic )
         {
            uu_exp( ic ) = cFE->value_at_IP( UU_EXP, L_UU_EXP, ic ) ;
         }
      }
      else if( TDISC == BDF2 && t_it->iteration_number() != 1 )
      {
         for( size_t ic=0 ; ic<nb_c ; ++ic )
         {
            uu_exp( ic ) = 
               2.0 * cFE->value_at_IP( UU_EXP, L_UU_EXP, ic )
             - 0.5 * cFE->value_at_IP( UU_EXP_EXP, L_UU_EXP_EXP, ic ) ;
         }
      }
      double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
      m_xx += alpha/dt ;
      if( TDISC == BDF2 && t_it->iteration_number() != 1 ) m_xx *= 1.5 ;
    
      for( size_t ic=0 ; ic<nb_c ; ++ic )
      {
         PEL_ASSERT( uu_exp( ic ) != PEL::bad_double() ) ;
         rhs( ic ) += alpha * uu_exp( ic ) / dt ;
      }
   }
   
   if( GAMMA != 0 )
   {
      m_xx += GAMMA->cell_value_at_IP( t_it, cFE ) ;
   }
   
   if( AA != 0 )
   {
      for( size_t d=0 ; d<aa.size() ; ++d )
      {
         aa( d ) = AA->cell_value_at_IP( t_it, cFE, d ) ;
      }

      norm_aa = 0.0 ;
      for( size_t d=0 ; d<aa.size() ; ++d )
      {
         norm_aa += aa( d ) * aa( d ) ;
      }
      norm_aa = PEL::sqrt( norm_aa ) ;
   }
}
