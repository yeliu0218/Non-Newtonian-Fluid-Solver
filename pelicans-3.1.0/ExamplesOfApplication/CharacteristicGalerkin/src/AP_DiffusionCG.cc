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

#include <AP_DiffusionCG.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>

using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_DiffusionCG const* AP_DiffusionCG::PROTOTYPE = new AP_DiffusionCG() ;

//----------------------------------------------------------------------
AP_DiffusionCG:: AP_DiffusionCG( void )
//----------------------------------------------------------------------
   : FE_Galerkin( "AP_DiffusionCG" )
{
}

//----------------------------------------------------------------------
AP_DiffusionCG*
AP_DiffusionCG:: create_replica( PEL_Object* a_owner, 
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_DiffusionCG* result = new AP_DiffusionCG( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_DiffusionCG:: AP_DiffusionCG( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_Galerkin( a_owner, dom, exp )
   , TT( dom->set_of_discrete_fields()->item( "temperature" ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , DENS( prms->item( exp->string_data( "param_density" ) ) )
   , COND( prms->item( exp->string_data( "param_conductivity" ) ) )
   , CP( prms->item( exp->string_data( "param_specific_heat" ) ) ) 
   , POW( prms->item( exp->string_data( "param_volumic_power" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( 0 )
   , bFE( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "AP_DiffusionCG:: AP_DiffusionCG" ) ;

   check_field_storage_depth( TT, PEL::max( L_UPDATE, L_EXPLICIT ) ) ;

   check_param_nb_components( DENS, "param_density", 1 ) ;
   check_param_nb_components( COND, "param_conductivity", 1 ) ;
   check_param_nb_components( CP, "param_specific_heat", 1 ) ;
   check_param_nb_components( POW, "param_volumic_power", 1 ) ;

   cFE = dom->create_LocalFEcell( this ) ;
   bFE = dom->create_LocalFEbound( this ) ;

   cFE->require_field_calculation( TT, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( TT, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( TT, PDE_LocalFE::N ) ;

   DENS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   COND->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   CP->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   POW->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   add_one_convected_field( TT, L_EXPLICIT ) ;

   PDE_LinkDOF2Unknown* tt_link =  PDE_LinkDOF2Unknown::create( 0, TT, true ) ;
   NMB = PDE_SystemNumbering::create( this, tt_link ) ;
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   F = A->create_vector( this ) ;
   U = A->create_vector( this ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
AP_DiffusionCG:: ~AP_DiffusionCG( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: transfer_calculation_requirements_for_material_derivative( 
         PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: transfer_calculation_requirements_for_material_derivative" ) ;
   PEL_CHECK_PRE( 
   transfer_calculation_requirements_for_material_derivative_PRE( fe ) ) ;

   fe->require_field_calculation( TT, PDE_LocalFE::N ) ;
   DENS->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
   CP->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: reset_discrete_problem( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: reset_discrete_problem" ) ;
   PEL_CHECK_PRE( reset_discrete_problem_PRE( t_it ) ) ;
   
   start_total_timer( "AP_DiffusionCG:: reset_discrete_problem" ) ;
   // --------------
   
   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;
   
   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   U->re_initialize( n_glob, n_loc ) ;
   
   U_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( U ) ;
   
   stop_total_timer() ;
   // -------------
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: terminate_discrete_problem( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: terminate_discrete_problem" ) ;
   PEL_CHECK_PRE( terminate_discrete_problem_PRE( t_it ) ) ;
   
   start_total_timer( "AP_DiffusionCG:: terminate_discrete_problem" ) ;
   start_assembling_timer() ;
   // ---------------------

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      build_cell_contribution_to_creation( t_it, cFE ) ;
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      build_bound_contribution_to_creation( t_it, bFE ) ;
   }
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------

   estimate_unknowns() ;
   
   stop_solving_timer() ;
   // -----------------

   update_fields() ;
   
   stop_total_timer() ;
   // ---------------
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   
   reset_discrete_problem( t_it ) ;

   start_total_timer( "AP_DiffusionCG:: do_one_inner_iteration" ) ;
   start_assembling_timer() ;
   // ---------------------

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      build_cell_contribution_to_material_derivative( t_it, cFE ) ;
      build_cell_contribution_to_creation( t_it, cFE ) ;
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      build_bound_contribution_to_creation( t_it, bFE ) ;
   }
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------
   
   estimate_unknowns() ;

   stop_solving_timer() ;
   // -----------------

   update_fields() ;

   stop_total_timer() ;
   // ---------------
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: build_cell_contribution_to_creation( 
                                            FE_TimeIterator const* t_it, 
                                            PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: build_cell_contribution_to_creation" ) ;
   
   fe->set_row_and_col_fields( TT, TT ) ;

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                                1,
                                fe->col_field_node_connectivity(), 
                                1 ) ;

   for( fe->start_IP_iterator( QRP ) ; fe->valid_IP() ; fe->go_next_IP() )
   {
      double coef = COND->cell_value_at_IP( t_it, fe  ) ;
      add_grad_row_dot_grad_col( ELEMENT_EQ, fe, coef ) ;
      coef = POW->cell_value_at_IP( t_it, fe  ) ;
      add_row( ELEMENT_EQ, fe, coef ) ;
   }

   PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
}

//----------------------------------------------------------------------
GE_QRprovider const*
AP_DiffusionCG:: QRprovider_for_material_derivative( void ) const
//----------------------------------------------------------------------
{
   GE_QRprovider const* result = QRP ;

   PEL_CHECK_POST( QRprovider_for_material_derivative_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: build_cell_contribution_to_material_derivative( 
                                              FE_TimeIterator const* t_it, 
                                              PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: build_cell_contribution_to_material_derivative" ) ;
   PEL_CHECK_PRE( build_cell_contribution_to_material_derivative_PRE( t_it, 
                                                                      fe ) ) ;

   double dt = t_it->time_step() ;

   fe->set_row_and_col_fields( TT, TT ) ;

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 1,
                           fe->col_field_node_connectivity(), 1 ) ;

   GE_QRprovider const* qrp = QRprovider_for_material_derivative() ;
   for( fe->start_IP_iterator( qrp ) ; fe->valid_IP() ; fe->go_next_IP() )
   {
      double coef = DENS->cell_value_at_IP( t_it, fe ) * 
                    CP->cell_value_at_IP( t_it, fe ) / dt ;
      add_row_col( ELEMENT_EQ, fe, coef ) ;
      add_row( ELEMENT_EQ, fe, coef*fe->value_at_IP( TT, L_EXPLICIT ) ) ;
   }
   
   PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: build_bound_contribution_to_creation(
               	                                FE_TimeIterator const* t_it,
                                                PDE_LocalFEbound* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: build_bound_contribution_to_creation" ) ;
   
   fe->set_row_and_col_fields( TT, TT ) ;

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                           1,
                           fe->col_field_node_connectivity(), 
                           1 ) ;

   if( BCs->has_BC( fe->color(), TT ) )
   {
      PEL_ModuleExplorer const* bc_exp =  BCs->BC_explorer( fe->color(), TT ) ;
      string bc_type = bc_exp->string_data( "type" ) ;
      if( bc_type=="convection" )
      {
         double h = bc_exp->double_data( "convection_coefficient" ) ;
         double Tinf = bc_exp->double_data( "far_field_temperature" ) ;

         for( fe->start_IP_iterator( QRP ) ; fe->valid_IP() ; fe->go_next_IP() )
         {
            add_row_col( ELEMENT_EQ, fe, h ) ;
            add_row( ELEMENT_EQ, fe, h*Tinf ) ;
         }
      }
   }
   
   PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: add_row_col( PDE_LocalEquation* leq,
                              PDE_LocalFE const* fe, 
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: add_row_col" ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;

   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=i ; j<nb_nodes ; ++j )
      {
         double xx = N_row( i ) * N_col( j ) ;
         xx *= c_w ;

         leq->add_to_matrix( xx, i, j ) ;

         if( i!=j ) leq->add_to_matrix( xx, j, i ) ;
      }
   }
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: add_grad_row_dot_grad_col( PDE_LocalEquation* leq,
                                            PDE_LocalFEcell const* fe,
                                            double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: add_grad_row_dot_grad_col" ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;

   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=i ; j<nb_nodes ; ++j )
      {
         double xx = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            xx += dN_row( i, d ) * dN_col( j, d ) ;
         }
         xx *= c_w ;
         leq->add_to_matrix( xx, i, j ) ;
         if( i != j ) leq->add_to_matrix( xx, j, i ) ;
      }
   }
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: add_row( PDE_LocalEquation* leq,
                          PDE_LocalFE const* fe,
                          double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: add_row" ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;

   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; i++ )
   {
      double xx = N_row( i ) ;
      xx *= c_w ;

      leq->add_to_vector( xx, i ) ;
   }
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: estimate_unknowns( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: estimate_unknowns" ) ;
   
   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, U ) ;
   SOLVER->unset_matrix() ;
}

//----------------------------------------------------------------------
void
AP_DiffusionCG:: update_fields( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_DiffusionCG:: update_fields" ) ;  
   
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( U, U_LOC ) ;
   TT->update_free_DOFs_value( L_UPDATE, U_LOC, NMB->link() ) ;
}
