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

#include <AP_StokesCG.hh>

#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_TwoBlocksMethod.hh>
#include <LA_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_StokesCG const* AP_StokesCG::PROTOTYPE = new AP_StokesCG() ;

//----------------------------------------------------------------------
AP_StokesCG:: AP_StokesCG( void )
//----------------------------------------------------------------------
   : FE_Galerkin( "AP_StokesCG" )
{
}

//----------------------------------------------------------------------
AP_StokesCG*
AP_StokesCG:: create_replica( PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_StokesCG* result = new AP_StokesCG( a_owner, dom, prms, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_StokesCG:: AP_StokesCG( PEL_Object* a_owner, 
                         PDE_DomainAndFields const* dom, 
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_Galerkin( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( "velocity" ) )
   , PP( dom->set_of_discrete_fields()->item( "pressure" ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , DENS( prms->item( exp->string_data( "param_density" ) ) )
   , VISC( prms->item( exp->string_data( "param_viscosity" ) ) )
   , RHSU( prms->item( exp->string_data( "param_force" ) ) )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( 0 )
   , S( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
   , P_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "AP_StokesCG:: AP_StokesCG" ) ;

   check_field_storage_depth( PP, L_UPDATE ) ;
   check_field_storage_depth( UU, PEL::max( L_UPDATE, L_EXPLICIT ) ) ;

   check_param_nb_components( DENS, "param_density", 1 ) ;
   check_param_nb_components( VISC, "param_viscosity", 1 ) ;
   check_param_nb_components( RHSU, "param_force", dom->nb_space_dimensions());

   if( PEL_Exec::communicator()->nb_ranks()>1 )
   {
      PEL_Error::object()->raise_plain(
         "\"AP_StokesCG\" does not support parallel mode") ;
   }
   
   cFE = dom->create_LocalFEcell( this ) ;

   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;

   add_one_convected_field( UU, L_EXPLICIT ) ;

   DENS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   VISC->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   RHSU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU, 
                                          "sequence_of_the_components",
                                          true ) ;
   PDE_LinkDOF2Unknown* pp_link = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
   
   NMB_U = PDE_SystemNumbering::create( this, uu_link ) ;
   NMB_P = PDE_SystemNumbering::create( this, pp_link ) ;

   PEL_ModuleExplorer const* ee = 
                      exp->create_subexplorer( 0, "LA_TwoBlocksMethod" ) ;
   SOLVER = LA_TwoBlocksMethod::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   B = A->create_matrix( this ) ;
   C = A->create_matrix( this ) ;
   
   U  = A->create_vector( this ) ;
   P  = A->create_vector( this ) ;
   if( SOLVER->S_is_required() ) S = A->create_vector( this ) ;
   
   F = A->create_vector( this ) ;
   G = A->create_vector( this ) ;
   
   SOLVER->set_matrix_prototype( A ) ;
}

//----------------------------------------------------------------------
AP_StokesCG:: ~AP_StokesCG( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
AP_StokesCG:: transfer_calculation_requirements_for_material_derivative( 
         PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: transfer_calculation_requirements_for_material_derivative" ) ;
   PEL_CHECK_PRE( 
   transfer_calculation_requirements_for_material_derivative_PRE( fe ) ) ;

   fe->require_field_calculation( UU, PDE_LocalFE::N ) ;
   DENS->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//----------------------------------------------------------------------
void
AP_StokesCG:: reset_discrete_problem( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: reset_discrete_problem" ) ;
   PEL_CHECK_PRE( reset_discrete_problem_PRE( t_it ) ) ;

   start_total_timer( "AP_StokesCG:: reset_discrete_problem" ) ;
   // --------------

   NMB_U->reset() ;
   NMB_P->reset() ;
   
   size_t nv_glob = NMB_U->nb_global_unknowns() ;
   size_t np_glob = NMB_P->nb_global_unknowns() ;

   size_t nv_loc = NMB_U->nb_unknowns_on_current_process() ;
   size_t np_loc = NMB_P->nb_unknowns_on_current_process() ;

   re_initialize_matrices_and_vectors( nv_glob, np_glob, nv_loc, np_loc ) ;
   
   SOLVER->re_initialize_internals( nv_glob, np_glob, nv_loc, np_loc ) ;

   U_LOC->re_initialize( NMB_U->link()->unknown_vector_size() ) ;
   P_LOC->re_initialize( NMB_P->link()->unknown_vector_size() ) ;
   
   NMB_U->define_scatters( U ) ;
   NMB_P->define_scatters( P ) ;
   
   stop_total_timer() ;
   // -------------
}

//----------------------------------------------------------------------
void
AP_StokesCG:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_StokesCG:: do_one_inner_iteration" ) ;
   // --------------

   reset_discrete_problem( t_it ) ;

   start_assembling_timer() ;
   // -------------------

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      build_cell_contribution_to_material_derivative( t_it, cFE ) ;
      build_cell_contribution_to_creation( t_it, cFE ) ;
   }

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------

   estimate_unknowns() ;
   
   stop_solving_timer() ;
   // -----------------

   if( !SOLUTION_ACHIEVED )
   {
      PEL_Error::object()->raise_plain( "solution failure" ) ;
   }   

   update_fields() ;
   
   stop_total_timer() ;
   // -------------
}

//----------------------------------------------------------------------
void
AP_StokesCG:: terminate_discrete_problem( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: terminate_discrete_problem" ) ;
   PEL_CHECK_PRE( terminate_discrete_problem_PRE( t_it ) ) ;

   start_total_timer( "AP_StokesCG:: terminate_discrete_problem" ) ;
   start_assembling_timer() ;
   // ---------------------

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      build_cell_contribution_to_creation( t_it, cFE ) ;
   }

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------

   estimate_unknowns() ;
   
   stop_solving_timer() ;
   // -----------------

   if( !SOLUTION_ACHIEVED )
   {
      PEL_Error::object()->raise_plain( "solution failure" ) ;
   }   

   update_fields() ;
   
   stop_total_timer() ;
   // ---------------
}

//----------------------------------------------------------------------
GE_QRprovider const*
AP_StokesCG:: QRprovider_for_material_derivative( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: QRprovider_for_material_derivative" ) ;

   GE_QRprovider const* result = QRP ;

   PEL_CHECK_POST( QRprovider_for_material_derivative_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_StokesCG:: build_cell_contribution_to_material_derivative( 
                                          FE_TimeIterator const* t_it, 
                                          PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: build_cell_contribution_to_material_derivative" ) ;
   PEL_CHECK_PRE( build_cell_contribution_to_material_derivative_PRE( t_it,
                                                                      fe ) ) ;

   static doubleVector aauu( UU->nb_components() ) ;

   fe->set_row_and_col_fields( UU, UU ) ;

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                           UU->nb_components(),
                           fe->col_field_node_connectivity(), 
                           UU->nb_components() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   GE_QRprovider const* qrp = QRprovider_for_material_derivative() ;
   for( fe->start_IP_iterator( qrp ) ; fe->valid_IP() ; fe->go_next_IP() ) 
   {
      double aa = DENS->cell_value_at_IP( t_it, fe )/t_it->time_step() ;
      add_row_dot_col( ELEMENT_EQ, fe, aa ) ;
      for( size_t d=0 ; d<nb_dims ; ++d ) 
      {
         aauu(d) = aa * fe->value_at_IP( UU, L_EXPLICIT, d ) ;
      }
      add_row( ELEMENT_EQ, fe, aauu ) ;
   }
   PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;
}

//----------------------------------------------------------------------
void
AP_StokesCG:: build_cell_contribution_to_creation(
                                          FE_TimeIterator const* t_it, 
                                          PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: build_cell_contribution_to_creation" ) ;

   fe->set_row_and_col_fields( UU, UU ) ;
   //  --------------------------------

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                           UU->nb_components(),
                           fe->col_field_node_connectivity(), 
                           UU->nb_components() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   doubleVector coef( nb_dims ) ;
   for( fe->start_IP_iterator( QRP ) ; fe->valid_IP() ; fe->go_next_IP() ) 
   {
      add_D_row_dot_grad_col( ELEMENT_EQ, fe, 
                              VISC->cell_value_at_IP( t_it, fe ) ) ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         coef( d ) = RHSU->cell_value_at_IP( t_it, fe, d ) ;
      }
      add_row( ELEMENT_EQ, fe, coef ) ;
   }
   PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;
 
   fe->set_row_and_col_fields( PP, UU ) ;
   //  --------------------------------

   ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                           1,
                           fe->col_field_node_connectivity(), 
                           UU->nb_components() ) ;

   for( fe->start_IP_iterator( QRP ) ; fe->valid_IP() ; fe->go_next_IP() ) 
   {
      add_row_div_col( ELEMENT_EQ, fe, -1.0 ) ;
   }
   PDE::assemble_in_matrix_vector_0( B, G, ELEMENT_EQ, NMB_P, NMB_U ) ;

   if( S != 0 )
   {
      fe->set_row_and_col_fields( PP, PP ) ;
      //  --------------------------------

      ELEMENT_EQ->initialize( fe->row_field_node_connectivity(), 
                              1,
                              fe->col_field_node_connectivity(), 
                              1 ) ;

      for( fe->start_IP_iterator( QRP ) ; fe->valid_IP() ; fe->go_next_IP() ) 
      {
         size_t nb_nodes = fe->nb_basis_functions( row ) ;
         double xx = fe->weight_of_IP() ;
         for( size_t i=0 ; i<nb_nodes ; ++i )
         {
            ELEMENT_EQ->add_to_vector( xx, i, 0 ) ;
         }
      }
      PDE::assemble_in_vector_1( S, ELEMENT_EQ, NMB_P ) ;
   }
}

//----------------------------------------------------------------------
void
AP_StokesCG:: add_D_row_dot_grad_col( PDE_LocalEquation* leq,
                                     PDE_LocalFEcell const* fe, 
                                     double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: add_D_row_dot_grad_col" ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = coef*fe->weight_of_IP() ;

   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; i++ )
   {
      for( size_t j=0 ; j<nb_nodes ; j++ )
      {
         for( size_t d2=0 ; d2<nb_dims ; d2++ )
         {
            double x = c_w * dN_row( i, d2 ) ;
            double y = x * dN_col( j, d2 ) ;
            for( size_t d1=0 ; d1<nb_dims ; d1++ )
            {
               leq->add_to_matrix( x*dN_col( j, d1 ), i, j, d1, d2 ) ;
               leq->add_to_matrix( y, i, j, d1, d1) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
AP_StokesCG:: add_row_div_col( PDE_LocalEquation* leq,
                              PDE_LocalFEcell const* fe,
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: add_row_div_col" ) ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->polyhedron()->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;

   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            double xx = N_row( i ) * dN_col( j, d ) ;
            xx *= c_w ;
            leq->add_to_matrix( xx, i, j, 0, d ) ; 
         }
      }
   }
}

//----------------------------------------------------------------------
void
AP_StokesCG:: add_row_dot_col( PDE_LocalEquation* leq,
                              PDE_LocalFEcell const* fe, 
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: add_row_dot_col" ) ;

   size_t nbc = fe->field( row )->nb_components() ;
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

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
            if( i != j ) leq->add_to_matrix( xx, j, i, ic, ic  ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
AP_StokesCG:: add_row( PDE_LocalEquation* leq,
                      PDE_LocalFEcell const* fe,
                      doubleVector const& coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: add_row" ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;

   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double yy = fe->weight_of_IP() * N_row( i ) ;
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = yy * coef( ic ) ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
AP_StokesCG:: estimate_unknowns( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: estimate_unknowns" ) ;
   
   //???  devrait etre L_EXPLICIT
   //???  PP->extract_unknown_DOFs_value( L_UPDATE, P_LOC, NMB_P->link() ) ;
   //???  NMB_P->scatter()->set( P_LOC, P ) ;

   //???  LA MATRICE C EST TOUJOURS NULLE ???
   
   A->synchronize() ;
   B->synchronize() ;
   C->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   
   if( S != 0 )
   {
      S->synchronize() ;
      S->set_as_reciprocal( S, 1.0e-100, 1.0 ) ;
      SOLVER->set_S( S ) ;
   }
   
   SOLVER->set_system( A, B, F, G, C ) ;
   
   SOLVER->estimate_unknowns( false, U, false, P ) ;
   
   //??? jamais testé
   SOLUTION_ACHIEVED = SOLVER->successful_estimation() ;
   
   SOLVER->unset_system() ;
}

//-----------------------------------------------------------------------
void
AP_StokesCG:: update_fields( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: update_fields" ) ;
   
   LA_Scatter const* sca = NMB_U->scatter() ;
   sca->get( U, U_LOC ) ;
   UU->update_free_DOFs_value( L_UPDATE, U_LOC, NMB_U->link() ) ;

   sca = NMB_P->scatter() ;
   sca->get( P, P_LOC ) ;
   PP->update_free_DOFs_value( L_UPDATE, P_LOC, NMB_P->link() ) ; 
}

//-----------------------------------------------------------------------
void
AP_StokesCG:: re_initialize_matrices_and_vectors( size_t nv_glob,
                                                 size_t np_glob,
                                                 size_t nv_loc,
                                                 size_t np_loc )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StokesCG:: re_initialize_matrices_and_vectors" ) ;

   A->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   B->re_initialize( np_glob, nv_glob, np_loc, nv_loc ) ;
   C->re_initialize( np_glob, np_glob, np_loc, np_loc ) ;

   U->re_initialize( nv_glob, nv_loc ) ;
   P->re_initialize( np_glob, np_loc ) ;
   S->re_initialize( np_glob, np_loc ) ; //????? si S non nul ???

   F->re_initialize( nv_glob, nv_loc ) ;
   G->re_initialize( np_glob, np_loc ) ;
}

