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

#include <AP_NavierStokes1G.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_TwoBlocksMethod.hh>
#include <LA_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>
#include <GE_Point.hh>

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

#include <FE.hh>
#include <FE_LocalBCsBuilder.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ; using std::ostringstream ; 

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_NavierStokes1G const* 
AP_NavierStokes1G:: PROTOTYPE = new AP_NavierStokes1G() ;

//---------------------------------------------------------------------------
AP_NavierStokes1G:: AP_NavierStokes1G( void )
//--------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_NavierStokes1G" )
{
}

//---------------------------------------------------------------------------
AP_NavierStokes1G*
AP_NavierStokes1G:: create_replica( PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_NavierStokes1G* result = new AP_NavierStokes1G( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes1G:: AP_NavierStokes1G( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , NB_DIMS( dom->nb_space_dimensions() )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data( "velocity") ) )
   , L_UU( exp->int_data( "level_of_velocity" ) )
   , UU_EXP( 0 )
   , L_UU_EXP( PEL::bad_index() )
   , PP( dom->set_of_discrete_fields()->item( exp->string_data( "pressure") ) )
   , L_PP( exp->int_data( "level_of_pressure" ) )
   , PP_EXP( 0 )
   , L_PP_EXP( PEL::bad_index() )
   , TDISC( NoTime )
   , L2_STAB( false )
   , ALPHA( 0 )
   , ALPHA_EXP( 0 )
   , MU( prms->item( exp->string_data( "param_viscous" ) ) )
   , RHSU( 0 )
   , BCs( dom->set_of_boundary_conditions() )
   , LOCAL_BC( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , L( 0 )
   , M( 0 )
   , S( 0 )
   , K( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
   , P_LOC( LA_SeqVector::create( this, 0 ) )
   , INIT_DISCRETE_P( exp->has_entry( "initialize_discrete_pressure" ) ?
                      exp->bool_data( "initialize_discrete_pressure" ) :
                      false )
   , INIT_DISCRETE_V( exp->has_entry( "initialize_discrete_velocity" ) ?
                      exp->bool_data( "initialize_discrete_velocity" ) :
                      false )
{
   PEL_LABEL( "AP_NavierStokes1G:: AP_NavierStokes1G" ) ;

   //????? separer en deux : une partie qui lit le jdd et une
   //????? partie qui fait les tests et fixe les requests
   check_field_nb_components( UU, NB_DIMS ) ;
   check_field_storage_depth( UU, L_UU ) ;
   check_field_nb_components( PP, 1 ) ;
   check_field_storage_depth( PP, L_PP ) ;
   
   PDE_SetOfDiscreteFields const* sdf = dom->set_of_discrete_fields() ;
   if( exp->has_module( "time_discretization" ) )
   {
      PEL_ModuleExplorer const* se = 
                         exp->create_subexplorer( 0, "time_discretization") ;
      UU_EXP = sdf->item( se->string_data( "velocity_explicit" ) ) ;
      L_UU_EXP = se->int_data( "level_of_velocity_explicit" ) ;
      check_field_storage_depth( UU_EXP, L_UU_EXP ) ;
      check_field_nb_components( UU_EXP, NB_DIMS ) ;
      
      cFE->require_field_calculation( UU_EXP, PDE_LocalFE::N )  ;
      bFE->require_field_calculation( UU_EXP, PDE_LocalFE::N )  ;
      
      PP_EXP = sdf->item( se->string_data( "pressure_explicit" ) ) ;
      L_PP_EXP = se->int_data( "level_of_pressure_explicit" ) ;
      check_field_storage_depth( PP_EXP, L_PP_EXP ) ;
      check_field_nb_components( PP_EXP, 1 ) ;
      
      cFE->require_field_calculation( PP_EXP, PDE_LocalFE::N ) ; //???? pas tjrs
      cFE->require_field_calculation( PP_EXP, PDE_LocalFE::dN ) ; //???? pas tjrs
      bFE->require_field_calculation( PP_EXP, PDE_LocalFE::N ) ; //????
      
      ALPHA = prms->item( se->string_data( "param_unsteady" ) ) ;
      check_param_nb_components( ALPHA, "param_unsteady", 1 ) ;
      ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      ALPHA->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
      
      string const& tt = se->string_data( "type" ) ;
      if( tt == "Euler" )
      {
         TDISC = Euler ;
         if( se->has_module( "kinetic_energy_control" ) )
         {
            L2_STAB = true ;
            PEL_ModuleExplorer* e =
                      se->create_subexplorer( 0, "kinetic_energy_control" ) ;
            ALPHA_EXP = prms->item( 
                        e->string_data( "param_unsteady_explicit" ) ) ;
            check_param_nb_components( 
                        ALPHA_EXP, "param_unsteady_explicit", 1 ) ;
            ALPHA_EXP->transfer_cell_calculation_requirements( 
                                                 cFE, FE_Parameter::Val ) ;
            ALPHA_EXP->transfer_bound_calculation_requirements( 
                                                  bFE, FE_Parameter::Val ) ;
            e->destroy() ; e = 0 ;
         }
      }
      else if( tt == "BDF2" )
      {
         TDISC = BDF2 ;
         UU_EXP_EXP = sdf->item( 
                           se->string_data( "velocity_explicit_explicit" ) ) ;
         L_UU_EXP_EXP = se->int_data( "level_of_velocity_explicit_explicit" ) ;
         
         check_field_nb_components( UU_EXP_EXP, NB_DIMS ) ;
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

   check_param_nb_components( MU, "param_viscous", 1 ) ;
   
   if( exp->has_entry( "param_source" ) )
   {
      RHSU = prms->item( exp->string_data( "param_source" ) ) ;
      check_param_nb_components( RHSU, "param_source", NB_DIMS ) ;
      RHSU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ); 
   }

   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;   
   cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ; //???? pas tjrs
   bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;

   check_param_nb_components( MU, "param_diffusion", 1 ) ;
   MU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   if( exp->string_data( "viscosity_term" ) == "mu_laplacian_uu" )
      LAPL_UU = true ;
   else if( exp->string_data( "viscosity_term" ) == "div_mu_D_uu" )
      LAPL_UU = false ;
   else
      PEL_Error::object()->raise_bad_data_value( exp, "viscosity_term", 
                                 "\"mu_laplacian_uu\" or \"div_mu_D_uu\"") ;

   if( exp->has_module( "advection" ) )
   {
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "advection" ) ;
      AA = prms->item( e->string_data( "param_advective_velocity" ) ) ;
      check_param_nb_components( AA, "param_advective_velocity", NB_DIMS ) ;
      AA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      e->destroy() ; e = 0 ;
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
   if( SOLVER->L_is_required() ) L = A->create_matrix( this ) ;
   if( SOLVER->MV_is_required() )
   {
      PEL_ASSERT( ALPHA != 0 ) ; //????????????
      M = A->create_matrix( this ) ;
   }
   
   F = A->create_vector( this ) ;
   G = A->create_vector( this ) ;
   if( SOLVER->S_is_required() ) S = A->create_vector( this ) ;
   if( SOLVER->K_is_required() ) K = A->create_vector( this ) ;
   
   U  = A->create_vector( this ) ;
   P  = A->create_vector( this ) ;
   
   SOLVER->set_matrix_prototype( A ) ;

   std::string nn = "L_multilevel_preconditioner" ;
   configure_multilevel_preconditioner( exp, nn, dom, NMB_P ) ;

   nn = "Mv_multilevel_preconditioner";
   configure_multilevel_preconditioner( exp, nn, dom, NMB_U ) ;

   nn = "A_multilevel_preconditioner" ;
   configure_multilevel_preconditioner( exp, nn, dom, NMB_U ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes1G:: ~AP_NavierStokes1G( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: ~AP_NavierStokes1G" ) ;
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   SOLVER->set_indent( indent() ) ;
   
   start_total_timer( "AP_NavierStokes1G:: do_one_inner_iteration" ) ;
   // ---------------------

   reset_discrete_problem() ;

   if( SOLVER->dtinv_is_required() )
   {
      double xx = 1.0/t_it->time_step() ;
      if( TDISC == BDF2 && t_it->iteration_number() != 1 ) xx *= 1.5 ;
      SOLVER->set_dtinv( xx ) ;
   }

   start_assembling_timer() ;
   // ---------------------

   loop_on_cells( t_it ) ;
   loop_on_bounds( t_it ) ;

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------
   
   if( S != 0 )
   {
      S->synchronize() ;
      S->set_as_reciprocal( S ) ;
      SOLVER->set_S( S ) ;
   }
   if( L != 0 )
   {
      L->synchronize() ;
      SOLVER->set_L( L ) ;
   }
   if( K != 0 )
   {
      K->synchronize() ;
      SOLVER->set_K( K ) ;
   }
   if( M != 0 ) 
   {
      M->synchronize() ;
      SOLVER->set_MV( M ) ;
   }
   A->synchronize() ;
   B->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   SOLVER->set_system( A, B, F, G ) ;
   
   // this initialization should be clearly understood, eg for projection 
   // methods: it means that field PP at level L_PP must contain an
   // approximation of the explicit pressure
   if( INIT_DISCRETE_P )
   {
      PP->extract_unknown_DOFs_value( L_PP, P_LOC, NMB_P->link() ) ;
      NMB_P->scatter()->set( P_LOC, P ) ;
   }
   
   if( INIT_DISCRETE_V )
   {
      UU->extract_unknown_DOFs_value( L_UU, U_LOC, NMB_U->link() ) ;
      NMB_U->scatter()->set( U_LOC, U ) ;
   }
   
   SOLVER->estimate_unknowns( INIT_DISCRETE_V, U, INIT_DISCRETE_P, P ) ;
   
   stop_solving_timer() ;

   if( !SOLVER->successful_estimation() )
   {
      PEL_Error::object()->display_info(
         "*** AP_NavierStokes1G error:\n"
         "    solution failure" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "   update of " << UU->name() 
                    << "(" << L_UU << ") and " << PP->name() 
                    << "(" << L_PP << ")" << endl ;
      }
      NMB_U->scatter()->get( U, U_LOC ) ;
      UU->update_free_DOFs_value( L_UU, U_LOC, NMB_U->link() ) ;

      NMB_P->scatter()->get( P, P_LOC ) ;
      PP->update_free_DOFs_value( L_PP, P_LOC, NMB_P->link() ) ;
   }
   
   SOLVER->unset_system() ;

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: print_additional_times( std::ostream& os,
                                            size_t indent_width) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: print_additional_times" ) ;

   SOLVER->print_times( indent_width ) ;
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: print" ) ;

   FE_OneStepIteration:: print( os, indent_width ) ;

   std::string const s( indent_width+3, ' ' ) ;
   os << s << "kinetic energy control: " ;
   if( L2_STAB )
   {
      os << "yes" << std::endl ;
   }
   else
   {
      os << "no" << std::endl ;
   }

   os << s << "initialize discrete pressure: " ;
   if( INIT_DISCRETE_P )
   {
      os << "yes" << std::endl ;
   }
   else
   {
      os << "no" << std::endl ;
   }

   os << s << "initialize discrete velocity: " ;
   if( INIT_DISCRETE_V )
   {
      os << "yes" << std::endl ;
   }
   else
   {
      os << "no" << std::endl ;
   }
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: loop_on_cells( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{ 
   PEL_LABEL( "AP_NavierStokes1G:: loop_on_cells" ) ;
   
   double m_xx, mu ;
   doubleVector aa( NB_DIMS ) ;
   doubleVector rhs( NB_DIMS ) ;
   
   //????? there are two many loops on integration points

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), NB_DIMS,
                              cFE->col_field_node_connectivity(), NB_DIMS ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         compute_coefs_at_IP( t_it, m_xx, aa, mu, rhs ) ;
                 
         if( AA != 0 )
         {
            if( L2_STAB )
            {
               FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa,  0.5 ) ;               
               FE::add_vvgrad_row_col( ELEMENT_EQ, cFE, aa, -0.5 ) ;
            }
            else
            {
               FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;               
            }
         }

         FE::add_row( ELEMENT_EQ, cFE, rhs ) ;
         FE::add_row_col_S( ELEMENT_EQ, cFE, m_xx ) ;
 
         if( LAPL_UU )
            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, mu ) ;
         else
            FE::add_grad_row_D_col_S( ELEMENT_EQ, cFE, mu ) ;

      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;

      cFE->set_row_and_col_fields( PP, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), NB_DIMS ) ;

      cFE->start_IP_iterator( QRP ) ;      
      for(  ; cFE->valid_IP() ; cFE->go_next_IP() ) 
      {
         FE::add_row_div_col( ELEMENT_EQ, cFE, -1.0 ) ;
      }
      PDE::assemble_in_matrix_vector_0( B, G, ELEMENT_EQ, NMB_P, NMB_U ) ;

      if( S != 0 )
      {
         cFE->set_row_and_col_fields( PP, PP ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                 cFE->col_field_node_connectivity(), 1 ) ;

         size_t nb_nodes = cFE->nb_basis_functions( row ) ;
         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            double c_w = cFE->weight_of_IP() ;
            if( FE::geometry() == FE::axisymmetrical )
            {
               c_w *= cFE->coordinates_of_IP()->coordinate(0) ;
            }
            for( size_t i=0 ; i<nb_nodes ; ++i )
            {  
               //????? on suppose ici que sum phi_j = 1
               //????? pas vrai pour les elements de Fanny
               double xx = c_w* cFE->N_at_IP( row, i ) ;
               ELEMENT_EQ->add_to_matrix( xx, i, i ) ;     
            }
         }
         
         //??????
         for( size_t il=0 ; il<ELEMENT_EQ->nb_rows() ; il++ )
         {
            size_t nn = ELEMENT_EQ->row_node( il ) ;
            if( NMB_P->link()->DOF_is_unknown( nn ) )
            {
               size_t ig = NMB_P->global_unknown_for_DOF( nn, 0, 0 ) ;
               double aii = ELEMENT_EQ->matrix_item( il, il ) ;
               S->add_to_item( ig, aii ) ;
            }
         }
      }

//???????? lenteur avec assemble_A ????
      if( M != 0  )
      {
         cFE->set_row_and_col_fields( UU, UU ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), NB_DIMS,
                                 cFE->col_field_node_connectivity(), NB_DIMS ) ;

         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
            FE::add_row_col_S( ELEMENT_EQ, cFE, alpha ) ;
         }
         PDE::assemble_in_matrix_0( M, ELEMENT_EQ, NMB_U ) ;
      }


//???????? lenteur avec assemble_MPl ????
      if( L != 0  )
      {
         cFE->set_row_and_col_fields( PP, PP ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                 cFE->col_field_node_connectivity(), 1 ) ;
         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            double xx = ALPHA->cell_value_at_IP( t_it, cFE ) ;
            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0/xx ) ;
            if( K != 0 )
            {
               if( L2_STAB )
               {
                  xx = PEL::sqrt( xx * 
                                  ALPHA_EXP->cell_value_at_IP( t_it, cFE ) ) ;
               }
               for( size_t d=0 ; d<NB_DIMS ; ++d )
               {
                  //????? on reutilise le tableau rhs
                  rhs( d ) = cFE->gradient_at_IP( PP_EXP, L_PP_EXP, d ) / xx ;
               }
               FE:: add_grad_row( ELEMENT_EQ, cFE, rhs ) ;
            }
         }
         if( K == 0 )
         {
            PDE::assemble_in_matrix_0( L, ELEMENT_EQ, NMB_P ) ;
         }
         else
         {        
            PDE::assemble_in_matrix_vector_0( L, K, ELEMENT_EQ, NMB_P ) ;
         }
      }     
   }
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: loop_on_bounds( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{ 
   PEL_LABEL( "AP_NavierStokes1G:: loop_on_bounds" ) ;

   size_t nbc = UU->nb_components() ;

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
               ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 
                                       nbc,
                                       bFE->col_field_node_connectivity(), 
                                       nbc ) ;

               LOCAL_BC->build_current_BC( ELEMENT_EQ, bFE, t_it, QRP ) ;

               PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;
            }
         }
      }
   }

   if( L2_STAB )
   {
      double aa_dot_n ;

      for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
      {
         GE_Color const* color = bFE->color() ;
         if( BCs->has_BC( color, UU ) )
         {
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;
            std::string const& bc_type = ee->string_data( "type" ) ;
            if( bc_type == "NeumannVelocity" )
            {
               bFE->set_row_and_col_fields( UU, UU ) ;
               ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(),
                                       nbc,
                                       bFE->col_field_node_connectivity(),
                                       nbc ) ;

               bFE->start_IP_iterator( QRP ) ;
               for( ; bFE->valid_IP() ; bFE->go_next_IP() )
               {
                  aa_dot_n = 0.0 ;
                  double xx = ALPHA->bound_value_at_IP( t_it, bFE ) ;
                  for( size_t d=0 ; d<NB_DIMS ; ++d )
                  {
                     aa_dot_n += xx * bFE->value_at_IP( UU_EXP, L_UU_EXP, d )
                                    * bFE->outward_normal()->component( d ) ;
                  }
                  FE::add_row_col_S( ELEMENT_EQ, bFE, 0.5*aa_dot_n ) ;
               }
               PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;
            }
         }
      }
   }

//???????????????????????????????
   if( L != 0 )
   {
      for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
      {
         if( BCs->has_BC( bFE->color(), PP ) )
         {
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( bFE->color(), PP ) ;
            std::string const& bc_type = ee->string_data( "type" ) ;
            if( bc_type == "Dirichlet_pressure_penalization" )
            {
               double penal = ee->double_data( "penalization_coefficient" ) ;
               bFE->set_row_and_col_fields( PP, PP ) ;
               ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 1,
                                       bFE->col_field_node_connectivity(), 1 ) ;

               bFE->start_IP_iterator( QRP ) ; 
               for( ; bFE->valid_IP() ; bFE->go_next_IP() )
               {
                  double xx = ALPHA->bound_value_at_IP( t_it, bFE ) ;
                  FE::add_row_col_S( ELEMENT_EQ, bFE, penal/xx ) ;
                  if( K != 0 )
                  {
                     if( L2_STAB )
                     {
                        xx = PEL::sqrt( xx * 
                             ALPHA_EXP->bound_value_at_IP( t_it, bFE ) ) ;
                     }
                     FE::add_row( ELEMENT_EQ, bFE, 
                         penal * bFE->value_at_IP( PP_EXP, L_PP_EXP ) / xx ) ;
                  }
               }
               if( K == 0 )
               {
                  PDE::assemble_in_matrix_0( L, ELEMENT_EQ, NMB_P ) ;
               }
               else
               {
                  PDE::assemble_in_matrix_vector_0( L, K, ELEMENT_EQ, NMB_P ) ;                  
               }
            }  
         }
      }
   }
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: reset_discrete_problem( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: reset_discrete_problem" ) ;
   
   NMB_U->reset() ;
   NMB_P->reset() ;
   
   size_t nv_glob = NMB_U->nb_global_unknowns() ;
   size_t np_glob = NMB_P->nb_global_unknowns() ;

   size_t nv_loc = NMB_U->nb_unknowns_on_current_process() ;
   size_t np_loc = NMB_P->nb_unknowns_on_current_process() ;

   A->re_initialize( nv_glob, nv_glob, nv_loc, np_loc ) ;
   B->re_initialize( np_glob, nv_glob, np_loc, nv_loc ) ;
   if( L != 0 ) L->re_initialize( np_glob, np_glob, np_loc, np_loc ) ;
   if( M != 0 ) M->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;

   U->re_initialize( nv_glob, nv_loc ) ;
   P->re_initialize( np_glob, np_loc ) ;

   F->re_initialize( nv_glob, nv_loc ) ;
   G->re_initialize( np_glob, np_loc ) ;
   
   if( S != 0 ) S->re_initialize( np_glob, np_loc ) ;
   if( K != 0 ) K->re_initialize( np_glob, np_loc ) ;
   
   SOLVER->re_initialize_internals( nv_glob, np_glob, nv_loc, np_loc ) ;

   U_LOC->re_initialize( NMB_U->link()->unknown_vector_size() ) ;
   P_LOC->re_initialize( NMB_P->link()->unknown_vector_size() ) ;
   
   NMB_U->define_scatters( U ) ;
   NMB_P->define_scatters( P ) ;   
}

//---------------------------------------------------------------------------
void
AP_NavierStokes1G:: compute_coefs_at_IP( FE_TimeIterator const* t_it,
                                         double& m_xx,
                                         doubleVector& aa,
                                         double& mu,
                                         doubleVector& rhs ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1G:: compute_coefs_at_IP" ) ;
   
   m_xx = 0.0 ;
   mu = MU->cell_value_at_IP( t_it, cFE ) ;
   rhs.set( 0.0 ) ;
   
   if( RHSU != 0 )
   {
      for( size_t d=0 ; d<NB_DIMS ; ++d )
      {
         rhs( d ) += RHSU->cell_value_at_IP( t_it, cFE, d ) ;
      }
   }
   
   if( ALPHA != 0 )
   {
      double dt = t_it->time_step() ;   
      doubleVector uu_exp( NB_DIMS, PEL::bad_double() ) ;
      if( TDISC == Euler || ( TDISC == BDF2 && t_it->iteration_number() == 1 ) )
      {
         for( size_t d=0 ; d<NB_DIMS ; ++d )
         {
            uu_exp( d ) = cFE->value_at_IP( UU_EXP, L_UU_EXP, d ) ;
         }
      }
      else if( TDISC == BDF2 && t_it->iteration_number() != 1 )
      {
         for( size_t d=0 ; d<NB_DIMS ; ++d )
         {
            uu_exp( d ) = 
               2.0 * cFE->value_at_IP( UU_EXP, L_UU_EXP, d )
             - 0.5 * cFE->value_at_IP( UU_EXP_EXP, L_UU_EXP_EXP, d ) ;
         }
      }
      
      double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
      m_xx = alpha/dt ;
      if( TDISC == BDF2 && t_it->iteration_number() != 1 ) m_xx *= 1.5 ;

      double rt_exp = ( L2_STAB ? 
            PEL::sqrt( alpha * ALPHA_EXP->cell_value_at_IP( t_it, cFE ) ) / dt : 
            alpha / dt ) ;
      for( size_t d=0 ; d<NB_DIMS ; ++d )
      {
         PEL_ASSERT( uu_exp( d ) != PEL::bad_double() ) ;
         rhs( d ) += rt_exp * uu_exp( d ) ;
      }
   }
   
   if( AA != 0 )
   {
      for( size_t d=0 ; d<aa.size() ; ++d )
      {
         aa( d ) = AA->cell_value_at_IP( t_it, cFE, d ) ;
      }
   }
}
