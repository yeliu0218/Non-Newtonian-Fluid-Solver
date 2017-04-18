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

#include <AP_NavierStokes2CFV.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <LA_SeqVector.hh>
#include <LA_Solver.hh>
#include <LA_SeqMatrix.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE.hh>
#include <PDE_CursorFEside.hh>
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
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <AP_NavierStokes1System.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;
using std::ofstream ;

AP_NavierStokes2CFV const* 
AP_NavierStokes2CFV::PROTOTYPE = new AP_NavierStokes2CFV() ;

struct AP_NavierStokes2CFV_ERROR
{
   static void n0( PDE_LocalFEbound const* fe, double int_v_nor ) ;
   static void n1( PEL_ModuleExplorer const* exp,
                   std::string const& keyword, 
                   size_t size) ;
   static void n2( void ) ;
} ;

//---------------------------------------------------------------------------
AP_NavierStokes2CFV:: AP_NavierStokes2CFV( void )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( "AP_NavierStokes2CFV" )
{
}

//---------------------------------------------------------------------------
AP_NavierStokes2CFV*
AP_NavierStokes2CFV:: create_replica( PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_NavierStokes2CFV* result = 
                      new AP_NavierStokes2CFV( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes2CFV:: AP_NavierStokes2CFV( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "velocity" ) ) )
   , L_UPDATE_UU( exp->int_data( "velocity_level_to_update" ) )
   , L_EXPLICIT_UU( PEL::bad_index() )
   , PP( dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "pressure" ) ) )
   , L_UPDATE_PP( exp->int_data( "pressure_level_to_update" ) )
   , L_EXPLICIT_PP( PEL::bad_index() )
   , IS_UNSTEADY( false )
   , C_DT( 0.0 )
   , C_VGRAD( exp->double_data( "coef_vgrad_v" ) )
   , MU( exp->double_data( "coef_viscous" ) )
   , PI( prms->item( exp->string_data( "param_source" ) ) )
   , NB_ITER_MAX( exp->int_data( "nb_iterations_max" ) )
   , RELAX( exp->has_entry( "relaxation_coefficient" ) ? 
            exp->double_data( "relaxation_coefficient" ) : 1.0 )
   , TOL( exp->double_data( "newton_tolerance" ) )
   , LAMBDA( exp->double_data( "lambda_in_pressure_stabilization" ) ) 
   , h_EXPONENT( exp->double_data( "h_exponent_in_infsup_stabilization" ) )
   , CLUSTER_BD( 0 )
   , BCs( dom->set_of_boundary_conditions() )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP_PI( GE_QRprovider::object( 
             exp->string_data( "quadrature_rule_provider_for_source" ) ) )
   , idx_UU( 0 )
   , idx_PP( 1 )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , CONTEXT( PEL_ContextSimple::create( this ) )
   , XX( 0 )
   , TT( 0 )
{
   PEL_LABEL( "AP_NavierStokes2CFV:: AP_NavierStokes2CFV" ) ;

   check_field_storage_depth( UU, L_UPDATE_UU ) ;
   check_field_storage_depth( PP, L_UPDATE_PP ) ;

   if( exp->has_module( "unsteady" ) )
   {
      IS_UNSTEADY = true ;
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "unsteady" ) ;
      L_EXPLICIT_UU = se->int_data( "level_of_explicit_velocity" ) ;
      C_DT = se->double_data( "coef_dt_v" ) ;
      se->destroy() ;

      check_field_storage_depth( UU, L_EXPLICIT_UU ) ;
      if( L_UPDATE_UU == L_EXPLICIT_UU )
         AP_NavierStokes2CFV_ERROR::n2() ;
   }

   check_param_nb_components( PI, "source", 
                              dom->nb_space_dimensions() ) ;

   sFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::node ) ;

   sFE->require_field_calculation( PP, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( PP, PDE_LocalFE::node ) ;

   PI->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
   PI->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
   PI->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   if( exp->has_entry( "color_of_cluster_boundaries" ) )
   {
      CLUSTER_BD = GE_Color::object(
                   exp->string_data( "color_of_cluster_boundaries" ) ) ;
   }

   //*****************
   // if the boundary conditions are "Dirichlet" on the whole boundary,
   // the pressure at node 0 is imposed
   //*****************
   PP->modify_DOF( true, PP->DOF_value( 0, 0, 0 ), 0 ) ;

   PEL_Vector* vec = PEL_Vector::create( 0, 2 ) ;
   vec->set_at( idx_UU, PDE_LinkDOF2Unknown::create( 0, UU,
                                             "sequence_of_the_components",
                                             true ) ) ;
   vec->set_at( idx_PP, PDE_LinkDOF2Unknown::create( 0, PP, true ) ) ;
   std::string ordering = "sequence_of_the_discrete_fields" ;
   NMB = PDE_SystemNumbering::create( this, vec, ordering ) ;  
   vec->destroy() ; vec = 0 ;

   //*****************
   // on pourrait eviter d'assembler toute la jacobienne en conservant
   // sur un niveau particulier les termes qui ne dependent pas des
   // iterations de Newton
   //*****************

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;
   
   size_t nb_global_unk = NMB->nb_global_unknowns() ;
   size_t nb_local_unk  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( nb_global_unk, nb_global_unk,
                     nb_local_unk, nb_local_unk ) ;
   F->re_initialize( nb_global_unk, nb_local_unk ) ;
   X->re_initialize( nb_global_unk, nb_local_unk ) ;

   NMB->define_scatters( X ) ;

   size_t nn = NMB->link( 0 )->unknown_vector_size() ;
   size_t mm = NMB->link( 1 )->unknown_vector_size() ;
   if( mm > nn ) nn = mm ;
   X_LOC->re_initialize( nn ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   XX = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object("DV_X"), XX ) ;

   TT = PEL_Double::create( CONTEXT, 0.0 ) ;
   CONTEXT->extend( PEL_Variable::object("DS_T"), TT ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes2CFV:: ~AP_NavierStokes2CFV( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
size_t
AP_NavierStokes2CFV:: nb_unknowns( void ) const
//---------------------------------------------------------------------------
{
   return( 2 ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField*
AP_NavierStokes2CFV:: field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: field" ) ;
   PEL_CHECK_PRE( field_PRE( i_unk ) ) ;

   PDE_DiscreteField* result = 0 ;
   if( i_unk == idx_UU ) result = UU ;
   if( i_unk == idx_PP ) result = PP ;

   PEL_CHECK_POST( field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
AP_NavierStokes2CFV:: level_of_field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: level_of_field" ) ;
   PEL_CHECK_PRE( level_of_field_PRE( i_unk ) ) ;

   size_t result = PEL::bad_index() ;
   if( i_unk == idx_UU ) result = L_UPDATE_UU ;
   if( i_unk == idx_PP ) result = L_UPDATE_PP ;

   PEL_CHECK_POST( level_of_field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
AP_NavierStokes2CFV:: link_DOF_2_unknown( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( link_DOF_2_unknown_PRE( i_unk ) ) ;

   PDE_LinkDOF2Unknown const* result = NMB->link( i_unk ) ;

   PEL_CHECK_POST( link_DOF_2_unknown_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_NavierStokes2CFV:: build_function_and_jacobian( 
                                           FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: build_function_and_jacobian" ) ;
   
   A->nullify() ;
   F->nullify() ;

   start_assembling_timer() ;
   // ---------------------

   loop_on_cells( t_it ) ;
   loop_on_sides( t_it ) ;
   loop_on_bounds( t_it ) ;

   A->synchronize() ;
   F->synchronize() ;
   
   stop_assembling_timer() ;
   // --------------------
}

//---------------------------------------------------------------------------
LA_SeqVector const*
AP_NavierStokes2CFV:: create_function( PEL_Object* a_owner,
                                            size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: create_function" ) ;
   PEL_CHECK_PRE( create_function_PRE( a_owner, i_unk ) ) ;

   size_t nn = NMB->link( i_unk )->unknown_vector_size() ;
   LA_SeqVector* result = LA_SeqVector::create( a_owner, nn ) ;
   NMB->scatter( i_unk )->get( F, result ) ;
   
   PEL_CHECK_POST( create_function_POST( result, a_owner, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SeqMatrix const*
AP_NavierStokes2CFV:: create_jacobian( PEL_Object* a_owner,
                                            size_t i_eq, 
                                            size_t j_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: create_jacobian" ) ;
   PEL_CHECK_PRE( create_jacobian_PRE( a_owner, i_eq, j_unk ) ) ;
   
   LA_SeqMatrix const* result =
                PDE::create_extracted_block( a_owner, A, NMB, i_eq, j_unk ) ;
   
   PEL_CHECK_POST( create_jacobian_POST( result, a_owner, i_eq, j_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
AP_NavierStokes2CFV:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_NavierStokes2CFV:: do_one_inner_iteration" ) ;
   // --------------

   size_t iter = 0 ;
   do
   {
      iter++ ;
      if( iter >= NB_ITER_MAX )
         PEL_Error::object()->raise_plain( "Newton convergence failure" ) ;

      build_function_and_jacobian( t_it ) ;

      start_solving_timer() ;
      // ------------------

      //??? Assuming that the matrix strusture doesn't change ?
      if( SOLVER->matrix_is_set() )
         SOLVER->reset_matrix() ;
      else
         SOLVER->set_matrix( A ) ;
      
      SOLVER->solve( F, X ) ;
      
      if( !SOLVER->solution_is_achieved() )
      {
         PEL_Error::object()->raise_plain(
            "AP_NavierStokes2CFV:: do_one_inner_iteration\n"
            "    Linear system resolution failure\n" ) ;
      }

      stop_solving_timer() ;
      // -----------------

      LA_Scatter const* sca = NMB->scatter( idx_UU ) ;
      PDE_LinkDOF2Unknown const* link = NMB->link( idx_UU ) ;
      PEL_ASSERT( UU == link->field() ) ;
      sca->get( X, X_LOC ) ;
      UU->add_to_free_DOFs_value( L_UPDATE_UU, X_LOC, link, RELAX ) ;
      if( verbose_level() >= 1 ) PRINTED_RES_UU = X_LOC->max_norm() ;

      sca = NMB->scatter( idx_PP ) ;
      link = NMB->link( idx_PP ) ;
      PEL_ASSERT( PP == link->field() ) ;
      sca->get( X, X_LOC ) ;
      PP->add_to_free_DOFs_value( L_UPDATE_PP, X_LOC, link, RELAX ) ;
      if( verbose_level() >= 1 ) PRINTED_RES_PP = X_LOC->max_norm() ;      
      
   } while( !convergence_achieved( iter ) ) ;

   stop_total_timer() ;
   // ---------------
}

//------------------------------------------------------------------------
bool
AP_NavierStokes2CFV:: convergence_achieved( size_t iter ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: convergence_achieved" ) ;
   
   double rhs = F->two_norm() ;
   bool result = ( rhs < TOL ) ;

   if( verbose_level() >= 1 )
   {
      ios_base::fmtflags original_flags = PEL::out().flags() ;
      PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

      if( iter==1 )
      {
         PEL::out() << indent()
                    << setw( 3+15 ) << "max(d_U)" 
                    << setw( 15 )   << "max(d_P)" 
                    << setw( 15 )   << "L2(rhs)" ;
         if( SOLVER->is_iterative() )
            PEL::out() << setw( 9 ) << "its" ;
         PEL::out() << endl ;
      }
      PEL::out() << indent() 
                 << setw( 3 ) << iter
                 << setprecision( 6 ) << setw( 15 )
                 << PRINTED_RES_UU
                 << setprecision( 6 ) << setw( 15 ) 
                 << PRINTED_RES_PP
                 << setprecision( 6 ) << setw( 15 ) 
                 << rhs ;
      if( SOLVER->is_iterative() )
      {
         PEL::out() << setw( 9 ) << SOLVER->nb_iterations_achieved() ;
      }
      PEL::out() << endl ;
      PEL::out().flags( original_flags ) ;
   }

   return( result ) ;
}
   
//************************************************************************
// L'implementation de cette fonction est trop longue et
// devrait etre splittee pour mieux faire apparaitre les différents termes
//************************************************************************
//------------------------------------------------------------------------
void
AP_NavierStokes2CFV:: loop_on_sides( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: loop_on_sides" ) ;

   size_t nb_dims = sFE->nb_space_dimensions() ;
   double xx  = PEL::bad_double() ;
   double rhs = PEL::bad_double() ;

   size_t nb_done = 0 ;
   size_t nb_skip = 0 ;

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Color const* color = sFE->color() ;
      bool do_stab = ( CLUSTER_BD == 0 ) ||
           ( (CLUSTER_BD != 0) && !color->is_matching( CLUSTER_BD ) ) ;

      size_t n_v_K = fe_K->global_node( UU, 0 ) ;
      size_t n_v_L = fe_L->global_node( UU, 0 ) ;
      size_t n_p_K = fe_K->global_node( PP, 0 ) ;
      size_t n_p_L = fe_L->global_node( PP, 0 ) ;

      size_t u_p_K, u_p_L ;
      double value_p_K, value_p_L ;

      if( NMB->link( idx_PP )->DOF_is_unknown( n_p_K, 0 ) )
      {
         u_p_K = NMB->global_unknown_for_DOF( n_p_K, 0, idx_PP ) ;
         value_p_K = PP->DOF_value( L_UPDATE_UU, n_p_K, 0 ) ;
      }
      else
      {
         u_p_K = PEL::bad_index() ;
         value_p_K = PP->DOF_imposed_value( n_p_K, 0 ) ;
      }

      if( NMB->link( idx_PP )->DOF_is_unknown( n_p_L, 0 ) )
      {
         u_p_L = NMB->global_unknown_for_DOF( n_p_L, 0, idx_PP ) ;
         value_p_L = PP->DOF_value( L_UPDATE_UU, n_p_L, 0 ) ;
      }
      else
      {
         u_p_L = PEL::bad_index() ;
         value_p_L = PP->DOF_imposed_value( n_p_L, 0 ) ;
      }

      GE_Vector const* normal = sFE->normal() ;

      double area = FE::side_measure( sFE ) ;

      double d_K = sFE->distance_to_adjacent_finite_volume_center( 0 ) ;
      double d_L = sFE->distance_to_adjacent_finite_volume_center( 1 ) ;
      double d_KL = d_K + d_L ;

      double rhs_stab = 0.0 ;
      double csta = 0.0 ;
      if( do_stab )
      {
         nb_done++ ;
         double h_K = fe_K->polyhedron()->inter_vertices_maximum_distance() ;
         double h_L = fe_L->polyhedron()->inter_vertices_maximum_distance() ;
         csta = LAMBDA * ( PEL::pow( h_K, h_EXPONENT ) + 
                           PEL::pow( h_L, h_EXPONENT ) ) * area ;

         rhs_stab = csta * ( value_p_K - value_p_L ) ;
         if( u_p_K != PEL::bad_index() ) add_to_G_item( u_p_K,  rhs_stab ) ;
         if( u_p_L != PEL::bad_index() ) add_to_G_item( u_p_L, -rhs_stab ) ;

         xx = -csta ;
         if( u_p_K != PEL::bad_index() ) add_to_C_item( u_p_K, u_p_K,  xx ) ;
         if( u_p_L != PEL::bad_index() && u_p_K != PEL::bad_index() ) 
            add_to_C_item( u_p_L, u_p_K, -xx ) ;

         xx = -xx ;
         if( u_p_K != PEL::bad_index() && u_p_K != PEL::bad_index() ) 
            add_to_C_item( u_p_K, u_p_L,  xx ) ;
         if( u_p_L != PEL::bad_index() ) add_to_C_item( u_p_L, u_p_L, -xx ) ;
      }
      else
      {
         nb_skip++ ;
      }

      double conv_flux = 0.0 ; 
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         conv_flux += ( UU->DOF_value( L_UPDATE_UU, n_v_K, ic ) * d_L + 
                        UU->DOF_value( L_UPDATE_UU, n_v_L, ic ) * d_K ) 
                      * normal->component( ic ) ;
      }
      conv_flux = conv_flux * area / d_KL ;

      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         size_t u_v_K = NMB->global_unknown_for_DOF( n_v_K, ic, idx_UU ) ;
         size_t u_v_L = NMB->global_unknown_for_DOF( n_v_L, ic, idx_UU ) ;

         double value_v_K = UU->DOF_value( L_UPDATE_UU, n_v_K, ic ) ;
         double value_v_L = UU->DOF_value( L_UPDATE_UU, n_v_L, ic ) ;

         // *** div( v v )

         rhs = - C_VGRAD * conv_flux * 0.5 * ( value_v_K + value_v_L ) 
               - C_VGRAD * rhs_stab  * 0.5 * value_v_K ;
         add_to_F_item( u_v_K,  rhs ) ;
         add_to_F_item( u_v_L, -rhs ) ;

         xx = C_VGRAD * ( conv_flux + rhs_stab ) * 0.5 ;
         add_to_A_item( u_v_K, u_v_K,  xx ) ;
         add_to_A_item( u_v_L, u_v_K, -xx ) ;

         xx = C_VGRAD * conv_flux * 0.5 ;
         add_to_A_item( u_v_K, u_v_L,  xx ) ;
         add_to_A_item( u_v_L, u_v_L, -xx ) ;

         for( size_t jc=0 ; jc<nb_dims ; ++jc )
         {
            size_t uKj = NMB->global_unknown_for_DOF( n_v_K, jc, idx_UU ) ;
            xx = C_VGRAD * 0.5 *( value_v_K + value_v_L ) * 
                              area * d_L / d_KL * normal->component( jc ) ;
            add_to_A_item( u_v_K, uKj,  xx ) ;
            add_to_A_item( u_v_L, uKj, -xx ) ;

            size_t uLj = NMB->global_unknown_for_DOF( n_v_L, jc, idx_UU ) ;
            xx = C_VGRAD * 0.5 *( value_v_K + value_v_L ) * 
                              area * d_K / d_KL * normal->component( jc ) ;
            add_to_A_item( u_v_K, uLj,  xx ) ;
            add_to_A_item( u_v_L, uLj, -xx ) ;
         }
         
         if( do_stab )
         {
            xx = C_VGRAD * csta * 0.5 * value_v_K ;
            add_to_D_item( u_v_K, u_p_K,  xx ) ;
            add_to_D_item( u_v_L, u_p_K, -xx ) ;

            xx = -xx ;
            add_to_D_item( u_v_K, u_p_L,  xx ) ;
            add_to_D_item( u_v_L, u_p_L, -xx ) ;
         }

         // *** - delta(v) 

         xx = MU * area / d_KL ;
         add_to_A_item( u_v_K, u_v_K,  xx ) ;
         add_to_A_item( u_v_L, u_v_K, -xx ) ;

         rhs = - xx * ( value_v_K - value_v_L ) ;
         add_to_F_item( u_v_K,  rhs ) ;
         add_to_F_item( u_v_L, -rhs ) ;

         xx = - MU * area / d_KL ;
         add_to_A_item( u_v_K, u_v_L,  xx ) ;
         add_to_A_item( u_v_L, u_v_L, -xx ) ;

         // *** -div(v)

         xx = - area * d_L / d_KL * normal->component( ic ) ;
         add_to_B_item( u_p_K, u_v_K,  xx ) ;
         add_to_D_item( u_v_K, u_p_K,  xx ) ;
         add_to_B_item( u_p_L, u_v_K, -xx ) ;
         add_to_D_item( u_v_K, u_p_L, -xx ) ;

         xx = - area * d_K / d_KL * normal->component( ic ) ;
         add_to_B_item( u_p_K, u_v_L,  xx ) ;
         add_to_D_item( u_v_L, u_p_K,  xx ) ;
         add_to_B_item( u_p_L, u_v_L, -xx ) ;
         add_to_D_item( u_v_L, u_p_L, -xx ) ;

         rhs = area * ( d_L * value_v_K + d_K * value_v_L ) / d_KL 
                    * normal->component( ic ) ;
         add_to_G_item( u_p_K,  rhs ) ;
         add_to_G_item( u_p_L, -rhs ) ;

         rhs = area * d_L * ( value_p_K - value_p_L ) / d_KL
                    * normal->component( ic ) ;
         add_to_F_item( u_v_K,  rhs ) ;
         rhs = area * d_K * ( value_p_K - value_p_L ) / d_KL
                    * normal->component( ic ) ;
         add_to_F_item( u_v_L,  rhs ) ;
      }
   }

   if( verbose_level() >= 2 )
   {
      PEL::out() << indent() << "   stabilization computed on " 
                 << nb_done << " sides (" << nb_skip << " skipped)" << endl ;
   }
}

//------------------------------------------------------------------------
void
AP_NavierStokes2CFV:: loop_on_bounds( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: loop_on_bounds" ) ;

   size_t nb_dims = bFE->nb_space_dimensions() ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Color const* color = bFE->color() ;

      PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;

      size_t n_v_K = bFE->global_node( UU, 0 ) ;
      size_t n_p_K = bFE->global_node( PP, 0 ) ;

      GE_Vector const* normal = bFE->outward_normal() ;
      double area = FE::bound_measure( bFE ) ;

      double conv_flux = 0.0 ;

      string const& type = ee->string_data( "type" ) ;
      if( type == "Dirichlet" )
      {
         GE_Point const* pt = bFE->polyhedron()->center() ;
         XX->set( pt->coordinate_vector() ) ;
         TT->set( t_it->time() ) ;
         doubleVector uimp = 
                      ee->doubleVector_data( "imposed_value", CONTEXT ) ;
         if( uimp.size() != nb_dims ) 
            AP_NavierStokes2CFV_ERROR::n1( ee, "imposed_value", nb_dims );

         double h = bFE->distance_to_adjacent_finite_volume_center() ;

         double xx_l = MU * area / h ;

         double uimp_nor = 0.0 ;
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            size_t u_v_K = NMB->global_unknown_for_DOF( n_v_K, ic, idx_UU ) ;
            double value_v_K = UU->DOF_value( L_UPDATE_UU, n_v_K, ic ) ;
            add_to_A_item( u_v_K, u_v_K, xx_l ) ;
            add_to_F_item( u_v_K, -xx_l*( value_v_K - uimp( ic ) ) ) ; 

            uimp_nor += uimp( ic ) * normal->component( ic ) ;
         }
         
         size_t u_p_K ;
         double value_p_K ;

         if( NMB->link( idx_PP )->DOF_is_unknown( n_p_K, 0 ) )
         {
            u_p_K = NMB->global_unknown_for_DOF( n_p_K, 0, idx_PP ) ;
            value_p_K = PP->DOF_value( L_UPDATE_UU, n_p_K, 0 ) ;
         }
         else
         {
            u_p_K = PEL::bad_index() ;
            value_p_K = PP->DOF_imposed_value( n_p_K, 0 ) ;
         }

         if( u_p_K != PEL::bad_index() ) add_to_G_item( u_p_K, uimp_nor*area ) ;

         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            conv_flux += uimp( ic ) * normal->component( ic ) ;
         }
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            size_t u_v_K = NMB->global_unknown_for_DOF( n_v_K, ic, idx_UU ) ;
            add_to_F_item( u_v_K, - C_VGRAD * area * conv_flux * uimp(ic) ) ; 
         }
      }
      else
         raise_bad_BC_type( type, "\"Dirichlet\"", UU ) ;
   }
}

//------------------------------------------------------------------------
void
AP_NavierStokes2CFV:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: loop_on_cells" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      size_t n_v_K = cFE->global_node( UU, 0 ) ;

      doubleVector int_source( nb_dims ) ;
      cFE->start_IP_iterator( QRP_PI ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            int_source( d ) += cFE->weight_of_IP() * 
                               PI->cell_value_at_IP( t_it, cFE, d ) ;
         }
      }

      double vol = FE::cell_measure( cFE ) ;
      double xx_l = vol * C_DT / t_it->time_step();

      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         size_t u_v_K = NMB->global_unknown_for_DOF( n_v_K, ic, idx_UU ) ;
         double xx_r = - int_source( ic ) ;
         if( IS_UNSTEADY )
         {
            xx_r += xx_l * ( UU->DOF_value( L_UPDATE_UU, n_v_K, ic ) - 
                             UU->DOF_value( L_EXPLICIT_UU, n_v_K, ic ) ) ;
            add_to_A_item( u_v_K, u_v_K, xx_l ) ;
         }
         add_to_F_item( u_v_K, -xx_r ) ; 
      }
   }
}

//************************************************************************
// Les fonctions add_to_*_item sont des survivances du passé
// ---> a supprimer ???
//************************************************************************

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_A_item( size_t i_row, 
                                          size_t j_col, 
                                          double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_A_item" ) ;

   // no imposed dof for velocity
   PEL_ASSERT( i_row != PEL::bad_index() ) ;
   PEL_ASSERT( j_col != PEL::bad_index() ) ;

   A->add_to_item( i_row, j_col, xx ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_F_item( size_t i_row, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_F_item" ) ;

   PEL_ASSERT( i_row != PEL::bad_index() ) ;

   F->add_to_item( i_row, xx ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_B_item( size_t i_row, 
                                          size_t j_col, 
                                          double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_B_item" ) ;

   PEL_ASSERT( j_col != PEL::bad_index() ) ;

   if( i_row != PEL::bad_index() )
   {
      A->add_to_item( i_row, j_col, xx ) ;
   }
}

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_D_item( size_t i_row, 
                                          size_t j_col, 
                                          double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_D_item" ) ;

   PEL_ASSERT( i_row != PEL::bad_index() ) ;

   if( j_col != PEL::bad_index() )
   {
      A->add_to_item( i_row, j_col, xx ) ;
   }
}

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_G_item( size_t i_row, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_G_item" ) ;

   if( i_row != PEL::bad_index() )
   {
      F->add_to_item( i_row, xx ) ;
   }
}

//----------------------------------------------------------------------
void
AP_NavierStokes2CFV:: add_to_C_item( size_t i_row, 
                                          size_t j_col, 
                                          double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes2CFV:: add_to_C_item" ) ;

   if( ( i_row != PEL::bad_index() ) && ( j_col != PEL::bad_index() ) )
   {
      A->add_to_item( i_row, j_col, xx ) ;
   }
}

//internal--------------------------------------------------------------
void
AP_NavierStokes2CFV_ERROR:: n0( PDE_LocalFEbound const* fe,
                                     double int_v_nor )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_NavierStokes2CFV:" 
       << std::endl << std::endl ;
   msg << "the bound : " << endl ;
   fe->print_current_mesh( msg , 6 ) ;
   msg << "is an inflow bound, but it does not have" << endl ;
   msg << "a Dirichlet boundary condition" << endl ;
   msg << "advective flux : " << int_v_nor ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
AP_NavierStokes2CFV_ERROR:: n1( PEL_ModuleExplorer const* exp,
                                     std::string const& keyword,
                                     size_t size )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_NavierStokes2CFV:" 
       << std::endl << std::endl ;
   msg << "In module: " << exp->absolute_path_name() << endl ;
   msg << "the data of keyword: " << keyword << endl ;
   msg << "should be a DoubleVector of size: " << size ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
AP_NavierStokes2CFV_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_NavierStokes2CFV:" 
       << std::endl << std::endl ;
   msg << "the data of keyword: \"velocity_level_to_update\"" << endl ;
   msg << "                and: \"level_of_explicit_velocity\"" << endl ;
   msg << "should be different" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
