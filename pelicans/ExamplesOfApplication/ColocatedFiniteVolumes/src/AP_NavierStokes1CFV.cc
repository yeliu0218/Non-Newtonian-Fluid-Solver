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

#include <AP_NavierStokes1CFV.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <LA_SeqMatrix.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <AP_NavierStokes1System.hh>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using std::cout ;
using std::string ;
using std::endl ;

AP_NavierStokes1CFV const* 
AP_NavierStokes1CFV::PROTOTYPE = new AP_NavierStokes1CFV() ;

struct AP_NavierStokes1CFV_ERROR
{
   static void n0( PDE_LocalFEbound const* fe, double int_v_nor ) ;
   static void n1( PEL_ModuleExplorer const* exp,
                   std::string const& keyword, 
                   size_t size) ;
} ;

//---------------------------------------------------------------------------
AP_NavierStokes1CFV:: AP_NavierStokes1CFV( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_NavierStokes1CFV" )
{
}

//---------------------------------------------------------------------------
AP_NavierStokes1CFV*
AP_NavierStokes1CFV:: create_replica( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_NavierStokes1CFV* result = new AP_NavierStokes1CFV( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes1CFV:: AP_NavierStokes1CFV( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( 
                                        exp->string_data( "velocity" ) ) )
   , L_UPDATE_UU( exp->int_data( "velocity_level_to_update" ) )
   , L_EXPLICIT_UU( exp->int_data( "level_of_explicit_velocity" ) )
   , PP( dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "pressure" ) ) )
   , L_UPDATE_PP( exp->int_data( "pressure_level_to_update" ) )
   , L_EXPLICIT_PP( exp->int_data( "level_of_explicit_pressure" ) )
   , ALPHA( exp->double_data( "coef_unsteady" ) )
   , MU( exp->double_data( "coef_viscous" ) )
   , PI( prms->item( exp->string_data( "param_source" ) ) )
   , LAMBDA( exp->double_data( "lambda_in_pressure_stabilization" ) ) 
   , h_EXPONENT( exp->double_data( "h_exponent_in_infsup_stabilization" ) )
   , CLUSTER_BD( 0 )
   , BCs( dom->set_of_boundary_conditions() )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP_PI( GE_QRprovider::object( 
             exp->string_data( "quadrature_rule_provider_for_source" ) ) )
   , GLOBAL_EQ( 0 )
   , CONTEXT( PEL_ContextSimple::create( this ) )
   , XX( 0 )
   , TT( 0 )
   , T_EXPLICIT( PEL::bad_double() )
   , T_UPDATE( PEL::bad_double() )
{
   PEL_LABEL( "AP_NavierStokes1CFV:: AP_NavierStokes1CFV" ) ;

   check_field_storage_depth( UU, L_UPDATE_UU ) ;
   check_field_storage_depth( UU, L_EXPLICIT_UU ) ;
   check_field_storage_depth( PP, L_UPDATE_PP ) ;
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

   PEL_ModuleExplorer* se =
                    exp->create_subexplorer( 0, "AP_NavierStokes1System" ) ;
   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
   PDE_LinkDOF2Unknown* pp_link = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
   GLOBAL_EQ = AP_NavierStokes1System::create( this, se, uu_link, pp_link ) ;
   se->destroy() ;

   XX = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object("DV_X"), XX ) ;

   TT = PEL_Double::create( CONTEXT, 0.0 ) ;
   CONTEXT->extend( PEL_Variable::object("DS_T"), TT ) ;
}

//---------------------------------------------------------------------------
AP_NavierStokes1CFV:: ~AP_NavierStokes1CFV( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_NavierStokes1CFV:: do_before_time_stepping( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "AP_NavierStokes1CFV:: do_before_time_stepping" ) ;
   // --------------

   T_EXPLICIT = t_it->initial_time() ;
   T_UPDATE   = t_it->initial_time() ;

   if( verbose_level() >= 2 )
   {
      PEL::out() << indent() << "   integral of pressure(" 
                             << L_UPDATE_PP << ") = " 
                             << pressure_integral( L_UPDATE_PP ) << endl ;
   }

   stop_total_timer() ;
   // -------------
}

//---------------------------------------------------------------------------
void 
AP_NavierStokes1CFV:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_NavierStokes1CFV:: do_one_inner_iteration" ) ;
   // --------------

   T_EXPLICIT = T_UPDATE ;
   T_UPDATE   = t_it->time() ;

   //??? a nullify may be enough ?
   GLOBAL_EQ->re_initialize() ;

   double xx = 1.0/t_it->time_step() ;
   GLOBAL_EQ->set_leading_BDF_over_dt( xx ) ;

   start_assembling_timer() ;
   // ---------------------

   loop_on_cells( t_it ) ;
   loop_on_sides( t_it ) ;
   loop_on_bounds( t_it ) ;

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------

   GLOBAL_EQ->set_indent( indent() ) ;
   GLOBAL_EQ->set_initial_guess_P( PP, L_EXPLICIT_PP ) ; 
   GLOBAL_EQ->set_initial_guess_U( UU, L_UPDATE_UU ) ; 
   GLOBAL_EQ->estimate_unknowns() ;

   stop_solving_timer() ;
   // ------------------

   if( !GLOBAL_EQ->unknowns_are_solution() )
   {
      PEL_Error::object()->display_info(
         "*** AP_NavierStokes1CFV error\n"
         "    unable to solve the discrete problem" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "   update of " << UU->name()
                    << "(" << L_UPDATE_UU << ") and " << PP->name()
                    << "(" << L_UPDATE_PP << ")" << endl ;
      }
      UU->update_free_DOFs_value( L_UPDATE_UU,
                                  GLOBAL_EQ->unknown_vector_U(),
                                  GLOBAL_EQ->linkDOF2Unknown_U() ) ;
      PP->update_free_DOFs_value( L_UPDATE_PP,
                                  GLOBAL_EQ->unknown_vector_P(), 
                                  GLOBAL_EQ->linkDOF2Unknown_P() ) ;

      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "   integral of pressure("
                    << L_UPDATE_PP << ") = "
                    << pressure_integral( L_UPDATE_PP ) << endl ;
      }
   }

   stop_total_timer() ;
   // ----------------
}

//------------------------------------------------------------------------
void
AP_NavierStokes1CFV:: loop_on_sides( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: loop_on_sides" ) ;

   size_t nb_dims = sFE->nb_space_dimensions() ;
   double xx = PEL::bad_double() ;

   size_t nb_done = 0 ;
   size_t nb_skip = 0 ;

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Color const* color = sFE->color() ;

      size_t n_v_K = fe_K->global_node( UU, 0 ) ;
      size_t n_v_L = fe_L->global_node( UU, 0 ) ;
      size_t n_p_K = fe_K->global_node( PP, 0 ) ;
      size_t n_p_L = fe_L->global_node( PP, 0 ) ;

      GE_Vector const* normal = sFE->normal() ;

      double area = FE::side_measure( sFE ) ;

      double d_K = sFE->distance_to_adjacent_finite_volume_center( 0 ) ;
      double d_L = sFE->distance_to_adjacent_finite_volume_center( 1 ) ;
      double d_KL = d_K + d_L ;

      size_t u_p_K = GLOBAL_EQ->global_unknown_for_DOF_of_P( n_p_K ) ;
      size_t u_p_L = GLOBAL_EQ->global_unknown_for_DOF_of_P( n_p_L ) ;

      // *** matrix pressure - pressure : stabilization

      double stab_explicit = 0.0 ; 
      if( ( CLUSTER_BD == 0) ||
          ( (CLUSTER_BD != 0) && !color->is_matching( CLUSTER_BD ) ) )
      {
         nb_done++ ;
         double h_K = fe_K->polyhedron()->inter_vertices_maximum_distance() ;
         double h_L = fe_L->polyhedron()->inter_vertices_maximum_distance() ;
         double csta = LAMBDA * ( PEL::pow( h_K, h_EXPONENT ) + 
                                  PEL::pow( h_L, h_EXPONENT ) ) * area ;
         xx = -csta ;
         GLOBAL_EQ->add_to_C_item( u_p_K, u_p_K,  xx ) ;
         GLOBAL_EQ->add_to_C_item( u_p_L, u_p_K, -xx ) ;

         xx = -xx ;
         GLOBAL_EQ->add_to_C_item( u_p_K, u_p_L,  xx ) ;
         GLOBAL_EQ->add_to_C_item( u_p_L, u_p_L, -xx ) ;

         stab_explicit = csta * ( PP->DOF_value( L_EXPLICIT_PP, n_p_K ) - 
                                  PP->DOF_value( L_EXPLICIT_PP, n_p_L ) ) ;
      }
      else
      {
         nb_skip++ ;
      }
      
      double conv_flux = 0 ; 
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         conv_flux += ( UU->DOF_value( L_EXPLICIT_UU, n_v_K, ic ) * d_L + 
                        UU->DOF_value( L_EXPLICIT_UU, n_v_L, ic ) * d_K ) 
                      * normal->component( ic ) ;
      }
      conv_flux = conv_flux * area / d_KL + stab_explicit ;
      
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         size_t u_v_K = GLOBAL_EQ->global_unknown_for_DOF_of_U( n_v_K, ic ) ;
         size_t u_v_L = GLOBAL_EQ->global_unknown_for_DOF_of_U( n_v_L, ic ) ;

         // *** matrix velocity - velocity : div(v* v) - delta(v)
        
         xx = MU * area / d_KL + ALPHA * conv_flux * 0.5 ;
         GLOBAL_EQ->add_to_A_item( u_v_K, u_v_K,  xx ) ;
         GLOBAL_EQ->add_to_A_item( u_v_L, u_v_K, -xx ) ;

         xx = - MU * area / d_KL + ALPHA * conv_flux * 0.5 ;
         GLOBAL_EQ->add_to_A_item( u_v_K, u_v_L,  xx ) ;
         GLOBAL_EQ->add_to_A_item( u_v_L, u_v_L, -xx ) ;

         // *** matrix pressure - velocity : -div(v)

         xx = - area * d_L / d_KL * normal->component( ic ) ;
         GLOBAL_EQ->add_to_B_item( u_p_K, u_v_K,  xx ) ;
         GLOBAL_EQ->add_to_B_item( u_p_L, u_v_K, -xx ) ;

         xx = - area * d_K / d_KL * normal->component( ic ) ;
         GLOBAL_EQ->add_to_B_item( u_p_K, u_v_L,  xx ) ;
         GLOBAL_EQ->add_to_B_item( u_p_L, u_v_L, -xx ) ;
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
AP_NavierStokes1CFV:: loop_on_bounds( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: loop_on_bounds" ) ;

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
         TT->set( T_UPDATE ) ;
         doubleVector uimp = 
                      ee->doubleVector_data( "imposed_value", CONTEXT ) ;
         if( uimp.size() != nb_dims ) 
            AP_NavierStokes1CFV_ERROR::n1( ee, "imposed_value", nb_dims ) ;

         double h = bFE->distance_to_adjacent_finite_volume_center() ;

         double xx_l = MU * area / h ;

         double uimp_nor = 0.0 ;
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            size_t u_v_K = GLOBAL_EQ->global_unknown_for_DOF_of_U( n_v_K, ic ) ;
            GLOBAL_EQ->add_to_A_item( u_v_K, u_v_K, xx_l ) ;
            GLOBAL_EQ->add_to_F_item( u_v_K, xx_l*uimp(ic) ) ; 

            uimp_nor += uimp( ic ) * normal->component( ic ) ;
         }
         
         size_t u_p_K = GLOBAL_EQ->global_unknown_for_DOF_of_P( n_p_K ) ;
         GLOBAL_EQ->add_to_G_item( u_p_K, uimp_nor*area ) ;

         TT->set( T_EXPLICIT ) ;
         doubleVector uimp_explicit = 
                      ee->doubleVector_data( "imposed_value", CONTEXT ) ;
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            conv_flux += uimp_explicit( ic ) * normal->component( ic ) ;
         }
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            size_t u_v_K = GLOBAL_EQ->global_unknown_for_DOF_of_U( n_v_K, ic ) ;
            GLOBAL_EQ->add_to_F_item( u_v_K, 
                                      -ALPHA * area * conv_flux * uimp(ic) ) ; 
         }
      }
      else
         raise_bad_BC_type( type, "\"Dirichlet\"", UU ) ;
   }
}

//------------------------------------------------------------------------
void
AP_NavierStokes1CFV:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: loop_on_cells" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      size_t n_v_K = cFE->global_node( UU, 0 ) ;
      size_t n_p_K = cFE->global_node( PP, 0 ) ;

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
      double xx_l = vol * ALPHA / t_it->time_step();

      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         size_t u_v_K = GLOBAL_EQ->global_unknown_for_DOF_of_U( n_v_K, ic ) ;
         double xx_r = int_source( ic ) 
                     + xx_l * UU->DOF_value( L_EXPLICIT_UU, n_v_K, ic ) ;
         GLOBAL_EQ->add_to_F_item( u_v_K, xx_r ) ; 
         GLOBAL_EQ->add_to_A_item( u_v_K, u_v_K, xx_l ) ;
      }

      if( GLOBAL_EQ->MPl_is_required() )
      {
         size_t u_p_K = GLOBAL_EQ->global_unknown_for_DOF_of_P( n_p_K ) ;
         GLOBAL_EQ->add_to_MPl_item( u_p_K, vol ) ;
      }
   }
}

//------------------------------------------------------------------------
double
AP_NavierStokes1CFV:: pressure_integral( size_t level )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1CFV:: pressure_integral" ) ;

   double result = 0.0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      size_t n_p_K = cFE->global_node( PP, 0 ) ;
      result += cFE->polyhedron()->measure() * 
                PP->DOF_value( level, n_p_K, 0 ) ;
   }
   return( result ) ;
}

//internal--------------------------------------------------------------
void
AP_NavierStokes1CFV_ERROR:: n0( PDE_LocalFEbound const* fe,
                             double int_v_nor )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_NavierStokes1CFV:" 
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
AP_NavierStokes1CFV_ERROR:: n1( PEL_ModuleExplorer const* exp,
                             std::string const& keyword,
                             size_t size )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_NavierStokes1CFV:" 
       << std::endl << std::endl ;
   msg << "In module: " << exp->absolute_path_name() << endl ;
   msg << "the data of keyword: " << keyword << endl ;
   msg << "should be a DoubleVector of size: " << size ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
