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

#include <AP_AdvectionDiffusion1CFV.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_MemoryTracer.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

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
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using std::string ;
using std::endl ;

AP_AdvectionDiffusion1CFV const* 
AP_AdvectionDiffusion1CFV::PROTOTYPE = new AP_AdvectionDiffusion1CFV() ;

struct AP_AdvectionDiffusion1CFV_ERROR
{
   static void n0( PDE_LocalFEbound const* fe, double int_v_nor ) ;
   static void n1( PEL_ModuleExplorer const* exp,
                   std::string const& keyword, 
                   size_t size) ;
} ;

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1CFV:: AP_AdvectionDiffusion1CFV( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_AdvectionDiffusion1CFV" )
{
}

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1CFV*
AP_AdvectionDiffusion1CFV:: create_replica( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_AdvectionDiffusion1CFV* result = 
                      new AP_AdvectionDiffusion1CFV( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1CFV:: AP_AdvectionDiffusion1CFV( PEL_Object* a_owner,
                                                 PDE_DomainAndFields const* dom,
                                                 FE_SetOfParameters const* prms,
                                                 PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "unknown_field" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , ADV( prms->item( exp->string_data( "param_advective_velocity" ) ) )
   , ALPHA( exp->double_data( "coef_unsteady" ) )
   , KAPPA( exp->double_data( "coef_diffusion" ) )
   , PI( prms->item( exp->string_data( "param_source" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP_ADV( GE_QRprovider::object( 
     exp->string_data( "quadrature_rule_provider_for_advective_velocity" ) ) )
   , QRP_PI( GE_QRprovider::object( 
             exp->string_data( "quadrature_rule_provider_for_source" ) ) )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: AP_AdvectionDiffusion1CFV" ) ;

   check_param_nb_components( ADV, "param_advective_velocity", 
                              dom->nb_space_dimensions() ) ;

   sFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::node ) ;

   ADV->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
   ADV->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
   PI->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
   PI->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
   PI->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU, true ) ;
   NMB = PDE_SystemNumbering::create( this, uu_link ) ;  

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//---------------------------------------------------------------------------
AP_AdvectionDiffusion1CFV:: ~AP_AdvectionDiffusion1CFV( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_AdvectionDiffusion1CFV:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_AdvectionDiffusion1CFV:: do_one_inner_iteration" ) ;
   // --------------

   PEL_MemoryTracer::object()->start_event(
      "AP_AdvectionDiffusion1CFV::do_one_inner_iteration \""+UU->name()+"\"" ) ;
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
   PEL_MemoryTracer::object()->start_event( "assemble" ) ;
   loop_on_cells( t_it ) ;
   loop_on_sides( t_it ) ;
   loop_on_bounds( t_it ) ;
   PEL_MemoryTracer::object()->stop_event() ;

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ----------------

   PEL_MemoryTracer::object()->start_event( "estimate_unknown" ) ;
   A->synchronize() ;
   F->synchronize() ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
   PEL_MemoryTracer::object()->stop_event() ;
   
   stop_solving_timer() ;
   // ---------------
   
   if( ! SOLVER->solution_is_achieved() )
   {
      PEL_Error::object()->display_info(
         "*** AP_AdvectionDiffusion1CFV error\n"
         "    convergence failure when solving the discrete problem" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "   update of " << UU->name()
                    << "(" << L_UPDATE << ")" << endl ;
      }
      LA_Scatter const* sca = NMB->scatter() ;
      sca->get( X, X_LOC ) ;
      UU->update_free_DOFs_value( L_UPDATE, X_LOC, NMB->link() ) ;
   }
   
   stop_total_timer() ;
   // -------------

   PEL_MemoryTracer::object()->stop_event() ;
}

//------------------------------------------------------------------------
void
AP_AdvectionDiffusion1CFV:: loop_on_sides( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: loop_on_sides" ) ;

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      size_t n_K = fe_K->global_node( UU, 0 ) ;
      size_t n_L = fe_L->global_node( UU, 0 ) ;

      double int_v_nor = side_advective_flux( t_it ) ;

      double area = ( sFE->is_periodic() ? 
                         FE::side_measure( sFE, 0 ) : 
                      FE::side_measure( sFE ) ) ;

      double d_KL = sFE->distance_to_adjacent_finite_volume_center( 0 ) +
                    sFE->distance_to_adjacent_finite_volume_center( 1 ) ;

      double xx_0 =   KAPPA * area / d_KL ;
      double xx_1 = - KAPPA * area / d_KL ;

      if( int_v_nor > 0.0 ) xx_0 += int_v_nor ;
      if( int_v_nor < 0.0 ) xx_1 += int_v_nor ;

      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;
      size_t u_L = NMB->global_unknown_for_DOF( n_L, 0 ) ;

      A->add_to_item( u_K, u_K,  xx_0 ) ;
      A->add_to_item( u_L, u_K, -xx_0 ) ;

      A->add_to_item( u_K, u_L,  xx_1 ) ;
      A->add_to_item( u_L, u_L, -xx_1 ) ;
   }
}

//------------------------------------------------------------------------
void
AP_AdvectionDiffusion1CFV:: loop_on_bounds( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: loop_on_bounds" ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Color const* color = bFE->color() ;

      PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;

      size_t n_K = bFE->global_node( UU, 0 ) ;
      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;
      
      double int_v_nor = bound_advective_flux( t_it ) ;

      double area = FE::bound_measure( bFE ) ;

      string const& type = ee->string_data( "type" ) ;
      if( type == "Dirichlet" )
      {
         doubleVector const& uimp = ee->doubleVector_data( "imposed_value" ) ;
         if( uimp.size() != 1 )
            AP_AdvectionDiffusion1CFV_ERROR::n1( ee, "imposed_value", 1 ) ;

         double h = bFE->distance_to_adjacent_finite_volume_center() ;

         double xx_l = KAPPA * area / h ;
         double xx_r = xx_l * uimp( 0 ) - int_v_nor * uimp( 0 ) ;

         A->add_to_item( u_K, u_K, xx_l ) ;
         F->add_to_item( u_K, xx_r ) ; 
      }
      else if( type == "NeumannScalarCFV" )
      {
         doubleVector const& fimp = ee->doubleVector_data( "flux_value" ) ;
         if( fimp.size() != 1 )
            AP_AdvectionDiffusion1CFV_ERROR::n1( ee, "flux_value", 1 ) ;
         F->add_to_item( u_K, area*fimp( 0 ) ) ;
      }
      else
         raise_bad_BC_type( type, "\"Dirichlet\"\n\"NeumannScalarCFV\"", UU ) ;

      if( type != "Dirichlet" )
      {
         if( int_v_nor < 0.0 ) 
            AP_AdvectionDiffusion1CFV_ERROR::n0( bFE, int_v_nor ) ;
      }
   }
}

//------------------------------------------------------------------------
void
AP_AdvectionDiffusion1CFV:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: loop_on_cells" ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      size_t n_K = cFE->global_node( UU, 0 ) ;

      double int_source = 0.0 ;
      cFE->start_IP_iterator( QRP_PI ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         int_source += cFE->weight_of_IP() * 
                       PI->cell_value_at_IP( t_it, cFE ) ;
      }

      double vol = FE::cell_measure( cFE ) ;

      double xx_l = vol * ALPHA / t_it->time_step() ;
      double xx_r = int_source + xx_l * UU->DOF_value( L_EXPLICIT, n_K ) ;

      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;

      A->add_to_item( u_K, u_K,  xx_l ) ;
      F->add_to_item( u_K, xx_r ) ;
   }
}

//------------------------------------------------------------------------
double
AP_AdvectionDiffusion1CFV:: side_advective_flux( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: side_advective_flux" ) ;
   
   size_t nb_dims = sFE->nb_space_dimensions() ;

   GE_Vector const* normal = ( sFE->is_periodic() ? 
                               sFE->normal( 0 ) : sFE->normal() ) ;
   double result = 0.0 ;
   sFE->start_IP_iterator( QRP_ADV ) ;
   for( ; sFE->valid_IP() ; sFE->go_next_IP() )
   {
      double ww = ( sFE->is_periodic() ? 
                    sFE->weight_of_IP( 0 ) : sFE->weight_of_IP() ) ;
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         result += ww * ADV->side_value_at_IP( t_it, sFE, ic ) * 
                        normal->component( ic ) ;
      }
   }
   return( result ) ;
}

//------------------------------------------------------------------------
double
AP_AdvectionDiffusion1CFV:: bound_advective_flux( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_AdvectionDiffusion1CFV:: bound_advective_flux" ) ;
   
   size_t nb_dims = bFE->nb_space_dimensions() ;

   GE_Vector const* normal = bFE->outward_normal() ;

   double result = 0.0 ;
   bFE->start_IP_iterator( QRP_ADV ) ;
   for( ; bFE->valid_IP() ; bFE->go_next_IP() )
   {
      double ww = bFE->weight_of_IP() ;
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         result += ww * ADV->bound_value_at_IP( t_it, bFE, ic ) * 
                        normal->component( ic ) ;
      }
   }
   return( result ) ;
}

//internal--------------------------------------------------------------
void
AP_AdvectionDiffusion1CFV_ERROR:: n0( PDE_LocalFEbound const* fe,
                                   double int_v_nor )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_AdvectionDiffusion1CFV:" 
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
AP_AdvectionDiffusion1CFV_ERROR:: n1( PEL_ModuleExplorer const* exp,
                                   std::string const& keyword,
                                   size_t size )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** AP_AdvectionDiffusion1CFV:" 
       << std::endl << std::endl ;
   msg << "In module: " << exp->absolute_path_name() << endl ;
   msg << "the data of keyword: " << keyword << endl ;
   msg << "should be a DoubleVector of size: " << size ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
