#include <AP_FSINavierStokes.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>
#include <GE_Point.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_LocalBCsBuilder.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <AP_FSINavierStokesSystem.hh>

#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ; 

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_FSINavierStokes const* AP_FSINavierStokes::PROTOTYPE = new AP_FSINavierStokes() ;

//---------------------------------------------------------------------------
AP_FSINavierStokes:: AP_FSINavierStokes( void )
//--------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_FSINavierStokes" )
{
}

//---------------------------------------------------------------------------
AP_FSINavierStokes*
AP_FSINavierStokes:: create_replica( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSINavierStokes:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_FSINavierStokes* result = 
                        new AP_FSINavierStokes( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_FSINavierStokes:: AP_FSINavierStokes( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data("velocity") ) )
   , L_UU( exp->int_data( "level_of_explicit_velocity" ) )
   , L_UPDATE_UU( exp->int_data( "velocity_level_to_update" ) )
   , PP( dom->set_of_discrete_fields()->item( exp->string_data("pressure") ) )
   , L_PP( exp->int_data( "level_of_explicit_pressure" ) )
   , L_UPDATE_PP( exp->int_data( "pressure_level_to_update" ) )
   , ADV( false )
   , ORDER( exp->int_data( "time_order" ) )
   , ALPHA( prms->item( exp->string_data( "param_unsteady" ) ) )
   , MU( prms->item( exp->string_data( "param_viscous" ) ) )
   , RHSU( prms->item( exp->string_data( "param_source" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , LOCAL_BC( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , GLOBAL_EQ( 0 )
{
   PEL_LABEL( "AP_FSINavierStokes:: AP_FSINavierStokes" ) ;

   PEL_ASSERT( ORDER == 1 ) ; //????

   check_field_storage_depth( PP, L_UPDATE_PP ) ;
   if( ORDER == 1 )
      check_field_storage_depth( UU, PEL::max( L_UU, L_UPDATE_UU ) ) ;
   else if( ORDER == 2 )
      check_field_storage_depth( UU, PEL::max( L_UU+1, L_UPDATE_UU ) ) ;
   else
   {
      cout << "time_order = " << ORDER << endl ;
      PEL_Error::object()->raise_bad_data_value( exp, "time_order", 
                                                 "   1 or 2" );
   }

   check_param_nb_components( ALPHA, "param_unsteady", 1 ) ;
   check_param_nb_components( MU, "param_viscous", 1 ) ;
   check_param_nb_components( RHSU, "param_source", 
                              dom->nb_space_dimensions() ) ;

   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;   
   cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ; //???? pas tjrs
   bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;

   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

   ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   MU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   RHSU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ); 
   if( exp->string_data( "viscosity_term" ) == "mu_laplacian_uu" )
      LAPL_UU = true ;
   else if( exp->string_data( "viscosity_term" ) == "div_mu_D_uu" )
      LAPL_UU = false ;
   else
      PEL_Error::object()->raise_bad_data_value( exp, "viscosity_term", 
                                 "\"mu_laplacian_uu\" or \"div_mu_D_uu\"") ;

   if( exp->has_module( "advection" ) )
   {
      ADV = true ;      
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "advection" ) ;

      PEL_ModuleExplorer* se = e->create_subexplorer( e, "advective_field" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
         PDE_DiscreteField const* ff = dfs->item( sse->string_data( "field" ) ) ;
         AAs.push_back( ff ) ;
         size_t ll = sse->int_data( "level" ) ;
         L_AAs.push_back( ll ) ;
         check_field_storage_depth( ff, ll ) ;
         FE_Parameter* prm = prms->item( sse->string_data( "param_coef" ) ) ;
         prm->transfer_cell_calculation_requirements( cFE, 
                                                      FE_Parameter::Val ) ;
         COEF_AAs.push_back( prm ) ;
         sse->destroy() ;

         cFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
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
   
   PEL_ModuleExplorer* se = 
                        exp->create_subexplorer( 0, "AP_FSINavierStokesSystem" ) ;
   PDE_LinkDOF2Unknown* ulink = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
   PDE_LinkDOF2Unknown* plink = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
   GLOBAL_EQ = AP_FSINavierStokesSystem::create( this, se, ulink, plink) ;
   se->destroy() ;
}

//---------------------------------------------------------------------------
AP_FSINavierStokes:: ~AP_FSINavierStokes( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSINavierStokes:: ~AP_FSINavierStokes" ) ;
}

//---------------------------------------------------------------------------
void
AP_FSINavierStokes:: do_before_inner_iterations_stage( 
                                          FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSINavierStokes:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "AP_FSINavierStokes:: do_before_inner_iterations_stage" ) ;

   //??? GLOBAL_EQ->nullify_LHS_and_RHS() ;
   GLOBAL_EQ->re_initialize() ; //??? peut-etre trop ????

   size_t nbc = UU->nb_components() ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   PEL_ASSERT( nbc == nb_dims ) ;
   doubleVector coef( nbc ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double rt = ALPHA->cell_value_at_IP( t_it, cFE )/t_it->time_step() ;
         PEL_ASSERT( ORDER == 1 ) ;
         if( ORDER == 1 || t_it->iteration_number() == 1 ) 
         {
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               coef( d ) = cFE->value_at_IP( UU, L_UU, d ) * rt ;
            }
         }
         FE::add_row( ELEMENT_EQ, cFE, coef ) ;
      }

      GLOBAL_EQ->assemble_A_F_explicit( ELEMENT_EQ ) ;
   }

   stop_total_timer() ;

   PEL_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
AP_FSINavierStokes:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FSINavierStokes:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_FSINavierStokes:: do_one_inner_iteration" ) ;
   // ---------------------

//??????? faire tout en un ?????
   GLOBAL_EQ->nullify_for_new_internal_iteration() ;
   GLOBAL_EQ->copy_explicit_terms() ;

   double xx = 1.0/t_it->time_step() ;
   if( ORDER == 2 && t_it->iteration_number() != 1 ) xx *= 1.5 ;
   GLOBAL_EQ->set_leading_BDF_over_dt( xx ) ;

   start_assembling_timer() ;
   // ---------------------

   loop_on_cells( t_it ) ;
   loop_on_bounds( t_it ) ;

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------

   GLOBAL_EQ->set_indent( indent() ) ;
   GLOBAL_EQ->set_initial_guess_P( PP, L_PP ) ; 
   GLOBAL_EQ->set_initial_guess_U( UU, L_UPDATE_UU ) ; 
   GLOBAL_EQ->estimate_unknowns() ;
   stop_solving_timer() ;

   if( !GLOBAL_EQ->unknowns_are_solution() )
   {
      PEL_Error::object()->raise_plain(
         "*** AP_FSINavierStokes error:\n"
         "    solution failure" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      if( verbose_level() >= 2 )
      {
         cout << indent() << "   update of " << UU->name() 
              << "(" << L_UPDATE_UU << ") and " << PP->name() 
              << "(" << L_UPDATE_PP << ")" << endl ;
      }
      UU->update_free_DOFs_value( L_UPDATE_UU,
                                  GLOBAL_EQ->unknown_vector_U(), 
                                  GLOBAL_EQ->system_numbering_U()->link() ) ;
      
      PP->update_free_DOFs_value( L_UPDATE_PP,
                                  GLOBAL_EQ->unknown_vector_P(), 
                                  GLOBAL_EQ->system_numbering_P()->link() ) ;
   }

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
AP_FSINavierStokes:: loop_on_cells( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{ 
   PEL_LABEL( "AP_FSINavierStokes:: loop_on_cells" ) ;
   size_t nbc = UU->nb_components() ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   PEL_ASSERT( nbc == nb_dims ) ;
   doubleVector aa( nb_dims ) ;
   doubleVector coef( nbc ) ;
   doubleVector tdmu( nb_dims ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double mu = MU->cell_value_at_IP( t_it, cFE ) ;
        
         if( ADV )
         {
            aa.set( 0.0 ) ;
            for( size_t i=0 ; i<AAs.size() ; ++i )
            {
               double xx = COEF_AAs[i]->cell_value_at_IP( t_it, cFE ) ;
               for( size_t d=0 ; d<nb_dims ; ++d )
               {
                  aa( d ) += xx * cFE->value_at_IP( AAs[i], L_AAs[i], d ) ;
               }
            }

            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;
         }
         double rt = ALPHA->cell_value_at_IP( t_it, cFE )/t_it->time_step() ;
         double bq = rt ;
         if( ORDER == 1 || t_it->iteration_number() == 1 ) 
         {
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
//	       coef( d ) = RHSU->cell_value_at_IP( t_it, cFE, d )
//                         + cFE->value_at_IP( UU, L_UU, d ) * rt ;
               coef( d ) = RHSU->cell_value_at_IP( t_it, cFE, d ) ;
            }
         }
         else if( ORDER == 2 && t_it->iteration_number() != 1 )
         {
//?????????????? on ne doit pas passer ici
            PEL_ASSERT( false ) ;
            bq *= 1.5 ;
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               coef( d ) = RHSU->cell_value_at_IP( t_it, cFE, d )
                         + 2.0 * cFE->value_at_IP( UU, L_UU  , d ) * rt
                         - 0.5 * cFE->value_at_IP( UU, L_UU+1, d ) * rt ;
            }
         }
         FE::add_row( ELEMENT_EQ, cFE, coef ) ;
         FE::add_row_col_S( ELEMENT_EQ, cFE, bq ) ;
 
         if( LAPL_UU )
            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, mu ) ;
         else
            FE::add_grad_row_D_col_S( ELEMENT_EQ, cFE, mu ) ;

      }
      GLOBAL_EQ->assemble_A_F( ELEMENT_EQ ) ;

      cFE->set_row_and_col_fields( PP, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(),  1,
                              cFE->col_field_node_connectivity(),  nbc ) ;

      cFE->start_IP_iterator( QRP ) ;      
      for(  ; cFE->valid_IP() ; cFE->go_next_IP() ) 
      {
         FE::add_row_div_col( ELEMENT_EQ, cFE, -1.0 ) ;
      }
      GLOBAL_EQ->assemble_B_G( ELEMENT_EQ ) ;

      if( GLOBAL_EQ->MPl_is_required() )
      {
         cFE->set_row_and_col_fields( PP, PP ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
	   		         cFE->col_field_node_connectivity(), 1 ) ;

         size_t nb_nodes = cFE->nb_basis_functions( row ) ;
         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            for( size_t i=0 ; i<nb_nodes ; ++i )
            {  
               double xx =  cFE->weight_of_IP()* cFE->N_at_IP( row, i ) ;
               ELEMENT_EQ->add_to_matrix( xx, i, i ) ;     
            }
         }
         GLOBAL_EQ->assemble_MPl( ELEMENT_EQ ) ;
      }

//???????? lenteur avec assemble_A ????
      if( GLOBAL_EQ->MV_is_required() )
      {
         cFE->set_row_and_col_fields( UU, UU ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                                 cFE->col_field_node_connectivity(), nbc ) ;

         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
            FE::add_row_col_S( ELEMENT_EQ, cFE, alpha ) ;
         }
         GLOBAL_EQ->assemble_MV( ELEMENT_EQ ) ;
      }


//???????? lenteur avec assemble_MPl ????
      if( GLOBAL_EQ->L_is_required() )
      {
         cFE->set_row_and_col_fields( PP, PP ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                 cFE->col_field_node_connectivity(), 1 ) ;
         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0/alpha ) ;
         }
         GLOBAL_EQ->assemble_L( ELEMENT_EQ )  ;
      }     
   }
}

//---------------------------------------------------------------------------
void
AP_FSINavierStokes:: loop_on_bounds( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{ 
   PEL_LABEL( "AP_FSINavierStokes:: loop_on_bounds" ) ;

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

               GLOBAL_EQ->assemble_A_F( ELEMENT_EQ ) ;
            }
         }
      }
   }

//???????????????????????????????
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
               FE::add_row_col_S( ELEMENT_EQ, bFE, penal ) ;
            }
            GLOBAL_EQ->assemble_L( ELEMENT_EQ ) ;
         }  
      }
   }
}
