#include <MI_NavierStokes.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Double.hh>
#include <PEL_Exec.hh>
#include <PEL_Variable.hh>

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

#include <FE_PrintVariables.hh>
#include <PEL_DataWithContextExp.hh>
#include <MI_NavierStokesSystem.hh>

#include <iostream>
#include <sstream>
#include <cmath>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ; 

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

MI_NavierStokes const* MI_NavierStokes::PROTOTYPE = new MI_NavierStokes() ;

//---------------------------------------------------------------------------
MI_NavierStokes:: MI_NavierStokes( void )
//--------------------------------------------------------------------------
        : FE_OneStepIteration( "MI_NavierStokes" )
        , COM(0)
        , COORDS(0)
        , TT(0)
{
}

//---------------------------------------------------------------------------
MI_NavierStokes*
                MI_NavierStokes:: create_replica( PEL_Object* a_owner,
                PDE_DomainAndFields const* dom,
                FE_SetOfParameters const* prms,
                PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
        PEL_LABEL( "MI_NavierStokes:: create_replica" ) ;
        PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

        MI_NavierStokes* result = 
                        new MI_NavierStokes( a_owner, dom, prms, exp ) ;

        PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
        return( result ) ;
}

//---------------------------------------------------------------------------
MI_NavierStokes:: MI_NavierStokes( PEL_Object* a_owner,
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
        , RE( prms->item( exp->string_data( "param_Reynolds" ) ) )
        , RHSU( prms->item( exp->string_data( "param_source" ) ) )
        , VISC( prms->item( exp->string_data( "param_viscous" ) ) )
        , BCs( dom->set_of_boundary_conditions() )
        , BC_STRESS( false )
        , LOCAL_BC( 0 )
        , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
        , QRP( GE_QRprovider::object(
                   exp->string_data( "quadrature_rule_provider" ) ) )
        , cFE( dom->create_LocalFEcell( this ) )
        , bFE( dom->create_LocalFEbound( this ) )
        , GLOBAL_EQ( 0 )
        , Nonlin( false )
        , NonlinType( "" )
        , NonlinMaxiter( 1 )
        , NonlinResidual( 0. )
        , NonlinIter( 0 )
        , L_NonlinU(1)
        , Flowrate(false)
        , FlowrateCP(0)
        , FlowrateMaxiter(1)
        , COM( PEL_Exec::communicator() )
	    , CONTEXT( PEL_ContextSimple::create( this ) )
	    , COORDS( PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) )
            , TT(PEL_Double::create( CONTEXT, 0.0 ) )
{
    PEL_LABEL( "MI_NavierStokes:: MI_NavierStokes" ) ;

    if( exp->has_module( "FE_PrintVariables" ) ){
        PEL_ModuleExplorer* ee =
                         exp->create_subexplorer( 0, "FE_PrintVariables" ) ;
        FE_PrintVariables::create( this, ee) ;
        ee->destroy() ; ee=0;
    }

    PEL_ASSERT( ORDER == 1 ) ; //????

    check_field_storage_depth( PP, L_UPDATE_PP ) ;
    if ( ORDER == 1 )
        check_field_storage_depth( UU, PEL::max( L_UU+1, L_UPDATE_UU ) ) ;
    else if ( ORDER == 2 )
        check_field_storage_depth( UU, PEL::max( L_UU+2, L_UPDATE_UU ) ) ;
    else
    {
        cout << "time_order = " << ORDER << endl ;
        PEL_Error::object()->raise_bad_data_value( exp, "time_order",
                "   1 or 2" );
    }

    check_param_nb_components( ALPHA, "param_unsteady", 1 ) ;
    check_param_nb_components( RE,    "param_Reynolds", 1 ) ;
    check_param_nb_components( VISC,  "param_viscous", 1 ) ;
    check_param_nb_components( RHSU,  "param_source",
                               dom->nb_space_dimensions() ) ;    

    cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ; //???? pas tjrs
    bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;	

    PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

    ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    RE->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    VISC->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    RHSU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val );
        
   if( exp->string_data( "viscosity_term" ) == "mu_laplacian_uu" )
      LAPL_UU = true ;
   else if( exp->string_data( "viscosity_term" ) == "div_mu_D_uu" )
      LAPL_UU = false ;
   else
      PEL_Error::object()->raise_bad_data_value( exp, "viscosity_term",
                                 "\"mu_laplacian_uu\" or \"div_mu_D_uu\"") ;

    if ( exp->has_module( "advection" ) )
    {
        ADV = true ;
        PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "advection" ) ;

        PEL_ModuleExplorer* se = e->create_subexplorer( e, "advective_field" ) ;
        se->start_module_iterator() ;
        for ( ; se->is_valid_module() ; se->go_next_module() )
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

    PEL_ModuleExplorer* se = 
                        exp->create_subexplorer( 0, "MI_NavierStokesSystem" ) ;
    PDE_LinkDOF2Unknown* ulink = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
    PDE_LinkDOF2Unknown* plink = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
    GLOBAL_EQ = MI_NavierStokesSystem::create( this, se, ulink, plink) ;
    se->destroy();

    if ( exp->has_entry( "boundary_conditions_types" ) )
    {
        // example:  boundary_conditions_types="FE_NormalVelocityBC"   
        stringVector const& bt =
            exp->stringVector_data( "boundary_conditions_types" ) ;   
        string const& fname = exp->string_data( "velocity" ) ;
        LOCAL_BC = FE_LocalBCsBuilder::create( this, dom, fname, bt, prms ) ;
        LOCAL_BC->transfer_calculation_requirements( bFE ) ;
    }
    if ( exp->has_entry( "param_bc_stress" ) )
    {
        BC_STRESS_VALUE = prms->item( exp->string_data( "param_bc_stress" ) ) ; 
        check_param_nb_components( BC_STRESS_VALUE,  "param_bc_stress", 1 ) ; 
        BC_STRESS_VALUE->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
        BC_STRESS = true ;
    }
    if ( exp->has_module( "NonlinearTerm" ) )
    {
        se = exp->create_subexplorer( 0, "NonlinearTerm" ) ;
        Nonlin = true;
        NonlinType 		= se->string_data("type");
        NonlinMaxiter 	= se->int_data("maxiter");
        NonlinResidual 	= se->double_data("residual");
        L_NonlinU 	= se->int_data("storage_level");
        se->destroy();
    }    

	if ( exp->has_module( "flowrate" ) )
	{
		se = exp->create_subexplorer( 0, "flowrate" ) ;
		Flowrate = true;
		FlowrateResidual = se->double_data("residual");
		FlowrateTarget =   se->double_data("target");
		FlowrateMaxiter	=  se->int_data("maxiter");
		FlowrateCP      =  se->int_data("component");
		FlowratedPdL.push_back(se->double_data("dPdL0"));
		FlowratedPdL.push_back(se->double_data("dPdL1"));
		FlowrateFr.push_back(se->double_data("Fl0"));
		FlowrateMethod  =  se->string_data("method");
		FlowrateScale   =  se->double_data("scale");
		se->destroy();
	}

	// Context :
	CONTEXT->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
        CONTEXT->extend( PEL_Variable::object( "DS_T" ), TT ) ;
}

//---------------------------------------------------------------------------
MI_NavierStokes:: ~MI_NavierStokes( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_NavierStokes:: ~MI_NavierStokes" ) ;
}

//---------------------------------------------------------------------------
void
MI_NavierStokes:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_NavierStokes:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
        
   start_total_timer( "MI_NavierStokes:: do_one_inner_iteration" ) ;
   // ---------------------

   GLOBAL_EQ->nullify_for_new_internal_iteration() ;
   GLOBAL_EQ->add_explicit_terms() ; // A_explicit and F_explicit
   GLOBAL_EQ->add_Outer_terms() ;    // A_Outer and F_Outer
      
   setup_Inner( t_it ); // contains GLOBAL_EQ->add_Inner_terms() and  
                        // GLOBAL_EQ->add_Flowrate_terms() ;  

// Loop for eventual flowrate steps
    size_t flowrateIter  = 0;
	double flowrateError = 0;
    do
	{
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
             "*** MI_NavierStokes error:\n"
             "    solution failure" ) ;
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
                                  GLOBAL_EQ->system_numbering_U()->link() ) ;
        
            PP->update_free_DOFs_value( L_UPDATE_PP,
                                  GLOBAL_EQ->unknown_vector_P(), 
                                  GLOBAL_EQ->system_numbering_P()->link() ) ;
        }
		
        if ( Flowrate && flowrateIter < FlowrateMaxiter )
        {
		    FlowrateFr.push_back( compute_flowrate( UU, L_UPDATE_UU ) ); //.add updated flowrate
	    	increase_indent();
		    flowrateError = fabs( *(FlowrateFr.end()-1) - FlowrateTarget );
			
	    	PEL::out().setf(std::ios::scientific);
	    	PEL::out() 
	    			<< indent() << "Current Step: dPdL = " << *( FlowratedPdL.end() - 1 )  
				    << ", Rate = " << *( FlowrateFr.end() -1 )
					<< ", Error = " << flowrateError << endl;
		
	    	if ( flowrateError > FlowrateResidual )
	    	{
	    		update_pressuredrop();
	    		PEL::out() << indent() << "New Step:    : dPdL = " << *( FlowratedPdL.end() -1 ) << endl;
	    	}
	    	else
	    	{
	    		FlowrateFr.pop_back();
	    	}
	    		
	    	decrease_indent();
	    	flowrateIter++;
	    }
	    	
    } while( Flowrate && flowrateIter < FlowrateMaxiter && flowrateError > FlowrateResidual);
	
	NonlinIter++;
	// ------------------------------------------------------------------------- 
       
    stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
MI_NavierStokes:: update_pressuredrop(  )
//---------------------------------------------------------------------------
{
	PEL_CHECK_PRE(FlowrateFr.size() >= 2);
	PEL_CHECK_PRE(FlowratedPdL.size() >= 2);
	PEL_CHECK_PRE(FlowratedPdL.size() == FlowrateFr.size() );
	
	size_t N = FlowratedPdL.size()-1;
	double Fn = FlowrateFr[N] - FlowrateTarget;
	double Fnp1 = FlowrateFr[N-1] - FlowrateTarget;
	double xn = FlowratedPdL[N];
	double xnp1 = FlowratedPdL[N-1];
	
	if (Fn != Fnp1)
	{
		double deltax = ( xn - xnp1 ) / ( Fn - Fnp1 ) * Fn;
		FlowratedPdL.push_back( xn - deltax);
	}
	else if (Fn < 0.)
	{
		// do one regula falsi step
		FlowratedPdL.push_back( xn * ( 1.+5.*FlowrateResidual ) );
	}
	else if (Fn > 0. )
	{
		FlowratedPdL.push_back( xn * ( 1.-5.*FlowrateResidual ) );
	}
	    
    GLOBAL_EQ->nullify_for_new_internal_iteration() ;
    GLOBAL_EQ->add_explicit_terms() ;
    GLOBAL_EQ->add_Outer_terms();
    GLOBAL_EQ->add_Inner_terms();
    GLOBAL_EQ->add_Flowrate_terms( *(FlowratedPdL.end()-1) );
}

//------------------------------------------------------------------------------
/**
 * @brief Wrapper to setup the outer system.
 * 
 * @param t_it Time iterator
 * 
 * Functions:
 * - Call the appropriate loops
 * 
 * @todo
 * - Check whether the linkDOF has to be reset (after e.g. a grid modification)
 * - What does set BDF do ?
 */
void
MI_NavierStokes:: setup_Outer( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	
	PEL_LABEL( "MI_NavierStokes:: setup_Outer" ) ;
	
	// ---------------------
	// Reset the links:
	// UU_link->reset() ;
    // PP_link->reset() ;

	// Set BDF ???
	double xx = 1.0/t_it->time_step() ;
	if( ORDER == 2 && t_it->iteration_number() != 1 ) xx *= 1.5 ;
	GLOBAL_EQ->set_leading_BDF_over_dt( xx ) ;

	// Assemble Discrete equation
	start_assembling_timer();
	
	// Remove old entries in the matrix:
	GLOBAL_EQ->re_initialize();  
	
	// Call the loops
	loop_on_cells_Outer( t_it );
	loop_on_bounds_Outer( t_it );
	
	// Do I need to do something else ?
	stop_assembling_timer();
}

//---------------------------------------------------------------------------
/**
 * @brief Wrapper for the setup of the inner iterations
 * 
 * @param t_it Time iterator
 * 
 * - Call loops
 * @todo
 * - Re-initialize -> Don't forget the initialguess ?
 * - remove contributions form the BlockAssembledSystem
 */
void
MI_NavierStokes:: setup_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: setup_Inner" ) ;
	
	double inner_rhs_scale = 0.;
	
	start_assembling_timer();

	// Remove old entries in the matrix:
    GLOBAL_EQ->nullify_A_F_Inner() ;

	// Setup coeff of the BDF
	double xx = 1.0/t_it->time_step() ;
	if ( ORDER == 2 && t_it->iteration_number() != 1 ) xx *= 1.5 ;
	GLOBAL_EQ->set_leading_BDF_over_dt( xx ) ;
	
	// Copy the value of the n+1 velocity to the nonlinear location level
	if ( Nonlin )
	{
		PEL_CHECK( UU->storage_depth() >= L_NonlinU );
		UU->copy_DOFs_value( L_UPDATE_UU, L_NonlinU );
	}
	
	if ( Flowrate )
		inner_rhs_scale = *(FlowratedPdL.end()-1);
	
	// Call the loops
	loop_on_cells_Inner( t_it );
	loop_on_bounds_Inner( t_it );
	

	// Assemble the main matrix system
    GLOBAL_EQ->add_Inner_terms() ;
	if ( Flowrate )
	    GLOBAL_EQ->add_Flowrate_terms( inner_rhs_scale ) ; // zero for (!Flowrate)
	
    GLOBAL_EQ->add_Stress_terms() ; 
	// Do I need to do something else ?
	stop_assembling_timer();
}

//---------------------------------------------------------------------------
/**
 * @brief Wrapper for the setup of the flowrate problem
 * 
 * @param t_it Time iterator
 * 
 */
void
MI_NavierStokes:: setup_Flowrate( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: setup_Flowrate" ) ;
	PEL_CHECK_PRE( Flowrate );
// NOT USED !!!	
	start_assembling_timer();
	
	double inner_rhs_scale = *(FlowrateFr.end()-1);
	
	// Assemble the main martix system
	//KWB GLOBAL_EQ->assemble_A_F( inner_rhs_scale ) ;
	GLOBAL_EQ->add_Flowrate_terms( inner_rhs_scale ) ;

	// Do I need to do something else ?
	stop_assembling_timer();
}

//---------------------------------------------------------------------------
/**
 * @brief Assemble finite element contributions on all cells for the outer iteration
 * @param t_it Time iterator
 * 
 * - System Contributions:
 *   - Mass form : \f$ \int_\Omega  Re \frac{\beta_q}{\delta t} u^{n+1} . v dV \f$
 *   - Advection term: \f$ \alpha \int_\Omega Re ( (u^{\star}\cdot \nabla) u^{n+1} \cdot v ) dv \f$
 *   - pressure-velocity: \f$ \int_\Omega p . (\nabla \cdot v) dV \f$
 * 
 * - RHS 
 *   - time discrete rhs:
 * 
 * - Algorithm dependent contributions:
 *   - Condensed pressure mass matrix MPL
 *   - L
 *   - MV
 * 
 * @todo
 * - Change advection term to advective parameter class
 */
void
MI_NavierStokes:: loop_on_cells_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
        PEL_LABEL( "MI_NavierStokes:: loop_on_cells_Outer" ) ;
        size_t nbc = UU->nb_components() ;
        size_t nb_dims = cFE->nb_space_dimensions() ;
        PEL_ASSERT( nbc == nb_dims ) ;

        doubleVector aa( nb_dims ) ;
        doubleVector coef( nbc ) ;
        doubleVector tdmu( nb_dims ) ;
        
	    double dt = t_it->time_step() ;
	    double Re = 1.;
	    double Alpha = 0.;

        for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
        {
		    // -------------------------------------------------------------------
		    // setup u dot v term
		    // - If flowrate problem is set the flowrate term will be added to
		    //   a different contribution in the rhs
		    // -------------------------------------------------------------------
            cFE->set_row_and_col_fields( UU, UU ) ;
            ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                                    cFE->col_field_node_connectivity(), nbc ) ;

            cFE->start_IP_iterator( QRP ) ;
            for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
            {
			    Re   = RE->cell_value_at_IP( t_it, cFE ) ;
			    Alpha= ALPHA->cell_value_at_IP( t_it, cFE ) ;
			    // Setup advection term
                if ( ADV  )
                {
                    aa.set( 0.0 ) ;
                    for ( size_t i=0 ; i<AAs.size() ; ++i )
                    {
                         double xx = COEF_AAs[i]->cell_value_at_IP( t_it, cFE ) ;
                         for ( size_t d=0 ; d<nb_dims ; ++d )
                         {
                             aa( d ) += xx * cFE->value_at_IP( AAs[i], L_AAs[i], d ) ;
                         }
                    }                        
                    FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, Re ) ;
                }

			    // Setup mass form associated with time and the rhs term
			    double rt = Alpha*Re/dt ;
			    double bq = rt ;
                if ( ORDER == 1 || t_it->iteration_number() == 1 )
                {
                    for ( size_t d=0 ; d<nb_dims ; ++d )
                    {
					    coef( d ) = (!Flowrate) * RHSU->cell_value_at_IP( t_it, cFE, d );
                    }
                }
                else if ( ORDER == 2 && t_it->iteration_number() != 1 )
                {
                    PEL_ASSERT( false ) ; // should not happen here!
                    bq *= 1.5 ;
                    for ( size_t d=0 ; d<nb_dims ; ++d )
                    {
                        coef( d ) = (!Flowrate) * RHSU->cell_value_at_IP( t_it, cFE, d )
							+ 2.0 * cFE->value_at_IP( UU, L_UU+1, d ) * rt
							- 0.5 * cFE->value_at_IP( UU, L_UU+2, d ) * rt ;
                    }
                }
                FE::add_row( ELEMENT_EQ, cFE, coef ) ;
                FE::add_row_col_S( ELEMENT_EQ, cFE, bq ) ;
             }
            GLOBAL_EQ->assemble_A_F_Outer( ELEMENT_EQ ) ; 

		// --------------------------------------------------------------------
		// Setup potential flowrate problem
		// --------------------------------------------------------------------
		    if ( Flowrate )
		    {
			ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
									cFE->col_field_node_connectivity(), nbc ) ;
	
			cFE->start_IP_iterator( QRP ) ;
			for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
			{	
				// Setup pressuredrop
				for ( size_t d=0 ; d<nb_dims ; ++d )
				{
					coef( d ) =  RHSU->cell_value_at_IP( t_it, cFE, d );
				}
				
				FE::add_row( ELEMENT_EQ, cFE, coef ) ;
			}
            GLOBAL_EQ->assemble_F_Flowrate( ELEMENT_EQ ) ; 
		    }
               
            // ---------------------------------------------------------------------
            // Velocity - pressure interactions
            // ---------------------------------------------------------------------
            cFE->set_row_and_col_fields( PP, UU ) ;
            ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(),  1,
                                    cFE->col_field_node_connectivity(),  nbc ) ;

            cFE->start_IP_iterator( QRP ) ;
            for (  ; cFE->valid_IP() ; cFE->go_next_IP() )
            {
                    FE::add_row_div_col( ELEMENT_EQ, cFE, -1.0 ) ;
            }
            GLOBAL_EQ->assemble_B_G( ELEMENT_EQ ) ;

            // ---------------------------------------------------------------------
            // MPL condensed pressure-pressure mass matrix
            // Required by: Augmented Lagrangian method OR Augmentation Projection
            // ---------------------------------------------------------------------
                if ( GLOBAL_EQ->MPl_is_required() )
                {
                        cFE->set_row_and_col_fields( PP, PP ) ;
                        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                                                        cFE->col_field_node_connectivity(), 1 ) ;

                        size_t nb_nodes = cFE->nb_basis_functions( row ) ;
                        cFE->start_IP_iterator( QRP ) ;
                        for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
                        {
				                double c_w = cFE->weight_of_IP() ;
  				                if ( FE::geometry() == FE::axisymmetrical )
				                {
                					c_w *= cFE->coordinates_of_IP()->coordinate(0) ;
				                }
                                for ( size_t i=0 ; i<nb_nodes ; ++i )
                                {
				                	double xx = c_w* cFE->N_at_IP( row, i ) ;
                                    ELEMENT_EQ->add_to_matrix( xx, i, i ) ;
                                }
                        }
                        GLOBAL_EQ->assemble_MPl( ELEMENT_EQ ) ;
                }

                // ---------------------------------------------------------------
                // MV 
                // Required by: Augmentation Projection
                // slow with assemble_A 
                // ---------------------------------------------------------------
                if ( GLOBAL_EQ->MV_is_required() )
                {
                        cFE->set_row_and_col_fields( UU, UU ) ;
                        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                                                                        cFE->col_field_node_connectivity(), nbc ) ;

                        cFE->start_IP_iterator( QRP ) ;
                        for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
                        {
                                double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
                                FE::add_row_col_S( ELEMENT_EQ, cFE, alpha ) ;
                        }
                        GLOBAL_EQ->assemble_MV( ELEMENT_EQ ) ;
                }

                // ------------------------------------------------------------
                // Matrix L: gradP gradP + alpha P P
                // Required by : Augmentation Projection, Yosida Projection 
                // slow with assemble_MPl 
                // ------------------------------------------------------------
                if ( GLOBAL_EQ->L_is_required() )
                {
                        cFE->set_row_and_col_fields( PP, PP ) ;
                        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                                cFE->col_field_node_connectivity(), 1 ) ;
                        cFE->start_IP_iterator( QRP ) ;
                        for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
                        {
                                double alpha = ALPHA->cell_value_at_IP( t_it, cFE ) ;
                                FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0/alpha ) ;
                        }
                        GLOBAL_EQ->assemble_L( ELEMENT_EQ )  ;
                }
        }
}

//---------------------------------------------------------------------------
/**
 * @brief Assemble finite element contributions on all cells in the inner iteration.
 * Add all contributions together.
 * 
 * @param t_it Time iterator
 * 
 * - System Contributions:
 *   - Viscosity form: \f$ \int_\Omega \mu (\nabla u + \nabla u^T):(\nabla v) dV \f$
 */
void
MI_NavierStokes:: loop_on_cells_Inner( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: loop_on_cells_Inner" ) ;
	size_t nbc = UU->nb_components() ;
	size_t nb_dims = cFE->nb_space_dimensions() ;
	PEL_ASSERT( nbc == nb_dims ) ;
	
	doubleVector aa( nb_dims ) ;
	doubleVector coef( nbc ) ;
	doubleVector tdmu( nb_dims ) ;
	
	double Visc = 1.;
	
	bool FirstIteration = false;
	
	if( NonlinIter == 0  )
	{
		if( t_it->is_started() )
		{
			if (t_it->iteration_number() <= 1)
				FirstIteration=true;
			else
				FirstIteration=false;
		}
	}

	for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
	{
		// setup u dot v term
		cFE->set_row_and_col_fields( UU, UU ) ;
		ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
								cFE->col_field_node_connectivity(), nbc ) ;

		cFE->start_IP_iterator( QRP ) ;
		for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
		{
			Visc = 1;
			if ( !FirstIteration )
				Visc = VISC->cell_value_at_IP( t_it, cFE ) ;

				FE::add_grad_row_D_col_S( ELEMENT_EQ, cFE, Visc ) ;
		}
		GLOBAL_EQ->assemble_A_F_Inner( ELEMENT_EQ ) ;
	}
}


//---------------------------------------------------------------------------
/**
 * @brief Setup finite element contributions on all bounds at the outer iterations stage.
 * @param t_it Time iterator
 * 
 * Assembled Forms:
 * - ???
 */
void
MI_NavierStokes:: loop_on_bounds_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_NavierStokes:: loop_on_bounds" ) ;

    size_t nbc = UU->nb_components() ;

    if ( LOCAL_BC != 0 ) // if boundary_conditions_types has a value
    {
        for ( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
        {
            GE_Color const* color = bFE->color() ;
            if ( BCs->has_BC( color, UU ) )
            {
                PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;
                LOCAL_BC->set_current_BC_type( ee->string_data( "type" ) ) ;
                if ( LOCAL_BC->current_BC_type_is_ok() )
                {
                    bFE->set_row_and_col_fields( UU, UU ) ;
                    ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(),
                                            nbc,
                                            bFE->col_field_node_connectivity(),
                                            nbc ) ;

                    LOCAL_BC->build_current_BC( ELEMENT_EQ, bFE, t_it, QRP ) ;

                    GLOBAL_EQ->assemble_A_F_Outer( ELEMENT_EQ ) ;
                }
                else	
			       PEL_Error::object()->display_info(
			    	    "*** MI_NavierStokes error:\n"
			           	"    boundary_conditions_types" );
            }
        }
    }

    // Boundary conditions
    for ( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
    {
	    bFE->set_calculation_point( bFE->polyhedron()->center() ) ;
        doubleVector value(0);
        if ( BCs->has_BC( bFE->color(), PP ) )
        {
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( bFE->color(), PP ) ;
            std::string const& bc_type = ee->string_data( "type" ) ;
            if ( bc_type == "Dirichlet_pressure_penalization" )
            {
                double penal = ee->double_data( "penalization_coefficient" ) ;
                bFE->set_row_and_col_fields( PP, PP ) ;
                ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 1,
                                        bFE->col_field_node_connectivity(), 1 ) ;

                bFE->start_IP_iterator( QRP ) ;
                for ( ; bFE->valid_IP() ; bFE->go_next_IP() )
                {
                    FE::add_row_col_S( ELEMENT_EQ, bFE, penal ) ;
                }
                GLOBAL_EQ->assemble_L( ELEMENT_EQ ) ;
            }
		else 
			raise_bad_BC_type( bc_type, "\"Dirichlet_pressure_penalization\"", PP);
        }
        if ( BCs->has_BC( bFE->color(), UU ) )
	    {           
            PEL_ModuleExplorer const* ee = BCs->BC_explorer( bFE->color(), UU ) ;
            std::string const& bc_type = ee->string_data( "type" ) ;		
	        if ( bc_type == "Stress" )
	        {
		        bFE->set_row_and_col_fields( UU, UU ) ; // (row,col) need only row
    		    ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), nbc,
					    bFE->col_field_node_connectivity(), nbc ) ;              

    		    PEL_DataWithContext const* value_bf 
                               = ee->abstract_data( 0, "value", CONTEXT ) ;
                                                            
                COORDS->set( bFE->calculation_point()->coordinate_vector() ) ;
                            
		        value = value_bf->to_double_vector() ;
		        value_bf->destroy() ;
		        bFE->start_IP_iterator( QRP ) ;
                
                if (BC_STRESS)
                {
                    double phiCC ;
                    doubleVector stress_value(UU->nb_components()) ;
    		        for( ; bFE->valid_IP() ; bFE->go_next_IP() )
	   	            {                                         
                        phiCC = BC_STRESS_VALUE->bound_value_at_IP( t_it, bFE ) ;
                        for (size_t i = 0 ; i< UU->nb_components() ; i++ )
                            stress_value(i) = phiCC*value(i) ;

                        FE::add_row( ELEMENT_EQ, bFE, stress_value ) ;
                    }
                }
                else 
                {
                    for( ; bFE->valid_IP() ; bFE->go_next_IP() )
                    {
			             FE::add_row( ELEMENT_EQ, bFE, value ) ;
                    }
		        }  
              
                GLOBAL_EQ->assemble_F_Stress( ELEMENT_EQ ) ; 
	        }
//	    Has to be commented out!! otherwise FE_BCupdate can't be used!
//	    else	
//		  raise_bad_BC_type ( bc_type, "\"Stress\"", UU);
	    }
    }
}

//---------------------------------------------------------------------------
/**
 * @brief Setup finite element contributions on all bounds at the inner 
 * iteration stage.
 * 
 * @param t_it Time iterator
 * 
 * Assembled Forms:
 * - None
 */
void 
MI_NavierStokes:: loop_on_bounds_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: loop_on_bounds_Inner" ) ;
}

//---------------------------------------------------------------------------
void
MI_NavierStokes:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: do_before_time_stepping" ) ;

	start_total_timer( "MI_NavierStokes:: do_before_time_stepping" ) ;

	FE_OneStepIteration::do_before_time_stepping( t_it ) ;

    GLOBAL_EQ->re_initialize() ;
    stop_total_timer() ;
}

//---------------------------------------------------------------------------
/**
 * @brief Do things before the inner iteration starts
 *
 * Jobs:
 * - Assemble the contributions contsant during the inner iteration
 * @param t_it
 */
void
MI_NavierStokes:: do_before_inner_iterations_stage( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_NavierStokes:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "MI_NavierStokes:: do_before_inner_iterations_stage" ) ;

   TT->set( t_it->time() ) ;

   NonlinIter = 0;
   //GLOBAL_EQ->re_initialize() ; // done within setup_Outer

   setup_Outer(t_it);   

   size_t nbc = UU->nb_components() ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   PEL_ASSERT( nbc == nb_dims ) ;
   doubleVector coef( nbc ) ;

   double Re = 1.0;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
        {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;
                
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
                {
         Re   = RE->cell_value_at_IP( t_it, cFE ) ;   
         double rt = Re*ALPHA->cell_value_at_IP( t_it, cFE )/t_it->time_step() ;
         PEL_ASSERT( ORDER == 1 ) ;
         if( ORDER == 1 || t_it->iteration_number() == 1 ) 
                {
                for( size_t d=0 ; d<nb_dims ; ++d )
                        {
           coef( d ) = cFE->value_at_IP( UU, L_UU+1, d ) * rt ;
           //     coef( d ) = cFE->value_at_IP( UU, L_UU, d ) * rt ;
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
/**
 * @brief Inner iterations stage:
 * @param t_it
 */
//----------------------------------------------------------------------------
void
MI_NavierStokes::do_inner_iterations_stage( FE_TimeIterator const* t_it )
{
	PEL_LABEL( "MI_NavierStokes:: do_inner_iterations_stage" ) ;
	do
    	do_one_inner_iteration( t_it ) ;
	while ( !inner_iterations_are_completed( t_it ) ) ;
}



//---------------------------------------------------------------------------
/**
 * @brief Convergence of the inner iteration stage
 *
 * Convergence of solution or maxiter reached
 * @todo
 *   - Implement norm comparator
 * @param t_it
 */
bool
MI_NavierStokes:: inner_iterations_are_completed( FE_TimeIterator const* t_it ) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: inner_iterations_are_completed" ) ;
	if ( Nonlin )
	{
		if ( NonlinIter > 1 )
			decrease_indent();
		double error = compute_L2_norm( UU, L_UU , L_NonlinU );
		PEL::out().setf(std::ios::scientific); 
		PEL::out() << indent() << "- nonlinear residual = " << error 
				<< ", iteration " << NonlinIter << std::endl;
		
		if ( NonlinIter > NonlinMaxiter  )
		{
			return true;
		}
		else if ( error < NonlinResidual )
		{
			return true;
		}
		else
		{
			increase_indent();
			return false;
		}
	}
	else
		return true;
}


//---------------------------------------------------------------------------
/**
 * @brief L2 norm for the velocity
 * @param FIELD 
 * @param level1 
 * @param level2 
 * @return L2 norm
 */
double
MI_NavierStokes::compute_L2_norm( PDE_DiscreteField const * FIELD,
								  size_t level1, size_t level2) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: compute_L2_square_norm" ) ;
	PEL_CHECK_INV( invariant() ) ;
	PEL_CHECK_PRE(FIELD->storage_depth() > level1);
	PEL_CHECK_PRE(FIELD->storage_depth() > level2);

	double error_norm = 0. ;

	size_t const nbc = FIELD->nb_components() ;


	for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
	{
		cFE->start_IP_iterator( QRP ) ;
		for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
		{
			doubleVector const& dv =
						cFE->coordinates_of_IP()->coordinate_vector() ;
			COORDS->set( dv ) ;
			double weight = cFE->weight_of_IP() ;
			if ( FE::geometry() == FE::axisymmetrical )
			{
				weight *= dv(0) ;
			}
			for ( size_t ic=0 ; ic<nbc ; ++ic   )
			{
				double const val1 =
							cFE->value_at_IP( FIELD, level1, ic );
				double const val2 =
							cFE->value_at_IP( FIELD, level2, ic );
				double const err = val1-val2;
				error_norm += weight*err*err ;
			}
		}
	}
	error_norm = COM->sum( error_norm ) ;
	return sqrt(error_norm);
}


//---------------------------------------------------------------------------
/**
 * @brief Calculate the flowrate.
 * @param FIELD 	link to velocity field
 * @param level1 	level of velocity field
 * @return value of the flowrate
 */
double
MI_NavierStokes::compute_flowrate( PDE_DiscreteField const * FIELD,
									size_t level1) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_NavierStokes:: compute_L2_square_norm" ) ;
	PEL_CHECK_INV( invariant() ) ;
	PEL_CHECK_PRE(FIELD->storage_depth() > level1);
	PEL_CHECK_PRE(FIELD->nb_components() > FlowrateCP);

	double flowrate = 0. ;

	for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
	{
		cFE->start_IP_iterator( QRP ) ;
		for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
		{
			doubleVector const& dv =
						cFE->coordinates_of_IP()->coordinate_vector() ;
			COORDS->set( dv ) ;
			double weight = cFE->weight_of_IP() ;
			if ( FE::geometry() == FE::axisymmetrical )
			{
				weight *= dv(0) ;
			}
			
			double const val1 =
						cFE->value_at_IP( FIELD, level1, FlowrateCP );
			flowrate += weight*val1 ;
			
		}
	}
	flowrate = COM->sum( flowrate  ) ;
	flowrate = flowrate * FlowrateScale;

	return flowrate;
}
