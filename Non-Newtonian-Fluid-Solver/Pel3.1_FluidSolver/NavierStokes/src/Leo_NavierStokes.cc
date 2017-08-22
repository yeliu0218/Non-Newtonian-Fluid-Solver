#include <Leo_NavierStokes.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
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
#include <FE_ScaledParameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <FE_PrintVariables.hh>
#include <FE_TensorFormAssembling.hh>
#include <PEL_DataWithContextExp.hh>
#include <ML_NavierStokesSystem.hh>
#include <Leo_Viscoplastic.hh>

#include <iostream>
#include <sstream>
#include <cmath>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

Leo_NavierStokes const* Leo_NavierStokes::PROTOTYPE = new Leo_NavierStokes() ;

//---------------------------------------------------------------------------
Leo_NavierStokes:: Leo_NavierStokes( void )
//--------------------------------------------------------------------------
        : FE_OneStepIteration( "Leo_NavierStokes" )
        , COM(0)
        , COORDS(0)
{
}

//---------------------------------------------------------------------------
Leo_NavierStokes*
                Leo_NavierStokes:: create_replica( PEL_Object* a_owner,
                PDE_DomainAndFields const* dom,
                FE_SetOfParameters const* prms,
                PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
        PEL_LABEL( "Leo_NavierStokes:: create_replica" ) ;
        PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

        Leo_NavierStokes* result =
                        new Leo_NavierStokes( a_owner, dom, prms, exp ) ;

        PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
        return( result ) ;
}

//---------------------------------------------------------------------------
Leo_NavierStokes:: Leo_NavierStokes( PEL_Object* a_owner,
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
        , STRESS( dom->set_of_discrete_fields()->item( exp->string_data("stress") ) )
        , GAMMADOT( dom->set_of_discrete_fields()->item( exp->string_data("gammadot") ) )
        , GAMMA( dom->set_of_discrete_fields()->item( exp->string_data("gamma") ) )
        , ADV( false )
        , ORDER( exp->int_data( "time_order" ) )
        , ALPHA( prms->item( exp->string_data( "param_unsteady" ) ) )
        , RE( prms->item( exp->string_data( "param_Reynolds" ) ) )
        , RHSU( prms->item( exp->string_data( "param_source" ) ) )
        , KAPPA( prms->item( exp->string_data( "param_consistency" ) ) )
        , BN( prms->item( exp->string_data( "param_Bingham" ) ) )
        , BCs( dom->set_of_boundary_conditions() )
        , BC_STRESS( false )
        , LOCAL_BC( 0 )
        , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
        , QRP( GE_QRprovider::object(
                   exp->string_data( "quadrature_rule_provider" ) ) )
        , cFE( dom->create_LocalFEcell( this ) )
        , bFE( dom->create_LocalFEbound( this ) )
        , GLOBAL_EQ( 0 )
        , Viscoplastic( true )
        , ViscoplasticType( "Bingham" )
        , ViscoplasticMaxiter( 1 )
        , ViscoplasticResidual( 0. )
        , ViscoplasticAug( 1. )
        , ViscoplasticIter( 0 )
        , FACTORIZE( true )
        , L_ViscoplasticU(1)
        , Flowrate(false)
        , FlowrateCP(0)
        , FlowrateMaxiter(1)
        , COM( PEL_Exec::communicator() )
	    , CONTEXT( PEL_ContextSimple::create( this ) )
	    , COORDS( PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) )
        , ElectricField(false)
        , FUNC_E(1.0)
{
    PEL_LABEL( "Leo_NavierStokes:: Leo_NavierStokes" ) ;

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
    check_param_nb_components( KAPPA, "param_consistency", 1 ) ;
    check_param_nb_components( BN,    "param_Bingham", 1 ) ;
    check_param_nb_components( RHSU,  "param_source",
                               dom->nb_space_dimensions() ) ;

    cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ; //???? pas tjrs
    bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;

    cFE->require_field_calculation( STRESS, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( GAMMADOT, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( GAMMA, PDE_LocalFE::N ) ;

    PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

    ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    RE->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    KAPPA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    BN->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
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
                        exp->create_subexplorer( 0, "ML_NavierStokesSystem" ) ;
    PDE_LinkDOF2Unknown* ulink = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
    PDE_LinkDOF2Unknown* plink = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
    PDE_LinkDOF2Unknown* stresslink = PDE_LinkDOF2Unknown::create( 0, STRESS,
                                          "sequence_of_the_components", true ) ;
    PDE_LinkDOF2Unknown* gammalink = PDE_LinkDOF2Unknown::create( 0, GAMMA,
                                          "sequence_of_the_components", true ) ;
//    PDE_LinkDOF2Unknown* gammadotlink = PDE_LinkDOF2Unknown::create( 0, GAMMADOT,
//                                          "sequence_of_the_components", true ) ;
//    GLOBAL_EQ = ML_NavierStokesSystem::create( this, se, ulink, plink, stresslink, gammalink, gammadotlink) ;

	GLOBAL_EQ = ML_NavierStokesSystem::create( this, se, ulink, plink, stresslink, gammalink) ;

    se->destroy();

    if ( exp->has_entry( "param_electric" ) )
    {
        ElectricField = true ;
        E_FIELD = prms->item( exp->string_data( "param_electric" ) ) ;
        check_param_nb_components( E_FIELD,"param_electric",
                                   dom->nb_space_dimensions() ) ;

        E_FIELD->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val );
    }

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
    if ( exp->has_module( "ViscoplasticTerm" ) )
    {
        se = exp->create_subexplorer( 0, "ViscoplasticTerm" ) ;
        Viscoplastic = true;
        ViscoplasticType        = se->string_data("type");
        ViscoplasticMaxiter     = se->int_data("maxiter");
        ViscoplasticResidual    = se->double_data("residual");
        L_ViscoplasticU     = se->int_data("storage_level");
        ViscoplasticAug     = se->double_data("augmentation");

        if (ViscoplasticType == "HerschelBulkley")
        {
        	HB_n = prms->item( exp->string_data( "param_Herschel_coef" ) ) ;
			check_param_nb_components( HB_n,    "param_Herschel_coef", 1 ) ;
			HB_n->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
        }
        if ( se->has_entry("factorize") ) FACTORIZE= se->bool_data("factorize");
        se->destroy();
    }
    else
    {
        PEL::out()
        << indent() << "The following module is missing:" << endl
        << indent() << "MODULE ViscoplasticTerm" << endl
        << indent() << "\ttype          = \"Bingham\" // Supported constitutive laws: Bingham, HerschelBulkley" << endl
        << indent() << "\tmaxiter       = 500   // Maximum number of augmented Lagrangian steps" << endl
        << indent() << "\tresidual      = 1.e-5 // Maximum H1 equivalent residual" << endl
        << indent() << "\tstorage_level = 3     // Storage level of the velocity for convergence " << endl
        << indent() << "\taugmentation  = 10.   // Augmentation parameter" << endl
        << indent() << "END MODULE ViscoplasticTerm" << endl;
        PEL_Error::object()->raise_missing_module(exp,"ViscoplasticTerm");
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

    PEL_CHECK_POST( L_ViscoplasticU > L_UU + ORDER );
    PEL_CHECK_POST( ViscoplasticType == "Bingham" || ViscoplasticType == "HerschelBulkley");
}

//---------------------------------------------------------------------------
Leo_NavierStokes:: ~Leo_NavierStokes( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "Leo_NavierStokes:: ~Leo_NavierStokes" ) ;
}

//---------------------------------------------------------------------------
void
Leo_NavierStokes:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "Leo_NavierStokes:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   //start_total_timer( "Leo_NavierStokes:: do_one_inner_iteration" ) ;

   // ======================================================================
   // = Step 1 (Navier Stokes Problem)
   // ---------------------------------------------------------------------

   PEL::out() << indent() << "Step 1: Navier Stokes problem" << endl;

   GLOBAL_EQ->nullify_for_new_internal_iteration() ;
   GLOBAL_EQ->add_explicit_terms() ; // A_explicit and F_explicit
   GLOBAL_EQ->add_Outer_terms() ;    // A_Outer and F_Outer

   GLOBAL_EQ->set_TminusRGAMMA( ViscoplasticAug ) ;

   setup_Inner( t_it ); // contains GLOBAL_EQ->add_Visco_terms() and
                        // GLOBAL_EQ->add_Flowrate_terms() and
                        // GLOBAL_EQ->add_BCstress_terms()

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
             "*** Leo_NavierStokes error:\n"
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

    // ======================================================================================
    // Step 2: Compute strain rate tensor GAMMA
    // ======================================================================================
    PEL::out() << indent() << "Step 2: Calculate gamma" << endl;
    Leo_Viscoplastic::compute_strain_rate_tensor( cFE, UU, GAMMADOT); // Compute GAMMADOT(u)
    if(ViscoplasticType == "Bingham")
    	Leo_Viscoplastic::update_strain_rate_tensor_d( cFE, STRESS,
             GAMMADOT, GAMMA, ViscoplasticAug, KAPPA, BN, FUNC_E, t_it );
    // Electrical field: FUNC_E =f(|E|) used to set B_new = B f(|E|)

    if(ViscoplasticType == "HerschelBulkley")
    {
    	Leo_Viscoplastic::update_strain_rate_tensor_d_HB( cFE, STRESS,
             GAMMADOT, GAMMA, ViscoplasticAug, KAPPA, BN, HB_n, FUNC_E, t_it );
    }


    // =====================================================================================
    // Step 3:  Update viscoplastic Lagrange multiplier Stress
    // =====================================================================================
    PEL::out() << indent() << "Step 3: Update Lagrange multiplier" << endl;
    Leo_Viscoplastic::update_Lagrange_multiplier( cFE, STRESS,
           GAMMADOT, GAMMA, ViscoplasticAug );

	ViscoplasticIter++;
	// -------------------------------------------------------------------------

    //stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
Leo_NavierStokes:: update_pressuredrop(  )
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
Leo_NavierStokes:: setup_Outer( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{

	PEL_LABEL( "Leo_NavierStokes:: setup_Outer" ) ;

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
Leo_NavierStokes:: setup_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: setup_Inner" ) ;

	double inner_rhs_scale = 0.;

	start_assembling_timer();

    // Remove old entries in the matrix:
    if ( FACTORIZE )  {
        if( ViscoplasticIter == 0 )
            GLOBAL_EQ->set_calculate_LHS(true);
        else
            GLOBAL_EQ->set_calculate_LHS(false);
    }
    GLOBAL_EQ->nullify_A_F_Inner() ;

	// Setup coeff of the BDF
	double xx = 1.0/t_it->time_step() ;
	if ( ORDER == 2 && t_it->iteration_number() != 1 ) xx *= 1.5 ;
	GLOBAL_EQ->set_leading_BDF_over_dt( xx ) ;

	// Copy the value of the n+1 velocity to the nonlinear location level
	if ( Viscoplastic )
	{
		PEL_CHECK( UU->storage_depth() >= L_ViscoplasticU );
		UU->copy_DOFs_value( L_UPDATE_UU, L_ViscoplasticU );
	}

	if ( Flowrate )
		inner_rhs_scale = *(FlowratedPdL.end()-1);

	// Call the loops
	//loop_on_cells_Inner( t_it );
	//loop_on_bounds_Inner( t_it );


	// Assemble the main matrix system
	//GLOBAL_EQ->add_Inner_terms() ;
	GLOBAL_EQ->add_Visco_terms() ;

	if ( Flowrate )
	    GLOBAL_EQ->add_Flowrate_terms( inner_rhs_scale ) ; // zero for (!Flowrate)

	GLOBAL_EQ->add_BCstress_terms() ;
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
 *   - Laplace part: \f$ \frac{r}{2} \int_\Omega {\dot \gamma}_{ij}(u) \ {\dot \gamma}_{ij}(v)  \ dV \f$
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
Leo_NavierStokes:: loop_on_cells_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
        PEL_LABEL( "Leo_NavierStokes:: loop_on_cells_Outer" ) ;
        size_t nbc = UU->nb_components() ;
        size_t nbc_gammadot = GAMMADOT->nb_components() ;
        size_t nb_dims = cFE->nb_space_dimensions() ;
        PEL_ASSERT( nbc == nb_dims ) ;

        doubleVector aa( nb_dims ) ;
        doubleVector coef( nbc ) ;
        doubleVector tdmu( nb_dims ) ;

	    double dt = t_it->time_step() ;
	    double Re = 1.;
	    double Alpha = 0.;

        PEL::out()
            << indent() << "Assemble the following matrices: " << endl
            << indent() << " - Advection term ( active = " << ADV << ")" << endl
            << indent() << " - Time dependent term ( LHS+RHS ) (active if ALPHA > 0)" << endl
            << indent() << " - Relaxed viscosity form" << endl;

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
                FE::add_row( ELEMENT_EQ, cFE, coef ) ; // right hand side: force + old timesteps
                FE::add_row_col_S( ELEMENT_EQ, cFE, bq ) ; // active time dependent term
                FE::add_grad_row_D_col_S( ELEMENT_EQ, cFE, ViscoplasticAug ) ; // relaxed "Viscosity" form
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

    PEL::out()
    << indent() << " - Assemble -2*(tau:grad(v)) = -tau:gammadot(v)";
    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        // Divergence matrix tau:grad(v)
        cFE->set_row_and_col_fields( UU, GAMMADOT ) ;
        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                                cFE->col_field_node_connectivity(), nbc_gammadot ) ;

        cFE->start_IP_iterator( QRP ) ;
        for (  ; cFE->valid_IP() ; cFE->go_next_IP() )
        {
            FE_TensorFormAssembling::add_grad_row_col_Leo( ELEMENT_EQ, cFE, -1. ) ;
        }
        GLOBAL_EQ->assemble_A_LagrangeMultDivMatrix( ELEMENT_EQ );
    }
}

//---------------------------------------------------------------------------
/**
 * @brief Assemble finite element contributions on all cells in the inner iteration.
 * Add all contributions together.
 *
 * @param t_it Time iterator
 *
 */
void
Leo_NavierStokes:: loop_on_cells_Inner( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: loop_on_cells_Inner" ) ;
	size_t nbc = UU->nb_components() ;
	size_t nb_dims = cFE->nb_space_dimensions() ;
	PEL_ASSERT( nbc == nb_dims ) ;
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
Leo_NavierStokes:: loop_on_bounds_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
    PEL_LABEL( "Leo_NavierStokes:: loop_on_bounds" ) ;

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
			    	    "*** Leo_NavierStokes error:\n"
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
Leo_NavierStokes:: loop_on_bounds_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: loop_on_bounds_Inner" ) ;
}

//---------------------------------------------------------------------------
void
Leo_NavierStokes:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: do_before_time_stepping" ) ;

	start_total_timer( "Leo_NavierStokes:: do_before_time_stepping" ) ;

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
Leo_NavierStokes:: do_before_inner_iterations_stage( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "Leo_NavierStokes:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "Leo_NavierStokes:: do_before_inner_iterations_stage" ) ;

   ViscoplasticIter = 0;

   //GLOBAL_EQ->re_initialize() ; // done within setup_Outer

   setup_Outer(t_it);

   start_assembling_timer();
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
                }
            }
            FE::add_row( ELEMENT_EQ, cFE, coef ) ;
         }

         GLOBAL_EQ->assemble_A_F_explicit( ELEMENT_EQ ) ;
   }

   stop_assembling_timer();
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
Leo_NavierStokes::do_inner_iterations_stage( FE_TimeIterator const* t_it )
{
	PEL_LABEL( "Leo_NavierStokes:: do_inner_iterations_stage" ) ;
    start_total_timer( "Leo_NavierStokes:: do_inner_iterations_stage" ) ;

    if( ElectricField )
    {

       double norm=param_euclidian_norm(t_it, E_FIELD);
       FUNC_E = param_yield_stress_function( norm ) ;

       if( verbose_level() >= 1 )
       {
            PEL::out() << indent() << "---------------------------------------------------------" << std::endl;
            PEL::out() << indent() << "E_FIELD: #comp = " << E_FIELD->nb_components() ;
            PEL::out() << ", |E| = " << norm  << ", func_E = " << FUNC_E << std::endl;

            PEL::out() << indent() << "---------------------------------------------------------" << std::endl;
       }
    }

	do
    	do_one_inner_iteration( t_it ) ;
	while ( !inner_iterations_are_completed( t_it ) ) ;
    stop_total_timer();
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
Leo_NavierStokes:: inner_iterations_are_completed( FE_TimeIterator const* t_it ) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: inner_iterations_are_completed" ) ;
	if ( Viscoplastic )
	{
		if ( ViscoplasticIter > 1 )
			decrease_indent();
        // Compute || u_current - u_old ||
		double error_velocity = compute_L2_norm( UU, L_UU , L_ViscoplasticU );
        // Compute ||gammadot-gamma||
        double error_strainrate=Leo_Viscoplastic::compute_norm_Dminusd( cFE, GAMMA, GAMMADOT );
        // Take the maximum as convergence error
        double error = max(error_velocity, error_strainrate);

		PEL::out().setf(std::ios::scientific);
		PEL::out() << indent()
        << "Iteration " << ViscoplasticIter << ": "
        << " || DeltaU || = " << error_velocity
        << ",  || Gamma - Gammadot || = " << error_strainrate
        << std::endl;

		if ( ViscoplasticIter > ViscoplasticMaxiter  )
		{
			return true;
		}
		else if ( error < ViscoplasticResidual )
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
Leo_NavierStokes::compute_L2_norm( PDE_DiscreteField const * FIELD,
								  size_t level1, size_t level2) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: compute_L2_square_norm" ) ;
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
Leo_NavierStokes::compute_flowrate( PDE_DiscreteField const * FIELD,
									size_t level1) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "Leo_NavierStokes:: compute_L2_square_norm" ) ;
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

//---------------------------------------------------------------------------
/**
 * @brief L2 norm of the electrical field
 * @param Parameter
 * @return L2 norm
 */
double
Leo_NavierStokes::param_euclidian_norm( FE_TimeIterator const* t_it , FE_Parameter const * param) const
//---------------------------------------------------------------------------
{
    PEL_LABEL( "Leo_NavierStokes:: param_euclidian_norm" ) ;

    double norm = 0. ;
    size_t const nbc = param->nb_components() ;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        cFE->set_calculation_point( cFE->polyhedron()->center() ) ;

        cFE->start_IP_iterator( QRP ) ;
        for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
        {
            double weight = cFE->weight_of_IP() ;

            for(size_t ic=0; ic< nbc; ic++)
            {
                double value = param->cell_value_at_IP( t_it, cFE, ic ) ;
                norm += weight * value*value ;
//                PEL::out() << ic << " cell_value_at_IP : " ;
//                PEL::out() << "   value=" << value << std::endl ;
            }
        }
    }

    norm = COM->sum( norm ) ;
    return sqrt(norm);
}

//---------------------------------------------------------------------------
/**
 * @brief computes f(|E|) where B_new = B f(|E|)
 *
 * \f$ \tau_Y = \alpha E_0^2 \left( \frac{\tanh(\sqrt{E_o/E_c})}{\sqrt{E_0/E_c}}\right)\f$
 * and
 * \f$ B = \frac{\tau_Y D}{\mu_0 U_0}\f$
 *
 * Choose B as \f$ B = \frac{\alpha E_c^2 D}{\mu_0 U_0}\f$, and \f$ E = E_0/E_c\f$. Then,
 *
 * \f[ f(|E|) = |E|^{3/2} \tanh \sqrt{|E|} \f]
 *
 * @param[in] param_norm |E|
 * @return function_value f(|E)
 */
double
Leo_NavierStokes::param_yield_stress_function( double const param_norm) const
//---------------------------------------------------------------------------
{
    PEL_LABEL( "Leo_NavierStokes:: param_euclidian_norm" ) ;

    return PEL::pow( param_norm,  3./2. ) * PEL::tanh( PEL::sqrt( param_norm) ) ;
}
