#include <MI_ViscoElastic.hh>

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
#include <FE_SetOfParameters.hh>
#include <FE_TauStab.hh>
#include <FE_TensorFormAssembling.hh>
#include <FE_TimeIterator.hh>

#include <FE_PrintVariables.hh>
#include <PEL_DataWithContextExp.hh>
#include <MI_ViscoElasticSystem.hh>

#include <UT_Viscoplastic.hh>
#include <UT_Viscoelastic.hh>

#include <LA_SymmetricMatrix.hh>
#include <LA_DenseMatrix.hh>
#include <LA_SeqVector.hh>
#include <UT_Viscoelastic.hh>

#include <iostream>
#include <sstream>
#include <cmath>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

MI_ViscoElastic const* MI_ViscoElastic::PROTOTYPE = new MI_ViscoElastic() ;

//---------------------------------------------------------------------------
MI_ViscoElastic:: MI_ViscoElastic( void )
//--------------------------------------------------------------------------
        : FE_OneStepIteration( "MI_ViscoElastic" )
        , COM(0)
        , COORDS(0)
{
}

//---------------------------------------------------------------------------
MI_ViscoElastic*
                MI_ViscoElastic:: create_replica( PEL_Object* a_owner,
                PDE_DomainAndFields const* dom,
                FE_SetOfParameters const* prms,
                PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
        PEL_LABEL( "MI_ViscoElastic:: create_replica" ) ;
        PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

        MI_ViscoElastic* result =
                        new MI_ViscoElastic( a_owner, dom, prms, exp ) ;

        PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
        return( result ) ;
}

//---------------------------------------------------------------------------
MI_ViscoElastic:: MI_ViscoElastic( PEL_Object* a_owner,
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
        , DD( dom->set_of_discrete_fields()->item( exp->string_data("dd") ) )
        , GAMMADOT( dom->set_of_discrete_fields()->item( exp->string_data("gammadot") ) )
        , SS( dom->set_of_discrete_fields()->item( exp->string_data("ss_field") ) )
        , L_SS( exp->int_data( "level_of_field_SS" ) )
        , SS_EXP( 0 )
        , L_SS_EXP( PEL::bad_index() )
        , SS_EXP_EXP( 0 )
        , L_SS_EXP_EXP( PEL::bad_index() )
        , TAU( dom->set_of_discrete_fields()->item( exp->string_data("stress_field") ) )
        , ADV( false )
        , ORDER( exp->int_data( "time_order" ) )
        , ALPHA( prms->item( exp->string_data( "param_unsteady" ) ) )
        , RE( prms->item( exp->string_data( "param_Reynolds" ) ) )
        , RHSU( prms->item( exp->string_data( "param_source" ) ) )
        , RHSU_DIV( prms->item( exp->string_data( "param_div" ) ) )
        , VISC( prms->item( exp->string_data( "param_viscous" ) ) )
        , KAPPA( prms->item( exp->string_data( "param_diffusion" ) ) )
        , PI( 0 )
        , ETA_P( prms->item( exp->string_data( "param_etap" ) ) )
        , LAMBDA( prms->item( exp->string_data( "param_lambda" ) ) )
        , ALPHA_SS( 0 )
        , TDISC( NoTime )
        , BCs( dom->set_of_boundary_conditions() )
        , LOCAL_BC( 0 )
        , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
        , TAUSTAB( 0 )
        , STAB( NoStab )
        , SGNTAU( 0.0 )
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
        , ElasticIter( 0 )
        , L_NonlinU(1)
        , Flowrate(false)
        , FlowrateCP(0)
        , FlowrateMaxiter(1)
        , COM( PEL_Exec::communicator() )
	    , CONTEXT( PEL_ContextSimple::create( this ) )
	    , COORDS( PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) )
{
    PEL_LABEL( "MI_ViscoElastic:: MI_ViscoElastic" ) ;

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
    if (dom->nb_space_dimensions() == 2)
        check_param_nb_components( RHSU_DIV,  "param_div", 3 ) ;
    else if (dom->nb_space_dimensions() == 3)
    {
        check_param_nb_components( RHSU_DIV,  "param_div", 6 ) ;
        PEL::out() << "Only 2D case is implemented!" << std::endl;
        PEL_Error:: exit() ;
    }
    check_field_storage_depth( SS, L_SS ) ;
    cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ; //???? pas tjrs
    bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
    bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;

    cFE->require_field_calculation( DD, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( GAMMADOT, PDE_LocalFE::N ) ;
    cFE->require_field_calculation( SS, PDE_LocalFE::N )  ;
    cFE->require_field_calculation( SS, PDE_LocalFE::dN ) ;
    bFE->require_field_calculation( SS, PDE_LocalFE::N )  ;
    cFE->require_field_calculation( TAU, PDE_LocalFE::N ) ;

    PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

    ALPHA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    RE->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    VISC->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    RHSU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val );
    RHSU_DIV->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val );

    check_param_nb_components( KAPPA, "param_diffusion", 1 ) ;
    KAPPA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    check_param_nb_components( ETA_P, "param_etap", 1 ) ;
    ETA_P->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    check_param_nb_components( LAMBDA, "param_lambda", 1 ) ;
    LAMBDA->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    size_t nbc_SS = SS->nb_components() ;

    if( exp->has_entry( "param_source_SS" ) )
    {
       PI = prms->item( exp->string_data( "param_source_SS" ) ) ;
       check_param_nb_components( PI, "param_source_SS", nbc_SS ) ;
       PI->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
    }
    PDE_SetOfDiscreteFields const* sdf = dom->set_of_discrete_fields() ;
    if( exp->has_module( "time_discretization_SS") )
    {
       PEL_ModuleExplorer const* se =
                         exp->create_subexplorer( 0, "time_discretization_SS") ;
       SS_EXP = sdf->item( se->string_data( "field_explicit_SS" ) ) ;
       L_SS_EXP = se->int_data( "level_of_field_explicit_SS" ) ;
       check_field_nb_components( SS_EXP, nbc_SS ) ;
       check_field_storage_depth( SS_EXP, L_SS_EXP ) ;

       cFE->require_field_calculation( SS_EXP, PDE_LocalFE::N )  ;

       ALPHA_SS = prms->item( se->string_data( "param_unsteady_SS" ) ) ;
       check_param_nb_components( ALPHA_SS, "param_unsteady_SS", 1 ) ;
       ALPHA_SS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

       string const& tt = se->string_data( "type" ) ;
       if( tt == "Euler" )
       {
          TDISC = Euler ;
       }
       else if( tt == "BDF2" )
       {
          TDISC = BDF2 ;
          SS_EXP_EXP = sdf->item( se->string_data( "field_explicit_explicit_SS" ) ) ;
          check_field_nb_components( SS_EXP_EXP, nbc_SS ) ;
          L_SS_EXP_EXP = se->int_data( "level_of_field_explicit_explicit_SS" ) ;
          check_field_storage_depth( SS_EXP_EXP, L_SS_EXP_EXP ) ;
          cFE->require_field_calculation( SS_EXP_EXP, PDE_LocalFE::N )  ;
       }
       else
       {
          PEL_Error::object()->raise_bad_data_value( se, "type", "\"Euler\" or \"BDF2\"" ) ;
       }
       se->destroy() ; se = 0 ;
    }

    if( exp->has_module("stabilization") )
    {
        PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "stabilization" ) ;
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
        TAUSTAB = FE_TauStab::create( this, sse ) ;
        cFE->require_field_calculation( SS, PDE_LocalFE::d2N ) ;
        se->destroy() ;
    }

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
                        exp->create_subexplorer( 0, "MI_ViscoElasticSystem" ) ;
    PDE_LinkDOF2Unknown* ulink = PDE_LinkDOF2Unknown::create( 0, UU,
                                          "sequence_of_the_components",
                                          true ) ;
    PDE_LinkDOF2Unknown* plink = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;
    PDE_LinkDOF2Unknown* dlink = PDE_LinkDOF2Unknown::create( 0, DD,
                                          "sequence_of_the_components", true ) ;
    PDE_LinkDOF2Unknown* slink = PDE_LinkDOF2Unknown::create( 0, SS,
                                          "sequence_of_the_components", true ) ;
    PDE_LinkDOF2Unknown* taulink = PDE_LinkDOF2Unknown::create( 0, TAU,
                                          "sequence_of_the_components", true ) ;

    GLOBAL_EQ = MI_ViscoElasticSystem::create( this, se, ulink, plink, dlink, slink, taulink) ;
    se->destroy();

    if ( exp->has_entry( "boundary_conditions_types" ) )
    {
        stringVector const& bt =
            exp->stringVector_data( "boundary_conditions_types" ) ;
        string const& fname = exp->string_data( "velocity" ) ;
        LOCAL_BC = FE_LocalBCsBuilder::create( this, dom, fname, bt, prms ) ;
        LOCAL_BC->transfer_calculation_requirements( bFE ) ;
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
}

//---------------------------------------------------------------------------
MI_ViscoElastic:: ~MI_ViscoElastic( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: ~MI_ViscoElastic" ) ;
}



//---------------------------------------------------------------------------
void
MI_ViscoElastic:: update_pressuredrop(  )
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
MI_ViscoElastic:: setup_Outer( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{

	PEL_LABEL( "MI_ViscoElastic:: setup_Outer" ) ;

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
MI_ViscoElastic:: setup_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: setup_Inner" ) ;

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
	GLOBAL_EQ->add_Div_terms() ;
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
MI_ViscoElastic:: setup_Flowrate( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: setup_Flowrate" ) ;
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
MI_ViscoElastic:: loop_on_cells_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
        PEL_LABEL( "MI_ViscoElastic:: loop_on_cells_Outer" ) ;
        size_t nbc = UU->nb_components() ;
        size_t nbc_dd = DD->nb_components() ;
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
					    coef( d ) = (!Flowrate) * ( RHSU->cell_value_at_IP( t_it, cFE, d ) +
                                                    RHSU_DIV->cell_value_at_IP( t_it, cFE, d ) );
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
            // Divergence vector -d:gammadot(v)
            // --------------------------------------------------------------------

            cFE->set_row_and_col_fields( UU, DD ) ;
            ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                                    cFE->col_field_node_connectivity(), nbc_dd ) ;
            cFE->start_IP_iterator( QRP ) ;
            for ( ; cFE->valid_IP() ; cFE->go_next_IP() )
            {
                 FE_TensorFormAssembling::add_grad_row_col( ELEMENT_EQ, cFE, -1. ) ;
            }
            GLOBAL_EQ->assemble_A_GammadotV( ELEMENT_EQ ) ;

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
					coef( d ) =  ( RHSU->cell_value_at_IP( t_it, cFE, d ) +
                                   RHSU_DIV->cell_value_at_IP( t_it, cFE, d ) );
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
MI_ViscoElastic:: loop_on_cells_Inner( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: loop_on_cells_Inner" ) ;
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
MI_ViscoElastic:: loop_on_bounds_Outer( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_ViscoElastic:: loop_on_bounds" ) ;

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
			    	    "*** MI_ViscoElastic error:\n"
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
		        for( ; bFE->valid_IP() ; bFE->go_next_IP() )
		        {
			        FE::add_row( ELEMENT_EQ, bFE, value ) ;
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
MI_ViscoElastic:: loop_on_bounds_Inner( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: loop_on_bounds_Inner" ) ;
}

//---------------------------------------------------------------------------
void
MI_ViscoElastic:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: do_before_time_stepping" ) ;

	start_total_timer( "MI_ViscoElastic:: do_before_time_stepping" ) ;

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
MI_ViscoElastic:: do_before_inner_iterations_stage( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "MI_ViscoElastic:: do_before_inner_iterations_stage" ) ;

   NonlinIter = 0;
   ElasticIter= 0;
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
MI_ViscoElastic::do_inner_iterations_stage( FE_TimeIterator const* t_it )
{
	PEL_LABEL( "MI_ViscoElastic:: do_inner_iterations_stage" ) ;

	do
    	do_one_inner_iteration( t_it ) ;
	while ( !inner_iterations_are_completed( t_it ) ) ;
}


//---------------------------------------------------------------------------
void
MI_ViscoElastic:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

//   start_total_timer( "MI_ViscoElastic:: do_one_inner_iteration" ) ;
   // ---------------------

    /*----------------------------------------------------------------------*/
    PEL::out() << "STEP 1: solve for (u,p) where tau, d are known " << std::endl;
    /*----------------------------------------------------------------------*/
    do
        do_one_inner_iteration_UP( t_it ) ;
    while ( !inner_iterations_are_completed_UP( t_it ) ) ;
    /*----------------------------------------------------------------------*/
    PEL::out() << "STEP 2: solve for d which is an L2 projection of gammadot(u) " << std::endl;
    /*----------------------------------------------------------------------*/
    do_one_inner_iteration_DEVSS( t_it ) ;

    /*----------------------------------------------------------------------*/
    PEL::out() << "STEP 3: solve for s using the velocity field from step 1 " << std::endl;
    /*----------------------------------------------------------------------*/

    do_one_inner_iteration_SS( t_it ) ;

    /*----------------------------------------------------------------------*/
    PEL::out() << "STEP 4: solve for tau which is an L2 projection " << std::endl;
    /*----------------------------------------------------------------------*/

    do_one_inner_iteration_TAU( t_it ) ;
 //  stop_total_timer() ;
   // -------------
   ElasticIter++;
}

//---------------------------------------------------------------------------
void
MI_ViscoElastic:: do_one_inner_iteration_UP( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_one_inner_iteration" ) ;
//   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "MI_ViscoElastic:: do_one_inner_iteration_UP" ) ;
   // ---------------------

   GLOBAL_EQ->nullify_for_new_internal_iteration() ;
   GLOBAL_EQ->add_explicit_terms() ; // A_explicit and F_explicit
   GLOBAL_EQ->add_Outer_terms() ;    // A_Outer and F_Outer
   GLOBAL_EQ->set_D() ;

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
             "*** MI_ViscoElastic error:\n"
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
//          PEL::out() << indent() << "L_UU = " << L_UU << std::endl;
//          PEL::out() << indent() << "L_UPDATE_UU= " << L_UPDATE_UU << std::endl;
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
MI_ViscoElastic:: do_one_inner_iteration_DEVSS( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_one_inner_iteration_DEVSS" ) ;
//   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "MI_MI_ViscoElastic:: do_one_inner_iteration_DEVSS" ) ;
   // --------------
   UT_Viscoplastic:: compute_strain_rate_tensor(cFE, UU, GAMMADOT) ;
   loop_on_cells_DEVSS();
   GLOBAL_EQ->estimate_unknowns_DEVSS() ;
   DD->update_free_DOFs_value( 0,
                               GLOBAL_EQ->unknown_vector_D(),
                               GLOBAL_EQ->system_numbering_D()->link() ) ;

   double error_DEVSS = UT_Viscoelastic::compute_norm_DminusGAMMADOT(cFE, DD, GAMMADOT);
   PEL::out() << "||DD-gammadot|| = " << error_DEVSS << std::endl;
   // -------------
   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
MI_ViscoElastic:: do_one_inner_iteration_SS( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_one_inner_iteration_SS" ) ;

   start_total_timer( "MI_ViscoElastic:: do_one_inner_iteration_SS" ) ;
   // --------------
   GLOBAL_EQ->re_initialize_SS();

   loop_on_cells_SS( t_it ) ;
   loop_on_bounds_SS( t_it ) ; // empty

   GLOBAL_EQ->estimate_unknowns_SS() ;
   SS->update_free_DOFs_value( 0,
                               GLOBAL_EQ->unknown_vector_S(),
                               GLOBAL_EQ->system_numbering_S()->link() ) ;
   SS->enforce_constraints_for_DOFs( L_SS ) ;  //JUST ADDED TO FIX THE IMPOSED DOFS AT INLET ?!

   // -------------
   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void
MI_ViscoElastic:: do_one_inner_iteration_TAU( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: do_one_inner_iteration_TAU" ) ;

   start_total_timer( "MI_ViscoElastic:: do_one_inner_iteration_TAU" ) ;
   // --------------
   GLOBAL_EQ->re_initialize_TAU();

   loop_on_cells_TAU( t_it ) ;
//   loop_on_bounds_TAU( t_it ) ; // empty

   GLOBAL_EQ->estimate_unknowns_TAU() ;
   TAU->update_free_DOFs_value( 0,
                               GLOBAL_EQ->unknown_vector_TAU(),
                               GLOBAL_EQ->system_numbering_TAU()->link() ) ;

   TAU->enforce_constraints_for_DOFs( 0 ) ;  //JUST ADDED TO FIX THE IMPOSED DOFS AT INLET ?!

   // -------------
   stop_total_timer() ;
}


//---------------------------------------------------------------------------
/**
 * @brief Convergence of the inner iteration stage For Explicit, there is no need to loop on the system of equations
 *
 */
bool
MI_ViscoElastic:: inner_iterations_are_completed( FE_TimeIterator const* t_it ) const
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_ViscoElastic:: inner_iterations_are_completed" ) ;

        return true;

}

// //---------------------------------------------------------------------------
// /**
//  * @brief Convergence of the inner iteration stage For implicit
//  *
//  */
// bool
// MI_ViscoElastic:: inner_iterations_are_completed( FE_TimeIterator const* t_it ) const
// //---------------------------------------------------------------------------
// {
//     PEL_LABEL( "MI_ViscoElastic:: inner_iterations_are_completed" ) ;
//     
//     double error_SS = compute_L2_norm( SS, L_SS , L_SS+1 ); // is not the correct norm!
//     PEL::out().setf(std::ios::scientific);
//     PEL::out() << indent() << "- SS residual = " << error_SS << std::endl;
//     PEL::out() << "***************************************" << std::endl;
//     PEL::out() << indent() << "ElasticIter = " <<  ElasticIter << std::endl;
//     if ( error_SS < 1.E-5 )
// {
//         return true;
// }
//     else if   ( ElasticIter >10)
// {
//        return true;
// }
//     else
// {
// }        return false;
// 
// }

//---------------------------------------------------------------------------
/**
 * @brief Convergence of the inner iteration stage (U,p)
 *
 * Convergence of solution or maxiter reached
 * @todo
 *   - Implement norm comparator
 * @param t_it
 */
bool
MI_ViscoElastic:: inner_iterations_are_completed_UP( FE_TimeIterator const* t_it ) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: inner_iterations_are_completed_UP" ) ;
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
MI_ViscoElastic::compute_L2_norm( PDE_DiscreteField const * FIELD,
								  size_t level1, size_t level2) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: compute_L2_square_norm" ) ;
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
//                          PEL::out() << indent() << "val1  = " << val1 << std::endl;
//                          PEL::out() << indent() << "val2 = " << val2 << std::endl;



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
MI_ViscoElastic::compute_flowrate( PDE_DiscreteField const * FIELD,
									size_t level1) const
//---------------------------------------------------------------------------
{
	PEL_LABEL( "MI_ViscoElastic:: compute_L2_square_norm" ) ;
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
 * @brief Calculates \f$
 */
void
MI_ViscoElastic::loop_on_cells_SS( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_ViscoElastic:: loop_on_cells_SS" ) ;
    double m_xx, norm_aa, kappa ;
    size_t nb_c = SS->nb_components() ;
    size_t NB_DIMS = UU->nb_components();
    doubleVector rhs( nb_c ) ;
    doubleVector rhs_ss( nb_c ) ;
    doubleVector aa( NB_DIMS ) ;
    
//    PEL::out()<< "# component of U ="<< UU->nb_components() << std::endl;
//    PEL::out()<< "# component of S ="<< nb_c << std::endl;
    
    for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
      cFE->set_row_and_col_fields( SS, SS ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_c,
                              cFE->col_field_node_connectivity(), nb_c ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         for( size_t i=0; i<UU->nb_components(); i++)
             aa(i) = cFE->value_at_IP( UU, L_UU, i ) ;

         compute_coefs_at_IP( t_it, m_xx, aa, norm_aa, kappa , rhs );
         compute_RHS_at_IP( t_it,rhs_ss );

         if( aa != 0 )    
            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;

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
         if( TAUSTAB != 0 )
         {
            double h = cFE->polyhedron()->inter_vertices_maximum_distance() ;
            taustab = TAUSTAB->tau( h, m_xx, norm_aa, kappa ) ;
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
      GLOBAL_EQ->assemble_A_F_SS( ELEMENT_EQ ) ;
   }
}

//---------------------------------------------------------------------------
void
MI_ViscoElastic::loop_on_bounds_SS( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: loop_on_bounds_SS" ) ;
}

//---------------------------------------------------------------------------
/**
 * @brief Calculates \f$ d-\dot\gamma = 0 \f$
 */
void
MI_ViscoElastic::loop_on_cells_DEVSS( void )
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_ViscoElastic:: loop_on_cells_DEVSS" ) ;
    size_t nbc_dd = DD->nb_components() ;
    doubleVector coef( nbc_dd ) ;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        cFE->set_row_and_col_fields( DD, DD ) ;
        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc_dd,
                                cFE->col_field_node_connectivity(), nbc_dd ) ;

        cFE->start_IP_iterator( QRP ) ;
        for (  ; cFE->valid_IP() ; cFE->go_next_IP() )
        {
           for ( size_t d=0 ; d<nbc_dd ; ++d ){
                coef( d ) =  cFE->value_at_IP( GAMMADOT, (size_t) 0, d ) ;
            }
            FE::add_row_col_S( ELEMENT_EQ, cFE, 1. ) ;
            FE::add_row( ELEMENT_EQ, cFE, coef, 1.0 ) ;
        }

       GLOBAL_EQ->assemble_A_F_DEVSS( ELEMENT_EQ );
    }
}

// coefs for the equation
//    m_xx u + aa.grad u - div( kappa grad u ) = rhs
//---------------------------------------------------------------------------
void
MI_ViscoElastic:: compute_coefs_at_IP( FE_TimeIterator const* t_it,
                                               double& m_xx,
                                               doubleVector& aa,
                                               double& norm_aa,
                                               double& kappa,
                                               doubleVector& rhs ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: compute_coefs_at_IP" ) ;

   size_t nb_c = SS->nb_components() ;

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

   if( ALPHA_SS != 0 )
   {
      double dt = t_it->time_step() ;
      doubleVector ss_exp( nb_c, PEL::bad_double() ) ;
      if( TDISC == Euler || ( TDISC == BDF2 && t_it->iteration_number() == 1 ) )
      {
         for( size_t ic=0 ; ic<nb_c ; ++ic )
         {
            ss_exp( ic ) = cFE->value_at_IP( SS_EXP, L_SS_EXP, ic ) ;
         }
      }
      else if( TDISC == BDF2 && t_it->iteration_number() != 1 )
      {
         for( size_t ic=0 ; ic<nb_c ; ++ic )
         {
            ss_exp( ic ) =
               2.0 * cFE->value_at_IP( SS_EXP, L_SS_EXP, ic )
             - 0.5 * cFE->value_at_IP( SS_EXP_EXP, L_SS_EXP_EXP, ic ) ;
         }
      }
      double alpha = ALPHA_SS->cell_value_at_IP( t_it, cFE ) ;
      m_xx += alpha/dt ;
      if( TDISC == BDF2 && t_it->iteration_number() != 1 ) m_xx *= 1.5 ;

      for( size_t ic=0 ; ic<nb_c ; ++ic )
      {
         PEL_ASSERT( ss_exp( ic ) != PEL::bad_double() ) ;
         rhs( ic ) += alpha * ss_exp( ic ) / dt ;
      }
   }

//    if( GAMMA != 0 )
//    {
//       m_xx += GAMMA->cell_value_at_IP( t_it, cFE ) ;
//    }

   norm_aa = 0.0 ;
   for( size_t d=0 ; d<aa.size() ; ++d )
   {
	   norm_aa += aa( d ) * aa( d ) ;
   }
   norm_aa = PEL::sqrt( norm_aa ) ;
}

// coefs for the equation
//    ds/dt + u.grad s = rhs_s
//---------------------------------------------------------------------------
void
MI_ViscoElastic:: compute_RHS_at_IP( FE_TimeIterator const* t_it,
                                      doubleVector& rhs_ss ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MI_ViscoElastic:: compute_RHS_at_IP" ) ;

   size_t dim= UU->nb_components() ;
   size_t SScompNumber = SS->nb_components() ; //here SS
   PEL_ASSERT( dim==2 );
   PEL_ASSERT( SScompNumber==3 );
   bool OK = false;
   double relax_time ;

   LA_SymmetricMatrix* SS_matrix = LA_SymmetricMatrix::create(0, dim);

   LA_SeqVector* eigenvals_log = LA_SeqVector::create(0, dim);
   LA_SeqVector* eigenvals = LA_SeqVector::create(0, dim);
   LA_DenseMatrix* R_matrix = LA_DenseMatrix::create(0, dim, dim); //R = eigenvectors

   LA_DenseMatrix* L_matrix = LA_DenseMatrix::create(0, dim, dim); // L = grad(velocity)
   LA_DenseMatrix* LR_matrix = LA_DenseMatrix::create(0, dim, dim); // LR
   LA_DenseMatrix* L_tilde_matrix = LA_DenseMatrix::create(0, dim, dim); // L_tilde

   LA_SeqVector* RHS_Diag = LA_SeqVector::create( 0, dim ) ; // lambda_dot lambda_inv (see equation (11))
   LA_SymmetricMatrix* RHS_matrix = LA_SymmetricMatrix::create( 0, dim ) ; // omega_tilde log(lambda) + log(lambda) omega_tilde_transpose (see equation (12 +13))
   LA_SymmetricMatrix* RHS_matrix_sum = LA_SymmetricMatrix::create( 0, dim ) ;
  // LA_DenseMatrix* RHS_matrix_Rt = LA_DenseMatrix::create( 0, dim, dim ) ;
  // LA_DenseMatrix* RHS_matrix_total = LA_DenseMatrix::create( 0, dim, dim ) ;

   LA_SeqVector* model_func = LA_SeqVector::create(0, dim) ;

   //step 1
   SS_matrix->set_item(0,0, cFE->value_at_IP(SS, 0, 0));
   SS_matrix->set_item(0,1, cFE->value_at_IP(SS, 0, 2));
   SS_matrix->set_item(1,1, cFE->value_at_IP(SS, 0, 1));

   SS_matrix->eigen_reduce(dim, eigenvals_log, R_matrix, OK);

   if (!OK) {
             PEL::out()<<"Error: Computation of eigenvalues/ eigenvectors in SS equation" << std::endl;
             PEL_Error:: exit();
   }

   //step 2
   for( size_t i=0 ; i<eigenvals->nb_rows(); i++)
   {
       double value = eigenvals_log->item(i);
       eigenvals->set_item(i, PEL::exp(value));
   }

   //step 3
   L_matrix->set_item(0,0, cFE->gradient_at_IP(UU, 0, 0,0));
   L_matrix->set_item(0,1, cFE->gradient_at_IP(UU, 0, 1,0));
   L_matrix->set_item(1,0, cFE->gradient_at_IP(UU, 0, 0,1));
   L_matrix->set_item(1,1, cFE->gradient_at_IP(UU, 0, 1,1));

   LR_matrix->nullify() ;
   L_tilde_matrix->nullify() ;
   LR_matrix->add_Mat_Mat(L_matrix, R_matrix, 1.0);
   L_tilde_matrix->add_tMat_Mat( R_matrix, LR_matrix, 1.0);

   relax_time =  LAMBDA->cell_value_at_IP( t_it, cFE ) ;
   double eps = 0.0 ;   // material constant related to the extensional viscosity
                    // we NEED it in the data_file, don't we?

   // step 4
   RHS_Diag->nullify();
   for (size_t i=0; i < dim; i++)
   {
      UT_Viscoelastic::elastic_model_PPT(model_func, eigenvals, relax_time, eps);
      double value = 2.*L_tilde_matrix->item(i,i) + model_func->item(i) / eigenvals->item(i) ;
      RHS_Diag->set_item(i,value);
   }

   // step 5
   RHS_matrix->nullify();
   double value1_RHS2 = 0.0;
   double value2_RHS2 = 0.0;
   double value3_RHS2 = 0.0;

   for (size_t i=0; i < dim; i++)
   {
        for (size_t j=i+1; j < dim; j++)
       {
                 if(eigenvals->item(i) != eigenvals->item(j))
            {
                value1_RHS2 = (PEL::log(eigenvals->item(j)) - PEL::log(eigenvals->item(i))) /
                                       (eigenvals->item(j) - eigenvals->item(i)) ;
                value2_RHS2 =    eigenvals->item(i)  * L_tilde_matrix->item(j,i)
                               + eigenvals->item(j)  * L_tilde_matrix->item(i,j) ;
                value3_RHS2 = value1_RHS2 * value2_RHS2 ;
                RHS_matrix->set_item(i,j,value3_RHS2) ;
            }
            else
            {
                value3_RHS2 =  L_tilde_matrix->item(i,j) + L_tilde_matrix->item(j,i) ;
                RHS_matrix->set_item(i,j,value3_RHS2) ;
            }
       }
    }

    RHS_matrix_sum->nullify();
    RHS_matrix_sum->add_to_diag(RHS_Diag) ;
    RHS_matrix_sum->add_Mat(RHS_matrix) ;

//     RHS_matrix_Rt->nullify();
//     RHS_matrix_Rt->add_Mat_tMat(RHS_matrix_sum, R_matrix, 1.0 ) ;
// 
//     RHS_matrix_total->nullify();
//     RHS_matrix_total->add_Mat_Mat(R_matrix, RHS_matrix_Rt, 1.0 ) ;
// 
//     rhs_ss( 0 ) = RHS_matrix_total->item(0,0);
//     rhs_ss( 1 ) = RHS_matrix_total->item(1,1);
//     rhs_ss( 2 ) = RHS_matrix_total->item(0,1);
    
    rhs_ss( 0 ) = RHS_matrix_sum->item(0,0);
    rhs_ss( 1 ) = RHS_matrix_sum->item(1,1);
    rhs_ss( 2 ) = RHS_matrix_sum->item(0,1);




    SS_matrix->destroy() ;
    eigenvals_log->destroy() ;
    eigenvals->destroy() ;
    R_matrix->destroy() ;
    L_matrix->destroy() ;
    LR_matrix->destroy() ;
    L_tilde_matrix->destroy() ;
    RHS_Diag->destroy() ;
    RHS_matrix->destroy() ;
    RHS_matrix_sum->destroy() ;
 //   RHS_matrix_Rt->destroy() ;
 //   RHS_matrix_total->destroy() ;
    model_func->destroy() ;

/*
    for (size_t i=0; i<SScompNumber;i++)
        PEL::out()<< "rhs_ss( " << i << ") = " <<   rhs_ss->item(i) << std::endl;
*/
}

//---------------------------------------------------------------------------
/**
 * @brief Calculates \f$ \tau-\frac{\eta_p}{\lambda} \left({\bf s - I} \right) = 0 \f$
 */
void 
MI_ViscoElastic::loop_on_cells_TAU( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
    PEL_LABEL( "MI_ViscoElastic:: loop_on_cells_STRESS" ) ;
    size_t nbc = SS->nb_components() ;
    doubleVector coef( nbc ) ;
    
    size_t nbc_stress = TAU->nb_components() ;
    doubleVector coef_stress( nbc_stress ) ;

    size_t dim= UU->nb_components(); 
    bool OK = false;
    double etap, lambda ;
    double  value; // eta/lambda;
    LA_SymmetricMatrix* SS_total = LA_SymmetricMatrix::create(0, dim);
    LA_SymmetricMatrix* SS_M = LA_SymmetricMatrix::create(0, dim);
    LA_DenseMatrix* SS_LHS= LA_DenseMatrix::create(0, dim , dim);
    LA_DenseMatrix* RHS_TAU = LA_DenseMatrix::create(0, dim, dim);
    LA_SeqVector* eigenvals_ss = LA_SeqVector::create(0, dim); // SS-> eigenvalues
    LA_DenseMatrix* eigenvects_ss = LA_DenseMatrix::create(0, dim, dim); // SS-> eigenvectors
    
    LA_DenseMatrix* RHS_TAU_GLOBAL_left = LA_DenseMatrix::create(0, dim, dim);
    LA_DenseMatrix* RHS_TAU_GLOBAL_total = LA_DenseMatrix::create(0, dim, dim);
    LA_SymmetricMatrix* SS_OLD = LA_SymmetricMatrix::create(0, dim);
    LA_SeqVector* eigenvals_ss_old = LA_SeqVector::create(0, dim); // SS-> eigenvalues
    LA_DenseMatrix* eigenvects_ss_old = LA_DenseMatrix::create(0, dim, dim); // SS-> eigenvectors




    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        cFE->set_row_and_col_fields( TAU, TAU ) ;
        ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc_stress,
                                cFE->col_field_node_connectivity(), nbc_stress ) ;

        cFE->start_IP_iterator( QRP ) ;
        for (  ; cFE->valid_IP() ; cFE->go_next_IP() )
        {
          SS_total->set_item(0,0, cFE->value_at_IP(SS, (size_t) 0, 0));
          SS_total->set_item(0,1, cFE->value_at_IP(SS, (size_t) 0, 2));
          SS_total->set_item(1,1, cFE->value_at_IP(SS, (size_t) 0, 1));
          SS_total->eigen_reduce(dim, eigenvals_ss, eigenvects_ss, OK);
    

          if (!OK) {
             PEL::out()<<"Error: Computation of eigenvalues/ eigenvectors in tau equation" << std::endl;
             PEL_Error:: exit();
          }

          SS_M->nullify();
          for (size_t d=0 ; d<dim ; ++d)
          { 
              double value2 = eigenvals_ss->item(d);
              SS_M->set_item(d,d, PEL::exp (value2));
          }

          SS_LHS->nullify();
          SS_LHS->add_Mat_Mat(eigenvects_ss, SS_M, 1.0);
          RHS_TAU->nullify() ;
          RHS_TAU->add_Mat_tMat( SS_LHS, eigenvects_ss, 1.0);


          SS_OLD->set_item(0,0, cFE->value_at_IP(SS,  L_SS+1, 0));
          SS_OLD->set_item(0,1, cFE->value_at_IP(SS,  L_SS+1, 2));
          SS_OLD->set_item(1,1, cFE->value_at_IP(SS,  L_SS+1, 1));
          SS_OLD->eigen_reduce(dim, eigenvals_ss_old, eigenvects_ss_old, OK);
    

          if (!OK) {
             PEL::out()<<"Error: Computation of eigenvalues/ eigenvectors in tau equation at time n" << std::endl;
             PEL_Error:: exit();
          }
         
          RHS_TAU_GLOBAL_left->nullify();
          RHS_TAU_GLOBAL_left->add_Mat_Mat( eigenvects_ss_old, RHS_TAU , 1.0);

          RHS_TAU_GLOBAL_total->nullify();
          RHS_TAU_GLOBAL_total->add_Mat_tMat( RHS_TAU_GLOBAL_left, eigenvects_ss_old, 1.0);

          etap  = ETA_P->cell_value_at_IP( t_it, cFE ) ;
          lambda  = LAMBDA->cell_value_at_IP( t_it, cFE ) ;   
          value = etap/lambda;
      
          coef( 0 ) = ( RHS_TAU_GLOBAL_total->item( 0,0 )-1.0)*value ;
          coef( 1 ) = ( RHS_TAU_GLOBAL_total->item( 1,1 )-1.0)*value ;
          coef( 2 ) = ( RHS_TAU_GLOBAL_total->item( 0,1 ))*value ;

          FE::add_row_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
          FE::add_row( ELEMENT_EQ, cFE, coef, 1.0 ) ;            
        }
        GLOBAL_EQ->assemble_A_F_TAU( ELEMENT_EQ );                                                                        
//        PDE::assemble_in_matrix_vector_0( A_STRESS, F_STRESS, ELEMENT_EQ, NMB_STRESS ) ;
    }
    SS_total->destroy() ;
    SS_M->destroy() ;
    SS_LHS->destroy() ;
    RHS_TAU->destroy() ;
    eigenvals_ss->destroy() ;
    eigenvects_ss->destroy() ;
    eigenvals_ss_old->destroy() ;
    eigenvects_ss_old->destroy() ;
    RHS_TAU_GLOBAL_total->destroy() ;
    RHS_TAU_GLOBAL_left->destroy() ;
    SS_OLD->destroy() ;
}
