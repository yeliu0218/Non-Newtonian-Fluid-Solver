/*
 *  Copyright : 
 *    "Institut de Radioprotection et de S�ret� Nucl�aire - IRSN" (1995-2008)
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

#include <FE_ViscosityHBParameter.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

#include <cmath>
#include <sstream>

using std::string ;

FE_ViscosityHBParameter const* 
FE_ViscosityHBParameter:: PROTOTYPE = new FE_ViscosityHBParameter() ;

//-------------------------------------------------------------------------
FE_ViscosityHBParameter:: FE_ViscosityHBParameter( void )
//-------------------------------------------------------------------------
	: FE_Parameter( "FE_ViscosityHBParameter" )
	, reg("Simple")
	, reg_param(0)
	, NBCS(0)
{
}

//-------------------------------------------------------------------------
FE_ViscosityHBParameter*
FE_ViscosityHBParameter:: create_replica( PEL_Object* a_owner,
        	                            PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ViscosityHBParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_ViscosityHBParameter* result = 
                              new FE_ViscosityHBParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_ViscosityHBParameter:: FE_ViscosityHBParameter( PEL_Object* a_owner,
        		                               PDE_DomainAndFields const* dom,
                                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
	: FE_Parameter( a_owner, exp->string_data( "name" ) )
	, FU( dom->set_of_discrete_fields()->item( exp->string_data( "vel_name" ) ) )
	, L_FU( exp->int_data( "vel_level" ) )
	, NBCS_FU( FU->nb_components() )
	, reg( exp->string_data( "reg") )
	, reg_param( exp->double_data( "reg_param") )
	, NBCS(0)
{
   // Check the link to the velocity field
   if( FU->nb_components() < 1 )
   {
	   std::string mesg =
			   "*** FE_ViscosityHBParameter : the associated field\n"
			   "    should have at least one component" ;
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
   
   
   //////////////////////////////////////////////////////////////////
   // Setup the list of parameters
   //////////////////////////////////////////////////////////////////
   
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "list_of_parameters" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
		// Load and check submodule
	   PEL_ModuleExplorer* ee = e->create_subexplorer( e ) ;
		
		// Define correct module name: -> erase leading param#
	   string module_name = ee->name();
	   size_t pos = module_name.find("#")+1;
	   if ( pos!= string::npos )
		   module_name.erase( 0, pos );
		
	   if ( module_name == "Bn" || module_name == "kappa" || module_name == "n" )
	   {
		   FE_Parameter* prm = 0 ;
		   string const& pmtype = ee->string_data( "type" ) ;
		   if( pmtype == "already_defined" )
		   {
			   FE_PARAMS_REALNAMES[ module_name ] = ee->string_data( "name" );
		   }
		   else if( pmtype == "to_be_defined" )
		   {
			   PEL_ModuleExplorer* eee = ee->create_subexplorer( ee, 
					   "FE_Parameter" ) ;
			   prm = FE_Parameter::make( this, dom, eee ) ;
			   FE_PARAMS_REALNAMES[ module_name ] = prm->name() ;
		   }
		   else
		   {
			   PEL_Error::object()->raise_bad_data_value( ee, "type",
								 "\"already_defined\" and \"to_be_defined\"" ) ;
		   }
	
		   if( prm != 0 && NBCS == 0 )
		   {
			   NBCS = prm->nb_components() ;
		   }
		   else if( prm != 0 && NBCS != prm->nb_components() )
		   {
			   std::string mesg =
					   "*** FE_ViscosityHBParameter error:\n"
					   "    should be defined from FE_Parameter instances\n"
					   "    with the same number of components\n";
			   std::ostringstream buffer;
			   buffer 	<< "    - NBCS = " << NBCS << "\n"
					    << "    - ParName = " << prm->name() << "\n" 
					    << "    - ParCmp  = " << prm->nb_components() << "\n";
			   mesg += buffer.str();
			   PEL_Error::object()->raise_plain( mesg ) ;
		   }
		   FE_PARAMS[ module_name ] = prm;
	   }
	   else
	   {
		   std::string mesg =
				   "*** FE_ConvexParameter:\n"
				   "    - accepted names of the modules:  Bn, kappa, n\n"
				   "      you supplied: ";
		   mesg += module_name;
		   PEL_Error::object()->raise_plain( mesg ) ;
	   }
   }
   e->destroy() ;

   if( FE_PARAMS.size() == 0 )
   {
	   std::string mesg =
			   "*** FE_ViscosityHBParameter error:\n"
			   "    empty list_of_parameters" ;
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
	
	// Postconditions:

}


//-------------------------------------------------------------------------
/**
 * @brief Execute the links. This function is responsible to initialize the
 * already existing parameters in FE_SetOfParameters. 
 * 
 * This is necessary as the collection of parameters is not available during 
 * creation of the oobject.
 * 
 * @param prms Collection of all system wide parameters.
 */
void 
FE_ViscosityHBParameter::do_the_links( FE_SetOfParameters const* prms ) 
//-------------------------------------------------------------------------
{
	PEL_LABEL( "FE_ViscosityHBParameter:: do_the_links" ) ;
	PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

	for( FE_PARAMS_MAP::iterator it = FE_PARAMS.begin(); it != FE_PARAMS.end(); ++it )
	{
		FE_Parameter* prm = it->second ;
		if( prm == 0 )
		{
			prm = prms->item( FE_PARAMS_REALNAMES[ it->first ] ) ;
			it->second = prm ;
			if( NBCS == 0 )
			{
				NBCS = prm->nb_components() ;
			}
			else if( NBCS != prm->nb_components() )
			{
				std::string mesg = "FE_ViscosityHBParameter should be defined from\n" ;
				mesg += "FE_Parameter instances with the same number of components" ;
				PEL_Error::object()->raise_plain( mesg ) ;
			}
		}
		prm->do_the_links( prms ) ;
	}
	
	// Set the parameters:
	ParBn= FE_PARAMS["Bn"];
	ParKappa = FE_PARAMS["kappa"];
	ParN = FE_PARAMS["n"];
	
	PEL_CHECK_POST( ParBn != 0 );
	PEL_CHECK_POST( ParKappa != 0 );
	PEL_CHECK_POST( ParN != 0 );
}



//-------------------------------------------------------------------------
FE_ViscosityHBParameter:: ~FE_ViscosityHBParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_ViscosityHBParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
FE_ViscosityHBParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
	PEL_LABEL( "FE_ViscosityHBParameter:: cell_value_at_pt" ) ;
	PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

	double gammadot=0.;
	FE::add_D_u_D_u_at_pt( gammadot, FU, L_FU, fe, 0.5 );
	gammadot = pow(gammadot,0.5);

	if (reg=="Simple")
	{
		return ( Simple( gammadot, 
				 		 ParBn->cell_value_at_pt( t_it, fe, 0 ),
						 ParKappa->cell_value_at_pt( t_it, fe, 0 ),
						 ParN->cell_value_at_pt( t_it, fe, 0 )
					   ) );
	}
	else
	{
		std::string mesg =
				"*** FE_ViscosityHBParameter : Incorrect Fluid model";
		PEL_Error::object()->raise_plain( mesg ) ;
	}
	return 2;
}

//-------------------------------------------------------------------------
double
FE_ViscosityHBParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
	PEL_LABEL( "FE_ViscosityHBParameter:: cell_value_at_IP" ) ;
	PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

	double gammadot=0.;
	FE::add_D_u_D_u_at_IP( gammadot, FU, L_FU, fe, 0.5 );
	gammadot = pow(gammadot,0.5);

	if (reg=="Simple")
	{
		return ( Simple( gammadot, 
				 		 ParBn->cell_value_at_IP( t_it, fe, 0 ),
						 ParKappa->cell_value_at_IP( t_it, fe, 0 ),
						 ParN->cell_value_at_IP( t_it, fe, 0 )
					   ) );
	}
	else
	{
		std::string mesg =
				"*** FE_ViscosityHBParameter : Incorrect Fluid model";
		PEL_Error::object()->raise_plain( mesg ) ;
	}
	return 2;
}

//-------------------------------------------------------------------------
double
FE_ViscosityHBParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ViscosityHBParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double gammadot=0.;
   FE::add_D_u_D_u_at_pt( gammadot, FU, L_FU, fe, 0.5 );
   gammadot = pow(gammadot,0.5);

   if (reg=="Simple")
   {
	   return ( Simple( gammadot, 
				 		 ParBn->bound_value_at_pt( t_it, fe, 0 ),
						 ParKappa->bound_value_at_pt( t_it, fe, 0 ),
						 ParN->bound_value_at_pt( t_it, fe, 0 )
					   ) );
   }
   else
   {
	   std::string mesg =
			   "*** FE_ViscosityHBParameter : Incorrect Fluid model";
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
   return 2;
}

//-------------------------------------------------------------------------
double
FE_ViscosityHBParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ViscosityHBParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double gammadot=0.;
   FE::add_D_u_D_u_at_IP( gammadot, FU, L_FU, fe, 0.5 );
   gammadot = pow(gammadot,0.5);

   if (reg=="Simple")
   {
	   return ( Simple( gammadot, 
				 		 ParBn->bound_value_at_IP( t_it, fe, 0 ),
						 ParKappa->bound_value_at_IP( t_it, fe, 0 ),
						 ParN->bound_value_at_IP( t_it, fe, 0 )
					   ) );
   }
   else
   {
	   std::string mesg =
			   "*** FE_ViscosityHBParameter : Incorrect Fluid model";
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
   return 2;
}

//-------------------------------------------------------------------------
void
FE_ViscosityHBParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ViscosityHBParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   
   fe->require_field_calculation( FU, PDE_LocalFE::dN ) ;
   
   ParBn->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
   ParKappa->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
   ParN->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
   
}

//-------------------------------------------------------------------------
void
FE_ViscosityHBParameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ViscosityHBParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

   fe->require_field_calculation( FU, PDE_LocalFE::dN ) ;
   
   ParBn->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   ParKappa->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   ParN->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
}




//-------------------------------------------------------------------------
/** @brief Calculate the viscosity of the Hershel-Bulkley law regularized using
 * the Simple approach.
 * 
 * \f[
 * \mu(\gamma(u)) := \kappa \ (\gamma(u)+\epsilon)^{n-1} + \frac{Bn}{\gamma(u)+\epsilon}
 * \f]
 * @param c concentration
 * @param gammadot
 * @return vsicosity of the Simple model
 * 
 * @todo
 * - Check preconditions: Bn, n and kappa available
 */
double 
FE_ViscosityHBParameter::Simple( double gammadot, 
								 double Bn, double kappa, double n ) const
//-------------------------------------------------------------------------
{
	return ( kappa*pow(gammadot+reg_param, n-1) + Bn/(gammadot+reg_param) );
}
