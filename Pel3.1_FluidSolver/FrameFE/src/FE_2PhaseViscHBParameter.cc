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

#include <FE_2PhaseViscHBParameter.hh>

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

FE_2PhaseViscHBParameter const* 
FE_2PhaseViscHBParameter:: PROTOTYPE = new FE_2PhaseViscHBParameter() ;

//-------------------------------------------------------------------------
FE_2PhaseViscHBParameter:: FE_2PhaseViscHBParameter( void )
//-------------------------------------------------------------------------
	: FE_Parameter( "FE_2PhaseViscHBParameter" )
	, Bn_1(1)
	, Bn_2(1)
	, n_1(2)
	, n_2(2)
	, kappa_r(2)
	, reg("Simple")
	, reg_param(0)
{
}

//-------------------------------------------------------------------------
FE_2PhaseViscHBParameter*
FE_2PhaseViscHBParameter:: create_replica( PEL_Object* a_owner,
        	                            PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_2PhaseViscHBParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_2PhaseViscHBParameter* result = 
                              new FE_2PhaseViscHBParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_2PhaseViscHBParameter:: FE_2PhaseViscHBParameter( PEL_Object* a_owner,
        		                               PDE_DomainAndFields const* dom,
                                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
	: FE_Parameter( a_owner, exp->string_data( "name" ) )
	, FC( dom->set_of_discrete_fields()->item( exp->string_data( "conc_name" ) ) )
	, L_FC( exp->int_data( "conc_level" ) )
	, NBCS_FC( 1 )
	, FU( dom->set_of_discrete_fields()->item( exp->string_data( "vel_name" ) ) )
	, L_FU( exp->int_data( "vel_level" ) )
	, NBCS_FU( 1 )

	, Bn_1( exp->double_data( "Bn_1") )
	, Bn_2( exp->double_data( "Bn_2") )
	, n_1( exp->double_data( "n_1") )
	, n_2( exp->double_data( "n_2"))
	, kappa_r( exp->double_data( "kappa_r") )
	, reg( exp->string_data( "reg"))
	, reg_param( exp->double_data( "reg_param") )
{
   if( FC->nb_components() != 1 )
   {
      std::string mesg =
         "*** FE_2PhaseViscHBParameter : the associated concentration field\n"
         "    should have exactly one component\n" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
   else if( FU->nb_components() < 1 )
   {
	   std::string mesg =
			   "*** FE_2PhaseViscHBParameter : the associated field\n"
			   "    should have at least one component" ;
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
}

//-------------------------------------------------------------------------
FE_2PhaseViscHBParameter:: ~FE_2PhaseViscHBParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_2PhaseViscHBParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS_FC ) ;
}

//-------------------------------------------------------------------------
double
FE_2PhaseViscHBParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
	PEL_LABEL( "FE_2PhaseViscHBParameter:: cell_value_at_pt" ) ;
	PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

	double c = fe->value_at_pt( FC, L_FC, ic ) ;
	double gammadot=0.;
	FE::add_D_u_D_u_at_pt( gammadot, FU, L_FU, fe, 0.5 );
	gammadot = pow(gammadot,0.5);

	if (reg=="Simple")
	{
		return ( Simple(c, gammadot) );
	}
	else
	{
		std::string mesg =
				"*** FE_2PhaseViscHBParameter : Incorrect Fluid model";
		PEL_Error::object()->raise_plain( mesg ) ;
	}
	return 2;
}

//-------------------------------------------------------------------------
double
FE_2PhaseViscHBParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
	PEL_LABEL( "FE_2PhaseViscHBParameter:: cell_value_at_IP" ) ;
	PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

	double c = fe->value_at_IP( FC, L_FC, ic ) ;
	double gammadot=0.;
	FE::add_D_u_D_u_at_IP( gammadot, FU, L_FU, fe, 0.5 );
	gammadot = pow(gammadot,0.5);

	if (reg=="Simple")
	{
		return ( Simple(c, gammadot) );
	}
	else
	{
		std::string mesg =
				"*** FE_2PhaseViscHBParameter : Incorrect Fluid model";
		PEL_Error::object()->raise_plain( mesg ) ;
	}
	return 2;
}

//-------------------------------------------------------------------------
double
FE_2PhaseViscHBParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_2PhaseViscHBParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double c = fe->value_at_pt( FC, L_FC, ic ) ;
   double gammadot=0.;
   FE::add_D_u_D_u_at_pt( gammadot, FU, L_FU, fe, 0.5 );
   gammadot = pow(gammadot,0.5);

   if (reg=="Simple")
   {
	   return ( Simple(c, gammadot) );
   }
   else
   {
	   std::string mesg =
			   "*** FE_2PhaseViscHBParameter : Incorrect Fluid model";
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
   return 2;
}

//-------------------------------------------------------------------------
double
FE_2PhaseViscHBParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_2PhaseViscHBParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double c = fe->value_at_IP( FC, L_FC, ic ) ;
   double gammadot=0.;
   FE::add_D_u_D_u_at_IP( gammadot, FU, L_FU, fe, 0.5 );
   gammadot = pow(gammadot,0.5);

   if (reg=="Simple")
   {
	   return ( Simple(c, gammadot) );
   }
   else
   {
	   std::string mesg =
			   "*** FE_2PhaseViscHBParameter : Incorrect Fluid model";
	   PEL_Error::object()->raise_plain( mesg ) ;
   }
   return 2;
}

//-------------------------------------------------------------------------
void
FE_2PhaseViscHBParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_2PhaseViscHBParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   
   fe->require_field_calculation( FC, PDE_LocalFE::N ) ;
   fe->require_field_calculation( FU, PDE_LocalFE::dN ) ;
}

//-------------------------------------------------------------------------
void
FE_2PhaseViscHBParameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_2PhaseViscHBParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

   fe->require_field_calculation( FC, PDE_LocalFE::N ) ;
   fe->require_field_calculation( FU, PDE_LocalFE::dN ) ;
}

//-------------------------------------------------------------------------
/** @brief Calculate the Bingham number depending on a concentration parameter
 * 
 * \f[
 * 	Bn:=c \ Bn_1 + (1-c) \ Bn_2
 * \f]
 * @param c concentration
 * @return local Bingham number
 */
double 
FE_2PhaseViscHBParameter::Bn( double c ) const
//-------------------------------------------------------------------------
{
	return ( c*Bn_1 + ( 1-c)*Bn_2 );
}

//-------------------------------------------------------------------------
/** @brief Calculate the consistency  depending on a concentration parameter
 * 
 * \f[
 * \kappa : =c  + (1-c) \ \frac{\hat \mu_1}{\hat \mu_2}
 * \f]
 * where \f$ \hat \mu_i:= \hat \kappa_i \left[ \frac{U_0}{D} \right]^{n_i-1} \f$
 * @param c concentration
 * @return local consitency
 */
double 
FE_2PhaseViscHBParameter::kappa( double c ) const
//-------------------------------------------------------------------------
{
	return ( c + ( 1-c)*kappa_r );
}


//-------------------------------------------------------------------------
/** @brief Calculate the power law index depending on a concentration parameter
 * 
 * \f[
 * n : = c n_1 + (1-c) n_2
 * \f]
 * where \f$ \hat \mu_i:= \hat \kappa_i \left[ \frac{U_0}{D} \right]^{n_i-1} \f$
 * @param c concentration
 * @return local consitency
 */
double 
FE_2PhaseViscHBParameter::n( double c ) const
//-------------------------------------------------------------------------
{
	return ( c*n_1 + ( 1.-c)*n_2 );
}

//-------------------------------------------------------------------------
/** @brief Calculate the viscosity of the Hershel-Bulkley law regularized using
 * the Simple approach.
 * 
 * \f[
 * \mu(c, \gamma(u)) := \kappa(c) \ (\gamma(u)+\epsilon)^{n(c)-1} + \frac{Bn(c)}{\gamma(u)+\epsilon}
 * \f]
 * @param c concentration
 * @param gammadot
 * @return vsicosity of the Simple model
 */
double 
FE_2PhaseViscHBParameter::Simple( double c, double gammadot ) const
//-------------------------------------------------------------------------
{
	return ( kappa(c)*pow(gammadot+reg_param, n(c)-1) + Bn(c)/(gammadot+reg_param) );
}
