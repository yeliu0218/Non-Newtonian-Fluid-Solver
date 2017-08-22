#include <AS_Viscoplastic.hh>

#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalEquation.hh>

#include <FE.hh>
#include <GE_Point.hh>
#include <doubleVector.hh>
#include <math.h>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;


/** @brief Compute the strain rate tensorial field \f$ {\dot \gamma}_{ij}(u) \f$ based on the
      velocity field
 *
 * Calculate \f$ \dot \gamma_{ij}(u) = u_{i,j}+u_{j,i} \f$
 * @param[in]  cFE finite elements
 * @param[in]  UU velocity field
 * @param[out] gammadot strain rate tensor field */
/*----------------------------------------------------------------------*/
void AS_Viscoplastic::compute_strain_rate_tensor( PDE_LocalFEcell* cFE,
        PDE_DiscreteField const* UU,
        PDE_DiscreteField* gammadot)
{
    PEL_LABEL( "AS_Viscoplastic:: compute_strain_rate_tensor" ) ;

    size_t dimension = UU->nb_components(), DcompNumber=0 ;
    size_t node=0;

    if (dimension == 2 ) DcompNumber=4;
    else if (dimension == 3 ) DcompNumber=6;
    else
    {
        cout << "!!! Dimension should be either 2 or 3" << endl;
        cout << "    Dimension " << dimension << " is not allowed in " << endl;
        cout << "    AS_Viscoplastic:: compute_strain_rate_tensor !!!" << endl;
        PEL_Error::exit();
    }

    if (DcompNumber != gammadot->nb_components())
    {
        cout << "!!! Strain rate tensor D(u) has the wrong dimension in" << endl;
        cout << "    AS_Viscoplastic:: compute_strain_rate_tensor !!!" << endl;
        PEL_Error::exit();
    }

    doubleVector StrainRateAtNode(DcompNumber);
    pair<size_t,double> CompAndFactor(0,0.);

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( gammadot ); ++nodeIdx )
        {
            // Nullify the strain rate vector
            for (size_t i=0; i < DcompNumber; ++i) StrainRateAtNode(i)=0.;

            // Set the calculation point at the strain rate node location
            cFE->set_calculation_point( cFE->local_node_location( gammadot, nodeIdx) );

            // Compute the strain rate at the node
            for ( size_t compIdx = 0; compIdx < dimension; ++compIdx )
            {
                for ( size_t direcIdx = 0; direcIdx < dimension; ++direcIdx )
                {
                    CompAndFactor=strain_rate_component_number(dimension,compIdx,
                                  direcIdx);
                    StrainRateAtNode(CompAndFactor.first)+=CompAndFactor.second*
                                                           cFE->gradient_at_pt( UU, 0, direcIdx, compIdx);
                }
            }
            doubleVector const& dv = 
                              cFE->calculation_point()->coordinate_vector() ;
            StrainRateAtNode(3)= 2.0*cFE->value_at_pt( UU, 0, 0)/dv(0);
            //StrainRateAtNode(3)=-StrainRateAtNode(0)-StrainRateAtNode(1);
            // Get the global node
            node=cFE->global_node( gammadot, nodeIdx);

            // Transfer the strain rate at the node into the strain rate field D
            for (size_t i=0; i < DcompNumber; ++i)
                gammadot->set_DOF_value(0, node, StrainRateAtNode(i), i);
        }
    }

}




/** @brief Return the component number of strain rate tensor stored as a
      vector and the pre-factor to compute it based on velocity derivatives
      given the velocity component, the derivative direction and the
      space dimension
      @param[in] dimension space dimension
      @param[in] compIdx velocity component
      @param[in] direcIdx direction index
	  @return (component of velocity gradient, prefactor)
 */
/*------------------------------------------------------------------------*/
pair<size_t,double> AS_Viscoplastic::strain_rate_component_number(
    size_t const &dimension,
    size_t const& compIdx, size_t const &direcIdx )
{
    PEL_LABEL( "AS_Viscoplastic:: strain_rate_component_number" ) ;

    size_t ipj = compIdx + direcIdx;
    pair<size_t,double> result(0,0.);

    if (dimension == 2)
    {
        switch (ipj)
        {
        case 0 :
            result.first = 0;
            result.second = 2.;
            break;
        case 1 :
            result.first = 2;
            result.second = 1.;
            break;
        case 2 :
            result.first = 1;
            result.second = 2.;
            break;
        default :
			PEL::out() << "Error in AS_Viscoplastic:: strain_rate_component_number " << std::endl;
        	PEL_Error:: exit() ;
        }
    }
    else if (dimension == 3)
    {
        switch (ipj)
        {
        case 0 :
            result.first = 0;
            result.second = 2.;
            break;
        case 1 :
            result.first = 3;
            result.second = 1.;
            break;
        case 2 :
            if (compIdx == direcIdx)
            {
                result.first = 1;
                result.second = 2.;
            }
            else
            {
                result.first = 4;
                result.second = 1.0;
            }
        case 3 :
            result.first = 5;
            result.second = 1.0;
            break;
        case 4 :
            result.first = 2;
            result.second = 2.;
            break;
        }
    }
    else
    {
        cout << "!!! Dimension should be either 2 or 3" << endl;
        cout << "    Dimension " << dimension << " is not allowed in " << endl;
        cout << "    AS_Viscoplastic:: strain_rate_component_number !!!" << endl;
        PEL_Error::exit();
    }

    return( result );

}




/** @brief Update the Lagrange multiplier field using the strain rate
      tensors D(u) and d
      @param[in] cFE finite elements
      @param[out] LAMBDA Lagrange multiplier field (Stress)
      @param[in] gammadot ( \f$ {\dot \gamma}_{ij} \f$ ) strain rate tensor field based on the velocity derivatives
      @param[in] gamma  ( \f$ {\gamma}_{ij} \f$ ) strain rate tensor field as an independent variable
      @param[in] aug_param (r) Lagrange augmentation parameter */
//-------------------------------------------------------------
void AS_Viscoplastic::update_Lagrange_multiplier( PDE_LocalFEcell* cFE,
        PDE_DiscreteField* LAMBDA,
        PDE_DiscreteField const* gammadot,
        PDE_DiscreteField const* gamma,
        double const& aug_param)
{
    PEL_LABEL( "AS_Viscoplastic:: update_Lagrange_multiplier" ) ;

    double gammadot_at_node=0.,gamma_at_node=0.,lambda_at_node=0.;
    size_t DcompNumber=gammadot->nb_components();
    size_t node=0;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( gammadot ); ++nodeIdx )
        {
            // Get the global node
            node=cFE->global_node( LAMBDA, nodeIdx);

            // Update each component of Lambda
            for (size_t i=0; i < DcompNumber; ++i)
            {
                gammadot_at_node=gammadot->DOF_value(0,node,i);
                gamma_at_node=gamma->DOF_value(0,node,i);
                lambda_at_node=LAMBDA->DOF_value(0,node,i);
                LAMBDA->set_DOF_value(0, node,
                                      lambda_at_node+aug_param*(gammadot_at_node-gamma_at_node), i);
            }
        }
    }

}




/** @brief Update the strain rate tensor field gamma using the Lagrange
      multiplier tensor and the strain rate tensor D(u)
 * Calculate the following quantity:
 * \f{eqnarray}
 * 	\gamma =
 * 	\begin{cases}
 * 		\left(
 * 			1-\frac{Bn}{\| T + r {\dot \gamma} (u)\|}
 * 		\right)
 * 		\frac{ T + r \dot \gamma (u)}{\kappa+r}
 * 		& \text{if } \| T + r \dot \gamma (u)\| > Bn  \\
 * 		0 & \text{if }  \|  T + r \dot \gamma (u) \| \leq Bn
 *	\end{cases}
 * \f}
      @param[in]  cFE finite elements
      @param[in]  LAMBDA Lagrange multiplier field
      @param[in]  gammadot strain rate tensor field based on the velocity derivatives
      @param[out] gamma strain rate tensor field as an independent variable
      @param[in]  aug_param Lagrange augmentation parameter
      @param[in]  kappa Consistency parameter
      @param[in]  Bn Bingham number
      @param[in]  func_E used with ElectricField, f(|E|), otherwise should be f(|E|) = 1.0
	  @param[in]  t_it Time iterator ( needed to access the Bingham number parameter )
	  @todo Treat the location of the Bingham number parameter in a better way
*/
/*-----------------------------------------------------------*/
void AS_Viscoplastic::update_strain_rate_tensor_d( PDE_LocalFEcell* cFE,
        PDE_DiscreteField const* LAMBDA,
        PDE_DiscreteField const* gammadot,
        PDE_DiscreteField* gamma,
        double const& aug_param,
		FE_Parameter const* kappa,
		FE_Parameter const* Bn,
		double func_E,
		FE_TimeIterator const* t_it)
{
    PEL_LABEL( "AS_Viscoplastic:: update_strain_rate_tensor_d" ) ;

    size_t DcompNumber=gammadot->nb_components();
    size_t node=0;
    doubleVector StrainRateLagMultAtNode(DcompNumber);
    double norm=0.,pre_factor=0.;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( gammadot ); ++nodeIdx )
        {
            // Get the global node
            node=cFE->global_node( LAMBDA, nodeIdx);

			// Get parameter value for the Bingham number
			cFE->set_calculation_point(
									   cFE->local_node_location( LAMBDA, nodeIdx ) ) ;
			double BnValue = func_E * Bn->cell_value_at_pt( t_it, cFE, 0 ) ;
			double kappaValue = kappa->cell_value_at_pt( t_it, cFE, 0 ) ;

            // Store lambda+r.D(u) in a vector
            for (size_t i=0; i < DcompNumber; ++i)
                StrainRateLagMultAtNode(i)=LAMBDA->DOF_value(0,node,i)+aug_param
                                           *gammadot->DOF_value(0,node,i);

            // Compute the Euclidian norm of lambda+r.D(u) at the node
            norm=Euclidian_norm(StrainRateLagMultAtNode);

            // Update each component of gamma
            if (norm > BnValue)
            {
		   //pre_factor=( 1. - BnValue/norm )/( aug_param + kappaValue );
                   pre_factor=( 1. - BnValue/norm )/( aug_param );
                for (size_t i=0; i < DcompNumber; ++i)
                    gamma->set_DOF_value(0, node,
                                     pre_factor*StrainRateLagMultAtNode(i), i);
            }
            else
                for (size_t i=0; i < DcompNumber; ++i)
                    gamma->set_DOF_value(0, node, 0., i);
        }
    }
}

/** @brief Update the strain rate tensor field gamma using the Lagrange
      multiplier tensor and the strain rate tensor D(u)
 * Calculate the following quantity:
 *
 * Solve
 * \f$ \kappa |\gamma|^n + r |\gamma| = \| T + r \dot \gamma (u)\| - Bn\f$
 *
 * \f{eqnarray}
 * 	\gamma =
 * 	\begin{cases}
 * 		\|\gamma\|
 * 		\frac{ T + r \dot \gamma (u)}{\| T + r \dot \gamma (u)\|}
 * 		& \text{if } \| T + r \dot \gamma (u)\| > Bn  \\
 * 		0 & \text{if }  \|  T + r \dot \gamma (u) \| \leq Bn
 *	\end{cases}
 * \f}
      @param[in]  cFE finite elements
      @param[in]  LAMBDA Lagrange multiplier field
      @param[in]  gammadot strain rate tensor field based on the velocity derivatives
      @param[out] gamma strain rate tensor field as an independent variable
      @param[in]  aug_param Lagrange augmentation parameter
      @param[in]  kappa Consistency parameter
      @param[in]  Bn Bingham number
      @param[in]  HB_n Herschel-Bulkley coefficient n
      @param[in]  func_E used with ElectricField, f(|E|), otherwise should be f(|E|) = 1.0
	  @param[in]  t_it Time iterator ( needed to access the Bingham number parameter )
	  @todo Treat the location of the Bingham number parameter in a better way
*/
/*-----------------------------------------------------------*/
void AS_Viscoplastic::update_strain_rate_tensor_d_HB( PDE_LocalFEcell* cFE,
        PDE_DiscreteField const* LAMBDA,
        PDE_DiscreteField const* gammadot,
        PDE_DiscreteField* gamma,
        double const& aug_param,
		FE_Parameter const* kappa,
		FE_Parameter const* Bn,
		FE_Parameter const* HB_n,
		double func_E,
		FE_TimeIterator const* t_it)
{
    PEL_LABEL( "AS_Viscoplastic:: update_strain_rate_tensor_d" ) ;

    size_t DcompNumber=gammadot->nb_components();
    size_t node=0;
    doubleVector StrainRateLagMultAtNode(DcompNumber);
    double norm=0.,pre_factor=0.;

    double BnValue, kappaValue, HBnValue;
    double funcHB, d_funcHB ;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( gammadot ); ++nodeIdx )
        {
            // Get the global node
            node=cFE->global_node( LAMBDA, nodeIdx);

			// Get parameter value for the Bingham number, consistency and HB_n
			cFE->set_calculation_point(cFE->local_node_location( LAMBDA, nodeIdx ) ) ;
			BnValue = func_E * Bn->cell_value_at_pt( t_it, cFE, 0 ) ;
			kappaValue = kappa->cell_value_at_pt( t_it, cFE, 0 ) ;
			HBnValue = HB_n->cell_value_at_pt( t_it, cFE, 0 ) ;

            // Store lambda+r.D(u) in a vector
            for (size_t i=0; i < DcompNumber; ++i)
                StrainRateLagMultAtNode(i)=LAMBDA->DOF_value(0,node,i)+aug_param
                                           *gammadot->DOF_value(0,node,i);

            // Compute the Euclidian norm of lambda+r.D(u) at the node
            norm=Euclidian_norm(StrainRateLagMultAtNode);

            if (norm > BnValue) //PEL_Error::exit();
            {
            	// use x = |gamma|
            	double x,x1, dx=1.;
            	int iter=1, maxiter=100;

            	x1 = (norm - BnValue)/(kappaValue+aug_param) ; // solution of n = 1
/*				double x0, x2;
            	x0 = (norm - BnValue - kappaValue) / aug_param ; // solution of n = 0
            	x2 = -aug_param/(2*kappaValue) // solution of n = 2
					 + PEL::sqrt( (norm - BnValue) / kappaValue + pow(aug_param,2) / pow(2.*kappaValue,2) );

            	PEL::out() << "--------Take solution of n=0 -------- " << x0 << std::endl;
            	PEL::out() << "--------Take solution of n=1 -------- " << x1 << std::endl;
        		PEL::out() << "--------Take solution of n=2 -------- " << x2 << std::endl;*/

           		x = x1; // initial value
           		if(HBnValue != 1. ) // if HBnValue = 1 then x = x1 is solution
           		{
           			if(x>=1){
           				while( fabs(dx) > 1.E-6 && iter < maxiter)
           				{
           					funcHB = kappaValue * pow(x,HBnValue) + aug_param * x - norm + BnValue ;
           					d_funcHB = HBnValue * kappaValue * pow(x,HBnValue-1.) + aug_param ;
           					dx = funcHB/d_funcHB ;
           					x -=  dx;
           					iter++;
           					if(x<0) PEL_Error::object()->raise_plain( "Newton method failed: x<0!" ) ;
           				}
           				if( iter >= maxiter){
           					PEL::out() << "END x = " << x << " iterations: " << iter << std::endl;
           					PEL_Error::object()->raise_plain( "x>=1: Newton method failed: Too many iterations >100!" ) ;
           				}
           			}
           			else if (x>0. && x<1.){
           				double left=0., right=1.;
           				double midpoint=0.5;
           				double lo, hi;
           				double funcHB_left, funcHB_mid, funcHB_right;
           				funcHB_left = kappaValue * pow(left,HBnValue) + aug_param * left - norm + BnValue ;
           				funcHB_right= kappaValue * pow(right,HBnValue) + aug_param * right - norm + BnValue ;

           				if(funcHB_left*funcHB_right > 0.)
           					PEL_Error::object()->raise_plain( "0<x<1: Bi-section method does NOT work!" ) ;

           				if(funcHB_left <= 0.0 ){
           					lo = left;
           					hi = right;
           				} else {
           					lo = right;
           					hi = left ;
           				}
           				midpoint = lo + (hi - lo)/ 2. ;
           				while (PEL::abs(hi - lo) > 1.E-6 && iter < maxiter) {
           					funcHB_mid  = kappaValue * pow(midpoint,HBnValue) + aug_param * midpoint - norm + BnValue ;
           					if( funcHB_mid <= 0.){
           						lo = midpoint ;
           					} else {
           						hi = midpoint ;
           					}
           					midpoint = lo + (hi - lo)/ 2. ;
           					iter++;
           				}
           				x=midpoint;
           				if( iter >= maxiter){
           					PEL::out() << "END x = " << x << " iterations: " << iter << std::endl;
           					PEL_Error::object()->raise_plain( "0<x<1: Bi-section method failed: Too many iterations >100!" ) ;
           				}
           			}
           			else PEL_Error::object()->raise_plain( "x = 0") ;
				}

//            	PEL::out() << "END x = " << x << " iterations: " << iter << std::endl;
            	pre_factor = x/norm;
            	for (size_t i=0; i < DcompNumber; ++i)
            		gamma->set_DOF_value(0, node,
            				pre_factor*StrainRateLagMultAtNode(i), i);
            }
            else
                for (size_t i=0; i < DcompNumber; ++i)
                    gamma->set_DOF_value(0, node, 0., i);
        }
    }
}


/** @brief Return the second invariant  of a symmetric tensor stored as a
      vector
      @param tensor tensor \f$ T \f$
      @return \f$ \|T\|_{II}=\sqrt{\frac{1}{2}T_{ij}T_{ij}} \f$
*/
/*---------------------------------------------------------------------*/
double AS_Viscoplastic::Euclidian_norm( doubleVector const& tensor )
{
    PEL_LABEL( "AS_Viscoplastic:: Euclidian_norm" ) ;
	PEL_CHECK_PRE( tensor.size()==4 || tensor.size()==6 );
    double norm=0;

    if (tensor.size() == 4)
    {
        norm=0.5*(pow(tensor(0),2.)+pow(tensor(1),2.)+pow(tensor(3),2.))+pow(tensor(2),2.);
        norm=sqrt(norm);
    }
    else if (tensor.size() == 6)
    {
        norm=0.5*(pow(tensor(0),2.)+pow(tensor(1),2.)+pow(tensor(2),2.))
             +pow(tensor(3),2.)+pow(tensor(4),2.)+pow(tensor(5),2.);
        norm=sqrt(norm);
    }
    else
    {
        cout << "!!! Number of components of tensor should be either 4 or 6"
        << endl;
        cout << "    Number of components " << tensor.size() << " is not allowed ";
        cout << " in" << endl;
        cout << "    AS_Viscoplastic:: Euclidian_norm !!!" << endl;
        PEL_Error::exit();
    }

    return( norm );

}

/** @brief Compute the norm of the difference between the strain rate
      tensors D(u) and d
      @param cFE finite elements
      @param DUU strain rate tensor field based on the velocity derivatives
      @param d strain rate tensor field as an independent variable */
/*------------------------------------------------------------*/
double AS_Viscoplastic::compute_norm_Dminusd( PDE_LocalFEcell* cFE,
        PDE_DiscreteField const* DUU,
        PDE_DiscreteField const* d)
{
    PEL_LABEL( "AS_Viscoplastic:: compute_norm_Dminusd" ) ;


    double duu_at_node=0.,d_at_node=0.,norm=0.,diff,norm_collective;
    size_t DcompNumber=DUU->nb_components();
    size_t node=0, NFE=0;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        NFE++;
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( DUU ); ++nodeIdx )
        {
            // Get the global node
            node=cFE->global_node( DUU, nodeIdx);

            // Update each component of Lambda
            for (size_t i=0; i < DcompNumber; ++i)
            {
                duu_at_node=DUU->DOF_value(0,node,i);
                d_at_node=d->DOF_value(0,node,i);
                diff = duu_at_node-d_at_node;
                if (fabs(duu_at_node)>1e-2) diff/=fabs(duu_at_node);
                norm += pow( diff, 2. );
            }
        }
    }

    PEL_Communicator const* pelCOMM=PEL_Exec::communicator();
    norm_collective = pelCOMM->sum( norm );

    norm_collective=sqrt(norm_collective/double(NFE));

    return( norm_collective );

}
