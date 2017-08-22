/*
 *  Copyright : 
 *    "Institut de Radioprotection et de Sûreté Nucléaire - IRSN" (1995-2008)
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

#ifndef FE_ViscosityHBParameter_HH
#define FE_ViscosityHBParameter_HH


#include <FE_Parameter.hh>

class GE_Vector ;
class PDE_DiscreteField ;
class FE;


/**
@brief PDE_Parameter describing a regularised Hershel-Bulkley fluid for two phase flow.

The viscosity due to the Herschley Bulkley model can be desribed as follows
 * \f[
 * 	 \mu(c) = \kappa \ \gamma^{n-1}(u) + \frac{Bn}{\gamma(u)}
 * \f]
 * where \f$ Bn, \kappa \f$ and \f$ n \f$ are the nondimensional parameters, realised as
 * FE_Parameter objects


Corresponding HDS for the Simple regularisation:
\code
MODULE FE_Parameter#ViscHB1
		concrete_name = "FE_ViscosityHBParameter"
		name = "Viscosity1"
		MODULE list_of_parameters
			MODULE param#Bn
				type = "to_be_defined"
				MODULE FE_Parameter
					concrete_name = "FE_UniformParameter"
					name = "Bn"
					value = $DV_Bn1
				END MODULE FE_Parameter
			END MODULE param#Bn
			MODULE param#kappa
				type = "to_be_defined"
				MODULE FE_Parameter
					concrete_name = "FE_UniformParameter"
					name = "kappa"
					value = $DV_kappa1
				END MODULE FE_Parameter
			END MODULE param#kappa
			MODULE param#n
				type = "to_be_defined"
				MODULE FE_Parameter
					concrete_name = "FE_UniformParameter"
					name = "n"
					value = $DV_n
				END MODULE FE_Parameter
			END MODULE param#n
		END MODULE list_of_parameters	
		// Link to velocity
		vel_name = "velocity" 	// Link to the velocity field
		vel_level = 0  	// Level of the velocity field 
		// Regularisation methods:
		reg = "Simple"
		reg_param = 1.E-3
	END MODULE FE_Parameter#ViscHB1
\endcode

PUBLISHED
 */
class FE_ViscosityHBParameter : public FE_Parameter
{
	public: //---------------------------------------------------------------

	//-- Instance delivery
		
		virtual void do_the_links( FE_SetOfParameters const* prms ) ;
		
   //-- Instance caracteristics

		virtual size_t nb_components( void ) const ;

   //-- Values on cells
      
		virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
										 PDE_LocalFEcell const* fe,
		   								 size_t ic = 0 ) const ;
      
		virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
										 PDE_LocalFEcell const* fe,
		   								 size_t ic = 0 ) const ;
      
   //-- Values on bounds
      
		virtual double bound_value_at_pt( FE_TimeIterator const* t_it,
										  PDE_LocalFEbound const* fe,
										  size_t ic = 0 ) const ;

		virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
										  PDE_LocalFEbound const* fe,
										  size_t ic = 0 ) const ;
      
	protected: //------------------------------------------------------------

	private: //--------------------------------------------------------------

		~FE_ViscosityHBParameter( void ) ;
		FE_ViscosityHBParameter( FE_ViscosityHBParameter const& other ) ;
		FE_ViscosityHBParameter& operator=( FE_ViscosityHBParameter const& other ) ;

		FE_ViscosityHBParameter( PEL_Object* a_owner,
								 PDE_DomainAndFields const* dom,
								 PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

		FE_ViscosityHBParameter( void ) ;

		virtual FE_ViscosityHBParameter* create_replica( 
				PEL_Object* a_owner,
				PDE_DomainAndFields const* dom,
				PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

		virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
		virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

   //-- Model implementation
		///@name Model implementation:
		//@{
		double Simple ( double gammadot, double Bn = 0, double kappa = 1, double n = 1 ) const;
		//@}
	
	
   //-- Class attributes

		static FE_ViscosityHBParameter const* PROTOTYPE ;

   //-- Attributes
	  	///@name Relations to fileds and parameters:
	  	//@{
		PDE_DiscreteField const* FU ; 	///< @brief Link to velocity field.
		size_t L_FU ; 					///< @brief Storage Level
		size_t NBCS_FU ; 				///< @brief number of components of the velocity field.
	  //@}
	  
		///@name Parameters from module explorer:
		//@{
		typedef std::map< std::string, FE_Parameter* > FE_PARAMS_MAP;
		typedef std::map< std::string, std::string > FE_PARAMS_STRING_MAP;
		FE_PARAMS_MAP FE_PARAMS ;
		FE_PARAMS_STRING_MAP FE_PARAMS_REALNAMES ;
		
		FE_Parameter* ParBn;	///<@brief Bingham number
		FE_Parameter* ParKappa;	///<@brief Consistency
		FE_Parameter* ParN;		///<@brief Power law exponent
		
		std::string reg;	///<@brief Regularisation Method
		double reg_param;	///<@brief Regularisation parameter
		//@}
		size_t NBCS;	///<@brief Number of components
} ;

#endif
