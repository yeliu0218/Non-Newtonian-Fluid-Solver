/*
 *  Copyright :
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#ifndef FE_ParisaParameter_HH
#define FE_ParisaParameter_HH


#include <FE_Parameter.hh>

class GE_Vector ;
class PDE_DiscreteField ;
class FE;


/**
@brief PDE_Parameter derivative to calculate the convec cobination of two quantities.

Let \f$ c_1 \f$, \f$ c_2 \f$ and \f$ \lambda \f$  be FE_Parameters, and \f$ c_s \f$, \f$ \mbox{pow}\f$ a
scalar.
The values of these parameters are:
\f[
	value(c):= \left[ (c_1)^{\text{pow}}\ \lambda +(1-\lambda) (c_2)^{\text{pow}} \right]^{1/\text{pow}} \ c_s
\f]

Corresponding HDS: The module param#c_pow is optional, if the module is missing \f$ \mbox{pow} = 1.0\f$
\code
MODULE FE_Parameter#Viscosity
	concrete_name = "FE_ParisaParameter"
	name = "Viscosity"
	MODULE list_of_parameters
		MODULE param#c_1
			type = "already_defined"
			name = "Viscosity1"
		END MODULE param#c_1
		MODULE param#c_2
			type = "already_defined"
			name = "Viscosity2"
		END MODULE param#c_2
		MODULE param#c_3
			type = "already_defined"
			name = "Viscosity3"
		END MODULE param#c_3
		MODULE param#c_s
			type = "to_be_defined"
			MODULE FE_Parameter
				concrete_name = "FE_UniformParameter"
				name = "c_s"
				value = < 1. >
			END MODULE FE_Parameter
		END MODULE param#c_s
		MODULE param#lambda
			type = "to_be_defined"
			MODULE FE_Parameter
				concrete_name = "FE_FieldParameter"
				name = "lambda"
				field_level = 0
				field_name = "CC"
			END MODULE FE_Parameter
		END MODULE param#lambda
                (optional)MODULE param#c_pow
			type = "to_be_defined"
			MODULE FE_Parameter
				concrete_name = "FE_UniformParameter"
				name = "c_pow"
				value = < 1. >
			END MODULE FE_Parameter
		END MODULE param#c_pow
	END MODULE list_of_parameters
END MODULE FE_Parameter#Viscosity
\endcode

 * @todo
 * - Check what happens if \f$ lambda >1 \f$ or \f$ lambda < 0 \f$
 * - Update documentation


PUBLISHED
*/
class FE_ParisaParameter : public FE_Parameter
{
	public: //---------------------------------------------------------------
   //-- Instance delivery and initialization
		///@name Instance delivery
		//@{
		virtual void do_the_links( FE_SetOfParameters const* prms ) ;
		//@}

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

			~FE_ParisaParameter( void ) ;
			FE_ParisaParameter( FE_ParisaParameter const& other ) ;
			FE_ParisaParameter& operator=( FE_ParisaParameter const& other ) ;

			FE_ParisaParameter( PEL_Object* a_owner,
								PDE_DomainAndFields const* dom,
								PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

			FE_ParisaParameter( void ) ;

			virtual FE_ParisaParameter* create_replica(
					PEL_Object* a_owner,
					PDE_DomainAndFields const* dom,
					PEL_ModuleExplorer const* exp  ) const ;



   //-- Prerequisites to calculations

			virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
			virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

   //-- Class attributes

			static FE_ParisaParameter const* PROTOTYPE ;

   //-- Attributes
	  	///@name Parameters from module explorer:
	  	//@{
		typedef std::map< std::string, FE_Parameter* > FE_PARAMS_MAP;
		typedef std::map< std::string, std::string > FE_PARAMS_STRING_MAP;
		FE_PARAMS_MAP FE_PARAMS ;
		FE_PARAMS_STRING_MAP FE_PARAMS_REALNAMES ;
		FE_Parameter*  ParFE_C1;		///<@brief First parameter
		FE_Parameter* ParFE_C2;	
                FE_Parameter* ParFE_C3;		///<@brief Second parameter
		FE_Parameter* ParFE_Cs;		///<@brief Global scale
		FE_Parameter* ParFE_Lambda;	///<@brief Convex combination parameter
        size_t NBCS;	///<@brief Number of components
        bool CPOW;
        
        FE_Parameter* ParFE_Cpow;		///<@brief Global scale
	  //@}
} ;

#endif
