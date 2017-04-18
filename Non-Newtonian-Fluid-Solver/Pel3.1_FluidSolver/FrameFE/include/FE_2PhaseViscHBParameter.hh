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

#ifndef FE_2PhaseViscHB_PARAMETER_HH
#define FE_2PhaseViscHB_PARAMETER_HH


#include <FE_Parameter.hh>

class GE_Vector ;
class PDE_DiscreteField ;
class FE;


/**
@brief PDE_Parameter describing a regularised Hershel-Bulkley fluid for two phase flow.

Let c be a given PDE_DiscreteField::' at a given level.
The corresponding viscosity of the Hershcel Bulkley fluid can be written as
\f[
	\mu(c, \gamma(u)) = \kappa(c) \ \gamma^{n(c)-1}(u) + \frac{Bn(c)}{\gamma(u)}
\f]
where \f$ Bn(c), \kappa(c) \f$ and \f$ n(c) \f$ depend on the concentraion in a convex manner.


Corresponding HDS:
\code
MODULE FE_Parameter#ViscHB
  concrete_name = "FE_2PhaseViscHBParameter"
  name = "viscosity"
  conc_name  = "c"
  conc_level = 0
  vel_name = "u"
  vel_level = 0
  Bn_1 = 5.
  Bn_2 = 3.
  kappa_r = 3.
  n_1 = 2.;
  n_2 = 3.;
  reg = "Simple"
  reg_param = 1e-3
END MODULE FE_Parameter#ViscHB
\endcode

PUBLISHED
*/
class FE_2PhaseViscHBParameter : public FE_Parameter
{
   public: //---------------------------------------------------------------

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

     ~FE_2PhaseViscHBParameter( void ) ;
      FE_2PhaseViscHBParameter( FE_2PhaseViscHBParameter const& other ) ;
      FE_2PhaseViscHBParameter& operator=( FE_2PhaseViscHBParameter const& other ) ;

      FE_2PhaseViscHBParameter( PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

      FE_2PhaseViscHBParameter( void ) ;

      virtual FE_2PhaseViscHBParameter* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

   //-- Model implementation
	///@name Model implementation:
	//@{
	double Bn( double c ) const;
	double kappa ( double c ) const;
	double n ( double c ) const;
	double Simple ( double c, double gammadot) const;
	//@}
	
	
   //-- Class attributes

      static FE_2PhaseViscHBParameter const* PROTOTYPE ;

   //-- Attributes
	  ///@name Relations to fileds and parameters:
	  //@{
      PDE_DiscreteField const* FC ; ///< @brief Link to concentration field.
      size_t L_FC ; ///< @brief Storage level
      size_t NBCS_FC ; ///< @brief number of components of the conentration field.
	  
	  PDE_DiscreteField const* FU ; ///< @brief Link to velocity field.
	  size_t L_FU ; ///< @brief Storage Level
	  size_t NBCS_FU ; ///< @brief number of components of the velocity field.
	  //@}
	  
	  ///@name Parameters from module explorer:
	  //@{
	  double Bn_1;///< @brief Global parameter
	  double Bn_2;
	  double n_1;
	  double n_2;
      double kappa_r;///< @brief Local in convex combination 
	  std::string reg;
	  double reg_param;
	  //@}
} ;

#endif
