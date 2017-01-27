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

#ifndef FE_2PhaseConvex_PARAMETER_HH
#define FE_2PhaseConvex_PARAMETER_HH


#include <FE_Parameter.hh>

class GE_Vector ;
class PDE_DiscreteField ;
class FE;


/**
@brief PDE_Parameter derivative to calculate the convec cobination of two quantities.

Let c (field_name in HDS) be a given PDE_DiscreteField::' at a given level.
The values of these parameters are:
\f[
	value(c):= \left[ c+(1-c) k_c \right] k_g
\f]
where \f$ k_g \f$ is a global constant and 
\f$ k_c \f$ is the convex combination constant.


Corresponding HDS:
\code
MODULE FE_Parameter#ConvexComb
  concrete_name = "FE_2PhaseConvexParameter"
  name = "Reynolds"
  field_name  = "c"
  field_level = 0
  global = 5.
  convex = 3.
END MODULE FE_Parameter#ConvexComb
\endcode

PUBLISHED
*/
class FE_2PhaseConvexParameter : public FE_Parameter
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

     ~FE_2PhaseConvexParameter( void ) ;
      FE_2PhaseConvexParameter( FE_2PhaseConvexParameter const& other ) ;
      FE_2PhaseConvexParameter& operator=( FE_2PhaseConvexParameter const& other ) ;

      FE_2PhaseConvexParameter( PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

      FE_2PhaseConvexParameter( void ) ;

      virtual FE_2PhaseConvexParameter* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

   //-- Class attributes

      static FE_2PhaseConvexParameter const* PROTOTYPE ;

   //-- Attributes
	  ///@name Relations to fileds and parameters:
      PDE_DiscreteField const* FF ; ///< @brief Link to concentration field.
      size_t L_FF ; ///< @brief active level of the discrete field.
      size_t NBCS ; ///< @brief number of components of the conentration field.
	  
	  ///@name Parameters from module explorer:
	  //@{
	  double global;///< @brief Global parameter
      double convex;///< @brief Local in convex combination 
	  //@}
} ;

#endif
