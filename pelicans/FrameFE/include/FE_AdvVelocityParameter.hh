/*
 *  Copyright 1995-2010 by IRSN
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

#ifndef FE_ADV_VELOCITY_PARAMETER_HH
#define FE_ADV_VELOCITY_PARAMETER_HH

/*

Advective velocity parameters.

In the next examples, we suppose that :
   0 : velocity current level (n+1)
   1 : velocity initial level (n)
   2 : velocity previous level (n-1)
   
Example 1 : implicit advective velocity

   MODULE PEL_Application
      MODULE PDE_DomainAndFields
         MODULE interior_fields
            MODULE velocity
               name = "velocity"
               storage_depth = 1
               ...
            END MODULE velocity
         END MODULE interior_fields
      END MODULE PDE_DomainAndFields
      MODULE PDE_SetOfParameters
         MODULE FE_Parameter#adv_velocity_param
            concrete_name = "FE_AdvVelocityParameter#order1"
            name = "advective_velocity"
            velocity_name = "velocity"
            velocity_level = 0
         END MODULE FE_Parameter#adv_velocity_param
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application
   
Example 2 : explicit advective velocity, order 1

   MODULE PEL_Application
      MODULE PDE_DomainAndFields
         MODULE interior_fields
            MODULE velocity
               name = "velocity"
               storage_depth = 2
               ...
            END MODULE velocity
         END MODULE interior_fields
      END MODULE PDE_DomainAndFields
      MODULE PDE_SetOfParameters
         MODULE FE_Parameter#adv_velocity_param
            concrete_name = "FE_AdvVelocityParameter#order1"
            name = "advective_velocity"
            velocity_name = "velocity"
            velocity_level = 1
         END MODULE FE_Parameter#adv_velocity_param
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application
   
Example 3 : explicit advective velocity, order 2

   MODULE PEL_Application
      MODULE PDE_DomainAndFields
         MODULE interior_fields
            MODULE velocity
               name = "velocity"
               storage_depth = 3
               ...
            END MODULE velocity
         END MODULE interior_fields
      END MODULE PDE_DomainAndFields
      MODULE PDE_SetOfParameters
         MODULE FE_Parameter#adv_velocity_param
            concrete_name = "FE_AdvVelocityParameter#order2"
            name = "advective_velocity"
            velocity_name = "velocity"
            initial_velocity_level = 1
            previous_velocity_level = 2 
         END MODULE FE_Parameter#adv_velocity_param
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application
   
Example 4 : explicit advective velocity, order 2 : user defined

   MODULE PEL_Application
      MODULE PDE_DomainAndFields
         MODULE interior_fields
            MODULE velocity
               name = "velocity"
               storage_depth = 3
               ...
            END MODULE velocity
         END MODULE interior_fields
      END MODULE PDE_DomainAndFields
      MODULE PDE_SetOfParameters
         MODULE FE_Parameter#adv_velocity_param
            concrete_name = "FE_AdvVelocityParameter#user_defined"
            name = "advective_velocity"
            velocity_name = "velocity"
            velocity_levels_table = < 1 2 >
            coefficients_table = < 2. -1. >
         END MODULE FE_Parameter#adv_velocity_param
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application

Example 5 : no advective velocity

   MODULE PEL_Application
      MODULE PDE_SetOfParameters
         advective_velocity = < 0. 0. >
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application
   
Example 6 : given advective velocity (function of the time and the coordinates)

   MODULE PEL_Application
      MODULE PDE_SetOfParameters
         MODULE FE_Parameter#advective_velocity
            concrete_name = "FE_SpaceTimeParameter"
            name = "advective_velocity"
            nb_components = 2
            $DS_x = component($DV_X,0)
            $DS_y = component($DV_X,1)
            value = vector( $DS_y, 1.0 - $DS_x )
         END MODULE FE_Parameter#advective_velocity
      END MODULE PDE_SetOfParameters
   END MODULE PEL_Application
PUBLISHED
*/

#include <FE_Parameter.hh>

#include <doubleVector.hh>
#include <size_t_vector.hh>

class PDE_DiscreteField ;

class PEL_EXPORT FE_AdvVelocityParameter : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance characteristics

      virtual size_t nb_components( void ) const ;

   //-- Values on cells

      virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;
      
      virtual double cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;
      
      virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;
      
   //-- Values on bounds

      virtual double bound_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;
      
      virtual double bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic = 0 ) const ;
      
      virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;

      virtual double bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic = 0 ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------
      
      enum FE_AdvVelocityLaw { order1, order2, user_defined } ;

   //-- Constructors

      FE_AdvVelocityParameter( void ) ;
     ~FE_AdvVelocityParameter( void ) ;
      FE_AdvVelocityParameter( FE_AdvVelocityParameter const& other ) ;
      FE_AdvVelocityParameter& operator=(
                               FE_AdvVelocityParameter const& other ) ;

      FE_AdvVelocityParameter( PEL_Object* a_owner,
                               FE_AdvVelocityLaw const id,
                               size_t_vector const& levels,
                               doubleVector const& coeffs,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in
      
      FE_AdvVelocityParameter( FE_AdvVelocityLaw const id,
                               std::string const& concrete_name ) ;

      virtual FE_AdvVelocityParameter* create_replica( 
                               PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp ) const ;
      
   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;

      virtual void prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const ;

      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

      virtual void prepare_for_gradient_on_bounds( 
                                                PDE_LocalFEbound* fe ) const ;
      
      virtual void prepare_for_value_on_sides( PDE_CursorFEside* fe ) const ;

      virtual void prepare_for_gradient_on_sides( PDE_CursorFEside* fe ) const ;

   //-- Class attributes
      
      static FE_AdvVelocityParameter const* PROTOTYPE_1 ;
      static FE_AdvVelocityParameter const* PROTOTYPE_2 ;     
      static FE_AdvVelocityParameter const* PROTOTYPE_U_DEF ;

   //-- Attributes
 
      FE_AdvVelocityLaw const ADV_VELO_LAW ;

      PDE_DiscreteField const* const VELOCITY ;
      size_t_vector const LEVELS ;
      doubleVector const COEFFS ;

} ;

#endif // FE_ADV_VELOCITY_PARAMETER_HH
