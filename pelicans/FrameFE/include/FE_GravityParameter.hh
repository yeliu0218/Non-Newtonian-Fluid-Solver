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

#ifndef FE_GRAVITY_PARAMETER_HH
#define FE_GRAVITY_PARAMETER_HH

/*

Parameters the value of which is (`rho'-`rho_shift')*`gravity'.

`rho' : is the fluid density (`FE_parameter::') (kg/m^3)
`rho_shift' : is a constant value (kg/m^3)
`gravity' : is the constant specific gravity (m/s^2)

PUBLISHED
*/

#include <FE_Parameter.hh>

#include <doubleVector.hh>

class PEL_EXPORT FE_GravityParameter : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual void do_the_links( FE_SetOfParameters const* prms ) ;

   //-- Instance characteristics

      virtual size_t nb_components( void ) const ;

   //-- Values on cells

      virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;
  
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

     ~FE_GravityParameter( void ) ;
      FE_GravityParameter( FE_GravityParameter const& other ) ;
      FE_GravityParameter& operator=( FE_GravityParameter const& other ) ;

      FE_GravityParameter( PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

      FE_GravityParameter( void ) ;

      virtual FE_GravityParameter* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
      virtual void prepare_for_value_on_sides( PDE_CursorFEside* fe ) const ;

   //-- Class attributes

      static FE_GravityParameter const* PROTOTYPE ;

   //-- Attributes

      doubleVector GRAV ;
      double RHOS ;
      FE_Parameter* RHO ;
      PEL_ModuleExplorer const* EXP ;
} ;

#endif
