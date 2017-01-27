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

#ifndef FG_STEFAN_1_VELOCITY_HH
#define FG_STEFAN_1_VELOCITY_HH

#include <FE_Parameter.hh>

class GE_Vector ;
class PDE_DiscreteField ;

/*
PUBLISHED
*/

class AP_Stefan1Velocity : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual void do_the_links( FE_SetOfParameters const* prms ) ;

   //-- Instance caracteristics

      virtual size_t nb_components( void ) const ;

   //-- Values on bounds

      virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

     ~AP_Stefan1Velocity( void ) ;
      AP_Stefan1Velocity( AP_Stefan1Velocity const& other ) ;
      AP_Stefan1Velocity& operator=( AP_Stefan1Velocity const& other ) ;

      AP_Stefan1Velocity( PEL_Object* a_owner,
                          PDE_DomainAndFields const* dom,
                          PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

      AP_Stefan1Velocity( void ) ;

      virtual AP_Stefan1Velocity* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

   //-- Class attributes

      static AP_Stefan1Velocity const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField const* TT ;
      size_t L_TT ;
      PEL_ModuleExplorer const* EXP ;
      FE_Parameter* KAPPA ;
      double S_DENS ;
      double LHEAT ;
      FE_Parameter* FLUX_EXT ;
      size_t NB_DIMS ;
      double TMEL ;
      bool CHECK_BOUNDING ;
      double TT_MIN ;
      double TT_MAX ;
} ;

#endif
