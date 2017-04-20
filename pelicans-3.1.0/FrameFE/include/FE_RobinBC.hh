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

#ifndef FE_ROBIN_BC_HH
#define FE_ROBIN_BC_HH

#include <FE_OneBCbuilder.hh>

#include <vector>

class FE_Parameter ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_RobinBC : public FE_OneBCbuilder
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual void read_boundary_condition( size_t idx,
		  		            PDE_DomainAndFields const* dom,
				            FE_SetOfParameters const* prms,
				            PEL_ModuleExplorer const* exp ) ;

   //-- Local discrete contribution of the boundary condition

      virtual void transfer_calculation_requirements( 
                                            PDE_LocalFEbound* fe ) const ;

      virtual void build( PDE_LocalEquation* leq,
		          PDE_LocalFEbound* fe,
		          FE_TimeIterator const* t_it,
                          GE_QRprovider const* qrp ) const ;


   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~FE_RobinBC( void ) ;
      FE_RobinBC( FE_RobinBC const& other ) ;
      FE_RobinBC& operator=( FE_RobinBC const& other ) ;

      FE_RobinBC( PEL_Object* a_owner,
	          PDE_DiscreteField const* ff,
	          GE_Color const* color ) ;

   //-- Plug in

      FE_RobinBC( void ) ;

      virtual FE_RobinBC* create_replica( PEL_Object* a_owner,
				          PDE_DiscreteField const* ff,
  				          GE_Color const* color ) const ;

   //-- Class attributes

      static FE_RobinBC const* PROTOTYPE ;

   //-- Attributes

      std::vector< FE_Parameter* > HH ;
      std::vector< FE_Parameter* > UU_OUT ;
} ;

#endif
