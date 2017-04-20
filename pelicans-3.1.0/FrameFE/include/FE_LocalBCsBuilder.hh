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

#ifndef FE_LOCAL_BCS_BUILDER_HH
#define FE_LOCAL_BCS_BUILDER_HH

#include <PEL_Object.hh>

#include <map>
#include <string>

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class FE_TimeIterator ;
class stringVector ;

class FE_OneBCbuilder ;
class FE_Parameter ;
class FE_SetOfParameters ;

/*
Builders of local discrete systems 
   - due to a given set of boundary conditions ;
   - related to a given field.

The given set of boundary conditions is extracted from the
"MODULE boudary_conditions" of the Hierarchical Data System
attached to a PDE_DomainAndFields object.

PUBLISHED
*/

class PEL_EXPORT FE_LocalBCsBuilder : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      /* 
      Create, initialize and return an instance, attached to the field
      called `field_name'. The handled boundary conditions will be those
         - defined in the "MODULE boundary_conditions" of the data hierarchy
           attached to `dom' ;
         - such that the value of the data of keyword "type" is an item
           of `bc_types'.
      Each item of `bc_types' should be the name of a class derived
      from `FE_OneBCbuilder::'.
      */ 
      static FE_LocalBCsBuilder* create( PEL_Object* a_owner,
			                 PDE_DomainAndFields const* dom,
                                         std::string const& field_name,
		    	                 stringVector const& bc_types,
			                 FE_SetOfParameters const* prms ) ;

   //-- Characteristics

      // the attached field
      PDE_DiscreteField const* field( void ) const ;

   //-- Loop on recorded BCs

      // Notify to `fe' of the spatial derivatives of the fields basis 
      // functions that will be requested during the iterations over the
      // bounds performed in `::build_current_BC'.
      void transfer_calculation_requirements( PDE_LocalFEbound* fe ) const ;

      // Notify that the subsequent call to `::build_current_BC' will handle
      // the boundary conditions defined in the "MODULE boundary_conditions" 
      // (Data Structure attached to the argument "dom" of `::create')
      // with `bc_type' as data of keyword "type".
      void set_current_BC_type( std::string const& bc_type ) ;

      // Does the previous call to `::set_current_BC_type' allows a
      // subsequent call to `::build_current_BC' ?
      bool current_BC_type_is_ok( void ) const ;

      // Add to `leq' the various terms of the discrete problem local to the
      // current bound of `fe', related to the boundary condition
      // identified by the revious call to `::set_current_BC_type'.
      void build_current_BC( PDE_LocalEquation* leq,
 		             PDE_LocalFEbound* fe,
		             FE_TimeIterator const* t_it,
                             GE_QRprovider const* qrp ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      FE_LocalBCsBuilder( void ) ;
     ~FE_LocalBCsBuilder( void ) ;
      FE_LocalBCsBuilder( FE_LocalBCsBuilder const& other ) ;
      FE_LocalBCsBuilder& operator=( FE_LocalBCsBuilder const& other ) ;

      FE_LocalBCsBuilder( PEL_Object* a_owner,
			  PDE_DomainAndFields const* dom,
                          std::string const& field_name,
			  stringVector const& bc_types,
			  FE_SetOfParameters const* prms ) ;

   //-- Attributes

      PDE_DiscreteField const* FF ;
      std::map< std::string, FE_OneBCbuilder* > BUILDERS ;
      std::map< std::string, FE_OneBCbuilder* >::const_iterator IT_BC_to_build ;
} ;

#endif
