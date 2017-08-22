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

#ifndef FE_ONE_BC_BUILDER_HH
#define FE_ONE_BC_BUILDER_HH

#include <PEL_Object.hh>

class PEL_ObjectRegister ;
class PEL_ModuleExplorer ;

class GE_Color ;
class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class FE_TimeIterator ;

class FE_SetOfParameters ;

#include <map>

/*
Builder of local discrete systems
   - do to a given type of boundary conditions ;
   - related to a given field.

FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT FE_OneBCbuilder : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance of a concrete subclass
      // of registration name `type', associated to `ff'.
      static FE_OneBCbuilder* make( PEL_Object* a_owner,
                                    std::string const& type,
				    PDE_DiscreteField const* ff,
                                    GE_Color const* color,
				    PDE_DomainAndFields const* dom,
				    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer const* exp ) ;

      void extend( GE_Color const* color,
		   PDE_DomainAndFields const* dom,
		   FE_SetOfParameters const* prms,
		   PEL_ModuleExplorer const* exp ) ;

      virtual void read_boundary_condition( 
                                         size_t idx,
		  		         PDE_DomainAndFields const* dom,
				         FE_SetOfParameters const* prms,
				         PEL_ModuleExplorer const* exp ) = 0 ;

   //-- Characteristics

      // the attached field
      PDE_DiscreteField const* field( void ) const ;

   //-- Local discrete contribution of the boundary condition

      // Notify to `fe' of the spatial derivatives of the fields basis 
      // functions that will be requested during the iterations over the
      // bounds performed in `::build'.
      virtual void transfer_calculation_requirements( 
                                            PDE_LocalFEbound* fe ) const = 0 ;

      virtual void build( PDE_LocalEquation* leq,
 		          PDE_LocalFEbound* fe,
		          FE_TimeIterator const* t_it,
                          GE_QRprovider const* qrp ) const = 0 ;


   protected: //--------------------------------------------------------------

   //-- Plug in

      virtual ~FE_OneBCbuilder( void ) ;

      // Registration of the concrete subclass under the name `a_type'.
      FE_OneBCbuilder( std::string const& a_type ) ;

      FE_OneBCbuilder( PEL_Object* a_owner, 
		       PDE_DiscreteField const* ff,
                       GE_Color const* color ) ;

      virtual FE_OneBCbuilder* create_replica(
                                           PEL_Object* a_owner,
				           PDE_DiscreteField const* ff,
  				           GE_Color const* color ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Colors

      size_t index_of_color( GE_Color const* color ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool read_boundary_condition_PRE( 
                                    size_t idx,
  		                    PDE_DomainAndFields const* dom,
				    FE_SetOfParameters const* prms,
				    PEL_ModuleExplorer const* exp ) const ;

      virtual bool transfer_calculation_requirements_PRE( 
                                    PDE_LocalFEbound* fe ) const ;

      virtual bool build_PRE( PDE_LocalEquation* leq,
 		              PDE_LocalFEbound* fe,
		              FE_TimeIterator const* t_it,
                              GE_QRprovider const* qrp ) const  ;

      virtual bool create_replica_PRE( PEL_Object* a_owner,
				       PDE_DiscreteField const* ff,
  				       GE_Color const* color ) const ;

      virtual bool create_replica_POST( FE_OneBCbuilder const* result,
                                        PEL_Object* a_owner,
				        PDE_DiscreteField const* ff,
  				        GE_Color const* color ) const ;

   private: //----------------------------------------------------------------

      FE_OneBCbuilder( void ) ;
      FE_OneBCbuilder( FE_OneBCbuilder const& other ) ;
      FE_OneBCbuilder& operator=( FE_OneBCbuilder const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      bool IS_PROTO ;
      PDE_DiscreteField const* FF ;
      size_t IDX ;
      std::map< GE_Color const*, size_t > COLOR_2_BCindex ;
} ;

#endif
