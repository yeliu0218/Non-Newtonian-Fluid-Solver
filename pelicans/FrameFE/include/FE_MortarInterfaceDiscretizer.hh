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

#ifndef FE_MORTAR_INTERFACE_DISCRETIZER_HH
#define FE_MORTAR_INTERFACE_DISCRETIZER_HH

#include <PEL_Object.hh>

#include <vector>

class PEL_ModuleExplorer ;
class PEL_Vector ;

class LA_Matrix ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_System ;
class PDE_DiscreteField ;
class PDE_InterfaceAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEmortarSide ;
class PDE_SetOfDomains ;
class PDE_SystemNumbering ;

class FE_OneStepIterationOpen ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_MortarInterfaceDiscretizer : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_MortarInterfaceDiscretizer* create( 
                         PEL_Object* a_owner,
                         PDE_InterfaceAndFields const* interf,
                         PEL_Vector const* dom_discs,
                         PEL_ModuleExplorer* exp ) ;

   //-- Instance characteristics
      
      PDE_InterfaceAndFields const* interface( void ) const ;
      
   //-- Unknowns

      virtual size_t nb_unknowns( void ) const ;

      virtual PDE_DiscreteField* field( size_t i_unk ) const ;
      
   //-- Discretization
         
      virtual PDE_LinkDOF2Unknown* 
              create_link_DOF_2_unknown( PEL_Object* a_owner, 
                                         size_t i_unk ) const ; 
      
      void assemble_contribution( LA_Matrix* matrix,
                                  LA_Vector* vector,
                                  PDE_SystemNumbering const* nmb,
                                  size_t interf_shift, 
                                  size_t domain_0_shift, 
                                  size_t domain_1_shift ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FE_MortarInterfaceDiscretizer( void ) ;
     ~FE_MortarInterfaceDiscretizer( void ) ;
      FE_MortarInterfaceDiscretizer( 
                               FE_MortarInterfaceDiscretizer const& other ) ;
      FE_MortarInterfaceDiscretizer& operator=( 
                               FE_MortarInterfaceDiscretizer const& other ) ;

      FE_MortarInterfaceDiscretizer( PEL_Object* a_owner,
                               PDE_InterfaceAndFields const* interf,
                               PEL_Vector const* dom_discs,
                               PEL_ModuleExplorer* exp ) ;

   //-- Internals

      static void check_field_consistency( PDE_DiscreteField const* ff,
                                           PDE_DiscreteField const* dom_ff,
                                           PEL_ModuleExplorer const* exp ) ;
                                          
   //-- Attributes

      PDE_InterfaceAndFields const* INTERF ;
      
      std::vector< PDE_DiscreteField* > LLs ;
      std::vector< PDE_LinkDOF2Unknown* > LL_links ;

      std::vector< size_t > IDs_0 ;
      std::vector< PDE_DiscreteField const* > UUs_0 ;

      std::vector< size_t > IDs_1 ;
      std::vector< PDE_DiscreteField const* > UUs_1 ;

      PDE_LocalEquation* ELEMENT_EQ ;
      PDE_LocalEquation* T_ELEMENT_EQ ;
      PDE_LocalFEmortarSide* msFE ;
      GE_QRprovider const* QRP ;
} ;

#endif
