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

#ifndef PDE_INTERFACE_BUILDER_HH
#define PDE_INTERFACE_BUILDER_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_Vector ;

class GE_Mpolyhedron ;
class GE_Point ;
class GE_SetOfPoints ;

class PDE_DiscOnMeshFE ;
class PDE_DomainBuilder ;
class PDE_SetOfDiscreteFields ;

class PEL_EXPORT PDE_InterfaceBuilder : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_InterfaceBuilder* create( PEL_Object* a_owner,
                                      PEL_ModuleExplorer* exp,
                                      PDE_DomainBuilder const* builder_0,
                                      PDE_DomainBuilder const* builder_1 ) ;

   //-- Access

      size_t nb_space_dimensions( void ) const ;

      PDE_SetOfDiscreteFields* set_of_discrete_fields( void ) const ;

      GE_SetOfPoints* set_of_vertices( void ) const ;

      PEL_Vector const* mortar_sides( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_InterfaceBuilder( void ) ;
     ~PDE_InterfaceBuilder( void ) ;
      PDE_InterfaceBuilder( PDE_InterfaceBuilder const& other ) ;
      PDE_InterfaceBuilder const& operator=( 
                            PDE_InterfaceBuilder const& other ) ;

      PDE_InterfaceBuilder( PEL_Object* a_owner,
                            PEL_ModuleExplorer* exp,
                            PDE_DomainBuilder const* builder_0,
                            PDE_DomainBuilder const* builder_1 ) ;

      void build_disc_ref_poly( PEL_ModuleExplorer const* exp ) ;

      void build_meshing( void ) ;

      void build_one_field( PEL_ModuleExplorer const* exp ) ;

      void connect_with_domain( size_t domain_id, 
                                PEL_Vector const* domain_bounds ) ;

      static size_t intersection_nature( 
                               GE_Point const* A0, GE_Point const* A1, 
                               GE_Point const* B0, GE_Point const* B1,
                               double& alpha, double& beta ) ;

   //-- Attributes

      size_t DIM ;
      PEL_ModuleExplorer* EXP ;
      PDE_SetOfDiscreteFields* FIELDS ;
      GE_SetOfPoints* VERTICES ;
      PEL_Vector* SIDES ;
      std::string D0_NAME ;
      std::string D1_NAME ;
      PEL_Vector* DISCS ;
} ;

#endif
