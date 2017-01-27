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

#ifndef PDE_INTERFACE_AND_FIELDS_HH
#define PDE_INTERFACE_AND_FIELDS_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;

class GE_SetOfPoints ;

class PDE_DomainAndFields ;
class PDE_InterfaceBuilder ;
class PDE_LocalFEmortarSide ;
class PDE_SetOfDiscreteFields ;

class PEL_EXPORT PDE_InterfaceAndFields : public PEL_Object
{
   public: //-------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_InterfaceAndFields* create( PEL_Object* a_owner,
                                             PEL_ModuleExplorer* exp ) ;
      
   //-- Characteristics

      // name 
      std::string const& name( void ) const ;
      
      // 'i_adj'-th adjacent domain
      PDE_DomainAndFields const* adjacent_domain( size_t i_adj ) const ;

   //-- Problem Definition

      PDE_SetOfDiscreteFields* set_of_discrete_fields( void ) const ;
      
      GE_SetOfPoints const* set_of_vertices( void ) const ;
      
   //-- Discretization

      PDE_LocalFEmortarSide* create_LocalFEmortarSide( 
                                                 PEL_Object* a_owner ) const ;

   //-- Persistence
      
      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;

   protected: //----------------------------------------------------------

   private: //------------------------------------------------------------

      PDE_InterfaceAndFields( void ) ;
     ~PDE_InterfaceAndFields( void ) ;
      PDE_InterfaceAndFields( PDE_InterfaceAndFields const& other ) ;
      PDE_InterfaceAndFields& operator=( 
                              PDE_InterfaceAndFields const& other ) ;
     
      PDE_InterfaceAndFields( PEL_Object* a_owner, 
                              PEL_ModuleExplorer* exp ) ;

   //-- Attributes

      std::string NAME ;
      PDE_DomainAndFields const* DOM0 ;
      PDE_DomainAndFields const* DOM1 ;
      PDE_InterfaceBuilder const* BUILDER ;
} ;

#endif
