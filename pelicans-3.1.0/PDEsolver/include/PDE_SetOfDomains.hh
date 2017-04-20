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

#ifndef PDE_SET_OF_DOMAINS_HH
#define PDE_SET_OF_DOMAINS_HH

#include <PEL_Object.hh>

#include <size_t_vector.hh>

class PEL_ModuleExplorer ;
class PEL_Vector ;
class PDE_DomainAndFields ;
class PDE_InterfaceAndFields ;

class PEL_EXPORT PDE_SetOfDomains : public PEL_Object
{
   public: //-------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_SetOfDomains* create( PEL_Object* a_owner,
                                       PEL_ModuleExplorer* exp ) ;

   //-- Volumic domains

      // number of `PDE_DomainAndFields::' objects that are recorded
      size_t nb_domains( void ) const ;

      // `i'-th recorded `PDE_DomainAndFields::' object
      PDE_DomainAndFields* domain( size_t i ) const ;

      // recorded `PDE_DomainAndFields::' object of name `name'
      PDE_DomainAndFields* domain( std::string const& name ) const ;

   //-- Interfaces between volumic domains

      // number of `PDE_InterfaceAndFields::' objects that are recorded
      size_t nb_interfaces( void ) const ;

      // `i'-th recorded `PDE_InterfaceAndFields::' object
      PDE_InterfaceAndFields* interface( size_t i ) const ;

      // recorded `PDE_InterfaceAndFields::' object of name `name'
      PDE_InterfaceAndFields* interface( std::string const& name  ) const ;

      // index of the `i_adj'-th adjacent domain of the `i'-th recorded
      // `PDE_InterfaceAndFields::'  object
      size_t index_of_interface_adjacent_domain( size_t i, 
                                                 size_t i_adj ) const ;

   protected: //----------------------------------------------------------

   private: //------------------------------------------------------------

      PDE_SetOfDomains( void ) ;
     ~PDE_SetOfDomains( void ) ;
      PDE_SetOfDomains( PDE_SetOfDomains const& other ) ;
      PDE_SetOfDomains& operator=( PDE_SetOfDomains const& other ) ;

      PDE_SetOfDomains( PEL_Object* a_owner,
                        PEL_ModuleExplorer* exp ) ;

      void resolve_conformal_adjacencies( PEL_ModuleExplorer const* exp ) ;

   //-- Attributes

      size_t VERB ;
      PEL_Vector* DOMAINS ;
      PEL_Vector* INTERFACES ;
      size_t_vector I_DOM_0 ;
      size_t_vector I_DOM_1 ;
} ;
#endif
