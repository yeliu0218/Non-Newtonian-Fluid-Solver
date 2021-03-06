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

#ifndef FE_SET_OF_PARAMETERS_HH
#define FE_SET_OF_PARAMETERS_HH

#include <PEL_Object.hh>

#include <map>
#include <string>

class PEL_ModuleExplorer ;

class PDE_DomainAndFields ;
class PDE_SetOfDomains ;

class FE_Parameter ;

/*
Sets of `FE_Parameter::' instances, uniquely determined by their name.  

PUBLISHED
*/

class PEL_EXPORT FE_SetOfParameters : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_SetOfParameters* create( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer* exp ) ;

      static FE_SetOfParameters* create( PEL_Object* a_owner,
                                         PDE_SetOfDomains const* sdoms,
                                         PEL_ModuleExplorer* exp ) ;
   //-- Access

      size_t nb_parameters( void ) const ;

      bool has( std::string const& name ) const ;

      FE_Parameter* item( std::string const& name ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------
      
      FE_SetOfParameters( void ) ;
     ~FE_SetOfParameters( void ) ;
      FE_SetOfParameters( FE_SetOfParameters const& other ) ;
      FE_SetOfParameters& operator=( FE_SetOfParameters const& other ) ;

      FE_SetOfParameters( PEL_Object* a_owner,
                          PDE_DomainAndFields const* dom,
                          PEL_ModuleExplorer* exp ) ;

      FE_SetOfParameters( PEL_Object* a_owner,
                          PDE_SetOfDomains const* sdoms,
                          PEL_ModuleExplorer* exp ) ;

   //-- Internals

      void build_uniform_parameters( PEL_ModuleExplorer* exp ) ;

      void add_one_parameter( FE_Parameter* prm ) ;
      
   //-- Attributes      

      std::map< std::string, FE_Parameter* > MAP ;
} ;

#endif
