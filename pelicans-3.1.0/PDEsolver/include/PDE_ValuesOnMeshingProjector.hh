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

#ifndef PDE_VALUES_ON_MESHING_PROJECTOR_HH
#define PDE_VALUES_ON_MESHING_PROJECTOR_HH

#include <PDE_ProjectorForDOFsSetting.hh>

class GE_Mpolyhedron ;

class PDE_ValuesOnMeshing ;

class PEL_EXPORT PDE_ValuesOnMeshingProjector : public PDE_ProjectorForDOFsSetting
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_ValuesOnMeshingProjector* create(
                          PEL_Object* a_owner,
                          PDE_DiscreteField* a_field,
                          size_t a_field_level,
                          PDE_ValuesOnMeshing* a_field_initializer,
                          PDE_DomainAndFields const* a_dom,
                          PEL_ModuleExplorer* a_exp ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      PDE_ValuesOnMeshingProjector( void ) ;
     ~PDE_ValuesOnMeshingProjector( void ) ;
      PDE_ValuesOnMeshingProjector(
                     PDE_ValuesOnMeshingProjector const& other ) ;
      PDE_ValuesOnMeshingProjector& operator=(
                     PDE_ValuesOnMeshingProjector const& other ) ;

      PDE_ValuesOnMeshingProjector(
                          PEL_Object* a_owner,
                          PDE_DiscreteField* a_field,
                          size_t a_field_level,
                          PDE_ValuesOnMeshing* a_field_initializer,
                          PDE_DomainAndFields const* a_dom,
                          PEL_ModuleExplorer* a_exp  ) ;
      
   //-- Projection

      virtual void compute_value_at_IP( PDE_LocalFEcell const* fe,
                                        doubleVector& result ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      mutable PDE_ValuesOnMeshing* INITIALIZER ;
      mutable GE_Mpolyhedron const* POLY ;
} ;

#endif
