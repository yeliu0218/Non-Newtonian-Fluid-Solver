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

#ifndef PDE_PROJECTOR_FOR_DOF_S_SETTING_HH
#define PDE_PROJECTOR_FOR_DOF_S_SETTING_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class doubleVector ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;
class PDE_SystemNumbering ;

/*
Objects that set the values of the DOFs of a `::PDE_DiscreteField' instance 
by performing a finite element L2 projection of a function whose
nature and definition is specific to each concrete subclass.

PUBLISHED
*/

class PEL_EXPORT PDE_ProjectorForDOFsSetting : public PEL_Object
{
   public: //----------------------------------------------------------------

   //-- Characteristics

      // instance whose DOFs values will be set
      PDE_DiscreteField const* field( void ) const ;

   //-- Projection
      
      // Perform the L2 projection, then set the values of the DOFs 
      // of `::field'.
      void project_and_update_field( void ) ;

   protected: //-------------------------------------------------------------

      PDE_ProjectorForDOFsSetting( PEL_Object* a_owner,
                                   PDE_DiscreteField* a_field,
                                   size_t a_field_level,
                                   PDE_DomainAndFields const* a_dom,
                                   PEL_ModuleExplorer* a_exp ) ;
      
      virtual ~PDE_ProjectorForDOFsSetting( void ) ;

   //-- Projection

      void add_field_requirement_on_cells( PDE_DiscreteField const* ff,
                                           int derivation_order ) ;

      virtual void compute_value_at_IP( PDE_LocalFEcell const* fe,
                                        doubleVector& result ) const = 0 ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool compute_value_at_IP_PRE( PDE_LocalFEcell const* fe,
                                            doubleVector& result ) const ;

      virtual bool invariant( void ) const ;

   private: //---------------------------------------------------------------

      PDE_ProjectorForDOFsSetting( void ) ;
      PDE_ProjectorForDOFsSetting( PDE_ProjectorForDOFsSetting const& other ) ;
      PDE_ProjectorForDOFsSetting& operator=( 
                                   PDE_ProjectorForDOFsSetting const& other ) ;

   //-- Attributes

      PDE_DiscreteField* const FIELD ;
      size_t const FIELD_LEVEL ;

      GE_QRprovider const* QRP ;
      PDE_LocalEquation* ELEMENT_EQ ;
      PDE_LocalFEcell* cFE ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;

      LA_Solver* SOLVER ;
} ;

#endif
