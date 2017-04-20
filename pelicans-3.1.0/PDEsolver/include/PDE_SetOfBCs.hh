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

#ifndef PDE_SET_OF_BCS_HH
#define PDE_SET_OF_BCS_HH

#include <PEL_Object.hh>
#include <stringVector.hh>

class PEL_List ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class GE_Color ;
class PDE_DiscreteField ;
class PDE_SetOfDiscreteFields ;

/*
Sets of BCs and/or macro_BCs.

Definition of BCs:
------------------

   Given a clusters of meshes identified by 
      - their color, and 
      - one or all components of a `PDE_DiscreteField::' object, 
   a BC (for Boundary Condition) is defined as a database of the 
   PELICANS  Hierarchical Data System, interrogeable thanks to a 
   `PEL_ModuleExplorer::' object. 

   The content of a BC, in terms of its data and their meaning, 
   is left to the user.

Definition of macro_BCs:
------------------------

   Given a clusters of meshes identified by their color, a macro_BC is 
   defined by an identifier. macro_BCs are used to define many related
   BCs by further specifying one or all components of 
   `PDE_DiscreteField::' objects together with additional specific data.

PUBLISHED
*/

class PEL_EXPORT PDE_SetOfBCs : public PEL_Object
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PDE_SetOfBCs const* create(
                              PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp,
                              PDE_SetOfDiscreteFields const* fields ) ;
   //-- Access
      
      // formal index that denotes all the components of a field
      static size_t const all_components ;

      // Is there a BC on the meshes identified by `color' for component `ic'
      // of `field' ? 
      bool has_BC( GE_Color const* color,
                   PDE_DiscreteField const* field,
                   size_t ic = all_components ) const ;

      // navigator on the data characterising the BC for component `ic'
      // of `field' on the meshes identified by `color'
      PEL_ModuleExplorer const* BC_explorer( 
                                         GE_Color const* color,
                                         PDE_DiscreteField const* field,
                                         size_t ic = all_components ) const ;

      // Is there a macro_BC on the meshes identified by `color'? 
      bool has_macro_BC( GE_Color const* color ) const ;

      // identifier of the macro_BC on the meshes identified by `color'
      std::string const& macro_BC( GE_Color const* color ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_SetOfBCs( void ) ;
     ~PDE_SetOfBCs( void ) ;
      PDE_SetOfBCs( PDE_SetOfBCs const& other ) ;
      PDE_SetOfBCs& operator=( PDE_SetOfBCs const& other ) ;

      PDE_SetOfBCs( PEL_Object* a_owner,
                    PEL_ModuleExplorer const* exp,
                    PDE_SetOfDiscreteFields const* fields ) ;

      void explore( PEL_ModuleExplorer* bcTree,
                    PDE_SetOfDiscreteFields const* fields,
                    GE_Color const* bcTree_color = 0 ) ;

      size_t macro_BC_index( GE_Color const* color ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      PEL_List* const BCS ;

      PEL_Vector* const MACRO_COLS ;
      stringVector MACRO_TYPE ;
} ;


#endif
