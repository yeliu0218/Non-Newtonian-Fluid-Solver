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

#ifndef PDE_BASIS_FUNCTION_HH
#define PDE_BASIS_FUNCTION_HH

#include <PEL_Object.hh>
#include <PDE_SetOfBasisFunctions.hh>

class PDE_BoundFE ;
class PDE_CellFE ;
class PDE_DiscreteField ;
class PDE_FaceFE ;
class PDE_MortarSideFE ;
class PDE_SetOfBasisFunctions ;

#include <vector>

class PEL_EXPORT PDE_BasisFunction : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Basic characteristics
      
      //Uniquely determine a basis functions in the data structure  
      //PDE_SetOfBasisFunctions (and equal to `PEL::bad_index()' 
      //if not included). This number is internally attributed 
      //when `PDE_SetOfBasisFunctions:: add' is called
      size_t id_number( void ) const ;
      
   //-- Activation and refinement

      virtual size_t refinement_level( void ) const = 0 ;

      bool is_active( void ) const ;

      void set_active( void ) ;

      void set_inactive( void ) ;

      bool is_refined( void ) const ;

      void set_refined( void ) ;

      void set_unrefined( void ) ;

   //-- Location

      virtual bool is_located_in_cell( PDE_CellFE const* a_cell ) const = 0 ;

      virtual bool is_located_on_bound( 
                                 PDE_BoundFE const* a_bound ) const = 0 ;

   //-- Discrete fields

      bool is_attached_to_valid_DOFs( void ) const ;

      void append_one_field( PDE_DiscreteField* ff, size_t n ) ;

      bool is_in_basis_function_set_of( PDE_DiscreteField const* ff ) const ;

      size_t node_of_DOF( PDE_DiscreteField const* ff ) const ;

      void start_field_iterator( void ) ;

      void go_next_field( void ) ;

      bool valid_field( void ) const ;

      PDE_DiscreteField* field( void ) const ;

      size_t node_of_field( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_BasisFunction( void ) ;

      PDE_BasisFunction( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
      
      bool is_located_in_cell_PRE( PDE_CellFE const* a_cell ) const ;

      bool is_located_on_bound_PRE( PDE_BoundFE const* a_bound ) const ;

   private: //----------------------------------------------------------

      PDE_BasisFunction( void ) ;
      PDE_BasisFunction( PDE_BasisFunction const& other ) ;
      PDE_BasisFunction& operator=( PDE_BasisFunction const& other ) ;

      void set_id_number( size_t id ) ;
      
      friend void PDE_SetOfBasisFunctions:: add( PDE_BasisFunctionCell* bf,
                                                 size_t e ) ; 

   //-- Attributes

      bool ACTIVE ;
      bool REFINED ;

      std::vector< size_t > FIELD_IDXS ;
      std::vector< PDE_DiscreteField* > FIELDS ;
      std::vector< size_t > FIELD_NODES ; 

      size_t I_FIELD ;
      size_t ID_NUMBER ;
} ;

#endif
