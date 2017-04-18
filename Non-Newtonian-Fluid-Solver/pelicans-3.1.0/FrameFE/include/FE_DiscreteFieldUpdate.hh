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

#ifndef FE_DISCRETE_FIELD_UPDATE_HH
#define FE_DISCRETE_FIELD_UPDATE_HH 

#include <FE_OneStepIteration.hh>

class FE_Parameter ;

class PDE_DiscreteField ;
class PDE_FieldComposition ;
class PDE_LocalFEcell ;

class PEL_Context ;
class PEL_Data ;

class doubleVector ;

/*
`FE_OneStepIteration::' objects performing the update of the
DOFs values `PDE_DiscreteField::' objects from data given in
a hierarchical data structure

The values should be defined as :

  - analytical values function of the time (variable $DS_T) and
    of the coordinates of the node of the field (variable $DV_X) :

    example :
       MODULE FE_OneStepIteration#set_with_analytical_value
          concrete_name = "FE_DiscreteFieldUpdate"
          field_name = "velocity"
          field_level = 0
          MODULE DOFs_values
             type = "from_analytic"
             value = vector( $DS_T*component($DV_X,0), $DS_T*component($DV_X,1) )
          END MODULE DOFs_values
       END MODULE FE_OneStepIteration#set_with_analytical_value

  - a `PDE_FieldComposition::' object defined in the module 
    field_compositions of the module PDE_DomainAndFields :

    remark :
       all the values of the fields needed for the expression of the
       composition are computed at a same field level
  
    example :
       MODULE FE_OneStepIteration#set_with_composition_value
          concrete_name = "FE_DiscreteFieldUpdate"
          field_name = "temperature"
          field_level = 0
          MODULE DOFs_values
             type = "from_field_composition"
             field_composition_name = "temperature"
             fields_level = 0
          END MODULE DOFs_values
       END MODULE FE_OneStepIteration#set_with_composition_value

  - a `FE_Parameter::' object defined in the module FE_SetOfParameters :
  
    example :
       MODULE FE_OneStepIteration#set_with_parameter_value
          concrete_name = "FE_DiscreteFieldUpdate"
          field_name = "temperature"
          field_level = 0
          MODULE DOFs_values
             type = "from_parameter"
             parameter_name = "temperature"
          END MODULE DOFs_values
       END MODULE FE_OneStepIteration#set_with_parameter_value

PUBLISHED
*/

class PEL_EXPORT FE_DiscreteFieldUpdate : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------

     ~FE_DiscreteFieldUpdate( void ) ;
      FE_DiscreteFieldUpdate( FE_DiscreteFieldUpdate const& other ) ;
      FE_DiscreteFieldUpdate& operator=(
                              FE_DiscreteFieldUpdate const& other ) ;

      FE_DiscreteFieldUpdate( PEL_Object* a_owner,
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer const* exp ) ;

      

   //-- Plug in

      FE_DiscreteFieldUpdate( void ) ;

      virtual FE_DiscreteFieldUpdate* create_replica( 
                              PEL_Object* a_owner,
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer* exp ) const ;
   //-- Field value

      enum FE_FieldUpdateType
      {
         from_analytic,
         from_composition,
         from_parameter
      } ;

      void set_field_values( size_t f_level,
                             FE_TimeIterator const* t_it ) ;
      void compute_field_value_at_pt( FE_TimeIterator const* t_it,
                                      PDE_LocalFEcell const* fe,
                                      doubleVector& result ) ;

      PEL_Context const* context_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Class attributes

      static FE_DiscreteFieldUpdate const* PROTOTYPE ;      

   //-- Attributes

      // Fields :
      PDE_DiscreteField* const FIELD ;
      size_t const FIELD_LEVEL ;

      // Local assembling :
      PDE_LocalFEcell* const cFE ;

      // Value :
      FE_FieldUpdateType VALUE_TYPE ;
      PEL_Data* VALUE ;
      PDE_FieldComposition* COMPO ;
      size_t COMPO_LEVEL ;
      FE_Parameter* PARAMETER ;
} ;

#endif
