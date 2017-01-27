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

#ifndef FE_PARAMETER_SAVER_HH
#define FE_PARAMETER_SAVER_HH

#include <FE_OneStepIteration.hh>

#include <size_t_vector.hh>
#include <stringVector.hh>

class PEL_Vector ;

class PDE_LocalFEcell ;

/*
  Saving of the values of `FE_Parameter::' objects for subsequent 
  post-processing.

  Saving types are:
     - "cell_values" :
          For each cell of the meshing, the result of 
          `FE_Parameter::cell_value' is saved, related to the cell.
     - "at_cell_centers" :
          For each cell of the meshing, `FE_Parameter::cell_value_at_pt' 
          is evaluated at the cell center and the result is saved, related
          to that cell.
     - "at_vertices" :
          For each vertex of the meshing and for all the cells sharing that
          vertex, `FE_Parameter::cell_value_at_pt' is evaluated. It is then
          checked whether the different values obtained (at the given vertex 
          from its cells) are identical. If not, a fatal error is raised,
          otherwise that value is saved, related to the considered vertex.
          
  Example:

     MODULE PEL_Application
        ...
        MODULE FE_SetOfParameters
           MODULE FE_Parameter#density
              ...
              name = "density_parameter"
           END MODULE FE_Parameter#density
           MODULE FE_Parameter#body_force
              ...
              name = "body_force_parameter"
           END MODULE FE_Parameter#body_force
        END MODULE FE_SetOfParameters
        MODULE FE_OneStepIteration
           concrete_name = "FE_SplitSystem"
           MODULE list_of_FE_OneStepIteration
              ...  
              MODULE FE_OneStepIteration#parameters_postprocessing
                 concrete_name = "FE_ParameterSaver"
                 MODULE parameters
                    MODULE param#0
                       type = "cell_values"
                       parameter_name = "density_parameter"
                       entry_name = "RHO"
                    END MODULE param#0
                    MODULE param#1
                       type = "at_cell_centers"
                       parameter_name = "body_force_parameter"
                       entry_name = "F"
                    END MODULE param#1
                 END MODULE parameters
              END MODULE FE_OneStepIteration#parameters_postprocessing
           END MODULE list_of_FE_OneStepIteration
        END MODULE FE_OneStepIteration
              
  PUBLISHED
*/

class PEL_EXPORT FE_ParameterSaver : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      
   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields(
                      FE_TimeIterator const* t_it, PDE_ResultSaver* rs ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
 
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~FE_ParameterSaver( void ) ;
      FE_ParameterSaver( FE_ParameterSaver const& other ) ;
      FE_ParameterSaver& operator=( FE_ParameterSaver const& other ) ;

      FE_ParameterSaver( PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_ParameterSaver( void ) ;

      virtual FE_OneStepIteration* create_replica(
                                          PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FE_ParameterSaver const* PROTOTYPE ;

   //-- Attributes

      PEL_Vector* const PARAMS ;
      stringVector NAMES ;
      size_t_vector TYPES ;
      
      PDE_LocalFEcell*  cFE ;
      
      double const EPS_DBL ;
      double const MIN_DBL ;
} ;

#endif
