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

#ifndef FE_FIELD_READER_HH
#define FE_FIELD_READER_HH

#include <FE_OneStepIteration.hh>

#include <size_t_vector.hh>
#include <stringVector.hh>

class PEL_Double ;
class PEL_Context ;
class PEL_Vector ;

/*
   Servers for restoring node values of `PDE_DiscreteField::' objects
   from files generated with `FE_FieldSaver::' server.
   This restoration involves the same discretization (same meshing, same
   reference elements) for the field and the field stored by `FE_FieldSaver::'.
   Hence, the node localisation is verified during the restoration step.

   Example:
   
         MODULE FE_OneStepIteration#field_restoring
            concrete_name = "FE_FieldReader"
            file_name = join( "..", "DoSaving", "field_saving_00001.pel" )
            MODULE discrete_fields
               MODULE df#1
                  name = "velocity"
               END MODULE df#1
               MODULE df#2
                  name = "pressure"
               END MODULE df#2
               MODULE df#3
                  name = "temperature"
               END MODULE df#3
            END MODULE discrete_fields
         END MODULE FE_OneStepIteration#field_restoring

         
   In this example, the node values of the discrete fields of name
   "velocity", "pressure" and "temperature" for all levels are
   initialized with the values stored in "../DoSaving/field_saving_00001.pel".
            
   By default, the field is initialized with the field of the same name.
   It is possible to specify different names:
   
         MODULE FE_OneStepIteration#field_restoring
            concrete_name = "FE_FieldReader"
            ...
            MODULE discrete_fields
               MODULE df#1
                  name = "fluid_velocity"
                  stored_field_name = "velocity"
               END MODULE df#1
               ...
            END MODULE discrete_fields
         END MODULE FE_OneStepIteration#field_restoring

   In this example, the discrete field of name "fluid_velocity" is initialized
   with the values of the field of name "velocity" in the saving file.

   The node values can be initialized with an expression involving
   the nodes values (in a predefined variable $DS_val) read in the file:
   
         MODULE FE_OneStepIteration#field_restoring
            concrete_name = "FE_FieldReader"
            ...
            MODULE discrete_fields
               MODULE df#1
                  name = "velocity"
                  expression = 2.*$DS_val
               END MODULE df#1
               ...
            END MODULE discrete_fields
         END MODULE FE_OneStepIteration#field_restoring

   In this example, the discrete field of name "velocity" is initialized
   with the values of the field read in the file times two.

   By default, all the levels of the field are initialized with the level 0
   of the field in the saving file.
   It is possible to specify level per level the intialization of the field:

         MODULE FE_OneStepIteration#field_restoring
            concrete_name = "FE_FieldReader"
            ...
            MODULE discrete_fields
               MODULE df#1
                  name = "velocity"
                  level = 0
                  stored_field_name = "v"
                  stored_field_level = 1
               END MODULE df#1
               MODULE df#2
                  name = "velocity"
                  level = 1
                  stored_field_name = "v"
                  stored_field_level = 2
               END MODULE df#2
               ...
            END MODULE discrete_fields
         END MODULE FE_OneStepIteration#field_restoring   

   In this example, the discrete field of name "velocity" is initialized as:
      - level 0 initialized with the level 1 of the field of name "v"
      - level 1 initialized with the level 2 of the field of name "v".
   
PUBLISHED
*/

class PEL_EXPORT FE_FieldReader : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------

     ~FE_FieldReader( void ) ;
      FE_FieldReader( FE_FieldReader const& other ) ;
      FE_FieldReader& operator=( FE_FieldReader const& other ) ;

      FE_FieldReader( PEL_Object* a_owner,
                      PDE_DomainAndFields const* dom,
                      FE_SetOfParameters const* prms,
                      PEL_ModuleExplorer const* exp ) ;

      void restore_field( size_t i_field, PEL_ModuleExplorer* exp ) ;
      
   //-- Plug in

      FE_FieldReader( void ) ;

      virtual FE_OneStepIteration* create_replica(
                                PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) const ;
      
   //-- Class attributes

      static FE_FieldReader const* PROTOTYPE ;

   //-- Attributes

      PDE_DomainAndFields const* const DOM ;
      double const DBLE_EPS ;
      double const DBLE_MIN ;
      
      std::string OFILENAME ;
      
      PEL_Vector* FIELDS ; // PDE_DiscreteFields*
      PEL_Vector* EXPRS ; // PEL_Data*
      size_t_vector FIELD_LEVELS ;
      stringVector  FIELD_RNAMES ;
      size_t_vector FIELD_RLEVELS ;

      PEL_Context* CTX ;
      PEL_Double*  VAL ;
} ;

#endif
