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

#ifndef FE_FIELD_VALUE_HH
#define FE_FIELD_VALUE_HH

#include <FE_OneStepIteration.hh>

#include <doubleVector.hh>
#include <size_t_vector.hh>
#include <vector>

class FE_Parameter ;

class GE_Point ;

class PDE_DomainAndFields ;
class PDE_DiscreteField ;
class PDE_LocalFEcell ;
class PDE_PointInGridFE ;

class PEL_Data ;
class PEL_Double ;
class PEL_DoubleVector ;

class doubleArray2D ;

/*
PUBLISHED

   Servers for determination of the values of several fields, parameters
   or expressions (or several components of fields, parameters or expressions) 
   at a succession of points of the domain.

   * Points definition

     1/ List of points
     
        MODULE FE_OneStepIteration
           concrete_name = "FE_FieldValue"
           MODULE points_definition
              type = "list_of_points"
              points = array( < 0. 0. >, < 1. 2. > )
           END MODULE points_definition
        END MODULE FE_OneStepIteration

        The values are computed at points < 0. 0. > and < 1. 2. >.

     2/ Regular cutline
     
        MODULE FE_OneStepIteration
           concrete_name = "FE_FieldValue"
           MODULE points_definition
              type = "regular_cutline"
              first_endpoint  = < 0.333 0.1 0. >
              second_endpoint = < 0.333 0.9 0. >
              number_of_points = 5
           END MODULE points_definition
        END MODULE FE_OneStepIteration

        The values are computed at 5 points regularly defined between
        < 0.333 0.1 0. > and < 0.333 0.9 0. > (points included)
   
     3/ General cutline
     
        MODULE FE_OneStepIteration
           concrete_name = "FE_FieldValue"
           MODULE points_definition
              type = "cutline"
              first_endpoint  = < 0.333 0.1 0. >
              second_endpoint = < 0.333 0.9 0. >
              curvilinear_abscissae = < -1. 0. 0.5 >
           END MODULE points_definition
        END MODULE FE_OneStepIteration

        The values are computed at the points defined with their
        given curvilinear abscissae in [first_endpoint,second_endpoint]
        local coordinates.

   * Remark on points that are not found within the grid:

        Some points may be located outside the grid. In order to remove them
        automatically from the cutline, the optional entry of keyword
        `ignore_exterior_points' may be used (its default data being `false'):

           MODULE FE_OneStepIteration
              concrete_name = "FE_FieldValue"
              MODULE points_definition
                 type = ...
                 ignore_exterior_points = true
              END MODULE points_definition
           END MODULE FE_OneStepIteration
  
   * Fields, parameters and expression definition :

      MODULE FE_OneStepIteration
         ...
         MODULE fields
            MODULE field#1
               name = "FF"
            END MODULE field#1
            MODULE field#2
               name = "UU"
               component = 1
            END MODULE field#2
        END MODULE fields
        MODULE parameters
           MODULE parameter#1
              name = "param1"
              type = "at_points"
           END MODULE parameter#1
           MODULE parameter#2
              name = "param2"
              type = "cell_values"
              component = 1
           END MODULE parameter#2
        END MODULE parameters
        MODULE field_compositions
           MODULE field_composition#1
              name = "compo1"
           END MODULE field_composition#1
           MODULE field_composition#2
              name = "compo2"
              component = 1
           END MODULE field_composition#2
        END MODULE field_compositions
        MODULE expressions
           MODULE expression#1
              value = vector( $DS_T*component( $DV_X, 0 ) )
           END MODULE expression#1
        END MODULE expressions
      END MODULE FE_OneStepIteration

      The fields of name "FF" and "UU", as well as the expression t*x,
      are evaluated at the different points.

      The parameter of name "param1" is evaluated with
      `FE_Parameter::cell_value_at_pt' at the different points
      (`type' "at_points") and the parameter of name "param2" is
      evaluated with `FE_Parameter::cell_value' for the cells containing
      the different points (`type' "cell_values").

      `FE_FieldCompositionParameter::' objects are built with the field
      compositions of name "compo1" and "compo2" and are evaluated at the
      different points (all the values of the discrete fields are computed
      at at given level 0).

   * Post processing :

      1/ In a file:
      
         MODULE FE_OneStepIteration
            ...
            MODULE post_processing
               type = "one_file"
               file_name = "values"
            END MODULE post_processing
         END MODULE FE_OneStepIteration

         All the values computed are stored in an unique file
         of name "values.txt", one line per storage cycle

         #
         # FE_FieldValue generated file
         #
         #    Field values computed at:
         #        pt0 = ( 0.000000000e+00 , 0.000000000e+00 )
         #        pt1 = ( 0.100000000e+01 , 0.200000000e+01 )
         #
         #        time ####          uu # ...
            0.0000e+00      0.00000e+00   ...
            1.0000e-00      3.06164e-01   ...
            2.0000e+00      6.16229e-01   ...

         Rem : the description of the file (lines beginning with #) can
               be omitted with the optional keyword :
      
                   MODULE post_processing
                      ...
                      banner = false
                   END MODULE post_processing

         Rem : by default, the saving times correspond to the
               `::save_other_than_time_and_fields' times.
               It is possible to specify different saving times :
               
                   MODULE post_processing
                      ...
                      saving_times = < 0. 1. 2. 3. >
                   END MODULE post_processing
                   
      2/ In separated files:
      
         MODULE FE_OneStepIteration
            ...
            MODULE post_processing
               type = "separated_files"
               file_basename = "values"
            END MODULE post_processing
         END MODULE FE_OneStepIteration

         Output will be done in one separate text file at each saving
         cycle :
         - File name will be formed by the output_file_basename, the
           saving cycle number and the ".txt" extension.
         - Output format will be the following one :
            - Each line corresponds to one point
            - First the point coordinates (2 or 3 depending on the
              space dimension), then the curvilinear abscissae for
              cutline computation, followed by the corresponding field
              values in the order defined in the data deck, finally
              followed by the corresponding expression values in the
              order defined in the data deck

         #
         # FE_FieldValue generated file
         #
         #    Field values computed at time: 6
         #
         #             coordinates #      curve #         uu #
          0.00000e+00  0.00000e+00  0.00000e+00  1.00000e+00
          0.10000e+01  0.20000e+02  0.10000e+01  9.88063e-01

         Rem : the description of the file (lines beginning with #) can
               be ommited with the optionnal keyword :
      
                   MODULE post_processing
                      ...
                      banner = false
                   END MODULE post_processing

         Rem : by default, the saving times correspond to the
               `::save_other_than_time_and_fields' times.
               It is possible to specify different saving times :
               
                   MODULE post_processing
                      ...
                      saving_times = < 0. 1. 2. 3. >
                   END MODULE post_processing
                   
      3/ Performed by the `PDE_ResultSaver' object transmitted to
         `::save_other_than_time_and_fields'
      
         MODULE FE_OneStepIteration
            ...
            MODULE post_processing
               type = "result_saver"
               variable_name = "VAL"
            END MODULE post_processing
         END MODULE FE_OneStepIteration

         The array "VAL" is stored in the result saver
               VAL(i,j) : j-th values at the i-th points

         For the previous example : at each storing cycle,
             VAL(1,1) : value of "F" at <0. 0.>
             VAL(1,2) : value of "UU(1)" at <0. 0.>
             VAL(1,3) : value of the expression t*x at <0. 0.>
             VAL(2,1) : value of "F" at <1. 2.>
             VAL(2,2) : value of "UU(1)" at <1. 2.>
             VAL(2,3) : value of the expression t*x at <1. 2.>
*/

class PEL_EXPORT FE_FieldValue : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      
      virtual void do_after_time_adaptation( FE_TimeIterator const* t_it ) ;
      
   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
                                            FE_TimeIterator const* t_it,
                                            PDE_ResultSaver* rs ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      enum SavingType { one_file, separated_files, result_saver } ;         

     ~FE_FieldValue( void ) ;
      FE_FieldValue( FE_FieldValue const& other ) ;
      FE_FieldValue& operator=( FE_FieldValue const& other ) ;

      FE_FieldValue( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom,
                     FE_SetOfParameters const* prms,
                     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_FieldValue( void ) ;

      virtual FE_FieldValue* create_replica(
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const ;

   //-- Initialization

      void set_fields_and_expressions( PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp ) ;
      void set_points( PDE_DomainAndFields const* dom,
                       PEL_ModuleExplorer const* exp ) ;
      void set_post_processing( PDE_DomainAndFields const* dom,
                                PEL_ModuleExplorer const* exp ) ;
      
   //-- Computation

      void save_field_value( FE_TimeIterator const* t_it,
                             PDE_ResultSaver* rs ) ;
      
      void set_context( GE_Point const* pt, double time ) ;
      
      void compute_values( FE_TimeIterator const* t_it,
                           doubleArray2D& values ) ;
      void compute_values_at_point(
                           size_t pt_index, FE_TimeIterator const* t_it,
                           bool& found, doubleArray2D& values ) ;

   //-- Post-processing
      
      void initialize_one_file( void ) const ;
      void save_values_in_one_file(
                        double time, doubleArray2D const& values ) const ;
      
      void save_values_in_separated_files(
                        size_t i_cycle,
                        double time, doubleArray2D const& values ) const ;

      void print_name( std::ostream& os, std::string name,
                       size_t nb_col ) const ;
      void print_value( std::ostream& os, double value ) const ;

   //-- Class attributes

      static FE_FieldValue const* PROTOTYPE ;

   //-- Attributes

      // Points
      size_t const NB_SP_DIMS ;
      bool REMOVE_PTS ;
      std::vector<GE_Point const*> POINTS ;
      doubleVector ABS_PTS ;
      PDE_PointInGridFE* const PIG ;
      PDE_LocalFEcell* const cFE ;
      size_t_vector PT_CELL ;
      size_t NB_VALUES ;
      bool REDO ;
      
      // Fields to be saved
      std::vector<PDE_DiscreteField const*> FIELDS ;
      size_t_vector FIELD_COMPS ;
      
      // Parameters to evaluate
      std::vector< FE_Parameter const* > PARAMS ;
      size_t_vector PARAM_TYPES ;
      size_t_vector PARAM_COMPS ;
      
      // Expressions to evaluate
      std::vector< PEL_Data* > EXPRS ;
      size_t_vector EXPR_COMPS ;
      size_t_vector EXPR_DIMS ;
      
      // For the evaluation of expressions
      PEL_DoubleVector* COORDS ;
      PEL_Double* TT ;

      // Post-processing
      FE_FieldValue::SavingType SAVING ;
      bool BANNER ;
      std::string OFILENAME ;
      std::string RSNAME ;
      
      size_t SAVING_NUMBER ;
      double NEXT_SAVING_TIME ;
      doubleVector SAVING_TIMES ;
} ;

#endif
