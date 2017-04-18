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

#ifndef FE_FIELD_SAVER_HH
#define FE_FIELD_SAVER_HH

#include <FE_OneStepIteration.hh>

#include <size_t_vector.hh>
#include <doubleVector.hh>

class PEL_Module ;
class PEL_Vector ;

class PDE_DiscreteField ;

/*
   Servers for saving node values of `PDE_DiscreteField::' objects.
   The generated files can be reload with `FE_FieldReader::' server in
   order to initialize fields for new computation.

   - each saved cycle is stored in a separate file (one cycle per file)
   
      Example:
   
         MODULE FE_OneStepIteration#field_saving
            concrete_name = "FE_FieldSaver"
            MODULE post_processing
               type = "cycles_in_separate_files"
               file_basename = "fields_saving"
            END MODULE post_processing
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
         END MODULE FE_OneStepIteration#field_saving

      In this example, the node values of the discrete fields of name
      "velocity", "pressure" and "temperature" are saved.
      A sequence of text files named "fields_saving.00001.pel", 
      "fields_saving.00002.pel",... is created to store respectively 
      the first cycle, the second cycle,...
      (a sequence of companion binary files named "fields_saving.00001.pel.bin", 
       fields_saving.00002.pel.bin",... is also created to stored the double values).

   - only the last two cycles are stored

     example :
     
         MODULE FE_OneStepIteration#field_saving
            concrete_name = "FE_FieldSaver"
            MODULE post_processing
               type = "last_two_cycles"
               file_name_0 = join( getcwd(), "saving_0.pel" )
               file_name_1 = join( getcwd(), "saving_1.pel" )
            END MODULE post_processing
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
         END MODULE FE_OneStepIteration#field_saving

      The text files "saving_0.pel" and "saving_1.pel" are created to store
      the last two cycles ; the time of last modification of these files 
      identifies that of the more recent saving.
      (the companion binary files named "saving_0.pel.bin" and 
       "saving_1.pel.bin" are also created to store with the double values).       

 
   By default, the saving times correspond to the
   `::save_other_than_time_and_fields' times.
   It is possible to specify different saving times:
   
         MODULE FE_OneStepIteration#field_saving
            concrete_name = "FE_FieldSaver"
            MODULE post_processing
               ...
               saving_times = < 0. 1. 2. 3. >
            END MODULE post_processing
            ...
         END MODULE FE_OneStepIteration#field_saving
              
   By default, the format of the saving file is the "hybrid" PELICANS format
   (binary for double data, text elsewhere).
   It is possible to specify "text" format:

         MODULE FE_OneStepIteration#field_saving
            concrete_name = "FE_FieldSaver"
            MODULE post_processing
               ...
               output_format = "text" // default: "hybrid"
            END MODULE post_processing
            ...
         END MODULE FE_OneStepIteration#field_saving

   By default, the field level used for saving is 0.
   It is possible to specify different field level:
   
         MODULE FE_OneStepIteration#field_saving
            concrete_name = "FE_FieldSaver"
            ...
            MODULE discrete_fields
               MODULE df#1
                  name = "velocity"
                  level = 1
               END MODULE df#1
               ...
            END MODULE discrete_fields
         END MODULE FE_OneStepIteration#field_saving
   
PUBLISHED
*/

class PEL_EXPORT FE_FieldSaver : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      
      virtual void do_after_time_adaptation( FE_TimeIterator const* t_it ) ;
      
   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
                                            FE_TimeIterator const* t_it,
                                            PDE_ResultSaver* rs ) ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------
      
      // per_one_cycle :   save all the cycles, but one per file
      // last_two_cycles : save only the two last cycles
      enum FE_FieldSaverType 
      {
         none,
         per_one_cycle,
         last_two_cycles
      } ;
      
     ~FE_FieldSaver( void ) ;
      FE_FieldSaver( FE_FieldSaver const& other ) ;
      FE_FieldSaver& operator=( FE_FieldSaver const& other ) ;

      FE_FieldSaver( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom,
                     FE_SetOfParameters const* prms,
                     PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      FE_FieldSaver( void ) ;

      virtual FE_OneStepIteration* create_replica(
                                PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) const ;
      
   //-- Post-processing

      void initialize_file( size_t i_cycle, double time ) ;
      void add_field( size_t i_field, PEL_Module* mod ) ;
      void save_field_value( size_t i_cycle, FE_TimeIterator const* t_it ) ;
      
   //-- Class attributes

      static FE_FieldSaver const* PROTOTYPE ;

   //-- Attributes
      
      PDE_DomainAndFields const* const DOM ;
      double const DBLE_EPS ;
      double const DBLE_MIN ;
      
      PEL_Vector* FIELDS ; // PDE_DiscreteFields const*
      size_t_vector FIELD_LEVELS ;
      
      FE_FieldSaverType OTYPE ;
      std::string OFILEFORMAT ;
      std::string OFILEBNAME ;
      std::string OFILENAME1 ;
      std::string OFILENAME2 ;
      
      std::string OFILENAME ;
      size_t SAVING_NUMBER ;
      double NEXT_SAVING_TIME ;
      doubleVector SAVING_TIMES ;
} ;

#endif
