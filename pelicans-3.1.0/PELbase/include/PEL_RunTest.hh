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

#ifndef PEL_RUN_TEST_HH
#define PEL_RUN_TEST_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;
class PEL_String ;

/*
  PELICANS application intended to run a list of reference tests.

  GENERAL CONFIGURATION:
  ---------------------
  
  Typical data file associated to PEL_RunTest utility looks like:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                test_directories = < "dir1" "dir2" ... "dirn" >
             END MODULE PEL_Application

  All data files of name "data.pel" recursively found in the directories
  defined by "test_directories" keyword are executed with this that has
  launch current "peltest" application.
  For each case, a new directory, whose name is linked with original one
  but prefixed with "PELICANS_TESTS", will be created for execution. The
  generated files are compared to original ones.

  1/ The default data file name "data.pel" can be changed:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                data_filename = "data2.pel"
             END MODULE PEL_Application

  2/ The default executable can be changed:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                executable = join( "home", "pelicans", "lib", "exe0" )
             END MODULE PEL_Application

  3/ The default root directory "PELICANS_TESTS" can be changed:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                test_directories = < "dir1" "dir2" ... "dirn" >
                output_directories = < "dir1" "dir2" ... "dirn" >
             END MODULE PEL_Application

  4/ Additional options can be defined for test executions:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                run_options = < "-Call" >
             END MODULE PEL_Application

  5/ Optional tools for pattern file managment ( cf. `PEL_ModulePattern::'):
  
      - pattern file is built with all data files
        (tests are executed, and pattern is build):

      
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                MODULE pattern
                   type = "build"
                   pattern_filename = "pattern.pel" // created
                END MODULE pattern
             END MODULE PEL_Application

      - all data files are checked with an existing pattern file
        (tests are NOT executed):
    
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                MODULE pattern
                   type = "verify"
                   pattern_filename = "pattern.pel" // should exists
                END MODULE pattern
             END MODULE PEL_Application
  
      - pattern file is built with all data files
        (tests are executed, and pattern is build)
        and then all data files are checked with this pattern file:
      
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                MODULE pattern
                   type = "build_then_verify"
                   pattern_filename = "pattern.pel" // created
                END MODULE pattern
             END MODULE PEL_Application

   6/ All generated files (except resu) are compared with the original ones

         1/ Native files (with the GENE, PEL or CSV format) are compared
            with the "pelcmp" application (cf. `PEL_Comparator::').

            By default, the format of the files are identified by
            a motif appearing in their name: ".gene" corresponds to GENE, 
            ".pel" corresponds to PEL and ".csv" corresponds to CSV.

            Other correspondances between the motifs and the formats 
            can be defined:
                    
               MODULE PEL_Application
                  concrete_name = "peltest"
                  ...
                  MODULE pel_compare
                     motifs  = < ".su" ".csv" ".tic" >
                     formats = < "PEL" "CSV"  "GENE" >
                  END MODULE pel_compare
               END MODULE PEL_Application

         2/ The other files are compared line per line

         3/ It is nevertheless possible to ignore some of the generated files 
            in the comparison:
          
               MODULE PEL_Application
                  concrete_name = "peltest"
                  ...
                  files_to_ignore = < "toto.txt" >
               END MODULE PEL_Application

  6/ Test execution diagnostic:

     - the execution of the test raises an error.

          the message: "Test failed" is printed

     - some generated file differ from the original one:
     
          the message: "Test success ?... to be analyzed..." is printed

     - test is successful
     
          the message: "Test is successful" is printed

  8/ Name of the output file:
  
             MODULE PEL_Application
                concrete_name = "peltest"
                ...
                output_file = "resu"
             END MODULE PEL_Application

  9/ The comparison of double can be customized.

        9.1/ Comparison without any tolerance:
  
                MODULE PEL_Application
                   concrete_name = "peltest"
                   ...
                   MODULE double_comparison
                      type = "exact"
                   END MODULE double_comparison
                END MODULE PEL_Application

        9.2/ Comparison with `PEL::double_equality':
  
                MODULE PEL_Application
                   concrete_name = "peltest"
                   ...
                   MODULE double_comparison
                      type = "PEL_double_equality"
                      dbl_eps = 1.e-5
                      dbl_min = 1.e-10
                   END MODULE double_comparison
                END MODULE PEL_Application
   
        These customization always overread any setting performed
        via individual "config.pel" files


  TEST CASES CONFIGURATION:
  ------------------------
  
  Each test can be tuned with a configuration file named "config.pel"
  (that name may be changed).
  
  An example of such a file is:
   
     MODULE test_config
        run_options = < "-Call" >                      // Extra options
        mpi_machinefile = vector( host_name() )        // mpi machines
        mpi_options = vector( "-np", "2" )             // mpi options
        files_to_ignore = < "test.gene" "test.bgene" > // Ignored files
        MODULE PEL_Comparator           // Specific comparison options
           MODULE xxx                   // Non significant name
              filename = "save.txt"     // file under consideration
              ignore_data = < "UU".. >  // Ignored entries
              format = "GENE"           // format (optional, GENE,PEL or CSV)
           END MODULE xxx
        END MODULE PEL_Comparator
     END MODULE test_config
          
  PUBLISHED
*/

class PEL_EXPORT PEL_RunTest : public PEL_Application
{
   public: //-----------------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~PEL_RunTest( void ) ;
      PEL_RunTest( PEL_RunTest const& other ) ;
      PEL_RunTest& operator=( PEL_RunTest const& other ) ;

      PEL_RunTest( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

      PEL_RunTest( PEL_Object* a_owner, stringVector& args ) ;

   //-- Plug in

      PEL_RunTest( void ) ;

      virtual PEL_RunTest* create_replica(
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp ) const ;

      virtual PEL_RunTest* create_replica_from_args(
                                        PEL_Object* a_owner,
                                        stringVector& args ) const ;
   //-- Command line

      virtual void print_usage( void ) const ;
      virtual void print_operands( void ) const ;
      virtual void print_exit_status( void ) const ;

   //-- Program core execution      

      void do_tests( void ) ;
      
      void do_test( std::string const& base,
                    std::string const& print_name,
                    std::string const& test_directory,
                    std::string const& filename ) ;

      void notify_error( std::string const& msg ) ;
      void notify_differ( std::string const& msg ) ;     

      PEL_ModuleExplorer* create_config( PEL_Object* a_owner,
                                         std::string const& base ) ;
      
   //-- Output

      void do_report( PEL_ModuleExplorer const* conf,
                      int err_code,
                      std::string const& base_dir,
                      std::string const& test_dir,
                      std::string const& filename ) ;
      void do_compare( std::string const& first,
                       std::string const& second,
                       PEL_ModuleExplorer const* options,
                       std::string const& format ) ;

      std::ostream& out( void ) const ;
      
      void notify_test_failure( std::string const& msg ) ;
      void notify_test_ambiguous( std::string const& msg ) ;

   //-- Internal tools

      static void set_resu_dirs( stringVector const& test_dirs,
                                 std::string& common_root,
                                 stringVector& resu_dirs ) ;
      
      static void root_directory( std::string const& dir,
                                  std::string& root_dir,
                                  std::string& base_dir ) ;
      
      PEL_ModuleExplorer const* specific_option( 
                                         PEL_ModuleExplorer const* conf,
                                         std::string const& filename ) ;

      void set_mpi_run( void ) ;
      
      void check_mpi_opts( stringVector const& machines,
                           stringVector const& mpi_options ) ;

   //-- Class attributes

      static PEL_RunTest const* PROTOTYPE ;

   //-- Attributes

      bool VERBOSE ;

      // Executable:
      std::string TESTED_EXE ;
      std::string PELCMP_EXE ;

      // Tests:
      std::string DATA_FILE ;
      std::string CONFIG_FILE ;
      std::string COMMON_ROOT ;
      stringVector LIST ;
      stringVector DIR_IGNORED ;
      stringVector NAME ;
      stringVector DIR ;
      stringVector OPTIONS ;

      // Pattern:
      enum Pattern
      {
         ignore,
         build,
         verify,
         build_then_verify
      } ;
      Pattern ACTION ;
      bool VERIFY ;
      std::string PATTERN_FILE ;

      enum TestResult
      {
         test_successful,
         test_failed,
         test_ambiguous
      } ;
      TestResult TEST_RESULT ;
      size_t NB_SUCCESSFUL ;
      size_t NB_AMBIGUOUS ;
      size_t NB_FAILED ;
      
      // Files ignored
      stringVector GLOB_IGNORED ;

      // Output file
      std::string RESU ;

      // Parallel test
      std::string MPI_RUN ;

      // Default parameters for numerical comparisons when enabled
      bool DBL_EXACT ;
      double MY_DBL_EPS ;
      double MY_DBL_MIN ;

      // For execution context
      PEL_String* MODE ;
} ;

#endif
