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

#ifndef PEL_TIC_IO_HH
#define PEL_TIC_IO_HH

#include <iosfwd>
#include <string>

#include <PEL.hh>

class PEL_Module ;
class PEL_Object ;

class doubleVector ;
class intArray2D ;
class intVector ;
class doubleArray2D ;

/*
Input output managment for the data postprocessor TIC
  
PUBLISHED
*/

class PEL_EXPORT PEL_TICio
{

   public: //-----------------------------------------------------------

      // File storage type :
      //     Text : ASCII
      //     Binary : Binary FORTRAN (in local system)
      //     CBinary : Binary C (in local system)
      //     Binary_no_local : Binary FORTRAN (in no local system)
      enum TIC_FORMAT
      {
         Unspecified,
         Text,
         CBinary,
         Binary,
         Binary_no_local
      } ;
      
   //-- Write

      static void create_file( std::string const& file_name,
                               TIC_FORMAT file_format ) ;
      
      static void write_new_cycle(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               size_t nb_vars, size_t i_cycle ) ;

      static void write_int_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               int val ) ;

      static void write_double_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               double val ) ;

      static void write_doubleVector_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               doubleVector const& val ) ;

      static void write_intVector_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               intVector const& val ) ;

      static void write_doubleArray2D_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               doubleArray2D const& val ) ;

      static void write_intArray2D_variable(
                               std::string const& file_name,
                               TIC_FORMAT file_format,
                               std::string const& name,
                               intArray2D const& val ) ;
      
   //-- Read

      // file format of the TIC file named `file_name'
      static TIC_FORMAT file_format( std::string const& file_name ) ;
      
      /* Create a module storing entries obtained by reading
         the file of name `file', that was previously created by `self'.
         Restore saving cycle defined by `i_cycle', restore all cycles if
         `i_cycle'= `PEL::bad_index()'. */
      static PEL_Module* create_from_gene_file(
                                   PEL_Object* a_owner,
                                   std::string const& mod_name,
                                   std::string const& file_name,
                                   size_t i_cycle = PEL::bad_index() ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PEL_TICio( void ) ;
     ~PEL_TICio( void ) ;
      PEL_TICio( PEL_TICio const& other ) ;
      PEL_TICio& operator=( PEL_TICio const& other ) ;

      // Convert the integer `val' stored in another plateform system
      // into local system (swap of the octet of `val').
      static void convert_to_local( TIC_FORMAT file_format, int* val ) ;

      // Convert the float `val' stored in another plateform system
      // into local system (swap of the octet of `val').
      static void convert_to_local( TIC_FORMAT file_format, float* val ) ;

   //-- Write

      // Check if `file' is ok for writing.
      static void check_file( std::ofstream const& file ) ;

      static void open_file( std::string const& file_name,
                             TIC_FORMAT file_format,
                             bool append,
                             std::ofstream& file ) ;

      // Check if `name' follows the FORTRAN convention.
      static void check_variable_name( std::string const& name,
                                       bool is_integer ) ;

      static void write_variable(
                           std::string const& file_name,
                           TIC_FORMAT file_format,
                           std::string const& name,
                           bool is_integer,
                           int dim, int end,
                           double* XVAR ) ;      

      static void write_ASCII_int( std::ofstream &file, int value ) ;
      static void write_ASCII_double( std::ofstream &file, double value ) ;

      // Write the lenght of the stored line (BINARY FORTRAN convention).
      static void write_length( std::ofstream &file, TIC_FORMAT file_format,
                                int length )  ;

      
      
   //-- Read

      // Check if `file' is ok for reading.
      static void check_file( std::ifstream const& file ) ;

      // Open file named `file_name' and find its format.
      static void open_file( std::string const& file_name,
                             TIC_FORMAT& file_format,
                             std::ifstream& file ) ;
      static void close_file( std::ifstream& file ) ;

      // Read a new cycle : read the number of the current cycle, and
      // the number of variables stored.
      static void read_new_cycle( std::ifstream& file,
                                  TIC_FORMAT file_format,
                                  size_t& nb_vars, size_t& i_cycle ) ;

      // Read a new variable.
      static void read_new_variable(
                              std::ifstream& file,
                              TIC_FORMAT file_format,
                              std::string& name,
                              size_t& type,
                              size_t& dim, size_t& end,
                              bool& is_integer ) ;

      // Read a table of floats.
      static void read_float( std::ifstream& file,
                              TIC_FORMAT file_format,
                              size_t dim, size_t end,
                              float* XVAR ) ;

      // Read a table of ints.
      static void read_int( std::ifstream& file,
                            TIC_FORMAT file_format,
                            size_t dim, size_t end,
                            int* IVAR ) ;

      // Read the lenght of the stored line (BINARY FORTRAN convention).
      static int read_length( std::ifstream &file, TIC_FORMAT file_format,
                              int expected_length )  ;
} ;

#endif
