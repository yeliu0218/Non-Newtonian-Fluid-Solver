/*
 *  Copyright :
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#ifndef PEL_MATLABio_HH
#define PEL_MATLABio_HH

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
Input output managment for the data postprocessor MATLAB

PUBLISHED
*/

class PEL_MATLABio
{

   public: //-----------------------------------------------------------

      // File storage type :
      //     Text : ASCII
      enum MATLAB_FORMAT
      {
         Unspecified,
         Text
      } ;

   //-- Write

      static void create_file( std::string const& file_name,
                               MATLAB_FORMAT file_format ) ;


      static void write_int_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               int val ) ;

      static void write_double_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               double val ) ;

      static void write_gridXY_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format, int dimension,
                               std::string const& name,
                               doubleVector const& valx, doubleVector const& valy, doubleVector const& valz ) ;

      static void write_doubleVector_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               doubleVector const& val ) ;

      static void write_doubleVector_variableUU(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               doubleVector const& val ) ;

      static void write_doubleVector_variableSTRESS(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               doubleVector const& val ) ;

      static void write_intVector_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               intVector const& val ) ;

      static void write_doubleArray2D_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               doubleArray2D const& val ) ;

      static void write_intArray2D_variable(
                               std::string const& file_name,
                               MATLAB_FORMAT file_format,
                               std::string const& name,
                               intArray2D const& val ) ;

   //-- Read

      // file format of the MATLAB file named `file_name'
      static MATLAB_FORMAT file_format( std::string const& file_name ) ;


   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PEL_MATLABio( void ) ;
     ~PEL_MATLABio( void ) ;
      PEL_MATLABio( PEL_MATLABio const& other ) ;
      PEL_MATLABio& operator=( PEL_MATLABio const& other ) ;

      // Convert the integer `val' stored in another plateform system
      // into local system (swap of the octet of `val').
      static void convert_to_local( MATLAB_FORMAT file_format, int* val ) ;

      // Convert the float `val' stored in another plateform system
      // into local system (swap of the octet of `val').
      static void convert_to_local( MATLAB_FORMAT file_format, float* val ) ;

   //-- Write

      // Check if `file' is ok for writing.
      static void check_file( std::ofstream const& file ) ;

      static void open_file( std::string const& file_name,
                             MATLAB_FORMAT file_format,
                             bool append,
                             std::ofstream& file ) ;

      static void write_variable(
                           std::string const& file_name,
                           MATLAB_FORMAT file_format,
                           std::string const& name,
                           bool is_integer,
                           int dim, int end,
                           double* XVAR ) ;

      static void write_variableUU(
                           std::string const& file_name,
                           MATLAB_FORMAT file_format,
                           std::string const& name,
                           bool is_integer,
                           int dim, int end,
                           double* XVAR ) ;

      static void write_variableSTRESS(
                           std::string const& file_name,
                           MATLAB_FORMAT file_format,
                           std::string const& name,
                           bool is_integer,
                           int dim, int end,
                           double* XVAR ) ;

      static void write_variableXYZ(
                           std::string const& file_name,
                           MATLAB_FORMAT file_format, int dimension,
                           std::string const& nameX, std::string const& nameY, std::string const& nameZ,
                           doubleVector const& valx, doubleVector const& valy, doubleVector const& valz) ;

      static void write_ASCII_int( std::ofstream &file, int value ) ;
      static void write_ASCII_double( std::ofstream &file, double value ) ;



   //-- Read

      // Check if `file' is ok for reading.
      static void check_file( std::ifstream const& file ) ;

      // Open file named `file_name' and find its format.
      static void open_file( std::string const& file_name,
                             MATLAB_FORMAT& file_format,
                             std::ifstream& file ) ;
      static void close_file( std::ifstream& file ) ;

} ;

#endif
