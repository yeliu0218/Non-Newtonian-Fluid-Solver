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

#ifndef PEL_1D_WRITER_HH
#define PEL_1D_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>

#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <fstream>

class size_t_array2D ;
class size_t_vector ;

/*
writers for data on 1D meshings and plotting programs operating with 
files where data are arranged in columns of numbers

For each cycle required to be saved, one or two new file are created with 
a name that follows the syntax :
   a basename given by the entry "files_basename" in the PDE_ResultSaver
              hierarchical data structure.
   + "_c.1d." for fields saved "at_cell_centers", or
     "_v.1d." for fields saved "at_vertices"
   + a five-figure number that is incremented at each new cycle saved.

a.e. with the entry PDE_ResultSaver/files_basename = "save"  the name
of the successive output files will be :
   save_c.1d.00001, save_v.1d.00001, save_c.1d.00002, save_v.1d.00002 ...

The first column of the files save_c.1d.* contains the coordinates of the 
cell centers, whereas the first column of the files save_v.1d.* contains 
the coordinates of the vertices. Each other column contains the associated
values of a component of a field.

The first lines of the files save_c.1d.* and save_v.1d.* begin with 
a # character and are intended to treated as comments by the 
plotting program. These comment lines give the meaning of the subsequent
columns.

PUBLISHED
*/

class PEL_EXPORT PEL_1Dwriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_1Dwriter( void ) ;
      PEL_1Dwriter( PEL_1Dwriter const& other ) ;
      PEL_1Dwriter& operator=( PEL_1Dwriter const& other ) ;
      
      PEL_1Dwriter( PEL_Object* a_owner,
   		     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_1Dwriter( void ) ;

      virtual PEL_1Dwriter* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;


   //-- Writing process main steps

      void determine_conditions( PEL_ModuleExplorer const* exp ) ;
      
   //-- Output file
      
      static void resize_and_add_cell_coords(  
                                           doubleArray2D& cell_matrix,
                                           doubleArray2D const& vertices,
                                           intVector const& cell_nb_verts,
                                           intArray2D const& cell2vert,
                                           size_t nb_columns ) ;

      static void resize_and_add_vert_coords( 
                                           doubleArray2D& vert_matrix,
                                           doubleArray2D const& vertices,
                                           size_t nb_columns ) ;

      static void add_one_field( doubleArray2D& matrix,
                                 size_t& current_column,
                                 doubleArray2D const& vals ) ;

      std::string output_file_name( size_t nb, std::string add_string ) ;

   //-- Data writing
      
      static void write_matrix( doubleArray2D const& matrix, 
                                double time,
                                stringVector const& var_names,
                                std::ofstream& os ) ;

   //-- Class attributes

      static PEL_1Dwriter const* PROTOTYPE ;

   //-- Attributes

      std::ios_base::openmode OPENMODE ;
      bool HAS_CELL_VARS ;
      stringVector CELL_VARS_NAMES ;
      bool HAS_VERT_VARS ;
      stringVector VERT_VARS_NAMES ;
      size_t DIM ;
      double TIME ;
      doubleArray2D CELL_MATRIX ;
      doubleArray2D VERT_MATRIX ;
      std::ofstream FILE_FOR_CELL_VARS ;
      std::ofstream FILE_FOR_VERT_VARS ;
      std::string FILEBASENAME ;
      std::string FILEXTENSION ;
      std::string FORMAT ;
} ;

#endif
