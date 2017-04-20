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

#ifndef PEL_GMV_WRITER_HH
#define PEL_GMV_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <fstream>
#include <doubleVector.hh>

class PEL_Data ;

/*
writers for the data postprocessor GMV

For each cycle required to be saved, a new file is created with 
a name that follows the syntax :
   a basename given by the entry "files_basename" in the PDE_ResultSaver
              hierarchical data structure.
   + ".gmv."
   + a five-figure number that is incremented at each new cycle saved.

a.e. with the entry PDE_ResultSaver/files_basename = "save"  the name
of the successive GMV output files will be :
   save.gmv.00001, save.gmv.00002 ...

In case of text output format and to save disk space, at each call to
save the grid a file named "basename"+"_grid.gmv."+a five-figure number 
that is incremented at each new saving of the grid will be created. Grid 
saving in the general GMV output files is replaced by a link to the 
corresponding grid saving file.

GMV treats specifically the time variable as its value appears 
automatically at the top right corner of the main viewer. If the
hierarchical data structure passed as an argument of method 'write_cycle' 
contains a variable named 'TIME' it will be considered by GMV as the 
problem time. This displaying can be cancelled from the GMV interface.

GMV treats specifically the velocity field as its value is automatically
considered as a vectorial field of components respectively named "U", "V" 
and "W". If a field is named "VELO" in the hierarchical data structure 
passed as an argument of method 'write_cycle' it  will be considered 
by GMV as the velocity field. This field can be defined "at_cell_centers" 
or "at_vertices". Anyway a velocity defined "at_cell_centers" will be 
averaged by GMV and turn to a field defined "at_vertices".

For fields that are saved "at_face_centers" in a 2D geometry, additional 
files are produced : at each call to save the grid a file named 
"basename"+"_grid_f.gmv."+a five-figure number that is incremented at 
each new saving of the grid will be created ; for each cycle required to 
be saved, a new file is created with a name that follows the syntax :
   a basename given by the entry "files_basename" in the PDE_ResultSaver
              hierarchical data structure.
   + "_f.gmv."
   + a five-figure number that is incremented at each new cycle saved.

PUBLISHED
*/

class PEL_EXPORT PEL_GMVwriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_GMVwriter( void ) ;
      PEL_GMVwriter( PEL_GMVwriter const& other ) ;
      PEL_GMVwriter& operator=( PEL_GMVwriter const& other ) ;
      
      PEL_GMVwriter( PEL_Object* a_owner,
   		     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_GMVwriter( void ) ;

      virtual PEL_GMVwriter* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;


   //-- Writing process main steps

      void determine_conditions( PEL_ModuleExplorer const* exp ) ;
      
      void write_cells( PEL_ModuleExplorer const* exp ) ;
      void write_faces( PEL_ModuleExplorer const* exp ) ;

      void write_fields( PEL_ModuleExplorer* exp ) ;
      void write_time_variable( double t, std::ofstream& file ) ;
      void write_velocity_field( PEL_ModuleExplorer const* fexp ) ;
      void write_one_field( std::string const& name,
			    PEL_ModuleExplorer const* fexp ) ;

   //-- Output file
      
      std::string output_file_name( size_t nb, std::string add_string ) ;

   //-- Data writing
      
      void write_data( std::string val, size_t length, std::ofstream& file ) ;
      void write_data( int val, std::ofstream& file ) ;
      void write_data( double val, std::ofstream& file) ;
      void write_data_endl( std::ofstream& file) ;

      void save_vertices( doubleArray2D const& vert,
                          std::ofstream& file ) ;
      void save_cells( intVector const& cell_vertices, 
                       intArray2D const& connec,
                       std::ofstream& file ) ;

   //-- Class attributes

      static PEL_GMVwriter const* PROTOTYPE ;

   //-- Attributes

      std::ios_base::openmode OPENMODE ;
      bool HAS_FACE_VARS ;
      size_t DIM ;
      std::ofstream FILE_FOR_PLAIN_VARS ;
      std::ofstream FILE_FOR_FACE2D_VARS ;
      std::string NAME_FILE_CELLS ;
      std::string NAME_FILE_FACES ;
      std::string FILEBASENAME ;
      std::string FILEXTENSION ;
      std::string FORMAT ;
      bool BINARY ;
      size_t ICYCLE ;
};

#endif
