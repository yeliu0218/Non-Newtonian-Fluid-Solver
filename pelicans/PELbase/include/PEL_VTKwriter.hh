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

#ifndef PEL_VTKwriter_HH
#define PEL_VTKwriter_HH

#include <PEL_DataOnMeshingWriter.hh>

#include <doubleArray2D.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <fstream>

class PEL_Data ;
class PEL_Communicator ;

/*
writers for data postprocessor using VTK XML file format as ParaView.

For each cycle required to be saved, a new file is created with 
a name that follows the syntax :
   a basename given by the entry "files_basename" in the PDE_ResultSaver
              hierarchical data structure.
   + T + cycle_number + ".vtu" ( sequential case )
     Ex : saveT1.vtu
   + T + cycle_number + "_" + rank() + ".vtu" ( parallel case )
     Ex : saveT1_2.vtu

Moreover, in the parallel case, a global ".pvtu" file is saved 
by the 0-ranked processor :
     Ex : saveT1.pvtu

The overall sequence of produced files is described under the 
ParaView native data file format in a single file with ".pvd"
extension. It is the file that should be opened with ParaView
(menu file/open) for subsequent postprocessing.

Notes :
o Text and binary formats are available.
o Used binary format is knowned in VTK jargon as
   "XML-embedded appended binary compressed non-encoded format".
o Compressed capability depends on linking with zlib library and
   compiling with ZLIB macro definition.

PUBLISHED
*/

class PEL_EXPORT PEL_VTKwriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;
      
      // IMPLEMENTATION : true
      virtual bool is_parallel_writer( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_VTKwriter( void ) ;
      PEL_VTKwriter( PEL_VTKwriter const& other ) ;
      PEL_VTKwriter& operator=( PEL_VTKwriter const& other ) ;
      
      PEL_VTKwriter( PEL_Object* a_owner,
   		     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_VTKwriter( void ) ;

      virtual PEL_VTKwriter* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;


   //-- Writing process main steps
      
      void write_pvd_file( PEL_ModuleExplorer const* exp,
                           std::string const& vtu_filename ) ;
      void store_meshing( PEL_ModuleExplorer const* exp ) ;
      void build_vtu( PEL_Module* vtk,
                      PEL_ModuleExplorer const* exp,
                      bool parallel )  ;
      void write_grid( PEL_Module* base,
                       bool parallel )  ;
      void write_fields( PEL_Module* base,
                         PEL_ModuleExplorer* exp,
                         bool parallel ) const ;
      void write_one_field( PEL_ModuleExplorer const* fexp,
                            PEL_Module* point_data,
                            PEL_Module* cell_data ) const ;
      
   //-- Output file
      
      std::string output_file_name( size_t nb,
                                    bool parallel,
                                    size_t rank ) ;
      
   //-- Data writing
      
      void write_vtk( PEL_ModuleExplorer* vtk,
                      std::ofstream& file,
                      size_t level,
                      bool parallel ) ;
      void start_output( size_t size, size_t number ) ;
      void write_double( std::ofstream& file, double val ) ;
      void write_int( std::ofstream& file, int val ) ;
      size_t store_int( int val ) ;
      void flush( std::ofstream& file )         ;
      void reserve_double( size_t size ) ;
      void check_allocated( size_t size ) ;

      void compress_segment( size_t seg ) ;
      

   //-- Class attributes

      static PEL_VTKwriter const* PROTOTYPE ;

   //-- Attributes

      size_t SPACE_DIM ;
      std::string BASE_FILENAME ;
      std::string PVD_FILENAME ;
      stringVector PVD_STRINGS ;
      bool BINARY ;
      size_t CYCLE_NUMBER ;

      static std::string UNAMED ;
      enum VTK_TYPE
      { VTK_UNDEF=0, VTK_LINE=3, VTK_TRIANGLE=5, VTK_QUAD=9, VTK_TETRA=10, VTK_HEXAHEDRON=12 } ;
         
      doubleArray2D VERTICES ;
      intVector LINEAR_CONNECTIVITY ;
      intVector OFFSET_CONNECTIVITY ;
      intVector CELL_TYPES ;
      size_t NB_CELLS ;
      char * BUFFER ;
      size_t ALLOCATED ;
      size_t OFFSET ;
      size_t CURRENT_LENGTH ;
      std::string ENCODING ;
      bool COMPRESS ;
      PEL_Communicator const* COM ;
      
};

#endif
