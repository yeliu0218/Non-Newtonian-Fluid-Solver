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

#ifndef PEL_OPEN_DX_WRITER_HH
#define PEL_OPEN_DX_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <fstream>
#include <doubleVector.hh>

class PEL_Data ;
class PEL_Map ;

/*
writers for the data postprocessor OpenDX
*/

class PEL_EXPORT PEL_OpenDXwriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_OpenDXwriter( void ) ;
      PEL_OpenDXwriter( PEL_OpenDXwriter const& other ) ;
      PEL_OpenDXwriter& operator=( PEL_OpenDXwriter const& other ) ;
      
      PEL_OpenDXwriter( PEL_Object* a_owner,
   		        PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_OpenDXwriter( void ) ;

      virtual PEL_OpenDXwriter* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Write

      void write_one_variable( std::string const& name,
			       PEL_Data const* val ) ;

      void write_grid( PEL_ModuleExplorer const* exp ) ;

      void write_field( PEL_ModuleExplorer const* exp ) ;

      void write_data_header( void ) ;
      void write_data( int val ) ;
      void write_data( double val ) ;
      void write_data_endl( void ) ;
      
      size_t save_vector( intVector const& I ) ;
      size_t save_array( intArray2D const& I ) ;
      size_t save_vector( doubleVector const& D ) ;
      size_t save_array( doubleArray2D const& D ) ;
      size_t save_array_with_vectors( doubleArray2D const& D ) ;
      size_t save_array_with_vectors( intArray2D const& I ) ;
      size_t save_scalar( double val ) ;
      size_t save_scalar( int val ) ;
      size_t new_object( void ) ;
      void add_serie( std::string const& name,
                      size_t obj,
                      int serie ) ;
      void finish( void ) ;

   //-- Class attributes

      static PEL_OpenDXwriter const* PROTOTYPE ;

   //-- Attributes

      std::ofstream file ;
      std::ofstream::pos_type file_pos ;
      std::ofstream binfile ;
      std::ofstream::pos_type binfile_pos ;
      size_t object_number ;
      size_t conn_object ;
      size_t pos_object ;
      PEL_Map* obj_map ;
      doubleVector t ;
      size_t icycle;
      bool binary ;
      std::string filename ;
      std::string bin_filename ;

};

#endif
