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

#ifndef PEL_FIELDVIEW_WRITER_HH
#define PEL_FIELDVIEW_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>

class doubleVector ;
class intVector ;
class size_t_vector ;
class stringVector ;

/*
Writers for the data postprocessor FieldView: http://www.ilight.com/

PUBLISHED
*/

class PEL_EXPORT PEL_FieldViewWriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_FieldViewWriter( void ) ;
      PEL_FieldViewWriter( PEL_FieldViewWriter const& other ) ;
      PEL_FieldViewWriter& operator=( PEL_FieldViewWriter const& other ) ;
      
      PEL_FieldViewWriter( PEL_Object* a_owner, 
                           PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_FieldViewWriter( void ) ;

      virtual PEL_FieldViewWriter* create_replica(
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;


   //-- Writing process main steps

      void write_meshing( PEL_ModuleExplorer const* exp ) ;

      void write_meshing_header( std::ofstream& file ) const ;

      void write_bound_colors( PEL_ModuleExplorer const* exp,
                               size_t& nb_bounds,
                               size_t_vector& db_color_idx,
                               std::ofstream& file ) const ;
      
      void write_vertices( PEL_ModuleExplorer const* exp,
                           std::ofstream& file ) const ;
      
      void write_boundaries( PEL_ModuleExplorer const* exp,
                             size_t nb_bounds,
                             size_t_vector const& db_color_idx,
                             std::ofstream& file ) const ;
      
      void write_cells( PEL_ModuleExplorer const* exp,
                        std::ofstream& file ) const ;

      void write_result_header( std::ofstream& file ) const ;

      void write_constants( PEL_ModuleExplorer const* exp,
                            std::ofstream& file ) const ;

      void write_fields( PEL_ModuleExplorer const* exp,
                         std::ofstream& file ) const ;
      
   //-- Output file
      
      std::string const& output_file_name(
         size_t cycle, std::string const& add_string ) const ;
         
   //-- Data writing

      void write_str_data( std::string const& val, 
                           std::ofstream& file ) const ;
      void write_data( int val, std::ofstream& file ) const  ;
      void write_data( intVector const& val, std::ofstream& file ) const  ;
      void write_data( double val, std::ofstream& file ) const  ;
      void write_data( doubleVector const& val, std::ofstream& file ) const  ;
      void write_data( stringVector const& val, std::ofstream& file ) const  ;
      void write_comment( std::string const& val, 
                          std::ofstream& file ) const  ;
      void write_data_endl( std::ofstream& file ) const  ;

   //-- Others
      
      int fv_encode_elem_header( int elem ) const ;

      void raise_internal_error( std::string const& mes ) const ;

   //-- Class attributes

      static PEL_FieldViewWriter const* PROTOTYPE ;

   //-- Attributes

      std::string const FILE_BASENAME ;
      bool const BINARY ;
      
      size_t DIM ;
      size_t ICYCLE ;
};

#endif
