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

#ifndef PEL_DATA_ON_MESHING_WRITER_HH
#define PEL_DATA_ON_MESHING_WRITER_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;
class doubleArray2D ;
class doubleVector ;
class intArray2D ;
class intVector ;
class stringVector ;

/*
writers which save data that are possibly associated to a meshing

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass.
   2. Create a static class pointer which defines the model of the concrete
      class (prototype). It is built calling the prototype constructor :
          `::PEL_DataOnMeshingWriter( std::string const& )'
   3. Implement the function :
          `::create_replica( PEL_Object*, PEL_ModuleExplorer const* ) const'
      which create a new instance of the concrete class, calling the constructor :
        `::PEL_DataOnMeshingWriter( PEL_Object* )'
   4. Implement a destructor
   5. Implement the virtual function :
          `::write_cycle( PEL_ModuleExplorer const*)'

PUBLISHED
*/

class PEL_EXPORT PEL_DataOnMeshingWriter : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_DataOnMeshingWriter* make( PEL_Object* a_owner,
                                            std::string const& name,
                                            PEL_ModuleExplorer const* exp ) ;

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) = 0 ;

      // Can `self' save distributed information in parallel context ?
      // IMPLEMENTATION : false
      virtual bool is_parallel_writer( void ) const ;

   //-- Facilities

      static int index_for_trash( void ) ;
      
      static double undefined_value( void ) ;

      static PEL_Module* create_meshing_module( PEL_Object* a_owner,
          std::string const& meshing_name,
          int nb_sp_dims, doubleArray2D const& vertices,
          intVector const& cell_nb_vertices, intArray2D const& cell2vertex,
          intVector const& cell_nb_faces,    intArray2D const& cell2face,
          intVector const& face_nb_vertices, intArray2D const& face2vertex,
          intArray2D const& face2cell,
          stringVector const& color_table, std::string const& halo_color_name,
          intArray2D const& macro_colors,
          intVector const& vertex_color,
          intVector const& cell_color,
          intVector const& face_color ) ;

      static PEL_Module* create_field_module( PEL_Object* a_owner,
                                              std::string const& field_name,
                                              std::string const& meshing, 
                                              std::string const& location,
                                              doubleArray2D const& value ) ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PEL_DataOnMeshingWriter( void ) ;

      PEL_DataOnMeshingWriter( std::string const& name ) ;

      PEL_DataOnMeshingWriter( PEL_Object* a_owner ) ;

      virtual PEL_DataOnMeshingWriter* create_replica( 
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Error management

      void raise_field_location_error(
                   std::string const& field_name,
                   std::string const& location,
                   stringVector const& allowed_locations ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool write_cycle_PRE( PEL_ModuleExplorer const* exp ) const ;

      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( PEL_DataOnMeshingWriter const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;
      
   private: //----------------------------------------------------------

      PEL_DataOnMeshingWriter( void ) ;
      PEL_DataOnMeshingWriter( PEL_DataOnMeshingWriter const& other ) ;
      PEL_DataOnMeshingWriter& operator=( 
                               PEL_DataOnMeshingWriter const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      bool IS_PROTO ;
};

#endif
