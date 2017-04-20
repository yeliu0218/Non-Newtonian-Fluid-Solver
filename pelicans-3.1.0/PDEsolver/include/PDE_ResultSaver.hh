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

#ifndef PDE_RESULT_SAVER_HH
#define PDE_RESULT_SAVER_HH

#include <PEL_Object.hh>

#include <PDE_DomainAndFields.hh>
#include <size_t_vector.hh>

class PEL_List ;
class PEL_ModuleExplorer ;
class intArray2D ;
class doubleArray2D ;
class doubleVector ;
class intVector ;

class GE_Mpolyhedron ;
class GE_Point ;
class GE_SetOfPoints ;

class PDE_DiscreteField ;
class PDE_SetOfDiscreteFields ;

/* 
Servers that save data in a file for subsequent post-processing.

The resulting file is organized in successive cycles.
Each cycle contains the values of a set of variables.
Variables have a name and may be of different types.

INSTANCE DELIVERY : `PDE_DomainAndFields::result_saver'
*/

class PEL_EXPORT PDE_ResultSaver : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Characteristics

      // the instance of  `PDE_DomainAndFields::' that is attached to `self'
      PDE_DomainAndFields const* attached_domain( void ) const ;
      
   //-- Cycles
      
      // Start a new cycle.
      void start_cycle( void ) ;

      // Terminate the current cycle.
      void terminate_cycle( void ) ;

      // Is there a cycle that is started and not terminated ?
      bool has_an_opened_cycle( void ) const ;

      // cycle number
      size_t cycle_number( void ) const ;
      
   //-- Variables(20.)

      // Has a value calling `name' already  been saved ?
      bool has_variable( std::string const& name ) const ;
      
      // Save `value' calling it `name'.
      void save_variable( double const& value,
                          std::string const& name ) ;

      // Save `value' calling it `name'.
      void save_variable( doubleVector const& value,
                          std::string const& name ) ;

      // Save `value' calling it `name'.
      void save_variable( doubleArray2D const& value,
                          std::string const& name ) ;

      // Save `value' calling it `name'.
      void save_variable( int const& value,
                          std::string const& name ) ;

      // Save `value' calling it `name'.
      void save_variable( intVector const& value,
                          std::string const& name ) ;

      // Save 'time' as a variable named "TIME".
      void save_time( double time ) ;
      
   //-- Grid: default management(30.)
      
      // Save the geometrical grid.
      void save_grid( void ) ;

      bool grid_is_saved( void ) const ;

      size_t nb_saved_cells( void ) const ;
      
      size_t nb_saved_faces( void ) const ;
      
      size_t nb_saved_vertices( void ) const ;
      
      size_t vertex_index( GE_Point const* vertex ) const ;
      
   //-- Grid: customized management(31.)
      
      static void update_active_vertices( GE_SetOfPoints const* vertices,
                                          GE_Mpolyhedron const* poly,
                                          size_t_vector& active_idx,
                                          size_t& nb_active ) ;
                                          
      // Modify `cell_nb_faces', `cell2face', `face2cell' to indicate 
      // that the il-th adjancent cell of `icell' is `iface'.
      void set_face_cell_link( size_t iface, size_t il, size_t icell,
                               intVector& cell_nb_faces, intArray2D& cell2face,
                               intArray2D& face2cell ) const ;

      void save_grid(
          size_t ndims, doubleArray2D const& vertices,
          intVector const& cell_nb_vertices, intArray2D const& cell2vertex,
          intVector const& cell_nb_faces,    intArray2D const& cell2face,
          intVector const& face_nb_vertices, intArray2D const& face2vertex,
          intArray2D const& face2cell,
          intVector const& vertex_color,
          intVector const& cell_color,
          intVector const& face_color ) ;

   //-- Fields: default management(35.)
      
      // For some items of the associated `PDE_SetOfDiscreteFields::' object, 
      // save a variable defined by the values of a reconstruction at 
      // some discrete locations in the geometrical grid (based on the 
      // `level'-th storage their DOFs). 
      // The considered items, the name of the associated variables, the 
      // discrete locations are given by the `PEL_ModuleExplorer::' object 
      // used when creating `self'.
      void save_fields( size_t level ) ;

   //-- Fields: customized management(35.)
         
      enum SavingLocation
      {
         AtVertices,
         AtCellCenters,
         InvalidLocation
      } ;
      
      void prepare_for_field_saving( SavingLocation location,
                                     std::string const& save_name,
                                     size_t nb_components,
                                     doubleArray2D& values,
                                     doubleVector& default_values ) ;
      
      bool field_saving_in_progress( void ) const ;
      
      SavingLocation saving_location( void ) const ;
      
      void save_field( doubleArray2D& values,
                       doubleVector const& default_values ) ;
      
   //-- Checks(40.)
      
      static double undefined_value( void ) ;

      /* 
      Check that the two values `old_val' and `new_val' at the same
      vertex `vertex' (but from two different cells sharing that vertex) 
      are identical.
      Remarks :
         - no check is performed if `old_val' = `::undefined_value' ;
         - the equality of `old_val' and `new_val' is checked calling
           `PEL::double_equality' with `dbl_eps' and `dbl_min' as 
            last two parameters.
      */
      static void check_value_consistency_at_vertex(
                       std::string const& banner,
                       GE_Point const* vertex,
                       double old_val, double new_val,
                       double dbl_eps, double dbl_min ) ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_ResultSaver( void ) ;
     ~PDE_ResultSaver( void ) ;
      PDE_ResultSaver( PDE_ResultSaver const& other ) ;
      PDE_ResultSaver& operator=( PDE_ResultSaver const& other ) ;

      friend PDE_ResultSaver* PDE_DomainAndFields::result_saver( void ) const ;
      friend PDE_ResultSaver* PDE_DomainAndFields::novel_result_saver( 
                                        PEL_ModuleExplorer const* exp ) const ;

      PDE_ResultSaver( PDE_DomainAndFields* dom,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Internals
      
      void build_reconstruction( PDE_DiscreteField const* field,
                                 size_t level,
                                 SavingLocation location,
                                 doubleArray2D& values,
                                 doubleVector& default_values ) ;

      void infer_active_grid_items( void ) ;
      
      void set_mesh_vertices( GE_Mpolyhedron const* poly,
                              size_t imesh,
                              intVector& mesh_nb_vertices,
                              intArray2D& mesh2vertex ) const ;
      
      std::string const& location2string( void ) const ;
      
   //-- Attributes

      PDE_DomainAndFields const* MY_DOM ;
      PEL_ModuleExplorer* MY_EXP ;
      bool OPENED_CYCLE ;
      size_t I_CYCLE ;
      GE_SetOfPoints const* VERTICES ;
      PDE_SetOfDiscreteFields const* FIELDS ;
      PEL_List* WRITERS ;
      PEL_Module* MOD ;
      PEL_Module* MOD_VARS ;
      PEL_Module* MOD_FF ;
      int NB_VARS ;
      int NB_FIELDS ;

      double const EPS_DBL ;
      double const MIN_DBL ;

      std::string MESHING_NAME ;
      size_t NB_ACTIVE_VERTS ;
      size_t NB_ACTIVE_CELLS ;
      size_t NB_ACTIVE_FACES ;
      size_t_vector ACTIVE_IDX_OF_VERT ;
      size_t_vector ACTIVE_IDX_OF_CELL ;
      bool PERIODIC ;
      bool DETECT_ACTIVE_ITEMS ;

      PDE_LocalFEcell* cFE ;
      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      
      SavingLocation FIELD_LOCATION ;
      std::string FIELD_NAME ;
} ;

#endif
