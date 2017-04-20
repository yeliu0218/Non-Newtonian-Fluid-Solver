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

#ifndef PDE_DOMAIN_VISU_HH
#define PDE_DOMAIN_VISU_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;
class boolVector ;
class size_t_vector ;

class GE_SetOfPoints ;

class PDE_CursorFEside ;
class PDE_LocalFE ;
class PDE_LocalFEcell ;
class PDE_LocalFEbound ;
class PDE_ResultSaver ;

#include <iosfwd>
#include <fstream>
#include <set>
#include <string>

/*
Application that creates an instance of `::PDE_DomainAndFields'
and saves, for subsequent post-processing:
   1. the meshing and the fields, as requested in the MODULE PDE_ResultSaver
      of the associated `::PDE_DomainAndFields' instance creational process ;
   2. a GMV file called "cells.gmv" that makes possible the visualization
      of the colors of the cells.
   3. a GMV file called "faces.gmv" that makes possible the visualization
      of the colors of the sides and of the bounds.
   4. a GMV file called "bounds.gmv" that makes possible the visualization
      of the colors of the bounds.

Example

   Consider the following data deck (only the relevant parts are shown):

   MODULE PEL_Application
      concrete_name = "PDE_DomainVisu"
      MODULE PDE_DomainAndFields
         ...
         MODULE GE_Meshing
            concrete_name = "GE_MefistoMeshing"
            mesh_polyhedron = < "GE_Triangle" "GE_Tetrahedron" >
            filename = "xyznpef.COQUE"
         END MODULE GE_Meshing
         MODULE interior_fields
            ...
         END MODULE interior_fields
         MODULE PDE_ResultSaver      
            writers = < "PEL_TICwriter" >
            files_basename = "save"
            ...
         END MODULE PDE_ResultSaver
      END MODULE PDE_DomainAndFields
   END MODULE PEL_Application
      
   Three files will be created:
      * "save.gene" for visualization of the fields initialization
      * "cells.gmv" for visualization of the colors of the cells
      * "faces.gmv" for visualization of the colors 
                    of the sides and of the bounds
      * "bound.gmv" for visualization of the colors of the bounds
                              
   In "faces.gmv", a surface (GMV notion) is created that contains 
   all the faces (either sides or bounds). Moreover, for each face color, a 
   surface field (GMV notion) is created, whose value is the identifier
   of that color (given by GE_Color::identifier) on the faces of that 
   color and -1 on all other faces. Finally, a surface field called
   all?colors (where ? is the number of face colors) is created, whose
   value is, on each face, the identifier of the color of that face.
   
   Typical steps to visualize faces of a given color, say "ma_jolie_couleur:
      1. Execute PDE_DomainVisu to produce the file "faces.gmv".
      2. Open the file "faces.gmv" with GMV.
      3. In the Display menu, choose Surfaces .
      4. In the Surfaces menu, choose Faces and Color by the Surf Field
         called ma_jolie_couleur. The faces whose color is "ma_jolie_couleur"
         will be displayed with the field value equal to the identifier of
         "ma_jolie_couleur", and the other faces with the field value -1.
         
   The file "bounds.gmv" is similar to "faces.gmv" except that it contains
   only informations relatives to the bounds (and is thus much smaller).

   In "cells.gmv",  for each cell color, a cell field (GMV notion) is created,
   whose value is the identifier of that color (given by GE_Color::identifier)
   on the cells of that color and -1 on all other cells. Finally, a cell field
   called all?colors (where ? is the number of cell colors) is created, whose
   value is, on each cell, the identifier of the color of that cell.
   
   Typical steps to visualize cells of a given color, say "ma_jolie_couleur:
      1. Execute PDE_DomainVisu to produce the file "cells.gmv".
      2. Open the file "cells.gmv" with GMV.
      3. In the Display menu, choose Cells .
      4. In the Cells menu, choose Color By: , then Cell Field
         and select ma_jolie_couleur. The cells whose color is 
         "ma_jolie_couleur" will be displayed with the field value equal to 
         the identifier of "ma_jolie_couleur", and the other cells with the 
         field value -1.

PUBLISHED
*/

class PEL_EXPORT PDE_DomainVisu : public PEL_Application
{
   
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_DomainVisu* create( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) ;

   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PDE_DomainVisu( void ) ;
      PDE_DomainVisu( PDE_DomainVisu const& other ) ;
      PDE_DomainVisu& operator=( PDE_DomainVisu const& other ) ;

      PDE_DomainVisu( PEL_Object* a_owner, 
		      PEL_ModuleExplorer const* exp ) ;

    //-- Plug in

      PDE_DomainVisu( void ) ;

      virtual PDE_DomainVisu* create_replica( 
	                         PEL_Object* a_owner,
	                         PEL_ModuleExplorer const* exp ) const ;

    //-- Internals

      void build_file_with_bounds( void ) ;
      void build_file_with_faces( void ) ;
      void build_file_with_cells( void ) ;
      
      void build_skin( GE_SetOfPoints* skin_vertices,
                       size_t& nb_skin_cells,
                       size_t_vector& glob_cell_2_skin_cell ) ;

      void write_gmv_nodev( std::ostream& ofile,
                            GE_SetOfPoints const* verts ) const ;
      void write_gmv_faces( std::ostream& ofile ) ;
      void write_gmv_faces( std::ostream& ofile, 
                            GE_SetOfPoints const* skin_vertices,
                            size_t nb_skin_cells,
                            size_t_vector const& glob_cell_2_skin_cell ) ;
      void write_gmv_cells( std::ostream& ofile ) ;
      void write_gmv_surface( std::ostream& ofile,
                              GE_SetOfPoints const* verts,
                              bool with_sides ) const ;
      void write_gmv_surfvars( std::ostream& ofile, bool with_sides ) const ;
      void write_gmv_variable( std::ostream& ofile ) const ;

      void display_color_info( std::string const& col_name,
                               size_t nb_s, size_t nb_b ) const ;
      void display_color_info( std::string const& col_name,
                               size_t nb_c ) const ;
      
      void append_matching_color_names( 
                                 PDE_LocalFE* fe,
                                 std::set< std::string >& colors ) const ;

      void append_matching_color_names( 
                                 PDE_CursorFEside* fe,
                                 std::set< std::string >& colors ) const ;

      void check_connexity( void ) ;

      void append_neigh( size_t icell, boolVector& ok_cell ) ;

    //-- Class attributes

      static PDE_DomainVisu const* PROTOTYPE ;
  
    //-- Attributes

      std::ofstream OFS ;
      bool TRACE_CELLS ;
      bool TRACE_SIDES ;
      bool TRACE_BOUNDS ;
      bool CHECK_CONN ;
      size_t NB_DIMS ;
      GE_SetOfPoints const* VERTICES ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;
      PDE_CursorFEside* sFE ;
      PDE_ResultSaver* SAVER ;
      size_t NBS ;
      bool GMV_FILES ;
} ;

#endif
