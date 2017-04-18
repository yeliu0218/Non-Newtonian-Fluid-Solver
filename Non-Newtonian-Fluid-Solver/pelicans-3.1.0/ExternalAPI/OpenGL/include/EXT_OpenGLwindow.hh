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

#ifndef EXT_OpenGLwindow_HH
#define EXT_OpenGLwindow_HH

#include <PEL_Object.hh>

#ifdef __APPLE__
   #include <GLUT/glut.h>  // for Max OS X
#else
   #include <GL/glut.h>    // for other systems
#endif

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>

class PEL_Vector ;
class size_t_vector ;
class size_t_array2D ;
class EXT_OpenGLviewer ;


// EXT_OpenGLwindow facilities.
// This toolbox class is designed to display fields and meshing on screen.
// Simple graphic user interface is provide usin GLUT facilities.

class EXT_OpenGLwindow : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      static EXT_OpenGLwindow* create( PEL_Object* a_owner,
                                       std::string const& name ) ;

      //-- Access

      // Is `self' current active window ?
      bool is_active( void) const ;

      void do_user_interface_loop( EXT_OpenGLviewer * viewer ) ;

      // Adjust `self' title.
      void set_title( std::string const& title ) ;

      // Declare color tablt index.
      void set_color_table( stringVector const& a_color_table,
                            intArray2D const& color_table_connectivity ) ;

      // Display set of points.
      void display_vertices( std::string const& name,
                             doubleArray2D const& vertices,
                             intVector const& vertex_color ) ;
      // Display meshing.
      void display_meshing( std::string const& name,
                            doubleArray2D const& vertices,
                            intArray2D const& connectivity,
                            intVector const& mesh_color,
                            size_t elem_space_dim  ) ;
      // Display discrete field.
      void display_field( std::string const& name,
                          std::string const& location,
                          doubleArray2D const& vertices,
                          intArray2D const& connectivity,
                          doubleArray2D const& values ) ;
   //-- Graphic primitives

      // Add simple polygon.
      void add_poly( GLenum mode,
                     doubleArray2D const& vertices,
                     size_t_vector const& poly,
                     size_t_array2D const& splitting,
                     doubleArray2D const& rgb,
                     double default_z_coordinate=0.0  ) ;
      // Reset defaults.
      void init_view( void ) ;
      // Fit view to window
      void auto_crop( void );

   //-- GLUT event handling

      void exit( void ) ;
      void display( void ) ;
      void idle( void ) ;
      void menu( int value ) ;
      void timer( int value ) ;
      void menuStatus(int status, int x, int y ) ;
      void mouse_click(int button, int state, int x, int y) ;
      void mouse_move(int x, int y ) ;
      void reshape(int width, int height ) ;
      void keyboard(unsigned char k, int x, int y) ;

      enum CycleMove{
             PREVIOUS_CYCLE=-6, NEXT_CYCLE=-7, FIRST_CYCLE=-8, LAST_CYCLE=-9,
             NO_CYCLE_CHANGE = -10 } ;
      enum { QUIT=-1, SURFACE=-2, HELP=-3, RESET=-4, CROP=-5,
             FIELDS=1000, MESHINGS=2000, COLORS=3000 } ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~EXT_OpenGLwindow(void ) ;
      EXT_OpenGLwindow( EXT_OpenGLwindow const& other ) ;
      EXT_OpenGLwindow& operator=( EXT_OpenGLwindow const& other ) ;

      EXT_OpenGLwindow( PEL_Object* a_owner,
                        std::string const& name ) ;

      static void split_polyhedron( size_t dim,
                                    size_t nb_vert,
                                    size_t_array2D& splitting ) ;
      void active( void ) ;
      void deal_with_selection_menu( void ) ;
      static void get_color( double &r, double &g, double &b,
                             double valn ) ;
      void compute_barycenter( doubleArray2D const& vertices ) ;

      static void display_help( std::ostream& os ) ;
      static void inverse_coordinate( double winx, double winy,
                                      double &x, double &y, double &z ) ;
      void  add_point( doubleArray2D const& vertices,
                       size_t idx ) ;
      void drawAxes(double length) ;
      void display_scalar( std::string const& location,
                           doubleArray2D const& vertices,
                           intArray2D const& connectivity,
                           doubleVector const& values ) ;
      void display_vector( std::string const& location,
                           doubleArray2D const& vertices,
                           intArray2D const& connectivity,
                           doubleArray2D const& values,
                           doubleVector const& valn ) ;
      void draw_arrow( doubleVector const& origin,
                       doubleVector const& vn,
                       doubleVector const& rgb ) ;
      void manage_sub_menu( int& submenu,
                            size_t& submenu_nb,
                            stringVector const& name,
                            size_t start_index,
                            bool & changed ) ;

      void set_projection( void ) ;

      void declare_meshing( std::string const& name ) ;

      void active_viewer( void ) ;

      //-- Window ID
      int ID ;
      static int MAIN ;

      //-- Viewpoint management
      GLdouble distance, elevation, azimuth;

      //-- Mouse events
      GLdouble x_0, y_0;
      bool moving ;
      bool scaling ;
      bool shift ;

      //-- View parameters
      GLdouble minWH;
      doubleVector barycenter ;
      doubleVector view ;
      double scale_x ;
      double scale_y ;
      double scale_z ;
      double aspect ;
      double a_distance ;
      double min_distance ;
      bool disp_3d ;
      size_t dim_space ;
      static const double z_upper ;
      doubleVector upper ;
      doubleVector lower ;
      int view_menu ;
      int cycle_menu ;

      //-- Window title
      std::string TITLE ;

      //-- Fields management
      stringVector fields_name ;
      int field_list ;
      intVector fields_list ;
      int field_sub_menu ;
      size_t field_sub_menu_nb ;
      bool surface_rendering ;
      bool surface_menu_status ;

      //-- Meshing management
      int meshing_sub_menu ;
      size_t meshing_sub_menu_nb ;
      intArray2D grid_list ;
      boolVector view_grid ;
      stringVector meshing_name ;

      //-- Color management
      boolVector used ;
      stringVector color_table ;
      intArray2D color_table_connectivity ;
      stringVector display_colors ;
      int color_sub_menu ;
      size_t color_sub_menu_nb ;
      int view_color ;

      //-- Miscellanous
      int selection_menu ;
      bool open_menu ;
      bool refresh ;
      std::string my_name ;
      EXT_OpenGLviewer * viewer ;

      CycleMove current_cycle ;
};

#endif
