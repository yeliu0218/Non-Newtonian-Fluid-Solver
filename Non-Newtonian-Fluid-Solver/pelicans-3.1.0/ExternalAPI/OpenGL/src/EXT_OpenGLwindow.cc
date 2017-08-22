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

#include <EXT_OpenGLwindow.hh>

#include <iostream>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>

#include <PEL.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>
#include <PEL_Vector.hh>
#include <EXT_GlutAPI.hh>
#include <EXT_OpenGLviewer.hh>

#include <iostream>
using std::endl ;

//----------------------------------------------------------------------
int EXT_OpenGLwindow:: MAIN = -1 ;
double const EXT_OpenGLwindow::z_upper = 1.0e-5 ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
static EXT_OpenGLwindow* Windows ;
//----------------------------------------------------------------------
extern "C"
{

//----------------------------------------------------------------------
EXT_OpenGLwindow*
my_object(void)
//----------------------------------------------------------------------
{
   return Windows  ;
}

//----------------------------------------------------------------------
static void
reshapeFunc( int width, int height )
//----------------------------------------------------------------------
{
   my_object()->reshape( width, height ) ;
}

//----------------------------------------------------------------------
static void
menuStatusFunc( int status, int x, int y )
//----------------------------------------------------------------------
{
   my_object()->menuStatus( status, x, y ) ;
}

//----------------------------------------------------------------------
static void
menuFunc( int value )
//----------------------------------------------------------------------
{
   my_object()->menu( value ) ;
}

//----------------------------------------------------------------------
static void
timerFunc( int value )
//----------------------------------------------------------------------
{
   my_object()->timer( value ) ;
   glutTimerFunc( value, timerFunc, value ) ;
}

//----------------------------------------------------------------------
static void
displayFunc(void)
//----------------------------------------------------------------------
{
   my_object()->display() ;
}

//----------------------------------------------------------------------
static void
mouseFunc(int button, int state, int x, int y)
//----------------------------------------------------------------------
{
   my_object()->mouse_click(button,state,x,y) ;
}

//----------------------------------------------------------------------
static void
mouseMotionFunc(int x, int y)
//----------------------------------------------------------------------
{
   my_object()->mouse_move(x,y) ;
}

//----------------------------------------------------------------------
static void
KeyboardFunc(unsigned char k, int x, int y)
//----------------------------------------------------------------------
{
   switch(k)
   {
      case ' ' :
         my_object()->keyboard(k,x,y) ;
         break ;
      case 'q' :
         my_object()->exit() ;
   }
}

//----------------------------------------------------------------------
} // extern "C"
//----------------------------------------------------------------------

//----------------------------------------------------------------------
EXT_OpenGLwindow*
EXT_OpenGLwindow:: create( PEL_Object* a_owner,
                           std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: create" ) ;

   EXT_OpenGLwindow* result = new EXT_OpenGLwindow( a_owner, name ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_OpenGLwindow:: EXT_OpenGLwindow( PEL_Object* a_owner,
                                     std::string const& name )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , barycenter( 3 )
   , view( 3 )
   , a_distance( 1.0 )
   , min_distance( 1.0 )
   , disp_3d( false )
   , dim_space( 0 )
   , upper(3)
   , lower(3)
   , view_menu( 0 )
   , cycle_menu( 0 )
   , fields_name(0)
   , fields_list(0)
   , field_sub_menu(0)
   , surface_rendering( false )
   , meshing_sub_menu(0)
   , grid_list(0,0)
   , view_grid(0)
   , meshing_name(0)
   , used( 0 )
   , color_table(0)
   , color_table_connectivity(0,0)
   , display_colors(0)
   , color_sub_menu(0)
   , view_color( -1 )
   , selection_menu(0)
   , open_menu( false )
   , my_name( name )
   , current_cycle( FIRST_CYCLE )
{
   EXT_GlutAPI::initialize_on_demand() ;
   lower.set( PEL::max_double() ) ;
   upper.set( -PEL::max_double() ) ;

   barycenter.set( 0.0 ) ;
   field_list = 0 ;

   ID = glutCreateWindow( name.c_str() ) ;
   EXT_GlutAPI::disable_window_destroy() ;

   glutDisplayFunc(displayFunc);
   glutReshapeFunc(reshapeFunc);
   glutMouseFunc(mouseFunc);
   glutMotionFunc(mouseMotionFunc);
   glutKeyboardFunc(KeyboardFunc);
   glutMenuStatusFunc(menuStatusFunc);
   glutTimerFunc(1000,timerFunc,1000);
   Windows = this ;
   MAIN = ID ;
   init_view() ;
   minWH = PEL::min( glutGet(GLUT_SCREEN_WIDTH),
                     glutGet(GLUT_SCREEN_HEIGHT) ) ;
   glEnable( GL_DEPTH_TEST ) ;
   glBlendFunc( GL_SRC_ALPHA, GL_DST_ALPHA ) ;

   glPointSize( 5.0 ) ;
}

//----------------------------------------------------------------------
EXT_OpenGLwindow:: ~EXT_OpenGLwindow( void )
//----------------------------------------------------------------------
{
   refresh=true ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: exit( void )
//----------------------------------------------------------------------
{
  std::cout << "Exiting ... " << std::endl;
  glutDestroyWindow(ID) ;
  PEL_System::exit(1);
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: reshape(int width, int height)
//----------------------------------------------------------------------
{
   aspect = (float) width / (float) height;
   minWH = PEL::min (width, height);

   glViewport(0, 0, (GLint) width, (GLint) height);
   set_projection() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: set_projection( void )
//----------------------------------------------------------------------
{
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   gluPerspective(60.0, aspect, 0.001, 1000.0);

   glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: mouse_click(int button, int state, int x, int y)
//----------------------------------------------------------------------
{
   moving = scaling = false ;
   shift = glutGetModifiers()==GLUT_ACTIVE_SHIFT ;

   if (button == GLUT_LEFT_BUTTON) {
      if (state == GLUT_DOWN) {
         moving = true;
      }
   }
   if (button == GLUT_MIDDLE_BUTTON) {
      if (state == GLUT_DOWN) {
         scaling = true;
      }
   }
   x_0 = x;
   y_0 = y;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: inverse_coordinate( double winx, double winy,
                             double &x, double &y, double &z )
//----------------------------------------------------------------------
{
   GLdouble P[16] ;
   GLdouble M[16] ;
   GLint V[4] ;
   glGetDoublev(GL_PROJECTION_MATRIX,P) ;
   glGetDoublev(GL_MODELVIEW_MATRIX,M) ;
   glGetIntegerv(GL_VIEWPORT,V) ;

   // Recovering origin projection Z coordinate
   x=y=z=0.0;
   double winz ;
   double dum ;

   gluProject( x, y, z, M, P, V, &dum, &dum, &winz ) ;
   gluUnProject( winx, winy, winz, M, P, V, &x, &y, &z ) ;

}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: mouse_move(int x, int y)
//----------------------------------------------------------------------
{
   if( scaling && ! shift )
   {
      // Uniform scaling
      double amount = ((double)(x-x_0))/minWH;
      distance *= (1.0+amount);
   }

   if( scaling && shift )
   {
      // No uniform scaling
      scale_x =  PEL::max( scale_x + ((double)(x-x_0))/100, 0.01 ) ;
      scale_y =  PEL::max( scale_y + ((double)(y-y_0))/100, 0.01 ) ;
      set_projection() ;
   }

   if( moving && shift )
   {
      // Rotation
      azimuth += (x-x_0)*PEL::pi()/minWH;
      if( disp_3d ) elevation += (y-y_0)*PEL::pi()/minWH;
   }

   if( moving && !shift )
   {
      // Translation
      GLdouble dx0, dy0, dz0 ;
      GLdouble dx, dy, dz ;
      int h = glutGet(GLUT_WINDOW_HEIGHT) ;

      double wx=x_0,wy=h-y_0 ;
      inverse_coordinate( wx,wy,dx0,dy0,dz0 ) ;
      wx=x ;wy=h-y ;
      inverse_coordinate( wx,wy,dx,dy,dz ) ;

      view(0) -= dx-dx0 ;
      view(1) -= dy-dy0 ;
      view(2) -= dz-dz0 ;
   }
   x_0 = x;
   y_0 = y;
   glutPostRedisplay() ;

}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow::do_user_interface_loop( EXT_OpenGLviewer* a_viewer )
//----------------------------------------------------------------------
{
   viewer = a_viewer ;

   active_viewer() ;
   display() ;

   glutMainLoop() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow::display( void )
//----------------------------------------------------------------------
{
   glutSetWindowTitle( TITLE.c_str() ) ;
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   glLoadIdentity();

   glScaled( scale_x, scale_y, scale_z ) ;


   double cosel=PEL::cos(elevation);
   // glTranslated( -view(0), -view(1), -view(2) ) ;

   gluLookAt( view(0)+distance*cosel*PEL::cos(azimuth),
              view(1)+distance*cosel*PEL::sin(azimuth),
              view(2)+distance*PEL::sin(elevation)
              , view(0),view(1),view(2)
              , PEL::cos(elevation+PEL::pi()/2.0)*PEL::cos(azimuth)
              , PEL::cos(elevation+PEL::pi()/2.0)*PEL::sin(azimuth)
              , PEL::sin(elevation+PEL::pi()/2.0) ) ;

//***************
//   11 Mai 2006 : les deux lignes suivantes ont été commentées pour éviter
//      l'affichage répétitif d'une erreur qui existe mais qui ne semble pas
//      nuire au déroulement de l'application. Affaire à instruire...
//***************
//   GLenum err  = glGetError() ;
//   if( err!=GL_NO_ERROR ) std::cout<< endl << "ERROR : "<<err<<std::endl ;

   glColor3f(1.0,1.0,1.0);
   drawAxes(1.0);
   if( field_list>0 ) glCallList(field_list) ;

   for( size_t color=0 ; color<color_table.size() ; color++ )
   {
      if( view_color==-1 || view_color==(int)color)
      {
         glColor3f(1.0,0.0,0.0);
         glLineWidth(2.0);
      }
      else
      {
         glColor3f(0.5,0.5,0.5);
         glLineWidth(1.0);
      }
      for( size_t i=0 ; i<view_grid.size() ; i++ )
         if( view_grid(i)  )
            glCallList(grid_list(i,color));
   }
   glLineWidth(1.0);

   glFlush();
   glutSwapBuffers();
   if( !open_menu )
   {
      deal_with_selection_menu();
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: drawAxes(double length)
//----------------------------------------------------------------------
{
   glBegin(GL_LINES);
   glVertex3f(0.0,0.0,0.0);
   glVertex3f(length,0.0,0.0);

   glVertex3f(0.0,0.0,0.0);
   glVertex3f(0.0,length,0.0);
   glVertex3f(0.0,0.0,0.0);
   glVertex3f(0.0,0.0,length);
   glEnd();
   glRasterPos3f(length,0.0,0.0);
   glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'X');
   glRasterPos3f(0.0,length,0.0);
   glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'Y');
   glRasterPos3f(0.0,0.0,length);
   glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'Z');
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: keyboard(unsigned char k, int x, int y)
//----------------------------------------------------------------------
{
   init_view() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: menuStatus(int status, int x, int y )
//----------------------------------------------------------------------
{
   open_menu = status==GLUT_MENU_IN_USE ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: menu( int value )
//----------------------------------------------------------------------
{
   if( value==QUIT )
   {
      exit() ;
   }
   else if( value==SURFACE )
   {
      surface_rendering = ! surface_rendering ;
      if( dim_space==3 )
      {
         if( surface_rendering)
            glEnable(GL_BLEND) ;
         else
            glDisable(GL_BLEND) ;
      }
   }
   else if( value==HELP )
   {
      display_help( std::cout ) ;
   }
   else if( value==RESET )
   {
      init_view() ;
   }
   else if( value==CROP )
   {
      auto_crop() ;
   }
   else if( value==FIRST_CYCLE || value==LAST_CYCLE ||
            value==PREVIOUS_CYCLE || value==NEXT_CYCLE )
   {
     current_cycle = (CycleMove) value ;
   }
   else
   {
      if( value>=FIELDS && value<(int)fields_list.size()+FIELDS )
      {
         int f = value-FIELDS ;
         field_list =
            ( field_list == fields_list(f) ? 0 : fields_list(f) ) ;
      }
      else if( value>=MESHINGS && value<MESHINGS+(int)meshing_name.size() )
      {
         int grid = value-MESHINGS;
         view_grid(grid) = !view_grid(grid) ;
      }
      else if( value>=COLORS && value<COLORS+(int)color_table.size()+1 )
      {
         size_t color = value-COLORS ;
         std::string sel = display_colors( color ) ;
         if( sel != "all" )
            view_color = color_table.index_of( sel ) ;
         else
            view_color = -1 ;
      }
   }
   refresh = true ;
   active_viewer() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: timer( int value )
//----------------------------------------------------------------------
{
   active_viewer() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: active_viewer( void )
//----------------------------------------------------------------------
{
   viewer->go_to( current_cycle ) ;
   current_cycle = NO_CYCLE_CHANGE ;
   if( refresh )
   {
      glutPostRedisplay() ;
      refresh=false ;
   }
}

//----------------------------------------------------------------------
bool
EXT_OpenGLwindow:: is_active( void) const
//----------------------------------------------------------------------
{
   return MAIN==ID ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: init_view( void )
//----------------------------------------------------------------------
{
   view = barycenter ;
   distance = 2.0 ;
   elevation = PEL::pi()/2.0 ;
   azimuth = -PEL::pi()/2.0 ;
   scale_x = scale_y = scale_z = 1.0 ;
   if( upper(0)>lower(0) && upper(1)>lower(1) )
   {
      double s = PEL::max_double() ;
      if(  upper(0)-lower(0) > 0.0 )
         s = 1.0/( upper(0)-lower(0)) ;
      if(  upper(1)-lower(1) > 0.0 )
         s = PEL::min( s, 1.0/( upper(1)-lower(1) ) ) ;
      if(  upper(2)-lower(2) > 0.0 )
         s = PEL::min( s, 1.0/( upper(2)-lower(2) ) ) ;
      scale_x = scale_y = scale_z = s ;
      distance = 2.0/s ;
   }

   set_projection() ;
   glutPostRedisplay() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: auto_crop( void )
//----------------------------------------------------------------------
{
   view = barycenter ;
   elevation = PEL::pi()/2.0 ;
   azimuth = -PEL::pi()/2.0 ;
   if( upper(0)>lower(0) && upper(1)>lower(1) )
   {
      if(  upper(0)-lower(0) > 0.0 )
         scale_x = 1.0/( upper(0)-lower(0)) ;
      if(  upper(1)-lower(1) > 0.0 )
         scale_y = 1.0/( upper(1)-lower(1)) ;
      if(  upper(2)-lower(2) > 0.0 )
         scale_z = 1.0/( upper(2)-lower(2)) ;
      distance = 2.0/scale_z ;
   }
   set_projection() ;
   glutPostRedisplay() ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: active( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: active" ) ;
   if( MAIN!=ID ) glutSetWindow( ID ) ;
   MAIN = ID ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: set_title( std::string const& title )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: set_title" ) ;

   TITLE = my_name + " - " + title ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: add_poly( GLenum mode,
                             doubleArray2D const& vertices,
                             size_t_vector const& poly,
                             size_t_array2D const& splitting,
                             doubleArray2D const& rgb,
                             double default_z_coordinate )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: add" ) ;
   PEL_CHECK_PRE( poly.size()==rgb.index_bound(1) ||
                  1==rgb.index_bound(1) ||
                  0==rgb.index_bound(1) ) ;

   static bool d3 = vertices.index_bound(0)==3 ;
   bool uniform = rgb.index_bound(1)==1 ;
   bool colored = rgb.index_bound(1)==poly.size() ;
   bool surf = rgb.index_bound(0)==4 && !d3 ;
   bool blend = rgb.index_bound(0)==4 && d3 ;

   double v = 0.0 ;
   static double att = 1.0 ;

   if( uniform )
   {
      if( blend )
         glColor4d( rgb(0,0), rgb(1,0), rgb(2,0), att*rgb(3,0) ) ;
      else
         glColor3d( rgb(0,0), rgb(1,0), rgb(2,0) ) ;
      if( surf ) v = rgb(3,0)*min_distance ;
   }

   size_t nb = splitting.index_bound(0) ;
   size_t nbs = splitting.index_bound(1) ;

   for( size_t j = 0 ; j<nb ; j++ )
   {
      glBegin( mode );
      for( size_t k = 0 ; k<nbs ; k++ )
      {
         size_t l = splitting(j,k) ;
         size_t idx = poly(l) ;
         if( colored )
         {
            if( blend )
               glColor4d( rgb(0,l), rgb(1,l), rgb(2,l), att*rgb(3,l) ) ;
            else
               glColor3d( rgb(0,l), rgb(1,l), rgb(2,l) ) ;
            if( surf ) v = rgb(3,l)*min_distance ;
         }
         if(d3)
            glVertex3d( vertices(0,idx), vertices(1,idx), vertices(2,idx) ) ;
         else if( surf )
            glVertex3d( vertices(0,idx), vertices(1,idx), v ) ;
         else
            glVertex3d( vertices(0,idx), vertices(1,idx),  default_z_coordinate) ;
      }
      glEnd();
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: draw_arrow( doubleVector const& origin,
                               doubleVector const& vn,
                               doubleVector const& rgb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: draw_arrow" ) ;
   PEL_CHECK( origin.size()==3 || origin.size()==2 ) ;
   PEL_CHECK( origin.size()==vn.size() ) ;
   PEL_CHECK( rgb.size()==3 || rgb.size()==4 ) ;

   static size_t dim = origin.size() ;
   static bool d3 = dim==3 ;
   doubleVector const* v ;
   doubleVector xt(dim) ;
   doubleVector xv(dim) ;
   doubleVector x1(dim) ;
   xv.set(0.0) ;
   xv(0) = -vn(1) ;
   xv(1) = vn(0) ;
   glBegin(GL_LINES);
   glLineWidth(2.0);
   glColor3d( rgb(0), rgb(1), rgb(2) ) ;
   bool surf = rgb.size()==4 ;

   for(size_t i=0 ; i<dim; i++ )
   {
      xt(i)=origin(i)+vn(i) ;
   }
   for( size_t l=0 ; l<3 ; l++ )
   {
      for( size_t t=0 ; t<2 ; t++ )
      {
         if( t==0 )
         {
            v = &xt ;
         }
         else if( l==1 )
         {
            v = &origin ;
         }
         else
         {
            v = &x1 ;
            int sign = l-1 ;
            for(size_t i=0 ; i<dim; i++ )
            {
               x1(i)=xt(i)+(sign*xv(i)-vn(i))*0.2 ;
            }
         }
         if(d3)
            glVertex3d( (*v)(0), (*v)(1), (*v)(2) ) ;
         else if( surf )
            glVertex3d( (*v)(0), (*v)(1), rgb(3) ) ;
         else
            glVertex3d( (*v)(0), (*v)(1), z_upper ) ;
      }
   }

   glLineWidth(1.0);
   glEnd();
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: add_point( doubleArray2D const& vertices,
                              size_t idx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: add_point" ) ;

   static bool d3 = vertices.index_bound(0)==3 ;
   if(d3)
      glVertex3d( vertices(0,idx), vertices(1,idx), vertices(2,idx) ) ;
   else
      glVertex3d( vertices(0,idx), vertices(1,idx), z_upper ) ;

}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: split_polyhedron( size_t dim,
                                     size_t nb_vert,
                                     size_t_array2D& splitting )
//----------------------------------------------------------------------
{
   if( nb_vert==8 && dim==3 ) // Cube
   {
      splitting.re_initialize(6,4) ;
      size_t i=0 ;
      splitting(i,0)=0;splitting(i,1)=3;splitting(i,2)=2;splitting(i,3)=1;i++;
      splitting(i,0)=4;splitting(i,1)=5;splitting(i,2)=6;splitting(i,3)=7;i++;
      splitting(i,0)=1;splitting(i,1)=2;splitting(i,2)=6;splitting(i,3)=5;i++;
      splitting(i,0)=2;splitting(i,1)=3;splitting(i,2)=7;splitting(i,3)=6;i++;
      splitting(i,0)=0;splitting(i,1)=4;splitting(i,2)=7;splitting(i,3)=3;i++;
      splitting(i,0)=0;splitting(i,1)=1;splitting(i,2)=5;splitting(i,3)=4;i++;

   }
   else if( nb_vert==4 && dim==3 ) // Cube
   {
      splitting.re_initialize(4,3) ;
      size_t i=0 ;
      splitting(i,0)=0;splitting(i,1)=1;splitting(i,2)=2;i++;
      splitting(i,0)=0;splitting(i,1)=1;splitting(i,2)=3;i++;
      splitting(i,0)=0;splitting(i,1)=2;splitting(i,2)=3;i++;
      splitting(i,0)=1;splitting(i,1)=2;splitting(i,2)=3;i++;

   }
   else
   {
      splitting.re_initialize(1,nb_vert) ;
      for( size_t i=0 ; i<nb_vert ; i++ ) splitting(0,i)=i ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: compute_barycenter( doubleArray2D const& vertices )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: compute_barycenter" ) ;
   // Barycenter
   size_t dim = vertices.index_bound(0) ;
   barycenter.set(0.0) ;
   size_t n = vertices.index_bound(1) ;
   for( size_t i=0 ; i<n ; i++ )
   {
      for( size_t j=0 ; j<dim ; j++ )
      {
         double c = vertices(j,i) ;
         barycenter(j) += c/n ;
         if( c<lower(j) ) lower(j) = c ;
         if( c>upper(j) ) upper(j) = c ;
      }
   }
   min_distance = PEL::max_double() ;
   a_distance = 0.0 ;
   for( size_t j=0 ; j<dim ; j++ )
   {
      double d = ( upper(j)-lower(j) ) ;
      a_distance += d / dim ;
      min_distance = PEL::min( min_distance, d ) ;
   }

   static bool first = true ;
   if( first )
   {
      init_view() ;
      first = false ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: set_color_table( stringVector const& a_color_table,
                                    intArray2D const& a_color_table_connectivity )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: set_color_table" ) ;
   color_table = a_color_table ;
   color_table_connectivity = a_color_table_connectivity ;
   used.re_initialize( a_color_table.size() ) ;
   used.set(false) ;
   refresh=true ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_vertices( std::string const& name,
                                     doubleArray2D const& vertices,
                                     intVector const& vertex_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: display_vertices" ) ;
   PEL_CHECK_PRE( vertices.index_bound(0)>=2 && vertices.index_bound(0)<=3 ) ;
   PEL_CHECK_PRE( vertex_color.size() == vertices.index_bound(1) ) ;

   dim_space = PEL::max( vertices.index_bound(0), dim_space ) ;

   declare_meshing( name ) ;

   size_t meshing = meshing_name.index_of( name ) ;

   disp_3d |= vertices.index_bound(0)==3 ;
   compute_barycenter( vertices ) ;

   size_t nbvert = vertices.index_bound(1) ;
   size_t nb_color = color_table.size() ;

   doubleArray2D null_color(3,0) ;

   for( size_t color =0 ; color<nb_color ; color++ )
   {
      glNewList(grid_list(meshing,color), GL_COMPILE) ;
      glBegin(GL_POINTS);
      for( size_t i=0 ; i<nbvert ; i++ )
         if( color_table_connectivity(color,vertex_color(i)) )
         {
            used(color)=true;
            add_point( vertices, i ) ;
         }

      glEnd();

      glEndList();
   }
   refresh=true ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: declare_meshing( std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: declare_meshing" ) ;

   if( !meshing_name.has( name ) )
   {
      meshing_name.append( name ) ;
      view_grid.append(false) ;
   }
   size_t nb_color = color_table.size() ;
   if( grid_list.index_bound(0)<meshing_name.size() ||
       grid_list.index_bound(1)<nb_color )
   {
      intArray2D copy(grid_list) ;
      grid_list.re_initialize( meshing_name.size(), nb_color ) ;
      for( size_t i=0 ; i<copy.index_bound(0) ; i++ )
         for( size_t j=0 ; j<copy.index_bound(1) ; j++ )
            grid_list(i,j) = copy(i,j) ;
      for( size_t i=0 ; i<grid_list.index_bound(0) ; i++ )
         for( size_t j=0 ; j<grid_list.index_bound(1) ; j++ )
            if( grid_list(i,j)==0) grid_list(i,j)=glGenLists(1) ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_meshing( std::string const& name,
                                    doubleArray2D const& vertices,
                                    intArray2D const& connectivity,
                                    intVector const& mesh_color,
                                    size_t elem_space_dim )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: display_meshing" ) ;
   PEL_CHECK_PRE( vertices.index_bound(0)>=2 && vertices.index_bound(0)<=3 ) ;
   PEL_CHECK_PRE( mesh_color.size() == connectivity.index_bound(1) ) ;

   dim_space = PEL::max( vertices.index_bound(0), dim_space ) ;
   disp_3d |= vertices.index_bound(0)==3 ;
   size_t const nbv = connectivity.index_bound(0) ;
   size_t_vector poly(nbv) ;
   size_t_array2D splitting(0,0) ;
   split_polyhedron( elem_space_dim, nbv, splitting ) ;
   size_t nb_poly = connectivity.index_bound(1) ;
   size_t nb_color = color_table.size() ;

   declare_meshing( name ) ;

   doubleArray2D null_color(3,0) ;
   size_t meshing = meshing_name.index_of( name ) ;

   for( size_t color =0 ; color<nb_color ; color++ )
   {
      glNewList( grid_list(meshing,color), GL_COMPILE );

      for( size_t i=0 ; i<nb_poly ; i++ )
      {
         if( color_table_connectivity( color,mesh_color(i) ) )
         {
            for( size_t j=0 ; j<nbv ; j++ )
            {
               poly(j)=connectivity(j,i) ;
            }
            add_poly( GL_LINE_LOOP, vertices, poly, splitting, null_color, z_upper ) ;
            used(color)=true ;
         }
      }
      glEndList();
   }

   refresh=true ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: get_color( double &r, double &g, double &b,
                              double valn )
//----------------------------------------------------------------------
{
   r = valn ;
   b = 1.0-r ;
   g = r*(1-r) ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_field( std::string const& name,
                                  std::string const& location,
                                  doubleArray2D const& vertices,
                                  intArray2D const& connectivity,
                                  doubleArray2D const& values )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: display_field" ) ;
   PEL_CHECK_PRE( values.index_bound(1)==connectivity.index_bound(1) ||
                  values.index_bound(1)==vertices.index_bound(1) ) ;

   dim_space = PEL::max( vertices.index_bound(0), dim_space ) ;
   disp_3d |= vertices.index_bound(0)==3 || surface_rendering ;
   if( !fields_name.has(name) )
   {
      fields_name.append( name ) ;
      fields_list.append( glGenLists(1) ) ;
   }

   size_t idx = fields_name.index_of( name ) ;
   int afield_list = fields_list(idx) ;
   glNewList(afield_list, GL_COMPILE);
   size_t const component = values.index_bound(0) ;
   double min = PEL::max_double() ;
   double max = -PEL::max_double() ;
   size_t nbval = values.index_bound(1) ;

   doubleVector valn(nbval) ;

   for( size_t i=0 ; i<nbval ; i++ )
   {
      double v = 0.0 ;
      if( component==1 )
      {
         v = values(0,i) ;
      }
      else
      {
         for( size_t j=0 ; j<component ; j++ )
         {
            v += PEL::sqr(values(j,i)) ;
         }
         v = PEL::sqrt(v) ;
      }
      valn(i) = v ;
      if( v<min ) min = v ;
      if( v>max ) max = v ;
   }

   double dx = max-min ;
   if( dx < 1.0e-6 ) dx=1.0 ;
   for( size_t i=0 ; i<nbval ; i++ )
   {
      valn(i) = (valn(i)-min)/dx ;
   }

   if( component!=vertices.index_bound(0) )
   {
      display_scalar( location, vertices, connectivity, valn ) ;
   }
   else
   {
      display_vector( location, vertices, connectivity, values, valn ) ;
   }
   glEndList();
   refresh=true ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_scalar( std::string const& location,
                                   doubleArray2D const& vertices,
                                   intArray2D const& connectivity,
                                   doubleVector const& values )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: display_scalar" ) ;
   bool uniform = location=="at_cell_centers" ;

   size_t const nbv = connectivity.index_bound(0) ;
   size_t_vector poly(nbv) ;
   size_t nb_poly = connectivity.index_bound(1) ;
   size_t_array2D splitting(0,0) ;
   split_polyhedron( vertices.index_bound(0), nbv, splitting ) ;

   size_t dim = ( uniform ? 1 : nbv ) ;
   doubleArray2D rgb((surface_rendering?4:3),dim) ;

   for( size_t i=0 ; i<nb_poly ; i++ )
   {
      if( uniform )
      {
         double valn = values(i) ;
         get_color( rgb(0,0), rgb(1,0), rgb(2,0), valn ) ;
         if( surface_rendering ) rgb(3,0)=valn ;
      }

      for( size_t j=0 ; j<nbv ; j++ )
      {
         poly(j)=connectivity(j,i) ;
         if(!uniform)
         {
            double valn = values(poly(j)) ;
            get_color( rgb(0,j), rgb(1,j), rgb(2,j), valn  ) ;
            if( surface_rendering )
            {
               rgb(3,j)=valn ;
            }
         }
      }
      size_t nb_sides = splitting.index_bound(1) ;

      if( nb_sides==4 )
         add_poly( GL_QUADS, vertices, poly, splitting, rgb ) ;
      else if(nb_sides ==3 )
         add_poly( GL_TRIANGLES, vertices, poly, splitting, rgb ) ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_vector( std::string const& location,
                                   doubleArray2D const& vertices,
                                   intArray2D const& connectivity,
                                   doubleArray2D const& values,
                                   doubleVector const& valn )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLwindow:: display_vector" ) ;
   bool uniform = location=="at_cell_centers" ;

   size_t nbval = values.index_bound(1) ;
   doubleVector rgb((surface_rendering?4:3)) ;
   size_t const dim = vertices.index_bound(0) ;
   doubleVector origin(dim) ;
   size_t nbv = connectivity.index_bound(0) ;
   doubleVector vn(dim) ;
   double denorm =
      a_distance/( dim==3 ? PEL::pow(nbval,0.33) : PEL::sqrt(nbval) ) ;

   for( size_t p=0 ; p<nbval ; p++ )
   {
      origin.set(0.0) ;
      double norm = 0.0 ;
      for( size_t j=0 ; j<dim ; j++ )
      {
         if( uniform )
         {
            for( size_t i=0 ; i<nbv ; i++ )
            {
               size_t v = connectivity(i,p) ;
               origin(j)+=vertices(j,v)/nbv ;
            }
         }
         else
            origin(j) = vertices(j,p) ;
         vn(j) = values(j,p) ;
         norm += PEL::sqr(vn(j)) ;
      }
      if( norm>1.0e-10 )
         norm = valn(p)*denorm/PEL::sqrt(norm) ;
      for( size_t j=0 ; j<dim ; j++ )
      {
         vn(j)*=norm ;
      }

      // get_color( rgb(0), rgb(1), rgb(2), valn(p)  ) ;
      // if(surface_rendering) rgb(3)=valn(p) ;
      rgb.set(1.0) ;
      if(surface_rendering) rgb(3)=0.0 ;
      draw_arrow(origin,vn,rgb) ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: manage_sub_menu( int& submenu,
                                    size_t& submenu_nb,
                                    stringVector const& name,
                                    size_t start_index,
                                    bool & changed )
//----------------------------------------------------------------------
{
   if( submenu!=0 && submenu_nb!=name.size() )
   {
      glutDestroyMenu( submenu ) ;
      submenu = submenu_nb = 0 ;
   }

   if( submenu==0 && name.size()>0 )
   {
      changed = true ;
      submenu = glutCreateMenu( menuFunc ) ;
      for( size_t i=0 ; i<name.size() ; i++ )
      {
         if( !name(i).empty() )
            glutAddMenuEntry( name(i).c_str(), start_index+i ) ;
      }
      submenu_nb = name.size() ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: deal_with_selection_menu( void )
//----------------------------------------------------------------------
{
   bool changed = false ;

   manage_sub_menu( field_sub_menu, field_sub_menu_nb, fields_name,
                    FIELDS, changed ) ;
   manage_sub_menu( meshing_sub_menu, meshing_sub_menu_nb, meshing_name,
                    MESHINGS, changed ) ;

   display_colors.re_initialize( 0 ) ;
   for( size_t i=0 ; i<color_table.size() ; i++ )
      if( used(i) ) display_colors.append(color_table(i)) ;
   display_colors.sort() ;
   if( color_table.size() > 1 ) display_colors.append("all") ;
   manage_sub_menu( color_sub_menu, color_sub_menu_nb, display_colors,
                    COLORS, changed ) ;

   changed = changed ||
      ( surface_rendering!=surface_menu_status ) ;
   if( changed && selection_menu!=0 )
   {
      glutDestroyMenu( selection_menu ) ;
      selection_menu = 0 ;
   }
   if( view_menu!=0 && changed)
   {
      glutDestroyMenu( view_menu ) ;
      view_menu = 0 ;
   }

   if( view_menu==0 )
   {
      changed = true ;
      view_menu = glutCreateMenu( menuFunc ) ;
      surface_menu_status = surface_rendering ;

      if( dim_space==2 && surface_rendering )
      {
         glutAddMenuEntry( "No surface", SURFACE ) ;
      }
      else if( dim_space==2 && !surface_rendering )
      {
         glutAddMenuEntry( "Surface", SURFACE ) ;
      }
      else if( dim_space==3 && surface_rendering )
      {
         glutAddMenuEntry( "No Transparency", SURFACE ) ;
      }
      else if( dim_space==3 && !surface_rendering )
      {
         glutAddMenuEntry( "Transparency", SURFACE ) ;
      }
      glutAddMenuEntry( "Reset", RESET ) ;
      glutAddMenuEntry( "Fit to window", CROP ) ;
   }

   if( cycle_menu==0 )
   {
      changed = true ;
      cycle_menu = glutCreateMenu( menuFunc ) ;
      glutAddMenuEntry( "First", FIRST_CYCLE ) ;
      glutAddMenuEntry( "Previous", PREVIOUS_CYCLE ) ;
      glutAddMenuEntry( "Next", NEXT_CYCLE ) ;
      glutAddMenuEntry( "Last", LAST_CYCLE ) ;
   }

   if( selection_menu==0 )
   {
      selection_menu=glutCreateMenu( menuFunc ) ;
      glutAttachMenu(GLUT_RIGHT_BUTTON);
      glutAddMenuEntry( "Quit", QUIT ) ;
      glutAddMenuEntry( "Help", HELP ) ;
      if( cycle_menu!=0 ) glutAddSubMenu( "Go to cycle...",
                                         cycle_menu ) ;
      if( view_menu!=0 ) glutAddSubMenu( "View",
                                         view_menu ) ;
      if( field_sub_menu!=0 ) glutAddSubMenu( "Fields",
                                              field_sub_menu ) ;
      if( meshing_sub_menu!=0 ) glutAddSubMenu( "Meshing",
                                                meshing_sub_menu ) ;
      if( color_sub_menu!=0 ) glutAddSubMenu( "Colors",
                                              color_sub_menu ) ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLwindow:: display_help( std::ostream& os )
//----------------------------------------------------------------------
{
   os<<endl<<endl<<"PELICANS quick viewer help"<<endl ;
   os<<"Mouse motion with left button pushed : translate scene"<<endl ;
   os<<"Mouse motion with left button and shift key pushed : rotate scene"<<endl ;
   os<<"Mouse motion with middle button pushed : zoom scene"<<endl ;
   os<<"Mouse motion with middle button and shift key pushed : modify X/Y ratio"<<endl ;
   os<<"<Space> key : reset view"<<endl ;
   os<<"<q> key : exit"<<endl ;
}
