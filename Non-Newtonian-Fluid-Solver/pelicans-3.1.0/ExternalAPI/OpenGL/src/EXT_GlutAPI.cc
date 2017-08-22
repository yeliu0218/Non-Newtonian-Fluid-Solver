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

#include <EXT_GlutAPI.hh>

#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>

#ifdef __APPLE__
   #include <GLUT/glut.h>  // for Max OS X
#else
   #include <GL/glut.h>    // for other systems
#endif


//----------------------------------------------------------------------
EXT_GlutAPI*
EXT_GlutAPI:: SINGLETON = new EXT_GlutAPI() ;
char** EXT_GlutAPI:: my_argv ;
int EXT_GlutAPI:: my_argc ;
bool EXT_GlutAPI::DISPLAY_DEFINED = false ;
//----------------------------------------------------------------------


//----------------------------------------------------------------------
EXT_GlutAPI:: EXT_GlutAPI( void )
//----------------------------------------------------------------------
   : PEL_ExternalAPI( "EXT_GlutAPI", 0 )
{
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_OpenGL" ), val ) ;
}


//----------------------------------------------------------------------
EXT_GlutAPI:: ~EXT_GlutAPI( void )
//----------------------------------------------------------------------
{
   delete [] my_argv ;
}

//----------------------------------------------------------------------
void
EXT_GlutAPI:: initialize( int& argc, char **& argv )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_GlutAPI:: mainAPI" ) ;
   PEL_CHECK( argc>0 ) ;

   my_argv = new char * [ argc ] ;
   my_argc=0 ;
   if( argc> 0 )
   {
      my_argv[my_argc++] = argv[0] ;
   }

   char** new_argv = new char * [ argc ] ;
   int new_argc = 0 ;
   DISPLAY_DEFINED = ! PEL_System::getenv( "DISPLAY" ).empty() ;

   for( size_t i=0 ; i<argc ; i++ )
   {
      std::string str = argv[i] ;
      if( str.find( "-display" )==0 )
      {
         DISPLAY_DEFINED = true ;
      }
      if( ( str.find( "-display" )==0 ) ||
          ( str.find( "-geometry" )==0 ) )
      {
         my_argv[my_argc++] = argv[i++] ;
         my_argv[my_argc++] = argv[i] ;
      }
      else if( ( str.find( "-iconic" )==0 ) ||
               ( str.find( "-indirect" )==0 ) ||
               ( str.find( "-direct" )==0 ) ||
               ( str.find( "-gldebug" )==0 ) ||
               ( str.find( "-sync" )==0 ) )
      {
         my_argv[my_argc++] = argv[i] ;
      }
      else
      {
         new_argv[new_argc++] = argv[i] ;
      }
   }
   if( my_argc>0 )
   {
      argc = new_argc ;
      argv = new_argv ;
   }
   else
   {
      delete [] new_argv ;
   }
}
#include <iostream>
void  glVerify(std::string const& what) {
  int glErr=glGetError();
//   if( glErr!=GL_NO_ERROR )
    std::cout<< std::endl << "ERROR : " << what
	     << " with error code " <<glErr<<std::endl ;
}

//----------------------------------------------------------------------
void
EXT_GlutAPI:: initialize_on_demand( void )
//----------------------------------------------------------------------
{
   static bool prem = true ;
   if( prem )
   {
      if( !DISPLAY_DEFINED )
      {
         PEL_Error::object()->raise_plain(
            "You must specify a valid DISPLAY environment variable to use OpenGL tools") ;
      }
      glutInit( &my_argc, my_argv ) ;
      glutInitWindowPosition( 0, 0 ) ;
      glutInitWindowSize( 400, 400 );
      //glutInitDisplayString("rgba double depth=24");
      glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
      prem = false ;
   }
}

//----------------------------------------------------------------------
void
EXT_GlutAPI:: disable_window_destroy( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_GlutAPI:: disable_window_destroy" ) ;
   FORMAL( "Former implementation no portable" ) ;
}
