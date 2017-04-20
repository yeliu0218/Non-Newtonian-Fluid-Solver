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

#include <EXT_OpenGLviewer.hh>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <EXT_OpenGLwindow.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_assertions.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>

using std::endl ;
using std::string ; using std::ostringstream ;

EXT_OpenGLviewer const* EXT_OpenGLviewer::PROTOTYPE = new EXT_OpenGLviewer() ;

//----------------------------------------------------------------------
EXT_OpenGLviewer:: EXT_OpenGLviewer( void )
//----------------------------------------------------------------------
   : PEL_Application( "viewer" )
   , ogl( 0 )
   , vertices( 0, 0 )
   , connectivity( 0, 0 )
{
}

//----------------------------------------------------------------------
EXT_OpenGLviewer*
EXT_OpenGLviewer:: create_replica( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: create_replica from exp" ) ;
   return 0 ;
}

//----------------------------------------------------------------------
EXT_OpenGLviewer*
EXT_OpenGLviewer:: create_replica_from_args( PEL_Object* a_owner,
                                             stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: create_replica" ) ;

   EXT_OpenGLviewer* result = new EXT_OpenGLviewer( a_owner, args ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
EXT_OpenGLviewer:: EXT_OpenGLviewer( PEL_Object* a_owner,
                                     stringVector& args )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , ogl( 0 )
   , vertices(0,0)
   , connectivity( 0, 0)
   , first( true )
{
   if(args.size()!=1) notify_error_in_arguments() ;

   filename = args(0) ;
   args.remove_at(0) ;

   ogl = EXT_OpenGLwindow::create( this, filename ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
EXT_OpenGLviewer:: ~EXT_OpenGLviewer( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void EXT_OpenGLviewer:: run( void )
//----------------------------------------------------------------------
{
   ogl->do_user_interface_loop( this ) ;
}

//----------------------------------------------------------------------
void EXT_OpenGLviewer:: go_to(
			EXT_OpenGLwindow::CycleMove relative_cycle_number )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: go_to" ) ;

   size_t requested_cycle = last_requested_cycle ;
   if( relative_cycle_number==EXT_OpenGLwindow::PREVIOUS_CYCLE &&
       last_displayed_cycle>1 )
      requested_cycle = last_displayed_cycle-1 ;
   else if( relative_cycle_number==EXT_OpenGLwindow::FIRST_CYCLE )
      requested_cycle = 1 ;
   else if( relative_cycle_number==EXT_OpenGLwindow::NEXT_CYCLE &&
            last_requested_cycle!=PEL::bad_index() )
      requested_cycle = last_displayed_cycle+1 ;
   else if( relative_cycle_number==EXT_OpenGLwindow::LAST_CYCLE )
      requested_cycle = PEL::bad_index() ;

//    std::cout << "EXT_OpenGLviewer : try to read : " ;
//    if( requested_cycle==PEL::bad_index() )
//       std::cout << " last cycle in file " ;
//    else
//       std::cout << " cycle number " << requested_cycle ;
//    std::cout << std::endl ;

   struct stat status ;
   if( stat(filename.c_str(),&status)!=0 )
   {
      PEL_Error::object()->raise_plain(
         "Unable to inquire for file "+filename ) ;
   }
   off_t last = status.st_size ;
   if( first ||
       ( last!=last_modif ) ||
       ( requested_cycle!=last_requested_cycle ) )
   {
      PEL_ModuleExplorer const* exp = read_cycle( 0, requested_cycle ) ;
      if( exp!=0 )
      {
         write_cycle(exp) ;
         exp->destroy() ;
         last_modif = last ;
         last_requested_cycle = requested_cycle ;
      }
   }
}

//----------------------------------------------------------------------
PEL_ModuleExplorer const*
EXT_OpenGLviewer:: read_cycle( PEL_Object * a_owner, size_t cycle_number )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: create_cycle" ) ;
   PEL_ModuleExplorer * result = 0 ;

   std::ifstream stream(filename.c_str()) ;
   if( !stream )
   {
      PEL_Error::object()->raise_plain(
         "Unable to open file "+filename ) ;
   }
   if( !first &&
       ( cycle_number>last_displayed_cycle || cycle_number==PEL::bad_index() ) )
   {
      stream.seekg( last_pos ) ;
   }

   std::string buff ;
   bool started = false ;
   std::string line ;
   std::ostringstream cycle_name ;
   cycle_name << "MODULE cycle_" ;
   if( cycle_number!=PEL::bad_index() ) cycle_name << cycle_number ;

   while( getline( stream, line ) && !stream.eof() )
   {
      bool has_cycle = line.find( cycle_name.str() ) < line.length() ;
      //if(has_cycle)
//      std::cout<<"line : "<<line<<std::endl ;
      if( has_cycle )
      {
         if( !started )
         {
            started = true ;
            buff = "" ;
         }
         else
         {
            started = false ;
            last_pos = stream.tellg() ;
            buff += line + "\n" ;
            if( cycle_number!=PEL::bad_index() ) break ;
         }
      }
      if( started ) buff += line + "\n" ;
   }
   if( !started && !buff.empty() )
   {
      first = false ;

      std::istringstream input_buff(buff) ;
      PEL_Module* mod = PEL_Module::create( 0, "save", input_buff ) ;
      PEL_ModuleExplorer * exp = PEL_ModuleExplorer::create( 0, mod ) ;
      exp->start_module_iterator() ;
      PEL_ASSERT( exp->is_valid_module() ) ;
      result = exp->create_subexplorer( a_owner ) ;
      exp->set_owner( result ) ;
      mod->set_owner( result ) ;

      last_displayed_cycle = result->int_data( "cycle_number" ) ;

   }
   PEL_CHECK_POST( result==0 || result->owner()==a_owner ) ;
   stream.close() ;

   return result ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLviewer:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: write_cycle" ) ;
   PEL_CHECK_PRE( exp!=0 ) ;
   PEL_CHECK_PRE( exp->has_entry( "cycle_number" ) ) ;

   int cycle_number = exp->int_data( "cycle_number" ) ;
   std::ostringstream title ;
   title << "cycle " << cycle_number ;
   
   if( exp->has_entry( "variables/TIME" ) )
   {
      title << " (t=" << exp->double_data( "variables/TIME" ) << ")" ;      
   }
   ogl->set_title( title.str() ) ;

   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      write_grid( se ) ;
      se->destroy() ;
   }
   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         PEL_ModuleExplorer* sse = se->create_subexplorer( 0 ) ;
         write_field( cycle_number, sse ) ;
         sse->destroy() ;
      }
      se->destroy() ;
   }
   if( exp->has_module( "integration_domain" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "integration_domain" ) ;
      write_integration_domain( se ) ;
      se->destroy() ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLviewer:: write_grid( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: write_grid" ) ;

   vertices = exp->doubleArray2D_data( "vertices" ) ;
   connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   intArray2D const& faces = exp->intArray2D_data( "face2vertex" ) ;
   size_t nbvert = vertices.index_bound(1) ;
   size_t nbmesh = connectivity.index_bound(1) ;
   size_t nb_dim_space = vertices.index_bound(0) ;
   if( exp->has_entry( "color_table" ) )
   {
      stringVector const& vec = exp->stringVector_data("color_table") ;

      intVector const& vertex_color = exp->intVector_data("vertex_color") ;
      intVector const& mesh_color = exp->intVector_data("cell_color") ;
      intVector const& face_color = exp->intVector_data("face_color") ;
      intArray2D const& macro_colors = exp->intArray2D_data("color_table_connectivity") ;
      ogl->set_color_table( vec, macro_colors ) ;

      ogl->display_vertices( "vertices",
                             vertices,
                             vertex_color ) ;
      ogl->display_meshing( "cells",
                            vertices,
                            connectivity,
                            mesh_color,
                            nb_dim_space ) ;
      ogl->display_meshing( "faces",
                            vertices,
                            faces,
                            face_color,
                            nb_dim_space-1 ) ;
   }
   else
   {
      intVector vertex_color(nbvert) ;
      vertex_color.set(0) ;
      intVector mesh_color(nbmesh) ;
      mesh_color.set(0) ;
      intVector face_color(faces.index_bound(1)) ;
      face_color.set(0) ;
      stringVector vec(1) ;
      vec(0) = "all" ;
      intArray2D macros_table(1,1) ;
      macros_table.set(1) ;

      ogl->set_color_table( vec, macros_table ) ;

      ogl->display_vertices( "vertices",
                             vertices,
                             vertex_color ) ;
      ogl->display_meshing( "interior",
                            vertices,
                            connectivity,
                            mesh_color,
                            nb_dim_space ) ;
      ogl->display_meshing( "boundary",
                            vertices,
                            faces,
                            face_color,
                            nb_dim_space-1 ) ;
   }
}

//----------------------------------------------------------------------
void
EXT_OpenGLviewer:: write_integration_domain( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: write_integration_domain" ) ;

   doubleArray2D const& verticesFS =
                           exp->doubleArray2D_data( "inner_boundary" ) ;
   size_t const n = verticesFS.index_bound(1) ;
   intArray2D connectivityFS( 2, n-1 ) ;
   for( size_t i=0 ; i<n-1 ; i++ )
   {
      connectivityFS(0,i)=i ;
      connectivityFS(1,i)=i+1 ;
   }
   intVector mesh_colorFS(n-1) ;
   mesh_colorFS.set(0) ;
   ogl->display_meshing( "free surface",
                         verticesFS,
                         connectivityFS,
                         mesh_colorFS,
                         1 ) ;
}

//----------------------------------------------------------------------
void
EXT_OpenGLviewer:: write_field( int cycle_number,
                                PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_OpenGLviewer:: write_field" ) ;

   doubleArray2D const& X = exp->doubleArray2D_data( "value" ) ;
   std::string const& name = exp->string_data( "name" ) ;
   std::string const& location = exp->string_data( "location" ) ;

   if( location != "at_cell_centers" && location != "at_vertices" )
   {
      ostringstream mesg ;
      mesg << "*** PELICANS viewer :" << endl ;
      mesg << "    unable to save fields with" << endl ;
      mesg << "    location = \"" << location << "\"" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   ogl->display_field( name, location, vertices, connectivity, X ) ;
}

//---------------------------------------------------------------------------
void
EXT_OpenGLviewer:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << usage_title( "viewer" )  ;
   PEL::out() << "filename" << endl << endl ;
   PEL::out() << "     Display in real-time results in PELICANS format" << endl ;
}

//---------------------------------------------------------------------------
void
EXT_OpenGLviewer:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << operands_title() ;
   PEL::out() << "     file" << endl
              << "          file name of the file to be displayed"
              << endl << endl ;
}

//---------------------------------------------------------------------------
void
EXT_OpenGLviewer:: print_exit_status( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << exit_status_title() ;
   PEL::out() << "     0    The file exists" << endl ;
   PEL::out() << "    >0    The file can't be read" << endl ;
}

