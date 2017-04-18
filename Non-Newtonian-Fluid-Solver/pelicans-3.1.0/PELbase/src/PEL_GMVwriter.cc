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

#include <PEL_GMVwriter.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Map.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::endl ; using std::cout ;

PEL_GMVwriter const* PEL_GMVwriter::PROTOTYPE = new PEL_GMVwriter() ;

//----------------------------------------------------------------------
PEL_GMVwriter:: PEL_GMVwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_GMVwriter" )
   , OPENMODE()
   , HAS_FACE_VARS( false )
   , DIM( PEL::bad_index() )
   , BINARY( false )
   , ICYCLE( 0 )
{
}

//----------------------------------------------------------------------
PEL_GMVwriter*
PEL_GMVwriter:: create_replica( PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_GMVwriter* result = new PEL_GMVwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_GMVwriter:: PEL_GMVwriter( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , OPENMODE()
   , HAS_FACE_VARS( false )
   , DIM( PEL::bad_index() )
   , BINARY( false )
   , ICYCLE( 0 )
{
   FILEBASENAME = exp->string_data( "files_basename" ) ;
   FILEXTENSION = ".gmv" ;
   OPENMODE =  std::ios::out ;

   FORMAT = "ascii" ;

   if( exp->string_data( "writing_mode" )=="binary" )
   {
      FORMAT = "ieeei4r4" ;
      BINARY = true ;
      OPENMODE = std::ios::out|std::ios::binary ;
      if( sizeof(int) != 4 )
      {
         std::ostringstream mesg ;
         mesg << "*** PEL_GMVwriter error:" << std::endl ;
         mesg << "    size_of(int) returns : " << sizeof(int) << std::endl ;
         mesg << "    and the size of the type \"int\" should be 4" ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_GMVwriter:: ~PEL_GMVwriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: determine_conditions( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: determine_conditions" ) ;

   if( DIM == PEL::bad_index() )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "meshing" ) ;
      DIM = se->int_data( "nb_sp_dims" ) ;
      se->destroy() ; se=0 ;
      if( exp->has_module( "fields" ) )
      {
         se = exp->create_subexplorer( 0, "fields" ) ;
         se->start_module_iterator() ;
         for( ; se->is_valid_module() ; se->go_next_module() )
         {
            PEL_ModuleExplorer const* sse = se->create_subexplorer( se ) ;
            std::string const& location = sse->string_data( "location" ) ;
            if( location == "at_face_centers" ) 
            {
               HAS_FACE_VARS = true ;
               break ;
            }
         }
         se->destroy() ; se=0 ;
      }
   }
   else
   {
      // *** peut etre à faire : vérifs de consistance.
   }
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   determine_conditions( exp ) ;

   // Current time
   ICYCLE = exp->int_data( "cycle_number" ) ;

   // Open file and test
   std::string fname = output_file_name( ICYCLE, "" ) ;
   FILE_FOR_PLAIN_VARS.open( fname.c_str(), OPENMODE ) ;
   if( FILE_FOR_PLAIN_VARS.fail() || !FILE_FOR_PLAIN_VARS.is_open() )
   {
      PEL_Error::object()->raise_plain(
	 "unable to create the GMV output file : " + fname ) ;
   }
   if( DIM==2 && HAS_FACE_VARS )
   {
      fname = output_file_name( ICYCLE, "_f" ) ;
      FILE_FOR_FACE2D_VARS.open( fname.c_str(), OPENMODE ) ;
      if( FILE_FOR_FACE2D_VARS.fail() || !FILE_FOR_FACE2D_VARS.is_open() )
      {
         PEL_Error::object()->raise_plain(
            "unable to create the GMV output file : " + fname ) ;
      }
   }

   write_data( "gmvinput ", 8, FILE_FOR_PLAIN_VARS ) ;
   write_data( FORMAT, 8, FILE_FOR_PLAIN_VARS ) ;
   write_data_endl( FILE_FOR_PLAIN_VARS ) ;
   write_data( "cycleno ", 8, FILE_FOR_PLAIN_VARS ) ;
   write_data( (int)ICYCLE, FILE_FOR_PLAIN_VARS ) ;
   write_data_endl( FILE_FOR_PLAIN_VARS ) ;

   if( DIM==2 && HAS_FACE_VARS )
   {
      write_data( "gmvinput ", 8, FILE_FOR_FACE2D_VARS ) ;
      write_data( FORMAT, 8, FILE_FOR_FACE2D_VARS ) ;
      write_data_endl( FILE_FOR_FACE2D_VARS ) ;
      write_data( "cycleno ", 8, FILE_FOR_FACE2D_VARS ) ;
      write_data( (int)ICYCLE, FILE_FOR_FACE2D_VARS ) ;
      write_data_endl( FILE_FOR_FACE2D_VARS ) ;
   }

   // Does a new grid must be saved?
   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      if( DIM==2 )
      {
         write_cells( se ) ;
         if( HAS_FACE_VARS ) write_faces( se ) ;
      }
      else
      {
         write_faces( se ) ;
      }
      se->destroy() ; se = 0 ;
   }

   if( DIM==2 )
   {
      write_data( "nodes ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( "fromfile ", 8, FILE_FOR_PLAIN_VARS ) ;
      std::string buf = "\""+PEL_System::basename( NAME_FILE_CELLS )+"\"" ;
      write_data( buf, buf.length(), FILE_FOR_PLAIN_VARS ) ;
      write_data_endl( FILE_FOR_PLAIN_VARS ) ;

      write_data( "cells ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( "fromfile ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( buf, buf.length(), FILE_FOR_PLAIN_VARS ) ;
      write_data_endl( FILE_FOR_PLAIN_VARS ) ;

      if( HAS_FACE_VARS )
      {
         write_data( "nodes ", 8, FILE_FOR_FACE2D_VARS ) ;
         write_data( "fromfile ", 8, FILE_FOR_FACE2D_VARS ) ;
         buf = "\""+PEL_System::basename( NAME_FILE_FACES )+"\"" ;
         write_data( buf, buf.length(), FILE_FOR_FACE2D_VARS ) ;
         write_data_endl( FILE_FOR_FACE2D_VARS ) ;

         write_data( "faces ", 8, FILE_FOR_FACE2D_VARS ) ;
         write_data( "fromfile ", 8, FILE_FOR_FACE2D_VARS ) ;
         write_data( buf, buf.length(), FILE_FOR_FACE2D_VARS ) ;
         write_data_endl( FILE_FOR_FACE2D_VARS ) ;
      }
   }
   else
   {
      write_data( "nodes ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( "fromfile ", 8, FILE_FOR_PLAIN_VARS ) ;
      std::string buf = "\""+PEL_System::basename( NAME_FILE_FACES )+"\"" ;
      write_data( buf, buf.length(), FILE_FOR_PLAIN_VARS ) ;
      write_data_endl( FILE_FOR_PLAIN_VARS ) ;

      write_data( "faces ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( "fromfile ", 8, FILE_FOR_PLAIN_VARS ) ;
      write_data( buf, buf.length(), FILE_FOR_PLAIN_VARS ) ;
      write_data_endl( FILE_FOR_PLAIN_VARS ) ;
   }

   // Saving of variables
   if( exp->has_module( "variables" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "variables" ) ;
      if( se->has_entry( "TIME" ) )
      {
         write_time_variable( se->double_data( "TIME" ), 
                              FILE_FOR_PLAIN_VARS ) ;
         if( DIM==2 && HAS_FACE_VARS )
            write_time_variable( se->double_data( "TIME" ), 
                                 FILE_FOR_FACE2D_VARS ) ;
      }
      se->destroy() ; se = 0 ;
   }

   // Saving of fields
   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
      write_fields( se ) ;
      se->destroy() ; se = 0 ;
   }

   write_data( "endgmv ", 8, FILE_FOR_PLAIN_VARS ) ;
   write_data_endl( FILE_FOR_PLAIN_VARS ) ;
   FILE_FOR_PLAIN_VARS.close() ;
   if( DIM==2 && HAS_FACE_VARS )
   {
      write_data( "endgmv ", 8, FILE_FOR_FACE2D_VARS ) ;
      write_data_endl( FILE_FOR_FACE2D_VARS ) ;
      FILE_FOR_FACE2D_VARS.close() ;
   }
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_cells( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL("PEL_GMVwriter:: write_cells") ;
   PEL_CHECK( exp!=0 ) ;

   // Open grid file and test
   NAME_FILE_CELLS = output_file_name( ICYCLE, "_grid" ) ;
   std::ofstream file( NAME_FILE_CELLS.c_str(), OPENMODE ) ;
   if( file.fail() || !file.is_open() )
   {
      PEL_Error::object()->raise_plain(
	 "write_cells : unable to create the GMV output file : " + NAME_FILE_CELLS ) ;
   }
   write_data( "gmvinput ", 8, file ) ;
   write_data( FORMAT, 8, file ) ;
   write_data_endl( file ) ;


   // Vertices
   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   write_data( "nodev ",8, file ) ;
   write_data( (int)vertices.index_bound( 1 ), file ) ;
   write_data_endl( file ) ;
   save_vertices( vertices, file ) ;

   // Cells
   intVector const& cell_nb_vertices = exp->intVector_data( "cell_nb_vertices" ) ;
   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   write_data( "cells ", 8, file ) ;
   write_data( (int)connectivity.index_bound( 1 ), file ) ;
   write_data_endl( file ) ;
 
   save_cells( cell_nb_vertices, connectivity, file ) ;

   write_data( "endgmv ", 8, file ) ;
   write_data_endl( file ) ;
   file.close() ;
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_faces( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL("PEL_GMVwriter:: write_faces") ;
   PEL_CHECK( exp!=0 ) ;

   // Open grid file and test
   std::string add_string = "_grid" ;
   if( DIM==2 ) 
   {
      PEL_ASSERT( HAS_FACE_VARS ) ;
      add_string += "_f" ;
   }
   NAME_FILE_FACES = output_file_name( ICYCLE, add_string ) ;
   std::ofstream file( NAME_FILE_FACES.c_str(), OPENMODE ) ;
   if( file.fail() || !file.is_open() )
   {
      PEL_Error::object()->raise_plain(
	 "write_faces : unable to create the GMV output file : " + NAME_FILE_FACES ) ;
   }
   write_data( "gmvinput ", 8, file ) ;
   write_data( FORMAT, 8, file ) ;
   write_data_endl( file ) ;

   // Vertices
   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   write_data( "nodev ",8, file ) ;
   write_data( (int)vertices.index_bound( 1 ), file ) ;
   write_data_endl( file ) ;
   save_vertices( vertices, file ) ;

   // Faces
   intArray2D const& cell2vertex = exp->intArray2D_data( "cell2vertex" ) ;
   intVector const& face_nb_verts = exp->intVector_data( "face_nb_vertices" ) ;
   intArray2D const& face2vertex = exp->intArray2D_data( "face2vertex" ) ;
   intArray2D const& face2cell = exp->intArray2D_data( "face2cell" ) ;
   write_data( "faces ", 8, file ) ;
   size_t nb_faces  = face2vertex.index_bound( 1 ) ;
   size_t nb_cells  = cell2vertex.index_bound( 1 ) ;
   write_data( (int)nb_faces, file ) ;
   write_data( (int)nb_cells, file ) ;
   write_data_endl( file ) ;
 
   for(  size_t m=0 ; m<nb_faces ; ++m )
   {
      write_data( face_nb_verts( m ), file ) ;
      for( size_t iv=0 ; iv<(size_t)face_nb_verts( m ) ; ++iv )
      {
         write_data( face2vertex( iv, m ) + 1, file ) ;
      }
      write_data( face2cell( 0, m ) + 1, file ) ;
      int icell = face2cell( 1, m ) ;
      if( icell == index_for_trash() )
      {
         icell = 0 ;
      }
      else
      {
         icell = icell+1 ;
      }
      write_data( icell, file ) ;
      write_data_endl( file ) ;
   }

   write_data( "endgmv ", 8, file ) ;
   write_data_endl( file ) ;
   file.close() ;
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_fields( PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: write_fields" ) ;
   PEL_CHECK( exp!=0 ) ;

   // Fields that are not the velocity field
   write_data( "variable ", 8, FILE_FOR_PLAIN_VARS ) ;
   write_data_endl( FILE_FOR_PLAIN_VARS ) ;
   if( DIM==2 && HAS_FACE_VARS )
   {
      write_data( "variable ", 8, FILE_FOR_FACE2D_VARS ) ;
      write_data_endl( FILE_FOR_FACE2D_VARS ) ;
   }
   exp->start_module_iterator() ;
   PEL_ModuleExplorer const* vexp = 0 ;
   for( ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer const* fexp = exp->create_subexplorer( 0 ) ;
      std::string const& name = fexp->string_data( "name" ) ;
      if( name!="VELO" )
      {
         write_one_field( name, fexp ) ;
      }
      else
      {
         vexp = exp->create_subexplorer( 0 ) ;
      }
      fexp->destroy() ; fexp = 0 ;
   }
   write_data( "endvars ", 8, FILE_FOR_PLAIN_VARS ) ;
   write_data_endl( FILE_FOR_PLAIN_VARS ) ;
   if( DIM==2 && HAS_FACE_VARS )
   {
      write_data( "endvars ", 8, FILE_FOR_FACE2D_VARS ) ;
      write_data_endl( FILE_FOR_FACE2D_VARS ) ;
   }

   // Velocity field
   if( vexp != 0 )
   {
      PEL_ASSERT( vexp->string_data( "name" )=="VELO" ) ;
      write_velocity_field( vexp ) ;
      vexp->destroy() ; vexp = 0 ;
   }
}

//-----------------------------------------------------------------------
void 
PEL_GMVwriter:: write_time_variable( double t, std::ofstream& file )
//-----------------------------------------------------------------------
{
   static const size_t sd = sizeof(float) ;
   static const size_t ss = sizeof(char) ;
   if( !BINARY )
   {
      file << "probtime " << t << std::endl ;
   }
   else
   {
      std::string tmp( "probtime" ) ;
      file.write(tmp.c_str(), 8*ss ) ;
      float fval = (float)t ;
      file.write((const char*)&fval,sd) ;
   }
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_velocity_field( PEL_ModuleExplorer const* fexp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: write_velocity_field" ) ;
   PEL_CHECK( fexp!=0 ) ;

   doubleArray2D const& X = fexp->doubleArray2D_data( "value" ) ;
   std::string const& location = fexp->string_data( "location" ) ;
   int pos = 0 ;
   if( location == "at_vertices" ) 
   {
      pos = 1 ;
   }
   else if( location == "at_face_centers" )
   {
      PEL_ASSERT( HAS_FACE_VARS ) ;
      pos = 2 ;
   }
   else if( location != "at_cell_centers" )
   {
      raise_field_location_error(
         "velocity", location,
         "at_vertices,at_face_centers,at_cell_centers" ) ;
   }

   std::ofstream& file = ( pos==2 && DIM==2 && HAS_FACE_VARS ) ? 
                            FILE_FOR_FACE2D_VARS : FILE_FOR_PLAIN_VARS ;

   write_data( "velocity ", 8, file ) ;
   write_data( pos, file ) ;
   write_data_endl( file ) ;

   // Field writing. If vectorial field create one new field by component 
   // adding "_" and the component number to the original name.
   size_t nb_components = X.index_bound( 0 ) ;
   size_t nb_values_by_component = X.index_bound(1) ;
   for( size_t j=0 ; j<nb_components ; j++ )
   {
      for( size_t i=0 ; i<nb_values_by_component ; i++ ) 
      {
         write_data( X(j,i), file ) ;
      }
      write_data_endl( file ) ;
   }
   for( size_t j=nb_components ; j<3 ; j++ )
   {
      for( size_t i=0 ; i<nb_values_by_component ; i++ ) 
      {
         write_data( 0., file ) ;
      }
      write_data_endl( file ) ;
   }
   write_data_endl( file ) ;
}

//----------------------------------------------------------------------
void
PEL_GMVwriter:: write_one_field( std::string const& name,
                                 PEL_ModuleExplorer const* fexp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: write_one_field" ) ;
   PEL_CHECK( fexp!=0 ) ;

   doubleArray2D const& X = fexp->doubleArray2D_data( "value" ) ;
   std::string const& location = fexp->string_data( "location" ) ;
   int pos = 0 ;
   if( location == "at_vertices" ) 
   {
      pos = 1 ;
   }
   else if( location == "at_face_centers" )
   {
      PEL_ASSERT( HAS_FACE_VARS ) ;
      pos = 2 ;
   }
   else if( location != "at_cell_centers" )
   {
      raise_field_location_error(
         name, location, "at_vertices,at_face_centers,at_cell_centers" ) ;
   }

   std::ofstream& file = ( pos==2 && DIM==2 && HAS_FACE_VARS ) ? 
                            FILE_FOR_FACE2D_VARS : FILE_FOR_PLAIN_VARS ;

   // Field writing. If vectorial field create one new field by component 
   // adding "_" and the component number to the original name.
   bool is_vectorial = ( X.index_bound(0)>1 ) ;
   for( size_t j=0 ; j<X.index_bound(0) ; j++ )
   {
      std::string rname = name ;
      std::ostringstream tmp ;
      if( is_vectorial )
      {
         tmp << j ;
         rname += ("_"+tmp.str()) ; 
      }
      write_data( rname, 8, file ) ;
      write_data( pos, file ) ;
      write_data_endl( file ) ;
      for( size_t i=0 ; i<X.index_bound(1) ; i++ ) 
      {
         write_data( X(j,i), file ) ;
      }
      write_data_endl( file ) ;
   }
   write_data_endl( file ) ;
}

//-----------------------------------------------------------------------
std::string
PEL_GMVwriter:: output_file_name( size_t nb, std::string add_string )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: output_file_name" ) ;
   PEL_CHECK( nb<99999 ) ;

   std::ostringstream tmp ;
   tmp << nb ;
   std::string nb_string = tmp.str() ;
   std::string result = FILEBASENAME+add_string+FILEXTENSION+".00000";
   result.replace( result.length()-nb_string.length(), 
                   nb_string.length(),
                   nb_string ) ;

   PEL_CHECK( result.length()==( FILEBASENAME.length()+add_string.length()
                                +FILEXTENSION.length()+6) ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: save_vertices( doubleArray2D const& vert,
                               std::ofstream& file )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: save_vertices" ) ;
   size_t dim     = vert.index_bound(0);
   size_t nb_vert = vert.index_bound(1);

   for( size_t i=0 ; i<nb_vert ; i++ ) 
   {
      for( size_t j=0 ; j<dim ; j++ ) 
      {
         write_data( vert(j,i), file ) ;
      }
      for( size_t j=dim ; j<(size_t)3 ; j++ ) 
      {
         write_data( 0., file ) ;
      }     
      write_data_endl( file ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: save_cells( intVector const& cell_vertices,
                            intArray2D const& connec,
                            std::ofstream& file )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: save_cells" ) ;
   
   size_t nb_cell = connec.index_bound(1) ;

   std::string elm_type = "" ;
   for( size_t i=0 ; i<nb_cell ; i++ ) 
   {
      size_t const vert_by_mesh = cell_vertices( i ) ;
      if( DIM==1 && vert_by_mesh==2 )
      {
         elm_type = "line " ;
      }
      else if( DIM==2 && vert_by_mesh==3 )
      {
         elm_type = "tri " ;
      }
      else if( DIM==2 && vert_by_mesh==4 )
      {
         elm_type = "quad " ;
      }
      else if( DIM==3 && vert_by_mesh==4 )
      {
         elm_type = "tet " ;
      }
      else if( DIM==3 && vert_by_mesh==8 )
      {
         // Instead of standard hex we use the Patran numbered hexahedron 
         elm_type = "phex8 " ;
      }
      else 
      {
         PEL_Error::object()->raise_plain( "Invalid polyhedron type for GMV output " ) ;
      }

      write_data( elm_type, 8, file ) ;
      write_data( (int) vert_by_mesh, file ) ;
      write_data_endl( file ) ;
      for( size_t j=0 ; j<vert_by_mesh ; j++ ) 
      {
         int nb = connec(j,i)+1 ;
         write_data( nb, file ) ;
      }
      write_data_endl( file ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: write_data( std::string val, size_t length, 
                            std::ofstream& file ) 
//-----------------------------------------------------------------------
{  
   static const size_t ss = sizeof(char) ;
   if( !BINARY )
   {
      file << val ;
   }
   else
   {
      std::string tmp( length, ' ' ) ;
      tmp.replace( 0, val.length(), val ) ;
      file.write( tmp.c_str(), length*ss ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: write_data( int val, std::ofstream& file ) 
//-----------------------------------------------------------------------
{
   static const size_t si = sizeof(int) ;
   if( !BINARY )
   {
      file << " " << val ;
   }
   else
   {
      file.write((const char*)&val,si) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: write_data( double val, std::ofstream& file ) 
//-----------------------------------------------------------------------
{
   static const size_t sd = sizeof(float) ;
   if( !BINARY )
   {
      file << " " << val ;
   }
   else
   {
      float fval = (float)val ;
      file.write((const char*)&fval,sd) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_GMVwriter:: write_data_endl( std::ofstream& file )
//-----------------------------------------------------------------------
{
   if( !BINARY )
   {
      file << std::endl ;
   }
}
