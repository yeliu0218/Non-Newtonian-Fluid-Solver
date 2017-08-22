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

#include <PEL_FieldViewWriter.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <iostream>
#include <fstream>
#include <sstream>

// *******************************************************************
// File: fv_reader_tags.h
//
/* Numeric tags (codes) for FIELDVIEW binary file format. */

#define FV_MAGIC	0x00010203	/* decimal 66051 */

/* Content of the file (grid only, results only or combined). */
#define FV_GRIDS_FILE           1
#define FV_RESULTS_FILE         2
#define FV_COMBINED_FILE        3

#define FV_NODES        	1001
#define FV_FACES        	1002
#define FV_ELEMENTS     	1003
#define FV_VARIABLES    	1004
#define FV_BNDRY_VARS   	1006
#define FV_ARB_POLY_FACES       1007
#define FV_ARB_POLY_ELEMENTS    1008
#define FV_ARB_POLY_BNDRY_VARS  1009

#define FV_TET_ELEM_ID          1
#define FV_HEX_ELEM_ID          2
#define FV_PRISM_ELEM_ID        3
#define FV_PYRA_ELEM_ID         4
#define FV_ARB_POLY_ELEM_ID     5

#define MAX_NUM_ELEM_FACES     6
#define BITS_PER_WALL  3
#define ELEM_TYPE_BIT_SHIFT    (MAX_NUM_ELEM_FACES*BITS_PER_WALL)

/* Values for "wall_info" array (see comments in fv_encode_elem_header). */
// #ifdef __STDC__
// #define A_WALL         (07u)
// #define NOT_A_WALL     (0u)
// #else
// #define A_WALL         (07)
// #define NOT_A_WALL     (0)
// #endif

// *******************************************************************

std::string const FILE_EXTENSION = ".uns" ;

// Dummy vectors:
intVector IBUFF2(2) ;
intVector IBUFF3(3) ;
intVector IBUFF4(4) ;
intVector IBUFF5(5) ;
doubleVector DBUFF4(4) ;

// Connectivity PEL -> FieldView
int const PEL2FV_HEXA[8]  = { 0, 1, 3, 2, 4, 5, 7, 6 } ;
int const PEL2FV_TETRA[4] = { 0, 1, 2, 3 } ;
int const PEL2FV_QUAD[4]  = { 0, 1, 2, 3 } ;
int const PEL2FV_TRIA[3]  = { 0, 1, 2 } ;

// FV_ELEMENTS:
int const TETRA_IDX = 1 ;
int const HEXA_IDX  = 2 ;
int const PRISM_IDX = 3 ;
int const PYRA_IDX  = 4 ;
   
PEL_FieldViewWriter const* 
PEL_FieldViewWriter::PROTOTYPE = new PEL_FieldViewWriter() ;

//----------------------------------------------------------------------
PEL_FieldViewWriter:: PEL_FieldViewWriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_FieldViewWriter" )
   , FILE_BASENAME()
   , BINARY( false )
   , DIM( PEL::bad_index() )
   , ICYCLE( 0 )
{
}

//----------------------------------------------------------------------
PEL_FieldViewWriter*
PEL_FieldViewWriter:: create_replica( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_FieldViewWriter* result = new PEL_FieldViewWriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_FieldViewWriter:: PEL_FieldViewWriter( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , FILE_BASENAME( exp->string_data( "files_basename" ) )
   , BINARY( exp->string_data( "writing_mode" )=="binary" )
   , DIM( PEL::bad_index() )
   , ICYCLE( 0 )
{
   PEL_CHECK_INV( invariant() ) ;

   if( BINARY )
   {
      if( sizeof(int) != 4 )
      {
         std::ostringstream mesg ;
         mesg << "*** PEL_FieldViewWriter error:" << std::endl ;
         mesg << "    size_of(int) returns : " << sizeof(int) << std::endl ;
         mesg << "    and the size of the type \"int\" should be 4" ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }
      if( sizeof(float) != 4 )
      {
         std::ostringstream mesg ;
         mesg << "*** PEL_FieldViewWriter error:" << std::endl ;
         mesg << "    size_of(float) returns : " << sizeof(float) << std::endl ;
         mesg << "    and the size of the type \"float\" should be 4" ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }
   }
}


//----------------------------------------------------------------------
PEL_FieldViewWriter:: ~PEL_FieldViewWriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;
   
   // Current time
   ICYCLE = exp->int_data( "cycle_number" ) ;

   // Grid file:
   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      if( DIM == PEL::bad_index() )
      {
         DIM = (size_t) se->int_data( "nb_sp_dims" ) ;
         if( DIM != 2 && DIM != 3 )
         {
            PEL_Error::object()->raise_plain(
               "*** PEL_FieldViewWriter error:\n"
               "    2D or 3D meshing is expected" ) ;
         }
      }
      write_meshing( se ) ;
      
      se->destroy() ; se = 0 ;
   }

   if( DIM == PEL::bad_index() )
   {
      PEL_Error::object()->raise_plain(
         "*** PEL_FieldViewWriter error:\n"
         "    meshing saving if expected at the first saving cycle" ) ;
   }

   // Result file:
   std::string const& file_name = output_file_name( ICYCLE, "" ) ;
   std::ios_base::openmode file_mode =
      ( BINARY ? std::ios::out|std::ios::binary : std::ios::out ) ;
   std::ofstream file( file_name.c_str(), file_mode ) ; 
   if( file.fail() || !file.is_open() )
   {
      PEL_Error::object()->raise_plain(
	 "*** PEL_FieldViewWriter error:\n"
         "    unable to create the FieldView output file: "+file_name ) ;
   }

   // Header:
   write_result_header( file ) ;

   // Contants:
   write_constants( exp, file ) ;

   // Fields:
   write_fields( exp, file ) ;

   file.close() ;
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_meshing( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter::  write_meshing" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( DIM !=  PEL::bad_index()) ;

   std::string const& file_name = output_file_name( ICYCLE, "_grid" ) ;
   std::ios_base::openmode file_mode =
      ( BINARY ? std::ios::out|std::ios::binary : std::ios::out ) ;
   std::ofstream file( file_name.c_str(), file_mode ) ; 
   if( file.fail() || !file.is_open() )
   {
      PEL_Error::object()->raise_plain(
	 "*** PEL_FieldViewWriter error:\n"
         "    unable to create the FieldView output file: "+file_name ) ;
   }

   // Header:
   write_meshing_header( file ) ;

   // Boundary colors:
   size_t nb_bounds = 0 ;
   size_t_vector db_color_idx(0) ;
   write_bound_colors( exp, nb_bounds, db_color_idx, file ) ;

   // Vertices:
   write_vertices( exp, file ) ;

   // Bounds:
   write_boundaries( exp, nb_bounds, db_color_idx, file ) ;

   // Cells:
   write_cells( exp, file ) ;
   
   file.close() ;
}

//----------------------------------------------------------------------
std::string const&
PEL_FieldViewWriter:: output_file_name( size_t cycle, 
                                        std::string const& add_string ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: output_file_name" ) ;

   static std::string result ;

   std::ostringstream tmp ;
   tmp << cycle  ;
   std::string nb_string = tmp.str() ;
   
   result = FILE_BASENAME+add_string+"_00000" ;
   result.replace( result.length()-nb_string.length(), 
                   nb_string.length(),
                   nb_string ) ;
   result += FILE_EXTENSION ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_meshing_header( std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_meshing_header" ) ;
   PEL_CHECK( file ) ;

   /* Output the magic number. */
   if( BINARY ) write_data( (int) FV_MAGIC, file ) ;
   
   /* Output file header and version number. */
   std::string const head = ( BINARY ? "FIELDVIEW" : "FIELDVIEW_Grids" ) ;
   write_str_data( head, file ) ;

   /* This version of the FIELDVIEW unstructured file is "3.0" */
   IBUFF2(0)=3 ; IBUFF2(1)=0 ;
   write_data( IBUFF2, file ) ;
   write_data_endl( file ) ;
   
   /* Comments */
   write_data_endl( file ) ;
   write_comment( "File generated by PELICANS for FieldView post-processing",
                  file ) ;
   write_comment( "Grid file", file ) ;
   write_data_endl( file ) ;
   
   
   // Nb grids
   int const nb_grids = 1 ;
   if( BINARY )
   {
      /* File type code - new in version 2.7 */
      write_data( FV_GRIDS_FILE, file ) ;
      
      /* Reserved field, always write a zero - new in version 2.6 */
      write_data( (int) 0, file ) ;
      
      /* Number of grids */
      write_data( nb_grids, file ) ;
   }
   else
   {
      write_str_data( "Grids", file ) ;
      write_data( nb_grids, file ) ;
      write_data_endl( file ) ;
      write_data_endl( file ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_bound_colors( PEL_ModuleExplorer const* exp,
                                          size_t& nb_bounds,
                                          size_t_vector& db_color_idx,
                                          std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_bound_colors" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( db_color_idx.size() == 0 ) ;
   PEL_CHECK( file ) ;
   
   intArray2D const& face2cell = exp->intArray2D_data( "face2cell" ) ;
   stringVector const& colors = exp->stringVector_data( "color_table" ) ;
   intVector const& face2col = exp->intVector_data( "face_color" ) ;
   size_t const nb_faces = face2col.size() ;

   nb_bounds = 0 ;
   size_t nb_bd_colors = 0 ;
   stringVector bd_colors(0) ;
   db_color_idx.resize( colors.size() ) ;
   db_color_idx.set( PEL::bad_index() ) ;
   for( size_t iface=0 ; iface<nb_faces ; ++iface )
   {
      if( face2cell(1,iface) == index_for_trash() )
      {
         size_t const bd_col = face2col(iface) ;
         if( db_color_idx(bd_col) == PEL::bad_index() )
         {
            db_color_idx(bd_col) = nb_bd_colors ;
            bd_colors.append( colors(bd_col) ) ;
            ++nb_bd_colors ;
         }
         ++nb_bounds ;
      }
   }

   if( !BINARY ) write_str_data( "Boundary Table", file ) ;
   write_data( (int) nb_bd_colors, file ) ;
   write_data_endl( file ) ;
   for( size_t icol=0 ; icol<nb_bd_colors ; ++icol )
   {
      if( BINARY )
      {
         IBUFF2(0) = 0 ; IBUFF2(1) = 0 ;
         write_data( IBUFF2, file ) ;
      }
      else
      {
         IBUFF3(0) = 0 ; IBUFF3(1) = 0 ; IBUFF3(2) = 0 ;
         write_data( IBUFF3, file ) ;
         write_str_data( " ", file ) ;
      }
      write_str_data( bd_colors(icol), file ) ;
      write_data_endl( file ) ;
   }
   write_data_endl( file ) ;
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_vertices( PEL_ModuleExplorer const* exp,
                                      std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_vertices" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( DIM == 2 || DIM == 3 ) ;

   double const epsilon = 1.e-5 ;

   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   int const nb_verts = vertices.index_bound( 1 ) ;
   int const nb_fv_verts = ( DIM == 3 ? nb_verts : 2*nb_verts ) ;
   if( BINARY )
   {
      IBUFF2(0) = FV_NODES ; IBUFF2(1) = nb_fv_verts ;
      write_data( IBUFF2, file ) ;
      for( size_t ic=0 ; ic<3 ; ic++ )
      {
         doubleVector coords( nb_fv_verts ) ;
         if( DIM == 3 )
         {
            vertices.extract_section( 0, ic, coords ) ;
         }
         else
         {
            if( ic == 2 )
            {
               for( int i=0 ; i<nb_verts ; ++i )
               {
                  coords(i) = 0. ;
                  coords(i+nb_verts) = epsilon ;
               }
            }
            else
            {
               for( int i=0 ; i<nb_verts ; ++i )
               {
                  double const xi = vertices( ic, i ) ;
                  coords(i) = xi ;
                  coords(i+nb_verts) = xi ;
               }
            }
         }
         write_data( coords, file ) ;
      }
   }
   else
   {
      write_str_data( "Nodes", file ) ;
      write_data( nb_fv_verts, file ) ;
      write_data_endl( file ) ;
      for( int i=0 ; i<nb_verts ; ++i )
      {
         write_data( vertices(0,i), file ) ;
         write_data( vertices(1,i), file ) ;
         double const z = ( DIM == 3 ? vertices(2,i) : 0. ) ;
         write_data( z, file ) ;
         write_data_endl( file ) ;
      }
      if( DIM == 2 )
      {
         for( int i=0 ; i<nb_verts ; ++i )
         {
            write_data( vertices(0,i), file ) ;
            write_data( vertices(1,i), file ) ;
            write_data( epsilon, file ) ;
            write_data_endl( file ) ;
         }
      }
      write_data_endl( file ) ;
   }   
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_boundaries( PEL_ModuleExplorer const* exp,
                                        size_t nb_bounds,
                                        size_t_vector const& db_color_idx,
                                        std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_boundaries" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( DIM == 2 || DIM == 3 ) ;

   size_t const nb_verts =
      exp->doubleArray2D_data( "vertices" ).index_bound( 1 ) ;
   intVector const& face_nb_verts = exp->intVector_data( "face_nb_vertices" ) ;
   intArray2D const& face2vertex = exp->intArray2D_data( "face2vertex" ) ;
   intArray2D const& face2cell = exp->intArray2D_data( "face2cell" ) ;
   intVector const& face2col = exp->intVector_data( "face_color" ) ;
   size_t const nb_faces = face2col.size() ;

   if( !BINARY )
   {
      write_str_data( "Boundary Faces", file ) ;
      write_data( (int) nb_bounds, file ) ;
      write_data_endl( file ) ;
   }

   for( size_t iface=0 ; iface<nb_faces ; ++iface )
   {
      if( face2cell(1,iface) == index_for_trash() )
      {
         int const bd_idx = (int) db_color_idx( face2col( iface ) )+1 ;
         if( BINARY )
         {
            IBUFF3(0) = FV_FACES ;
            IBUFF3(1) = bd_idx ;
            IBUFF3(2) = 1 ;
            write_data( IBUFF3, file ) ;
         }
         else
         {
            write_data( bd_idx, file ) ;
         }
         size_t const f_nb_verts = face_nb_verts(iface) ;
         size_t const fv_f_nb_verts = ( DIM == 3 ? f_nb_verts : 4 ) ;
         if( DIM == 3 )
         {
            const int* pel2fv = 0 ;
            if( f_nb_verts==3 )
            {
               pel2fv = PEL2FV_TRIA ;
            }
            else if( f_nb_verts==4 )
            {
               pel2fv = PEL2FV_QUAD ;
            }
            else
            {
               raise_internal_error(
                  "In 3D, triangle or quadrangle is expected as boundary" ) ;
            }
            IBUFF4(3) = 0 ;
            for( size_t iv=0 ; iv<f_nb_verts ; ++iv )
            {
               IBUFF4(iv) = face2vertex( pel2fv[iv], iface )+1 ;
            }
         }
         else
         {
            if( f_nb_verts != 2 )
            {
               raise_internal_error(
                  "In 2D, segment is expected as boundary" ) ;
            }
            IBUFF4(0) = face2vertex( 0, iface )+1 ;
            IBUFF4(1) = face2vertex( 1, iface )+1 ;
            IBUFF4(2) = face2vertex( 1, iface )+1+nb_verts ;
            IBUFF4(3) = face2vertex( 0, iface )+1+nb_verts ;
         }
         if( BINARY )
         {
            write_data( IBUFF4, file ) ;
         }
         else
         {
            write_data( (int) fv_f_nb_verts, file ) ;
            for( size_t iv=0 ; iv<fv_f_nb_verts ; ++iv )
            {
               write_data( IBUFF4( iv), file ) ;
            }
         }
         write_data_endl( file ) ;
      }
   }
   write_data_endl( file ) ;
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_cells( PEL_ModuleExplorer const* exp,
                                   std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_cells" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( DIM == 2 || DIM == 3 ) ;
   
   size_t const nb_verts =
                       exp->doubleArray2D_data( "vertices" ).index_bound( 1 ) ;
   intVector const& cell_nb_verts = exp->intVector_data( "cell_nb_vertices" ) ;
   intArray2D const& cell2vertex = exp->intArray2D_data( "cell2vertex" ) ;
   size_t const nb_cells = cell_nb_verts.size() ;
   
   if( BINARY )
   {
      IBUFF5.set(0) ;
      IBUFF5(0) = FV_ELEMENTS ;
      for( size_t icell=0 ; icell<nb_cells ; ++icell )
      {
         size_t const c_nb_verts = cell_nb_verts(icell) ;
         if( DIM == 3 && c_nb_verts == 4 ) /* tet count */
         {
            IBUFF5(TETRA_IDX) += 1 ;
         }
         else if( DIM == 3 && c_nb_verts == 8 ) /* hex count */
         {
            IBUFF5(HEXA_IDX) += 1 ;
         }
         else if( DIM == 2 && c_nb_verts == 3 ) /* prism count */
         {
            IBUFF5(PRISM_IDX) += 1 ;
         }
         else if( DIM == 2 && c_nb_verts == 4 ) /* hex count */
         {
            IBUFF5(HEXA_IDX) += 1 ;
         }
         else
         {
            raise_internal_error( "invalid cell" ) ;
         }
      }
      write_data( IBUFF5, file ) ;
   }
   else
   {
      write_str_data( "Elements", file ) ;
      write_data_endl( file ) ;
   }
   
   for( size_t icell=0 ; icell<nb_cells ; ++icell )
   {
      size_t const c_nb_verts = cell_nb_verts(icell) ;
      size_t const fv_c_nb_verts = ( DIM == 3 ? c_nb_verts : 2*c_nb_verts ) ;
      int elem = 0 ;
      const int* pel2fv = 0 ;
      if( DIM == 2 && c_nb_verts == 3 )
      {
         elem = PRISM_IDX ;
      }
      else if( DIM == 2 && c_nb_verts == 4 )
      {
         elem = HEXA_IDX ;
      }
      else if( DIM == 3 && c_nb_verts == 4 )
      {
         elem = TETRA_IDX ;
         pel2fv = PEL2FV_TETRA ;
      }
      else if( DIM == 3 && c_nb_verts == 8 )
      {
         elem = HEXA_IDX ;
         pel2fv = PEL2FV_HEXA ;
      }
      else
      {
         raise_internal_error( "Unexpected cell" ) ;
      }
      if( BINARY )
      {
         unsigned int header = (unsigned int) fv_encode_elem_header( elem ) ;
         file.write( (const char*)&header, sizeof(header) ) ;
      }
      else
      {
         write_data( elem, file ) ;
         write_data( (int) 1, file ) ;
         write_data_endl( file ) ;
      }
      intVector conn( fv_c_nb_verts ) ;
      if( DIM == 3 )
      {
         for( size_t iv=0 ; iv<fv_c_nb_verts ; ++iv )
         {
            conn( iv ) = cell2vertex( pel2fv[iv], icell )+1 ;
         }
      }
      else if( c_nb_verts == 3 )
      {
         conn(0) = cell2vertex( 0, icell )+1 ;
         conn(1) = cell2vertex( 0, icell )+1+nb_verts ;
         conn(2) = cell2vertex( 1, icell )+1+nb_verts ;
         conn(3) = cell2vertex( 1, icell )+1 ;
         conn(4) = cell2vertex( 2, icell )+1+nb_verts ;
         conn(5) = cell2vertex( 2, icell )+1 ;
      }
      else if( c_nb_verts == 4 )
      {
         conn(0) = cell2vertex( 0, icell )+1 ;
         conn(1) = cell2vertex( 1, icell )+1 ;
         conn(2) = cell2vertex( 3, icell )+1 ;
         conn(3) = cell2vertex( 2, icell )+1 ;
         conn(4) = cell2vertex( 0, icell )+1+nb_verts ;
         conn(5) = cell2vertex( 1, icell )+1+nb_verts ;
         conn(6) = cell2vertex( 3, icell )+1+nb_verts ;
         conn(7) = cell2vertex( 2, icell )+1+nb_verts ;
      }      
      write_data( conn, file ) ;
      write_data_endl( file ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_result_header( std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_result_header" ) ;
   PEL_CHECK( file ) ;

   /* Output the magic number. */
   if( BINARY ) write_data( (int) FV_MAGIC, file ) ;
   
   /* Output file header and version number. */
   std::string const head = ( BINARY ? "FIELDVIEW" : "FIELDVIEW_Results" ) ;
   write_str_data( head, file ) ;

   /* This version of the FIELDVIEW unstructured file is "3.0" */
   IBUFF2(0)=3 ; IBUFF2(1)=0 ;
   write_data( IBUFF2, file ) ;
   write_data_endl( file ) ;

   /* Comments */
   write_data_endl( file ) ;
   write_comment( "File generated by PELICANS for FieldView post-processing",
                  file ) ;
   write_comment( "Result file", file ) ;
   write_data_endl( file ) ;
   
   
   if( BINARY )
   {
      /* File type code - new in version 2.7 */
      write_data( FV_RESULTS_FILE, file ) ;

      /* Reserved field, always write a zero - new in version 2.6 */
      write_data( (int) 0, file ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_constants( PEL_ModuleExplorer const* exp,
                                       std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_constants" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( file ) ;

   if( !BINARY )
   {
      write_str_data("Constants", file ) ;
      write_data_endl( file ) ;
   }

   DBUFF4.set(0.) ; // time, fsmach, alpha and re.
   if( exp->has_module( "variables" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "variables" ) ;
      if( se->has_entry( "TIME" ) )
      {
         DBUFF4(0) = se->double_data( "TIME" ) ;
      }
      se->destroy() ; se = 0 ;
   }
   write_data( DBUFF4, file ) ;
   write_data_endl( file ) ;
   write_data_endl( file ) ;

   // Nb grids:
   int const nb_grids = 1 ;
   if( BINARY )
   {
      write_data( nb_grids, file ) ;
   }
   else
   {
      write_str_data( "Grids", file ) ;
      write_data( nb_grids, file ) ;
      write_data_endl( file ) ;
      write_data_endl( file ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_fields( PEL_ModuleExplorer const* exp,
                                    std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_fields" ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( file ) ;

   PEL_ModuleExplorer* se = 0 ;
   if( exp->has_module( "fields" ) )
      se = exp->create_subexplorer( 0, "fields" ) ;

   int nb_vars = 0 ;
   int nb_verts = 0 ;
   stringVector var_names( 0 ) ;
   if( se != 0 )
   {
      for( se->start_module_iterator() ;
           se->is_valid_module() ;
           se->go_next_module() )
      {
         PEL_ModuleExplorer const* fexp = se->create_subexplorer( 0 ) ;
         std::string const& name = fexp->string_data( "name" ) ;
         std::string const& location = fexp->string_data( "location" ) ;
         if( location != "at_vertices" )
         {
            raise_field_location_error( name, location, "at_vertices" ) ;
         }
         doubleArray2D const& X = fexp->doubleArray2D_data( "value" ) ;
         size_t const nb_comps = X.index_bound(0) ;
         nb_verts = X.index_bound(1) ;
         if( DIM == 2 ) nb_verts *= 2 ;
         if( nb_comps == 1 )
         {
            nb_vars += 1 ;
            var_names.append( name ) ;
         }
         else if( nb_comps == DIM )
         {
            nb_vars += 3 ;
            var_names.append( name+"_0; "+name ) ;
            var_names.append( name+"_1" ) ;
            var_names.append( name+"_2" ) ;
         }
         else
         {
            nb_vars += nb_comps ;
            for( size_t ic=0 ; ic<nb_comps ; ++ic )
            {
               std::ostringstream tmp ;
               tmp << ic ;
               var_names.append( name+"_"+tmp.str() ) ;
            }
         }
         fexp->destroy() ; fexp = 0 ;
      }
   }
   
   if( !BINARY ) write_str_data( "Variable Names", file ) ;
   write_data( nb_vars, file ) ;
   write_data_endl( file ) ;
   write_data( var_names, file ) ;
   write_data_endl( file ) ;
   
   if( !BINARY ) write_str_data( "Boundary Variable Names", file ) ;
   write_data( (int) 0, file ) ;
   write_data_endl( file ) ;
   write_data_endl( file ) ;

   if( se != 0 )
   {
      if( BINARY )
      {
         IBUFF2(0) = FV_NODES ; IBUFF2(1) = nb_verts ;
         write_data( IBUFF2, file ) ;
      }
      else
      {
         write_str_data( "Nodes", file ) ;
         write_data( nb_verts, file ) ;
         write_data_endl( file ) ;
         write_data_endl( file ) ;
      }

      if( BINARY )
      {
         write_data( FV_VARIABLES, file ) ;
      }
      else
      {
         write_str_data( "Variables", file ) ;
         write_data_endl( file ) ;
      }
      
      doubleVector x( nb_verts ) ;
      size_t ivar = 0 ;
      for( se->start_module_iterator() ;
           se->is_valid_module() ;
           se->go_next_module() )
      {
         PEL_ModuleExplorer const* fexp = se->create_subexplorer( 0 ) ;
         doubleArray2D const& X = fexp->doubleArray2D_data( "value" ) ;
         size_t const nb_comps = X.index_bound(0) ;
         size_t const nb_vs = X.index_bound(1) ;
         for( size_t ic=0 ; ic<nb_comps ; ++ic )
         {
            for( size_t i=0 ; i<nb_vs ; ++i )
            {
               x(i) = X(ic,i) ;
               if( DIM == 2 ) x(i+nb_vs) = X(ic,i) ;
            }
            write_comment( var_names(ivar++), file ) ;
            write_data( x, file ) ;
            write_data_endl( file ) ;
            write_data_endl( file ) ;
         }
         if( nb_comps == DIM && DIM == 2 )
         {
            x.set( 0. ) ;
            write_comment( var_names(ivar++), file ) ;
            write_data( x, file ) ;
            write_data_endl( file ) ;
            write_data_endl( file ) ;
         }
         fexp->destroy() ; fexp = 0 ;
      }
      if( BINARY )
      {
         write_data( FV_BNDRY_VARS, file ) ;
      }
   }
   
   if( se != 0 ) se->destroy() ;
}

//-----------------------------------------------------------------------
int
PEL_FieldViewWriter:: fv_encode_elem_header( int elem ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: fv_encode_elem_header" ) ;
   PEL_CHECK( elem == TETRA_IDX ||
              elem == HEXA_IDX ||
              elem == PRISM_IDX ||
              elem == PYRA_IDX ) ;

   int result = 0 ;
   switch( elem )
   {
      case TETRA_IDX:
         result = ( 1 << ELEM_TYPE_BIT_SHIFT ) ;
         break ;
      case PYRA_IDX:
         result = ( 2 << ELEM_TYPE_BIT_SHIFT ) ;
         break ;
      case PRISM_IDX:
         result = ( 3 << ELEM_TYPE_BIT_SHIFT ) ;
         break ;
      case HEXA_IDX:
         result = ( 4 << ELEM_TYPE_BIT_SHIFT ) ;
         break ;
      default:
         raise_internal_error( "Unexpected element" ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_str_data( std::string const& val,
                                      std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_str_data(string)" ) ;

   static const size_t ss = sizeof(char) ;
   static char cbuf[80] ;
   
   if( !BINARY )
   {
      file << val ;
   }
   else
   {
      size_t const len = val.size() ;
      for( size_t i = 0 ; i < ( len < 80 ? len : 80 ) ; ++i )
         cbuf[i] = val[i] ;
      for( size_t i = len ; i < 80 ; ++i )
        cbuf[i] = '\0';  /* pad with zeros */
      file.write( (const char*) cbuf, 80*ss ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data( int val, std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data(int)" ) ;

   static const size_t si = sizeof(int) ;
   if( !BINARY )
   {
      file << " " << val ;
   }
   else
   {
      file.write( (const char*) &val, si ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data( double val, std::ofstream& file )  const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data(double)" ) ;
   
   static const size_t sd = sizeof(float) ;
   float fval = (float)val ;
   if( !BINARY )
   {
      file << " " << fval ;
      if( (int) fval == fval ) file << ".0" ;
   }
   else
   {
      file.write( (const char*) &fval, sd ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data( intVector const& val, 
                                  std::ofstream& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data(intVector)" ) ;

   static const size_t si = sizeof(int) ;
   if( !BINARY )
   {
      size_t j=1 ;
      for( size_t i=0 ; i<val.size() ; ++i, ++j )
      {
         write_data( val(i), file ) ;
         if( j == 10 )
         {
            file << std::endl ;
            j = 0 ;
         }
      }
   }
   else
   {
      size_t const nb = val.size() ;
      int *dval = new int[nb] ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dval[i] = val(i) ;
      }      
      file.write( (const char*) dval, nb*si ) ;
      delete [] dval ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data( doubleVector const& val, 
                                  std::ofstream& file ) const 
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data(doubleVector)" ) ;
   
   static const size_t sd = sizeof(float) ;
   if( !BINARY )
   {
      size_t j=1 ;
      for( size_t i=0 ; i<val.size() ; ++i, ++j )
      {
         write_data( val(i), file ) ;
         if( j == 10 )
         {
            file << std::endl ;
            j = 0 ;
         }
      }
   }
   else
   {
      size_t const nb = val.size() ;
      float *dval = new float[nb] ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dval[i] = (float) val(i) ;
      }
      file.write( (const char*) dval, nb*sd ) ;
      delete [] dval ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data( stringVector const& val, 
                                  std::ofstream& file ) const 
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data(stringVector)" ) ;

   for( size_t i=0 ; i<val.size() ; ++i )
   {
      write_str_data( val(i), file ) ;
      write_data_endl( file ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_comment( std::string const& val, 
                                     std::ofstream& file ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_comment" ) ;   
   
   if( !BINARY )
   {
      file << "! " << val << std::endl ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: write_data_endl( std::ofstream& file ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FieldViewWriter:: write_data_endl" ) ;   
   
   if( !BINARY )
   {
      file << std::endl ;
   }
}

//-----------------------------------------------------------------------
void
PEL_FieldViewWriter:: raise_internal_error( std::string const& mes ) const 
//-----------------------------------------------------------------------
{
   PEL_Error::object()->raise_internal( "*** PEL_FieldViewWriter error\n"+mes ) ;
}
