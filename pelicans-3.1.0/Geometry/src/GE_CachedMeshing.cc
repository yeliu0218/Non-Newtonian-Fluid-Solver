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

#include <GE_CachedMeshing.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <fstream>
#include <iostream>
#include <string>

GE_CachedMeshing const* 
GE_CachedMeshing:: PROTOTYPE = new GE_CachedMeshing() ;

//-------------------------------------------------------------------------
GE_CachedMeshing:: GE_CachedMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_CachedMeshing" )
   , VERTICES( 0, 0 )
   , VERTEX_COLORS( 0 )
   , CELL_FACES( 0, 0 )
   , CELL_VERTICES( 0, 0 )
   , CELL_COLORS( 0 )
   , FACE_VERTICES( 0, 0 )
   , FACE_COLORS( 0 )
{
}

//-------------------------------------------------------------------------
GE_CachedMeshing*
GE_CachedMeshing:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp,
                                     size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_CachedMeshing* result = new GE_CachedMeshing(  a_owner, exp, dim_space ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_CachedMeshing:: GE_CachedMeshing( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp,
                                     size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , VERTICES( 0, 0 )
   , VERTEX_COLORS( 0 )
   , VERT_INDEX( PEL::bad_index() )
   , CELL_FACES( 0, 0 )
   , CELL_VERTICES( 0, 0 )
   , CELL_COLORS( 0 )
   , CELL_INDEX( PEL::bad_index() )
   , FACE_VERTICES( 0, 0 )
   , FACE_COLORS( 0 )
   , FACE_INDEX( PEL::bad_index() )
{
   PEL_LABEL( "GE_CachedMeshing:: GE_CachedMeshing" ) ;
   bool read = false ;
   bool has_meshing = exp->has_module( "GE_Meshing" ) ;
   std::string filename = ( exp->has_entry( "filename" ) ?
                            exp->string_data( "filename" ) : "" ) ;
   if( !filename.empty() )
   {
      std::ifstream in( filename.c_str() ) ;
      if( in )
      {
         in.close() ;
         read_meshing( filename ) ;
         read = true ;
      }
   }
   
   if( !read )
   {
      if( !has_meshing ) 
      {
         PEL_Error::object()->raise_plain(
            "GE_CachedMeshing : can't read cached meshing and no \n internal one is provided in missing GE_Meshing module" ) ;
      }
      
      PEL_ModuleExplorer* se = 
         exp->create_subexplorer( 0, "GE_Meshing" ) ;
      GE_Meshing* meshing = GE_Meshing::create( 0, se, dim_space ) ;

      // Vertices :
      build_vertices(meshing) ;

      // Cells :
      build_cells(meshing) ;

      // Faces :
      build_faces(meshing) ;
   
      se->destroy() ;
      meshing->destroy() ;
   }

   if( !filename.empty() && !read )
   {
      save_meshing( filename ) ;
   }
   
}

//------------------------------------------------------------------------
GE_CachedMeshing:: ~GE_CachedMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_CachedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
  
   return( VERTICES.index_bound(1)) ;
}

//------------------------------------------------------------------------
size_t
GE_CachedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: nb_cells" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( CELL_VERTICES.index_bound(1) ) ;
}

//------------------------------------------------------------------------
size_t
GE_CachedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: nb_faces" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( FACE_VERTICES.index_bound(1) ) ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   VERT_INDEX = 0 ;
}

//------------------------------------------------------------------------
bool
GE_CachedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( VERT_INDEX<nb_vertices() ) ;   
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   VERT_INDEX++ ;   
}

//------------------------------------------------------------------------
doubleVector const& 
GE_CachedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   static doubleVector result(nb_space_dimensions()) ;
   VERTICES.extract_section( 1, VERT_INDEX, result ) ;
   
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   CELL_INDEX = 0 ;
   
   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   CELL_INDEX++ ;
}

//------------------------------------------------------------------------
bool
GE_CachedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: valid_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return CELL_INDEX < nb_cells() ;
   
}

//------------------------------------------------------------------------
std::string const&
GE_CachedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLYHEDRON ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_CachedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   static size_t_vector result(CELL_VERTICES.index_bound(0)) ;

   CELL_VERTICES.extract_section( 1, CELL_INDEX, result ) ;
   
   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result );
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_CachedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   static size_t_vector result(CELL_FACES.index_bound(0)) ;

   CELL_FACES.extract_section( 1, CELL_INDEX, result ) ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   FACE_INDEX = 0 ;
   
   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   FACE_INDEX++ ;
}

//------------------------------------------------------------------------
bool
GE_CachedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: valid_face" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = FACE_INDEX<nb_faces() ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_CachedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = FACE_POLYHEDRON ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_CachedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   static size_t_vector result(FACE_VERTICES.index_bound(0)) ;

   FACE_VERTICES.extract_section( 1, FACE_INDEX, result ) ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_CachedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Color const* result  = static_cast<GE_Color const*>(
      VERTEX_COLORS->at( VERT_INDEX ) ) ;

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_CachedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result  = static_cast<GE_Color const*>(
      CELL_COLORS->at( CELL_INDEX ) ) ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_CachedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result  = static_cast<GE_Color const*>(
      FACE_COLORS->at( FACE_INDEX ) ) ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: build_vertices( GE_Meshing* meshing )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: build_vertices" ) ;

   size_t nb = meshing->nb_vertices() ;
   VERTICES.re_initialize(nb_space_dimensions(), nb ) ;
   VERTEX_COLORS = PEL_Vector::create( this, nb ) ;
   
   size_t i=0 ;
   for( meshing->start_vertex_iterator() ;
        meshing->valid_vertex() ;
        ++i, meshing->go_next_vertex() )
   {
      VERTICES.set_section( 1, i, meshing->vertex_coordinates() ) ;
      VERTEX_COLORS->set_at( i,
                             const_cast<GE_Color*>( meshing->vertex_color() ) ) ;
   }
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: build_cells( GE_Meshing* meshing )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: build_cells" ) ;

   size_t nb = meshing->nb_cells() ;

   meshing->start_cell_iterator() ;
   PEL_ASSERT( meshing->valid_cell() ) ;
   
   CELL_POLYHEDRON = meshing->cell_polyhedron_name() ;
   size_t nbs = meshing->cell_faces().size() ;
   size_t nbv = meshing->cell_vertices().size() ;

   CELL_FACES.re_initialize( nbs, nb ) ;
   CELL_VERTICES.re_initialize( nbv, nb ) ;
   CELL_COLORS = PEL_Vector::create( this, nb ) ;
   size_t i=0 ;
   
   for( meshing->start_cell_iterator() ;
        meshing->valid_cell() ;
        ++i, meshing->go_next_cell() )
   {
      PEL_ASSERT( meshing->cell_polyhedron_name() == CELL_POLYHEDRON ) ;
      
      size_t_vector const& connect = meshing->cell_vertices() ;
      size_t_vector const& faces = meshing->cell_faces() ;
      GE_Color* col = const_cast<GE_Color*>( meshing->cell_color() ) ;
      CELL_FACES.set_section( 1, i, faces ) ;
      CELL_VERTICES.set_section( 1, i, connect ) ;
      CELL_COLORS->set_at( i, col ) ;
   } 
}


//------------------------------------------------------------------------
void
GE_CachedMeshing:: build_faces( GE_Meshing* meshing )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: build_faces" ) ;

   size_t nb = meshing->nb_faces() ;

   meshing->start_face_iterator() ;
   PEL_ASSERT( meshing->valid_face() ) ;
   
   FACE_POLYHEDRON = meshing->face_polyhedron_name() ;
   size_t nbv = meshing->face_vertices().size() ;

   FACE_VERTICES.re_initialize( nbv, nb ) ;
   FACE_COLORS = PEL_Vector::create( this, nb ) ;
   size_t i=0 ;
   
   for( meshing->start_face_iterator() ;
        meshing->valid_face() ;
        ++i, meshing->go_next_face() )
   {
      PEL_ASSERT( meshing->face_polyhedron_name() == FACE_POLYHEDRON ) ;
      
      size_t_vector const& connect = meshing->face_vertices() ;
      GE_Color* col = const_cast<GE_Color*>( meshing->face_color() ) ;
      FACE_VERTICES.set_section( 1, i, connect ) ;
      FACE_COLORS->set_at( i, col ) ;
   } 
}


//------------------------------------------------------------------------
void
GE_CachedMeshing:: read_meshing( std::string const& filename )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: read_meshing" ) ;

   // Prevent to duplicate null color
   GE_Color::null_color() ;
   
   PEL_Module* read = PEL_Module::create( 0, "Meshing", filename ) ;
   PEL_ASSERT( read!=0 ) ;
   PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( read, read ) ;
   exp = exp->create_subexplorer( exp, "Meshing" ) ;
   
   stringVector const& col = exp->stringVector_data( "COLOR_TABLE" ) ;
   for( size_t i=0 ; i<col.size() ; i++ )
      GE_Color::extend( col(i) ) ;
   
   VERTICES = exp->doubleArray2D_data( "VERTICES" ) ;
   intVector const& colv = exp->intVector_data( "VERTEX_COLORS" ) ;
   
   VERTEX_COLORS = PEL_Vector::create( this, colv.size() ) ;
   for( size_t i=0 ; i<colv.size() ; i++ )
      VERTEX_COLORS->set_at( i, const_cast<GE_Color*>(
                              GE_Color::object( col( colv( i ) ) ) ) ) ;

   CELL_POLYHEDRON = exp->string_data( "CELL_POLYHEDRON" ) ;
   CELL_FACES.set( exp->intArray2D_data( "CELL_FACES" ) ) ;
   CELL_VERTICES.set( exp->intArray2D_data( "CELL_VERTICES" ) ) ;
   intVector const& colc = exp->intVector_data( "CELL_COLORS" ) ;
   CELL_COLORS = PEL_Vector::create( this, colc.size() ) ;
   for( size_t i=0 ; i<colc.size() ; i++ )
      CELL_COLORS->set_at( i, const_cast<GE_Color*>(
                              GE_Color::object( col( colc( i ) ) ) ) ) ;
   
   FACE_POLYHEDRON = exp->string_data( "FACE_POLYHEDRON" ) ;
   FACE_VERTICES.set( exp->intArray2D_data( "FACE_VERTICES" ) ) ;
   intVector const& cols = exp->intVector_data( "FACE_COLORS" ) ;
   FACE_COLORS = PEL_Vector::create( this, cols.size() ) ;
   for( size_t i=0 ; i<cols.size() ; i++ )
      FACE_COLORS->set_at( i, const_cast<GE_Color*>(
                              GE_Color::object( col( cols( i ) ) ) ) ) ;

   read->destroy() ;
   
}

//------------------------------------------------------------------------
void
GE_CachedMeshing:: save_meshing( std::string const& filename )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CachedMeshing:: save_meshing" ) ;

   PEL_Module* mod = PEL_Module::create( 0, "Meshing" ) ;
   
   stringVector const& colt = GE_Color::color_table() ;
   mod->add_entry( "COLOR_TABLE", PEL_StringVector::create( mod, colt ) ) ;
   mod->add_entry( "VERTICES", PEL_DoubleArray2D::create( mod, VERTICES ) ) ;
   
   intVector col(VERTEX_COLORS->index_limit()) ;
   for( size_t i=0 ; i<col.size() ; i++ )
      col(i) = static_cast<GE_Color const*>(
         VERTEX_COLORS->at( i ) )->identifier() ;
   mod->add_entry( "VERTEX_COLORS", PEL_IntVector::create( mod, col ) ) ;

   mod->add_entry( "CELL_POLYHEDRON",
                   PEL_String::create( mod, CELL_POLYHEDRON ) ) ;
   mod->add_entry( "CELL_FACES",
                   PEL_IntArray2D::create( mod, CELL_FACES ) ) ;
   mod->add_entry( "CELL_VERTICES",
                   PEL_IntArray2D::create( mod, CELL_VERTICES ) ) ;
   
   col.re_initialize( CELL_COLORS->index_limit() ) ;
   for( size_t i=0 ; i<col.size() ; i++ )
      col(i) = static_cast<GE_Color const*>(
         CELL_COLORS->at( i ) )->identifier() ;
   mod->add_entry( "CELL_COLORS", PEL_IntVector::create( mod, col ) ) ;
   
   mod->add_entry( "FACE_POLYHEDRON",
                   PEL_String::create( mod, FACE_POLYHEDRON ) ) ;
   mod->add_entry( "FACE_VERTICES",
                   PEL_IntArray2D::create( mod, FACE_VERTICES ) ) ;
   
   col.re_initialize( FACE_COLORS->index_limit() ) ;
   for( size_t i=0 ; i<col.size() ; i++ )
      col(i) = static_cast<GE_Color const*>(
         FACE_COLORS->at( i ) )->identifier() ;
   mod->add_entry( "FACE_COLORS", PEL_IntVector::create( mod, col ) ) ;
   
   mod->write( filename, "hybrid" ) ;
   
   mod->destroy() ;
}

