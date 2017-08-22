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

#include <GE_TransformedMeshing.hh>

#include <PEL_Bool.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;

struct GE_TransformedMeshing_ERROR
{
   static void n0( stringVector const& undef_vars,
                   stringVector const& predef_vars ) ;
   static void n1( size_t nb_sp_dims ) ;
} ;

GE_TransformedMeshing const* 
GE_TransformedMeshing:: PROTOTYPE = new GE_TransformedMeshing() ;

//-----------------------------------------------------------------------------
GE_TransformedMeshing:: GE_TransformedMeshing( void )
//-----------------------------------------------------------------------------
   : GE_Meshing( "GE_TransformedMeshing" )
   , INITIAL( 0 )
   , CTX( 0 )
   , COORDS( 0 )
   , HMIN( 0 )
   , IS_BOUNDARY( 0 )
   , VERT_IDX( 0 )
   , TRANSFO( 0 )
   , IVERT( PEL::bad_index() )
   , TRANSFORMED_COORDINATES( 0, 0 )
{
}

//-----------------------------------------------------------------------------
GE_TransformedMeshing*
GE_TransformedMeshing:: create_replica( PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp,
                                        size_t dim_space ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_TransformedMeshing* result =
                         new GE_TransformedMeshing( a_owner, exp, dim_space ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-----------------------------------------------------------------------------
GE_TransformedMeshing:: GE_TransformedMeshing( PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp,
                                               size_t dim_space )
//-----------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , INITIAL( 0 )
   , CTX( 0 )
   , COORDS( 0 )
   , HMIN( 0 )
   , IS_BOUNDARY( 0 )
   , VERT_IDX( 0 )
   , TRANSFO( 0 )
   , IVERT( PEL::bad_index() )
   , TRANSFORMED_COORDINATES( 0, 0 )
{
   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "GE_Meshing" ) ;
   INITIAL = GE_Meshing::create( this, se, dim_space ) ;
   se->destroy() ; se = 0 ;
   
   // Context for formula :
   {
      CTX = PEL_ContextSimple::create( this ) ;
      COORDS = PEL_DoubleVector::create( CTX, doubleVector(0) ) ;
      CTX->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;

      PEL_DataWithContext* data = exp->abstract_data( 0, "transformation" ) ;
      stringVector const& vars = data->undefined_variables( 0 ) ;
      if( vars.has( "IS_VERT_IDX" ) )
      {
         VERT_IDX = PEL_Int::create( CTX, PEL::bad_int() ) ;
         CTX->extend( PEL_Variable::object( "IS_VERT_IDX" ), VERT_IDX ) ;
      }
      if( vars.has( "DS_VERT_HMIN" ) )
      {
         HMIN = PEL_Double::create( CTX, PEL::bad_double() ) ;
         CTX->extend( PEL_Variable::object( "DS_VERT_HMIN" ), HMIN ) ;
      }
      if( vars.has( "BS_VERT_ON_BOUND" ) )
      {
         IS_BOUNDARY = PEL_Bool::create( CTX, false ) ;
         CTX->extend( PEL_Variable::object( "BS_VERT_ON_BOUND" ), IS_BOUNDARY ) ;
      }
      data->destroy() ; data = 0 ;
   }

   // Transformation expression:
   TRANSFO = exp->abstract_data( CTX, "transformation", CTX ) ;
   if( TRANSFO->data_type() != PEL_Data::DoubleVector )
   {
      PEL_Error::object()->raise_bad_data_type( exp, "transformation", 
                                                PEL_Data::DoubleVector ) ;
   }
   if( !TRANSFO->value_can_be_evaluated() )
   {
      stringVector const predef_vars =
                        "DV_X,IS_VERT_IDX,DS_VERT_HMIN,BS_VERT_ON_BOUND" ;
      GE_TransformedMeshing_ERROR:: n0( TRANSFO->undefined_variables(),
                                        predef_vars ) ;
   }

   // Meshing characteristics:
   set_meshing_data() ;
}

//------------------------------------------------------------------------
GE_TransformedMeshing:: ~GE_TransformedMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_TransformedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: nb_vertices" ) ;

   return( INITIAL->nb_vertices() ) ;
}

//------------------------------------------------------------------------
size_t
GE_TransformedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: nb_cells" ) ;

   return( INITIAL->nb_cells() ) ;
}

//------------------------------------------------------------------------
size_t
GE_TransformedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: nb_faces" ) ;

   return( INITIAL->nb_faces() ) ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   IVERT = 0 ;
   INITIAL->start_vertex_iterator() ;
}

//------------------------------------------------------------------------
bool
GE_TransformedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( INITIAL->valid_vertex() ) ;   
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++IVERT ;
   INITIAL->go_next_vertex() ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_TransformedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   static doubleVector result( nb_space_dimensions() ) ;
   if( result.size() != nb_space_dimensions() ) 
   {
      result.re_initialize( nb_space_dimensions() ) ;
   }
   
   TRANSFORMED_COORDINATES.extract_section( 0, IVERT, result ) ;
   
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   INITIAL->start_cell_iterator() ;
   
   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_TransformedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: valid_cell" ) ;

   bool result = ( INITIAL->valid_cell() )  ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   INITIAL->go_next_cell() ;
}

//------------------------------------------------------------------------
std::string const&
GE_TransformedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->cell_polyhedron_name() ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_TransformedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = INITIAL->cell_vertices() ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_TransformedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = INITIAL->cell_faces() ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   INITIAL->start_face_iterator() ;
   
   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_TransformedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: valid_face" ) ;

   bool result = ( INITIAL->valid_face() ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   INITIAL->go_next_face() ;
}

//------------------------------------------------------------------------
std::string const&
GE_TransformedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->face_polyhedron_name() ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_TransformedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = INITIAL->face_vertices() ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   INITIAL->print( os, indent_width+3 ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_TransformedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = INITIAL->vertex_color() ;

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_TransformedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = INITIAL->cell_color() ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_TransformedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = INITIAL->face_color() ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_TransformedMeshing:: set_meshing_data( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_TransformedMeshing:: set_meshing_data" ) ;
   
   size_t const nb_verts = INITIAL->nb_vertices() ;
   size_t const dim_space = nb_space_dimensions() ;
   doubleArray2D coords(0,0) ;
   coords.re_initialize( nb_verts, dim_space ) ;
   doubleVector hmin_table( nb_verts, PEL::max_double() ) ;
   boolVector is_boundary_table( 0 ) ;
   INITIAL->start_vertex_iterator() ;
   for( size_t iv = 0 ;
        INITIAL->valid_vertex() ;
        ++iv, INITIAL->go_next_vertex() )
   {
      coords.set_section( 0, iv, INITIAL->vertex_coordinates() ) ;
   }
   
   if( HMIN != 0 || IS_BOUNDARY != 0 || VERT_IDX != 0 )
   {
      if( HMIN != 0 || IS_BOUNDARY != 0 )
      {
         size_t_vector faces_to_cells( 0 ) ;
         if( IS_BOUNDARY != 0 )
         {
            faces_to_cells.re_initialize( INITIAL->nb_faces() ) ;
            faces_to_cells.set( 0 ) ;
         }
         for( INITIAL->start_cell_iterator() ;
              INITIAL->valid_cell() ;
              INITIAL->go_next_cell() )
         {
            if( IS_BOUNDARY != 0 )
            {
               size_t_vector const& cf = INITIAL->cell_faces() ;
               for( size_t i=0 ; i<cf.size() ; ++i )
               {
                  ++faces_to_cells( cf(i) ) ;
               }
            }
            if( HMIN != 0 )
            {
               size_t_vector const& vf = INITIAL->cell_vertices() ;
               for( size_t i=0 ; i<vf.size() ; ++i )
               {
                  size_t ii = vf(i) ;
                  for( size_t j=i+1 ; j<vf.size() ; ++j )
                  {
                     size_t jj = vf(j) ;
                     double d = 0. ;
                     for( size_t ic=0 ; ic<dim_space ; ++ic )
                     {
                        double const xx = coords(ii,ic)-coords(jj,ic) ;
                        d += xx*xx ;
                     }
                     d = PEL::sqrt( d ) ;
                     hmin_table(ii) = PEL::min( hmin_table(ii), d ) ;
                     hmin_table(jj) = PEL::min( hmin_table(jj), d ) ;
                  }
               }
            }
         }
         if( IS_BOUNDARY != 0 )
         {
            is_boundary_table.re_initialize( nb_verts ) ;
            is_boundary_table.set( false ) ;
            INITIAL->start_face_iterator() ;
            for( size_t i = 0 ;
                 INITIAL->valid_face() ;
                 INITIAL->go_next_face(), ++i )
            {
               if( faces_to_cells(i) == 1 )
               {
                  size_t_vector const& sv = INITIAL->face_vertices() ;
                  for( size_t iv=0 ; iv<sv.size() ; ++iv )
                  {
                     is_boundary_table( sv(iv) ) = true ;
                  }
               }
            }
	 }
      }
   }
   TRANSFORMED_COORDINATES.re_initialize( nb_verts, dim_space ) ;
   
   doubleVector coord( dim_space ) ;
 
   // Compute transformed coordinates
   for( size_t iv=0 ; iv<nb_verts ; iv++ )
   {
      coords.extract_section( 0, iv, coord ) ;
      COORDS->set( coord ) ;
      if( VERT_IDX ) VERT_IDX->set( (int) iv ) ;
      if( HMIN != 0 ) HMIN->set( hmin_table( iv ) ) ;
      if( IS_BOUNDARY != 0 ) IS_BOUNDARY->set( is_boundary_table( iv ) ) ;
      
      doubleVector const& result = TRANSFO->to_double_vector() ;
      if( result.size() != nb_space_dimensions() )
         GE_TransformedMeshing_ERROR:: n1( nb_space_dimensions() ) ;
      TRANSFORMED_COORDINATES.set_section( 0, iv, result ) ;
   }
   
}

//internal--------------------------------------------------------------
void 
GE_TransformedMeshing_ERROR:: n0( stringVector const& undef_vars,
                                  stringVector const& predef_vars )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** GE_TransformedMeshing error:" << endl
        << "    the data of keyword \"transformation\"" << endl
        << "    cannot be evaluated" << endl ;
   if( undef_vars.size() > 0 )
   {
      mesg << "    undefined variable(s): " << endl ;
      for( size_t i=0 ; i<undef_vars.size() ; ++i )
      {
         mesg << "       - \"" << undef_vars(i) << "\"" << endl ;
      }
   }
   if( predef_vars.size() > 0 )
   {
      mesg << "    you can used predefined variable(s): " << endl ;
      for( size_t i=0 ; i<predef_vars.size() ; ++i )
      {
         mesg << "       - \"" << predef_vars(i) << "\"" << endl ;
      }
   }
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
GE_TransformedMeshing_ERROR:: n1( size_t nb_sp_dims )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** GE_TransformedMeshing error:" << endl
        << "    the data of keyword \"transformation\"" << endl
        << "    should have " << nb_sp_dims << " components since" << endl
        << "    the number of space dimensions is " << nb_sp_dims ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}


