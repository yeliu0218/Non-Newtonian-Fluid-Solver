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

#include <GE_CurveWithSegments.hh>
    
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Data.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <size_t_vector.hh>
    
#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
    
using std::string ;
    
//-------------------------------------------------------------------------
GE_CurveWithSegments const* 
GE_CurveWithSegments::prototype = new GE_CurveWithSegments() ;
//-------------------------------------------------------------------------
    
//-------------------------------------------------------------------------
GE_CurveWithSegments:: GE_CurveWithSegments( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_CurveWithSegments" )
   , NB_SP_DIMS( PEL::bad_index() )
   , BENDS( 0, 0 )
   , NB_SEG_VERT( 0 )
   , X_SEG_VERT( 0, 0 )
   , VERT_COORD( 0 )
   , mesh2verts( 0 )
   , mesh2sides( 0 )
{
} 
    
//-------------------------------------------------------------------------
GE_CurveWithSegments*
GE_CurveWithSegments:: create_replica( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp,
                                      size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_CurveWithSegments* result = new GE_CurveWithSegments( a_owner,
                                                    exp,
                                                    dim_space ) ;
    
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;
}
    
//-------------------------------------------------------------------------
GE_CurveWithSegments:: GE_CurveWithSegments( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp,
                                           size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , NB_SP_DIMS( dim_space )
   , CLOSED_CURVE( exp->bool_data( "closed_curve" ) )
   , BENDS( exp->doubleArray2D_data( "bends" ) )
   , NB_SEG_VERT( 0 )
   , X_SEG_VERT( 0, 0 )
   , CELL_IT( false )
   , SIDE_IT( false )
   , VERT_COORD( dim_space )
   , mesh2verts( 2 )
   , mesh2sides( 2 )
{
   if( BENDS.index_bound( 1 ) != dim_space )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_CurveWithSegments error :\n"
         "    \"bends\" : invalid size for the points" ) ;
   }
   
   NB_SEG = BENDS.index_bound(0)-1 ;
   if( CLOSED_CURVE ) ++NB_SEG ;
   NB_SEG_VERT.re_initialize( NB_SEG ) ;
    
   doubleVector const& sbd = exp->doubleVector_data( "subdivisions" ) ;
   size_t old = 0 ;
   i_seg = 0 ;
   size_t nb_seg_vert_max = 0 ;
   bool ok = ( sbd(0) == 0.0 ) ;
   NB_VERTICES = 1 ;
   for( size_t i=1 ; i<sbd.size() ; ++i )
   {
      if( sbd(i) < sbd(i-1) )
      {
         ok = ok && ( sbd(i)==0.0 && sbd(i-1)==1.0 ) ;
         NB_SEG_VERT( i_seg ) = (i - old) ;
    	 if( NB_SEG_VERT(i_seg) > nb_seg_vert_max ) 
    	 {
            nb_seg_vert_max = NB_SEG_VERT( i_seg ) ;
    	 }
         old = i ;
         ++i_seg ;
      }
      else if( i == (sbd.size()-1) )
      {
         ok = ok && ( sbd(i) == 1.0 ) ;
         NB_SEG_VERT( i_seg ) = i-old+1 ;
         if( NB_SEG_VERT(i_seg) > nb_seg_vert_max ) 
         {
            nb_seg_vert_max = NB_SEG_VERT( i_seg ) ;
         }
         ++NB_VERTICES ;
      }
      else
      {
         ++NB_VERTICES ;
      }
   }
   if( CLOSED_CURVE ) --NB_VERTICES ; // original and final vertices identical

   ok = ok && ( i_seg == (NB_SEG-1) ) ;
   if( !ok )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_CurveWithSegments error :\n"
         "    \"subdivisions\" : invalid argument" ) ;
   }
    
   X_SEG_VERT.re_initialize( NB_SEG, nb_seg_vert_max ) ;
   size_t ii=0 ;
   for( i_seg=0 ; i_seg<NB_SEG ; ++i_seg)
   {
      for( i_seg_vert=0 ; i_seg_vert<NB_SEG_VERT( i_seg ) ; ++i_seg_vert )
      {
         X_SEG_VERT( i_seg, i_seg_vert ) = sbd( ii ) ;
         ++ii ;
      }
   }
    
   read_mesh_polyhedron( exp, 1, SIDE_POLY_NAME, CELL_POLY_NAME ) ;
} 
    
//------------------------------------------------------------------------
GE_CurveWithSegments:: ~GE_CurveWithSegments( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
} 
    
//------------------------------------------------------------------------
size_t
GE_CurveWithSegments:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
    
   return( NB_VERTICES ) ;
}
    
//------------------------------------------------------------------------
size_t
GE_CurveWithSegments:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: nb_cells" ) ;
   size_t result = ( CLOSED_CURVE ? nb_vertices() : nb_vertices()-1 ) ;
   return( result ) ;
}    
    
//------------------------------------------------------------------------
size_t
GE_CurveWithSegments:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: nb_faces" ) ;
   size_t result = nb_vertices() ;
   return( result ) ;
}   
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   i_vert=0 ;
   i_seg_vert= 0 ;
   i_seg=0 ;
}   
    
//------------------------------------------------------------------------
bool
GE_CurveWithSegments:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( i_vert<nb_vertices() ) ;
}   
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   ++i_vert ;
   ++i_seg_vert ;
   if( i_seg_vert >= NB_SEG_VERT(i_seg) )
   {
      ++i_seg ;
      i_seg_vert=1 ;
   }
}
    
//------------------------------------------------------------------------
doubleVector const& 
GE_CurveWithSegments:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
    
   if( i_seg_vert == 0 )
   {
      for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
      {
         VERT_COORD(d) = BENDS( i_seg, d ) ;
      }
   }
   else if( i_seg_vert == ( NB_SEG_VERT(i_seg)-1 ) )
   {
      if( CLOSED_CURVE && ( i_seg==(NB_SEG-1) ) )
      {
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
            VERT_COORD(d) = BENDS( 0, d ) ;
      }
      else
      {
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
            VERT_COORD(d) = BENDS( i_seg+1, d ) ;
      }
   }
   else
   {
      double xx = X_SEG_VERT( i_seg, i_seg_vert ) ;
      if( CLOSED_CURVE && ( i_seg==(NB_SEG-1) ) )
      {
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
            VERT_COORD(d) = BENDS( i_seg, d ) + 
                            xx * ( BENDS( 0, d ) - BENDS( i_seg, d ) ) ;
      }
      else
      {
         for( size_t d=0 ; d<NB_SP_DIMS ; ++d )
            VERT_COORD(d) = BENDS( i_seg, d ) + 
                            xx * ( BENDS( i_seg+1, d ) - BENDS( i_seg, d ) ) ;
      }
   }
   doubleVector const& result = VERT_COORD ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
} 
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;
    
   CELL_IT = true ;
   SIDE_IT = false ;
   i_mesh = 0 ;
    
   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}  
    
//------------------------------------------------------------------------
bool
GE_CurveWithSegments:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: valid_cell" ) ;
    
   bool result = CELL_IT && ( i_mesh < nb_cells() ) ;
    
   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}  
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;
    
   ++i_mesh ;
} 
    
//------------------------------------------------------------------------
std::string const&
GE_CurveWithSegments:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;
    
   std::string const& result = CELL_POLY_NAME ;
    
   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}
    
//------------------------------------------------------------------------
size_t_vector const& 
GE_CurveWithSegments:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;
    
   mesh2verts.re_initialize(2) ;
   mesh2verts(0) = i_mesh ;
   if( CLOSED_CURVE && ( i_mesh==nb_cells()-1 ) )
   {
      mesh2verts(1) = 0 ;
   }
   else
   {
      mesh2verts(1) = i_mesh+1 ;
   }
    
   size_t_vector const& result = mesh2verts ;
    
   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}
    
//????? pareil que cell_vertices ????? les sides et les vertices ont le
//????? meme numero car sont confondus
//------------------------------------------------------------------------
size_t_vector const& 
GE_CurveWithSegments:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;
    
   mesh2sides(0) = i_mesh ;
   if( CLOSED_CURVE && ( i_mesh==nb_cells()-1 ) )
   {
      mesh2sides(1) = 0 ;
   }
   else
   {
      mesh2sides(1) = i_mesh+1 ;
   }
   size_t_vector const& result = mesh2sides ;
    
   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result ) ;
}
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;
    
   CELL_IT = false ;
   SIDE_IT = true ;
   i_mesh = 0 ;
    
   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}  
    
//------------------------------------------------------------------------
bool
GE_CurveWithSegments:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: valid_face" ) ;
    
   bool result = SIDE_IT && ( i_mesh < nb_faces() ) ;
    
   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}   
    
//------------------------------------------------------------------------
void
GE_CurveWithSegments:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;
    
   ++i_mesh ;
} 
    
//------------------------------------------------------------------------
std::string const&
GE_CurveWithSegments:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;
    
   std::string const& result = SIDE_POLY_NAME ;
    
   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}   
    
//------------------------------------------------------------------------
size_t_vector const& 
GE_CurveWithSegments:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;
    
   mesh2verts.re_initialize(1) ;
   mesh2verts(0) = i_mesh ;
    
   size_t_vector const& result = mesh2verts ;
    
   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}   
    
//------------------------------------------------------------------------
GE_Color const*
GE_CurveWithSegments:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
    
   GE_Color const* result =  GE_Color::null_color() ;
    
   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
} 
    
//------------------------------------------------------------------------
GE_Color const*
GE_CurveWithSegments:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;
    
   GE_Color const* result = GE_Color::null_color() ;
    
   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}
    
//------------------------------------------------------------------------
GE_Color const*
GE_CurveWithSegments:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;
    
   GE_Color const* result = GE_Color::null_color() ;
    
   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}
    
//-------------------------------------------------------------------------
void
GE_CurveWithSegments:: check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CurveWithSegments:: check_mesh_polyhedron" ) ;
   PEL_CHECK( check_mesh_polyhedron_PRE(
                            dim_space, face_poly_name, cell_poly_name ) ) ;

   // Check the cell reference polyhedron :
   GE_ReferencePolyhedron const* cell_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   if( cell_poly_ref->nb_vertices() != 2 )
   {
      raise_invalid_cell_polyhedron( cell_poly_name ) ;
   }

   // Check the side reference polyhedron :
   GE_ReferencePolyhedron const* side_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( face_poly_name ) ;
   if( side_poly_ref->nb_vertices() != 1 )
   {
      raise_invalid_face_polyhedron( face_poly_name ) ;
   }
}
    
    
    
