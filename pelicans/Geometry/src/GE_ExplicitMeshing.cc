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

#include <GE_ExplicitMeshing.hh>

#include <PEL_Error.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Data.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>
#include <intVector.hh>
#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

#include <sstream>

using std::string ;

GE_ExplicitMeshing const* 
GE_ExplicitMeshing::PROTOTYPE = new GE_ExplicitMeshing() ;

//-------------------------------------------------------------------------
GE_ExplicitMeshing:: GE_ExplicitMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_ExplicitMeshing" )
   , EXP( 0 )
   , MESHES_EXP( 0 )
   , VERT_COORDS( doubleVector( 0 ) )
   , VERT_COLORS( 0 )
   , IVERT( PEL::bad_index() )
   , A_VERT_COORDS( 0 )
   , MESH_COLORS( 0 )
   , MESH_POLY_NAME( "" )
   , CELL_IT( false )
   , FACE_IT( false )
   , I_MESH( PEL::bad_index() )
   , I_GLOB_MESH( PEL::bad_index() )
   , MESH_IDX( 0 )
   , MESH_2_VERTS( 0 )
   , HAS_MESH_2_FACES( false )
   , MESH_2_FACES( 0 )
   , VALID_MESH( false )
   , DECAL_MESH( PEL::bad_index() )
{
}

//-------------------------------------------------------------------------
GE_ExplicitMeshing*
GE_ExplicitMeshing:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp,
                                     size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_ExplicitMeshing* result =
                       new GE_ExplicitMeshing( a_owner, exp, dim_space ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_ExplicitMeshing:: GE_ExplicitMeshing( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp,
                                         size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , EXP( exp->create_clone( this ) )
   , MESHES_EXP( 0 )
   , VERT_COORDS( exp->doubleVector_data( "vertices/coordinates" ) )
   , VERT_COLORS( PEL_Vector::create( this, 0 ) )
   , IVERT( PEL::bad_index() )
   , A_VERT_COORDS( dim_space )
   , MESH_COLORS( PEL_Vector::create( this, 0 ) )
   , MESH_POLY_NAME( "" )
   , CELL_IT( false )
   , FACE_IT( false )
   , I_MESH( PEL::bad_index() )
   , I_GLOB_MESH( PEL::bad_index() )
   , MESH_IDX( 0 )
   , MESH_2_VERTS( 0 )
   , HAS_MESH_2_FACES( false )
   , MESH_2_FACES( 0 )
   , VALID_MESH( false )
   , DECAL_MESH( PEL::bad_index() )
 
{
   PEL_LABEL( "GE_ExplicitMeshing:: GE_ExplicitMeshing" ) ;

   // Vertices checks:
   if( VERT_COORDS.size() != dim_space*nb_vertices() )
   {
      PEL_Error::object()->raise_data_error(
         exp, "vertices",
         "   the size of the table of coordinates\n"
         "   and the number of vertices are mismatched" ) ;
   }

   // Cell checks:
   size_t nb_elms = 0 ;
   for( start_cell_iterator() ; valid_cell() ; go_next_cell() )
   {
      ++nb_elms ;
   }
   if( nb_elms != nb_cells() )
   {
      PEL_Error::object()->raise_data_error(
         exp, "cells",
         "   the size of the table of connectivities\n"
         "   and the number of cells are mismatched" ) ;
   }
   
   
   // Face checks:
   nb_elms = 0 ;
   size_t nb_f_verts_max = 0 ;
   for( start_face_iterator() ; valid_face() ; go_next_face(), ++nb_elms )
   {
      nb_f_verts_max = PEL::max( nb_f_verts_max, face_vertices().size() ) ;
   }
   if( nb_elms != nb_faces() )
   {
      PEL_Error::object()->raise_data_error(
         exp, "faces",
         "   the size of the table of connectivities\n"
         "   and the number of cells are mismatched" ) ;
   }

   // Face cells numbering:
   size_t_vector nb_f_verts( nb_faces(), 0 ) ;
   size_t_array2D f_verts( nb_faces(), nb_f_verts_max, PEL::bad_index() ) ;
   nb_elms = 0 ;
   for( start_face_iterator() ; valid_face() ; go_next_face(), ++nb_elms )
   {
      size_t_vector const& face_v = face_vertices() ;
      nb_f_verts( nb_elms ) = face_v.size() ;
      for( size_t i=0 ; i<face_v.size() ; ++i )
      {
         f_verts( nb_elms, i ) = face_v(i) ;
      }
   }
   
   size_t_vector ref_poly_faces(0) ;
   size_t_vector poly_faces(0) ;
   nb_elms = 0 ;
   for( start_cell_iterator() ; valid_cell() ; go_next_cell(), ++nb_elms )
   {
      GE_ReferencePolyhedron const* ref_poly = cell_reference_polyhedron() ;
      size_t_vector const& cell_v = cell_vertices() ;
      size_t_vector const& cell_f = cell_faces() ;
      if( ref_poly->nb_vertices() != cell_v.size() ||
          ref_poly->nb_faces() != cell_f.size() )
      {
         PEL_Error::object()->raise_data_error(
            exp, "cells",
            "   bad size for connectivity table" ) ;
      }
      for( size_t i=0 ; i<ref_poly->nb_faces() ; ++i )
      {
         // Reference polyhedron face:
         ref_poly_faces.resize( ref_poly->nb_face_vertices( i ) ) ;
         for( size_t j=0 ; j<ref_poly->nb_face_vertices( i ) ; ++j )
         {
            ref_poly_faces( j ) = cell_v( ref_poly->face_vertex( i, j ) ) ;
         }
         ref_poly_faces.sort_increasingly() ;

         // Polyhedron face:
         size_t ii = cell_f( i ) ;
         poly_faces.resize( nb_f_verts( ii ) ) ;
         for( size_t j=0 ; j<nb_f_verts( ii ) ; ++j )
         {
            poly_faces( j ) = f_verts( ii, j ) ;
         }
         poly_faces.sort_increasingly() ;
         if( poly_faces != ref_poly_faces )
         {
            std::ostringstream msg ;
            msg << "   the face connectivity does not respect numbering"
                << std::endl
                << "   imposed by reference polyhedron" << std::endl
                << "      cell: " << nb_elms
                << " (" << cell_polyhedron_name() << ")" << std::endl
                << "      face: " << i << std::endl
                << "      face connectivity: " << poly_faces << std::endl
                << "      expected connectivity: " << ref_poly_faces ;
            PEL_Error::object()->raise_data_error( exp, "cells", msg.str() ) ;
         }
      }
   }
}

//------------------------------------------------------------------------
GE_ExplicitMeshing:: ~GE_ExplicitMeshing( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: ~GE_ExplicitMeshing" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( is_a_prototype() ) PROTOTYPE = 0 ;
}

//------------------------------------------------------------------------
size_t
GE_ExplicitMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: nb_vertices" ) ;
   size_t result = EXP->int_data( "vertices/number" ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
GE_ExplicitMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: nb_cells" ) ;
   size_t result = EXP->int_data( "cells/number" ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
GE_ExplicitMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: nb_faces" ) ;
   size_t result = EXP->int_data( "faces/number" ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   IVERT = 0 ;

   // Recover vertices color:
   PEL_ModuleExplorer* color_identification = 0 ;
   if( EXP->has_module( "vertices/colors" ) )
   {
      color_identification =
                 EXP->create_subexplorer( 0, "vertices/colors" ) ;
   }
   PEL_ModuleExplorer* special_color = 0 ;
   if( EXP->has_module( "vertices/special_colors" ) )
   {
      special_color =
         EXP->create_subexplorer( 0, "vertices/special_colors" ) ;
   }
   recover_color_table( color_identification, special_color,
                        nb_vertices(),
                        VERT_COLORS ) ;
   if( color_identification!=0 )
   {
      color_identification->destroy() ; color_identification = 0 ;
   }
   if( special_color!=0 )
   {
      special_color->destroy() ; special_color = 0 ;
   }
}

//------------------------------------------------------------------------
bool
GE_ExplicitMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( IVERT<nb_vertices() ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   IVERT++ ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_ExplicitMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   size_t nbc = nb_space_dimensions() ;
   for( size_t ic=0 ; ic<nbc ; ic++ )
   {
      A_VERT_COORDS( ic ) = VERT_COORDS( IVERT*nbc+ic )  ;
   }
   doubleVector const& result = A_VERT_COORDS ;

   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   CELL_IT = true ;
   FACE_IT = false ;

   if( MESHES_EXP!=0 )
   {
      destroy_possession( MESHES_EXP ) ;
   }
   MESHES_EXP = EXP->create_subexplorer(
                           this, "cells/polyhedra_and_connectivities" ) ;
   HAS_MESH_2_FACES = true ;
   
   PEL_ModuleExplorer* color_identification = 0 ;
   if( EXP->has_module( "cells/colors" ) )
   {
      color_identification =
                      EXP->create_subexplorer( 0, "cells/colors" ) ;
   }
   PEL_ModuleExplorer* special_color = 0 ;
   if( EXP->has_module( "cells/special_colors" ) )
   {
      special_color =
              EXP->create_subexplorer( 0, "cells/special_colors" ) ;
   }
   recover_color_table( color_identification, special_color, nb_cells(),
                        MESH_COLORS ) ;

   MESHES_EXP->start_entry_iterator() ;
   PEL_ASSERT( MESHES_EXP->is_valid_entry() ) ;
   
   take_new_mesh_partition() ;

   if( color_identification!=0 )
   {
      color_identification->destroy() ; color_identification = 0 ;
   }
   if( special_color!=0 )
   {
      special_color->destroy() ; special_color = 0 ;
   }

   I_GLOB_MESH = 0 ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_ExplicitMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: valid_cell" ) ;

   bool result = CELL_IT && VALID_MESH ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;
   I_MESH++ ;
   I_GLOB_MESH++ ;
   size_t index_mesh = DECAL_MESH*I_MESH ;
   if( index_mesh >= MESH_IDX->size() )
   {
      MESHES_EXP->go_next_entry() ;
      take_new_mesh_partition() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_ExplicitMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = MESH_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExplicitMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   for( size_t i=0 ; i<MESH_2_VERTS.size() ; i++ )
   {
      int ivert = (*MESH_IDX)( I_MESH*DECAL_MESH+i )  ;
      if( ivert<0 || ivert>=(int)nb_vertices() )
      {
         PEL_Error::object()->raise_plain(
                                      "Bad index of vertices for cell" ) ;
      }
      MESH_2_VERTS( i ) = (size_t) ivert ;
   }
   size_t_vector const& result = MESH_2_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExplicitMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t decal = I_MESH*DECAL_MESH + MESH_2_VERTS.size() ;
   for( size_t i=0 ; i<MESH_2_FACES.size() ; i++ )
   {
      int iside = (*MESH_IDX)( decal+i ) ;
      if( iside<0 || iside>=(int)nb_faces() )
      {
         PEL_Error::object()->raise_plain(
                                      "Bad index of faces for cell" ) ;
      }
      MESH_2_FACES( i ) = (size_t) iside  ;
   }
   size_t_vector const& result = MESH_2_FACES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   CELL_IT = false ;
   FACE_IT = true ;

   if( MESHES_EXP!=0 )
   {
      destroy_possession( MESHES_EXP ) ;
   }
   MESHES_EXP = EXP->create_subexplorer(
                             this, "faces/polyhedra_and_connectivities" ) ;
   HAS_MESH_2_FACES = false ;

   PEL_ModuleExplorer* color_identification = 0 ;
   if( EXP->has_module( "faces/colors" ) )
   {
      color_identification =
                        EXP->create_subexplorer( 0, "faces/colors" ) ;
   }
   PEL_ModuleExplorer* special_color = 0 ;
   if( EXP->has_module( "faces/special_colors" ) )
   {
      special_color =
                EXP->create_subexplorer( 0, "faces/special_colors" ) ;
   }
   recover_color_table( color_identification, special_color, nb_faces(),
                        MESH_COLORS ) ;
   if( color_identification!=0 )
   {
      color_identification->destroy() ; color_identification = 0 ;
   }
   if( special_color!=0 )
   {
      special_color->destroy() ; special_color = 0 ;
   }

   if( EXP->has_module( "bounds") )
   {
      intVector const& global_sides =
                    EXP->intVector_data( "bounds/adjacent_faces" ) ;
      size_t nb_bds = global_sides.size() ;
      
      PEL_ModuleExplorer* bd_color_identification =
                     EXP->create_subexplorer( 0, "bounds/colors" ) ;
      PEL_Vector* bd_colors = PEL_Vector::create( 0, 0 ) ;
      PEL_ModuleExplorer* bd_special_color = 0 ;
      recover_color_table( bd_color_identification, bd_special_color,
                           nb_bds,
                           bd_colors ) ;
      
      for( size_t i_bd=0 ; i_bd<nb_bds ; ++i_bd )
      {
         if( global_sides(i_bd)<0 || global_sides(i_bd)>=(int)nb_faces() )
         {
            PEL_Error::object()->raise_plain( "Bad index of bounds" ) ;
         }
         size_t i_side = (size_t) global_sides(i_bd) ;
         GE_Color* col = static_cast<GE_Color*>( bd_colors->at(i_bd) ) ;
         if( col!=0 ) MESH_COLORS->set_at( i_side, col ) ;
      }
      bd_colors->destroy() ; bd_colors = 0 ;
      bd_color_identification->destroy() ; bd_color_identification = 0 ;
   }
   MESHES_EXP->start_entry_iterator() ;
   PEL_ASSERT( MESHES_EXP->is_valid_entry() ) ;
   
   take_new_mesh_partition() ;

   I_GLOB_MESH = 0 ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_ExplicitMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: valid_face" ) ;

   bool result = FACE_IT && VALID_MESH ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;
   I_MESH++ ;
   I_GLOB_MESH++ ;
   size_t index_mesh = DECAL_MESH*I_MESH ;
   if( index_mesh >= MESH_IDX->size() )
   {
      MESHES_EXP->go_next_entry() ;
      take_new_mesh_partition() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_ExplicitMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = MESH_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExplicitMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   for( size_t i=0 ; i<MESH_2_VERTS.size() ; i++ )
   {
      int ivert = (*MESH_IDX)( I_MESH*DECAL_MESH+i ) ;
      if( ivert<0 || ivert>=(int)nb_vertices() )
      {
         PEL_Error::object()->raise_plain(
                                      "Bad index of vertices for face" ) ;
      }
      MESH_2_VERTS( i ) = (size_t) ivert  ;
   }
   size_t_vector const& result = MESH_2_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExplicitMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = GE_Color::null_color() ;
   if( VERT_COLORS->index_limit()>IVERT )
   {
      GE_Color const* in_list_color =
         static_cast<GE_Color const* >( VERT_COLORS->at( IVERT ) ) ;
      if( in_list_color!=0 )
      {
         result = in_list_color ;
      }
   }

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExplicitMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = GE_Color::null_color() ;
   if( MESH_COLORS->index_limit()>I_MESH )
   {
      GE_Color const* in_list_color =
         static_cast<GE_Color const* >( MESH_COLORS->at( I_GLOB_MESH ) ) ;
      if( in_list_color!=0 )
      {
         result = in_list_color ;
      }
   }

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExplicitMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = GE_Color::null_color() ;
   if( MESH_COLORS->index_limit()>I_MESH )
   {
      GE_Color const* in_list_color =
         static_cast<GE_Color const* >( MESH_COLORS->at( I_GLOB_MESH ) ) ;
      if( in_list_color!=0 )
      {
         result = in_list_color ;
      }
   }

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: recover_color_table(
                           PEL_ModuleExplorer* color_identification,
                           PEL_ModuleExplorer* special_color,
                           size_t nb_elements,
                           PEL_Vector* color_table ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: recover_color_table" ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   color_table->re_initialize( nb_elements ) ;

   // Colors :
   if( color_identification!=0 )
   {
      for( color_identification->start_entry_iterator() ;
           color_identification->is_valid_entry() ;
           color_identification->go_next_entry() )
      {
         string const& id_color = color_identification->keyword() ;
         GE_Color::extend( id_color ) ;
         GE_Color const* color = GE_Color::object( id_color ) ;
         intVector const& color2elem =
            color_identification->intVector_data( id_color ) ;
         for( size_t i=0 ; i<color2elem.size() ; i++ )
         {
            if( color2elem(i)<0 || color2elem(i)>=(int)nb_elements )
            {
               PEL_Error::object()->raise_plain(
                  "Bad index of a geometrical element on color \""
                  +color->name()+"\"" ) ;
            }
            size_t id_elem = (size_t) color2elem(i) ;
            PEL_Object const* old_color = color_table->at( id_elem ) ;
            if( old_color!=0 )
            {
               GE_Color const* c = static_cast<GE_Color const*>( old_color ) ;
               PEL_Error::object()->raise_plain(
                  "\""+c->name()+"\" and \""+color->name()
                              +"\" colorize the same geometrical element." ) ;
            }
            color_table->set_at( id_elem, const_cast<GE_Color*>( color ) ) ;
         }
      }
   }

   // Halo elements :
   if( special_color!=0 )
   {
      if( special_color->has_entry( "null" ) )
      {
         GE_Color const* color = GE_Color::null_color() ;
         intVector const& color2elem =
                                    special_color->intVector_data( "null" ) ;
         for( size_t i=0 ; i<color2elem.size() ; i++ )
         {
            if( color2elem(i)<0 || color2elem(i)>=(int)nb_elements )
            {
               PEL_Error::object()->raise_plain(
                  "Bad index of a geometrical element on color \""
                  +color->name()+"\"" ) ;
            }
            size_t id_elem = (size_t) color2elem(i) ;
            PEL_Object const* old_color = color_table->at( id_elem ) ;
            if( old_color!=0 )
            {
               GE_Color const* c = static_cast<GE_Color const*>( old_color ) ;
               PEL_Error::object()->raise_plain(
                  "\""+c->name()+"\" and \""+color->name()
                              +"\" colorize the same geometrical element." ) ;
            }
            color_table->set_at( id_elem, const_cast<GE_Color*>( color ) ) ;
         }
      }
      if( special_color->has_entry( "halo" ) )
      {
         GE_Color const* color = GE_Color::halo_color() ;
         intVector const& color2elem =
                                   special_color->intVector_data( "halo" ) ;
         for( size_t i=0 ; i<color2elem.size() ; i++ )
         {
            if( color2elem(i)<0 || color2elem(i)>=(int)nb_elements )
            {
               PEL_Error::object()->raise_plain(
                  "Bad index of a geometrical element on color \""
                  +color->name()+"\"" ) ;
            }
            size_t id_elem = (size_t) color2elem(i) ;
            PEL_Object const* old_color = color_table->at( id_elem ) ;
            if( old_color!=0 )
            {
               GE_Color const* c = static_cast<GE_Color const*>( old_color ) ;
               PEL_Error::object()->raise_plain(
                  "\""+c->name()+"\" and \""+color->name()
                              +"\" colorize the same geometrical element." ) ;
            }
            color_table->set_at( id_elem, const_cast<GE_Color*>( color ) ) ;
         }
      }
   }
}

//------------------------------------------------------------------------
void
GE_ExplicitMeshing:: take_new_mesh_partition( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExplicitMeshing:: take_new_mesh_partition" ) ;
   
   PEL_CHECK( MESHES_EXP != 0 ) ;
   VALID_MESH = MESHES_EXP->is_valid_entry() ;
   
   if( VALID_MESH )
   {
      MESH_POLY_NAME = MESHES_EXP->keyword() ;
      GE_ReferencePolyhedron const* ref_poly =
                    GE_Mpolyhedron::reference_polyhedron( MESH_POLY_NAME ) ;
      MESH_IDX = &MESHES_EXP->intVector_data( MESH_POLY_NAME ) ;
      I_MESH = 0 ;
      size_t nb_vert = ref_poly->nb_vertices() ;
      MESH_2_VERTS.re_initialize( nb_vert ) ;
      DECAL_MESH = nb_vert ;
      if( HAS_MESH_2_FACES )
      {
         size_t nb_conn = ref_poly->nb_faces() ;
         MESH_2_FACES.re_initialize( nb_conn ) ;
         DECAL_MESH += nb_conn ;
      }
      if( MESH_IDX->size()%DECAL_MESH != 0 )
      {
         PEL_Error::object()->raise_data_error(
            MESHES_EXP, MESH_POLY_NAME,
            "   bad size of the table of connectivity" ) ;
      }
   }
   else
   {
      MESHES_EXP = 0 ;
   }
}

//------------------------------------------------------------------------
bool
GE_ExplicitMeshing:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;
   PEL_ASSERT( IMPLIES( valid_cell() || valid_face(), MESH_IDX!=0 ) ) ;
   return( true ) ;
}


   
