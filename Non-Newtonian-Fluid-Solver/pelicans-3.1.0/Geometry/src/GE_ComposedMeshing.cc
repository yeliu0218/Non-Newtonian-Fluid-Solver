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

#include <GE_ComposedMeshing.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <iostream>
#include <sstream>
#include <string>

GE_ComposedMeshing const* 
GE_ComposedMeshing:: PROTOTYPE = new GE_ComposedMeshing() ;

//-------------------------------------------------------------------------
GE_ComposedMeshing:: GE_ComposedMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_ComposedMeshing" )
   , MESHINGS( 0 )
   , MESHING_COLS( 0 )
   , VERTICES( 0 )
   , VERT_LOC_TO_GLOB( 0, 0 )
   , VERT_INDEX( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , CELL_MESHING( 0 )
   , CELL_MESHING_INDEX( PEL::bad_index() )
   , NB_SIDES( PEL::bad_index() )
   , SIDE_LOC_TO_GLOB( 0, 0 )
   , SIDE_OK( 0 )
   , SIDE_COLORS( 0 )
   , SIDE_MESHING( 0 )
   , SIDE_MESHING_INDEX( PEL::bad_index() )
   , SIDE_INDEX( PEL::bad_index() )
{
}

//-------------------------------------------------------------------------
GE_ComposedMeshing*
GE_ComposedMeshing:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp,
                                     size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_ComposedMeshing* result = new GE_ComposedMeshing(  a_owner, exp, dim_space ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_ComposedMeshing:: GE_ComposedMeshing( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp,
                                         size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , MESHINGS( PEL_Vector::create( this, 0 ) )
   , MESHING_COLS( 0 )
   , VERTICES( GE_SetOfPoints::create( this, dim_space ) )
   , VERT_LOC_TO_GLOB( 0, 0 )
   , VERT_INDEX( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , CELL_MESHING( 0 )
   , CELL_MESHING_INDEX( PEL::bad_index() )
   , NB_SIDES( PEL::bad_index() )
   , SIDE_LOC_TO_GLOB( 0, 0 )
   , SIDE_OK( 0 )
   , SIDE_COLORS( PEL_Vector::create( this, 0 ) )
   , SIDE_MESHING( 0 )
   , SIDE_MESHING_INDEX( PEL::bad_index() )
   , SIDE_INDEX( PEL::bad_index() )
{
   PEL_LABEL( "GE_ComposedMeshing:: GE_ComposedMeshing:: GE_ComposedMeshing" ) ;   
   // Meshings :
   {
      PEL_ModuleExplorer* se = 
                          exp->create_subexplorer( 0, "list_of_GE_Meshing" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         PEL_ModuleExplorer* me = se->create_subexplorer( 0 ) ;
         GE_Meshing* meshing = GE_Meshing::create( MESHINGS, me, dim_space ) ;
         MESHINGS->append( meshing ) ;
         me->destroy() ; me = 0 ;
      }
      se->destroy() ; se = 0 ;
   }

   if( exp->has_entry( "meshing_names" ) )
   {
      MESHING_COLS = exp->stringVector_data( "meshing_names" ) ;
      if( MESHING_COLS.size() != MESHINGS->index_limit() )
      {
         PEL_Error::object()->raise_data_error(
                            exp, "meshing_names", "bad size" ) ;
      }
   }
   else
   {
      MESHING_COLS.re_initialize( MESHINGS->index_limit() ) ;
      for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
      {
         std::ostringstream os ;
         os << "m" << i ;
         MESHING_COLS(i) = os.str() ;
      }
   }
   
   // Vertices :
   build_vertices() ;

   // Cells :
   build_cells() ;

   // Sides :
   build_sides() ;
   
   // Check :
   // Check_meshing in o(nb_vertices()*nb_cells)
   if( exp->bool_data( "check_meshing" ) )
   {
      check_meshing() ;
   }
}

//------------------------------------------------------------------------
GE_ComposedMeshing:: ~GE_ComposedMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_ComposedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
  
   return( VERTICES->nb_points()) ;
}

//------------------------------------------------------------------------
size_t
GE_ComposedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: nb_cells" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( NB_CELLS ) ;
}

//------------------------------------------------------------------------
size_t
GE_ComposedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: nb_faces" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( NB_SIDES ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   VERT_INDEX = 0 ;
}

//------------------------------------------------------------------------
bool
GE_ComposedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( VERT_INDEX<nb_vertices() ) ;   
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   VERT_INDEX++ ;   
}

//------------------------------------------------------------------------
doubleVector const& 
GE_ComposedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   GE_Point const* pt = VERTICES->point( VERT_INDEX ) ;
   static doubleVector result(0) ;
   result.re_initialize( nb_space_dimensions() ) ;
   for( size_t i=0 ; i<nb_space_dimensions() ; ++i )
   {
      result( i ) =  pt->coordinate( i ) ;
   }
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   if( MESHINGS->index_limit() > 0 )
   {
      CELL_MESHING_INDEX = 0 ;
      CELL_MESHING =
         static_cast<GE_Meshing*>( MESHINGS->at( 0 ) ) ;
      CELL_MESHING->start_cell_iterator() ;
   }
   
   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   CELL_MESHING->go_next_cell() ;
   while( CELL_MESHING!=0 &&
          !CELL_MESHING->valid_cell() )
   {
      ++CELL_MESHING_INDEX ;
      if( CELL_MESHING_INDEX<MESHINGS->index_limit() )
      {
         CELL_MESHING =
            static_cast<GE_Meshing*>( MESHINGS->at( CELL_MESHING_INDEX ) ) ;
      CELL_MESHING->start_cell_iterator() ;
      }
      else
      {
         CELL_MESHING = 0 ;
      }
   }
   if(CELL_MESHING==0)
   {
      CELL_MESHING_INDEX = PEL::bad_index() ;
   }
}

//------------------------------------------------------------------------
bool
GE_ComposedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: valid_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( CELL_MESHING_INDEX!=PEL::bad_index() ) ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_ComposedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_MESHING->cell_polyhedron_name() ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ComposedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& c_verts = CELL_MESHING->cell_vertices() ;
   static size_t_vector result(0) ;
   result.re_initialize( c_verts.size() ) ;
   for( size_t i=0 ; i<c_verts.size() ; ++i )
   {
      result(i) = VERT_LOC_TO_GLOB( CELL_MESHING_INDEX, c_verts(i) ) ;
   }

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result );
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ComposedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& c_sides = CELL_MESHING->cell_faces() ;
   static size_t_vector result(0) ;
   result.re_initialize( c_sides.size() ) ;
   for( size_t i=0 ; i<c_sides.size() ; ++i )
   {
      result(i) = SIDE_LOC_TO_GLOB( CELL_MESHING_INDEX, c_sides(i) ) ;
   }

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   if( MESHINGS->index_limit() > 0 )
   {
      SIDE_MESHING_INDEX = 0 ;
      SIDE_INDEX = 0 ;
      SIDE_MESHING =
         static_cast<GE_Meshing*>( MESHINGS->at( SIDE_MESHING_INDEX ) ) ;
      SIDE_MESHING->start_face_iterator() ;
   }
   
   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   SIDE_MESHING->go_next_face() ;
   ++SIDE_INDEX ;
   while( SIDE_MESHING->valid_face() && !SIDE_OK(SIDE_INDEX) )
   {
      SIDE_MESHING->go_next_face() ;
      ++SIDE_INDEX ;
   }
   if( !SIDE_MESHING->valid_face() )
   {
      --SIDE_INDEX ;
      ++SIDE_MESHING_INDEX ;
      if( SIDE_MESHING_INDEX<MESHINGS->index_limit() )
      {
         SIDE_MESHING =
            static_cast<GE_Meshing*>( MESHINGS->at( SIDE_MESHING_INDEX ) ) ;
         SIDE_MESHING->start_face_iterator() ;
         ++SIDE_INDEX ;
         while( SIDE_MESHING->valid_face() && !SIDE_OK(SIDE_INDEX) )
         {
            SIDE_MESHING->go_next_face() ;
            ++SIDE_INDEX ;
         }
         PEL_ASSERT( SIDE_OK(SIDE_INDEX) ) ;
      }
      else
      {
         SIDE_MESHING_INDEX = PEL::bad_index() ;
         SIDE_INDEX = PEL::bad_index() ;
         SIDE_MESHING = 0 ;
      }
   }
}

//------------------------------------------------------------------------
bool
GE_ComposedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: valid_face" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = ( SIDE_MESHING_INDEX!=PEL::bad_index() ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_ComposedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = SIDE_MESHING->face_polyhedron_name() ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ComposedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& s_verts = SIDE_MESHING->face_vertices() ;
   static size_t_vector result(0) ;
   result.re_initialize( s_verts.size() ) ;
   for( size_t i=0 ; i<s_verts.size() ; ++i )
   {
      result(i) = VERT_LOC_TO_GLOB( SIDE_MESHING_INDEX, s_verts(i) ) ;
   }

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing const* m = static_cast<GE_Meshing const*>( MESHINGS->at(i) ) ;
      m->print( os, indent_width+3 ) ;
   }
}

//------------------------------------------------------------------------
GE_Color const*
GE_ComposedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Color const* result  = VERTICES->color( VERT_INDEX ) ;

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ComposedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const m_col = MESHING_COLS(CELL_MESHING_INDEX) ;
   std::string const c_col = m_col+"_"+CELL_MESHING->cell_color()->name() ;
   GE_Color::extend( c_col ) ;
   GE_Color const* result  = GE_Color::object( c_col ) ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ComposedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = static_cast<GE_Color const*>(
      SIDE_COLORS->at( SIDE_INDEX ) ) ;
   
   if( result == 0 ) result = GE_Color::null_color() ;
   
   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: build_vertices( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: build_vertices" ) ;

   size_t nb_max = 0 ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing const* m = static_cast<GE_Meshing const*>( MESHINGS->at(i) ) ;
      nb_max = PEL::max( nb_max, m->nb_vertices() ) ;
   }
   VERT_LOC_TO_GLOB.re_initialize( MESHINGS->index_limit(), nb_max ) ;
   VERT_LOC_TO_GLOB.set( PEL::bad_index() ) ;

   GE_Point* dummy_pt = GE_Point::create( 0, nb_space_dimensions() ) ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing* m = static_cast<GE_Meshing*>( MESHINGS->at(i) ) ;
      std::string const meshing_col = MESHING_COLS(i)+"_" ;

      // Build vertices :
      size_t j = 0 ;
      for( m->start_vertex_iterator() ; m->valid_vertex() ; ++j, m->go_next_vertex() )
      {
         dummy_pt->set_coordinates( m->vertex_coordinates() ) ;
         std::string pt_color ;
         if( m->vertex_color()!=GE_Color::null_color() )
         {
            pt_color = meshing_col+m->vertex_color()->name() ;
         }
         size_t index = PEL::bad_index() ;
         if( !VERTICES->has( dummy_pt ) )
         {
            index = VERTICES->nb_points() ;
            VERTICES->append( dummy_pt ) ;
         }
         else
         {
            index = VERTICES->index( dummy_pt ) ;
            GE_Color const* old_col = VERTICES->color( index ) ;
            if( old_col!=GE_Color::null_color() )
            {
               if( !pt_color.empty() )
               {
                  pt_color = old_col->name()+"_"+pt_color ;
               }
               else
               {
                  pt_color = old_col->name() ;
               }
            }
         }
         VERT_LOC_TO_GLOB(i,j) = index ;
         if( !pt_color.empty() )
         {
            GE_Color::extend( pt_color ) ;
            VERTICES->modify_color( index, GE_Color::object( pt_color ) ) ;
         }
      }
   }
   dummy_pt->destroy() ; dummy_pt = 0 ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: build_cells( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: build_cells" ) ;

   NB_CELLS = 0 ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing const* m = static_cast<GE_Meshing const*>( MESHINGS->at(i) ) ;
      NB_CELLS += m->nb_cells() ;
   }
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: build_sides( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: build_sides" ) ;

   size_t nb_max = 0 ;
   size_t nb_sid = 0 ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing* m = static_cast<GE_Meshing *>( MESHINGS->at(i) ) ;
      nb_max = PEL::max( nb_max, m->nb_faces() ) ;
      nb_sid += m->nb_faces() ;
   }
   
   SIDE_LOC_TO_GLOB.re_initialize( MESHINGS->index_limit(), nb_max ) ;
   SIDE_LOC_TO_GLOB.set( PEL::bad_index() ) ;

   SIDE_OK.re_initialize( nb_sid ) ;
   SIDE_OK.set( false ) ;

   SIDE_COLORS->re_initialize( nb_sid ) ;
   
   // Side data base :
   PEL_BalancedBinaryTree * side_centers =
      PEL_BalancedBinaryTree:: create( 0 ) ;

   size_t_vector bd_centers_to_glob( nb_sid ) ;
   size_t_vector bd_centers_to_side_glob( nb_sid ) ;
   
   size_t i_side = 0 ;
   NB_SIDES = 0 ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing* m = static_cast<GE_Meshing*>( MESHINGS->at(i) ) ;

      // Bounds :
      size_t_vector nb_adj_cells( m->nb_faces() ) ;
      {
         nb_adj_cells.set(0) ;
         for( m->start_cell_iterator() ; m->valid_cell() ; m->go_next_cell() )
         {
            size_t_vector const& c_sides = m->cell_faces() ;
            for( size_t is=0 ; is< c_sides.size() ; ++is )
            {
               ++nb_adj_cells( c_sides(is) ) ;
            }
         }
      }

      // Side colors :
      {
         std::string const meshing_col = MESHING_COLS(i)+"_" ;
         size_t side_index = i_side ;
         m->start_face_iterator() ;
         for( ; m->valid_face() ; ++side_index, m->go_next_face() )
         {
            if( m->face_color()!=GE_Color::null_color() )
            {
               std::string const ncol = meshing_col+m->face_color()->name() ;
               GE_Color::extend( ncol ) ;
               SIDE_COLORS->set_at(side_index,
                                   const_cast<GE_Color *>(GE_Color::object(ncol)) ) ;
            }
         }
      }

      // Sides :
      size_t i_loc = 0 ;
      m->start_face_iterator() ;
      for( ; m->valid_face() ; ++i_loc, ++i_side, m->go_next_face() )
      {
         bool new_side = true ;
         size_t i_glob = NB_SIDES ;
         size_t i_side_glob = i_side ;
         
         if( nb_adj_cells(i_loc)==1 ) // Bound
         {
            
            size_t_vector const& s_verts = side_global_vertices( i ) ;
            PEL_IndexSet* idx = PEL_IndexSet::create( 0,
                                                      s_verts,
                                                      side_centers->count() ) ;
            
            if( !side_centers->has( idx ) )
            {
               size_t ii = side_centers->count() ;
               bd_centers_to_glob( ii ) = i_glob ;
               bd_centers_to_side_glob( ii ) = i_side_glob ;
               side_centers->extend( idx ) ;
               idx->set_owner( side_centers ) ;
            }
            else
            {
               PEL_IndexSet const* l =
                  static_cast<PEL_IndexSet const*>( side_centers->item( idx ) ) ;
               size_t ii = l->id() ;
               idx->destroy() ;
               i_glob = bd_centers_to_glob( ii ) ;
               i_side_glob = bd_centers_to_side_glob( ii ) ;
               new_side = false ;
            }
         }

         if( new_side )
         {
            ++NB_SIDES ;
            SIDE_OK( i_side ) = true ;
         }
         else
         {
            GE_Color * s_old = static_cast<GE_Color *>( SIDE_COLORS->at( i_side_glob ) ) ;
            GE_Color * s_new = static_cast<GE_Color *>( SIDE_COLORS->at( i_side ) ) ;
            if( s_old!=0 || s_new!=0 )
            {
               if( s_new!=0 )
               {
                  if( s_old==0 )
                  {
                     SIDE_COLORS->set_at( i_side_glob, s_new ) ;
                  }
                  else
                  {
                     std::string new_color = s_old->name()+"_"+s_new->name() ;
                     GE_Color::extend( new_color ) ;
                     GE_Color* ge_new_color =
                        const_cast<GE_Color *>(GE_Color::object( new_color ))  ;
                     
                     SIDE_COLORS->set_at( i_side_glob, ge_new_color ) ;
                     SIDE_COLORS->set_at( i_side, ge_new_color ) ;
                  }
               }
               else
               {
                  SIDE_COLORS->set_at( i_side, s_old ) ;
               }
            }
         }
         SIDE_LOC_TO_GLOB( i, i_loc ) = i_glob ;
      }
   }
   side_centers->destroy() ; side_centers = 0 ;
}

//------------------------------------------------------------------------
size_t_vector const&
GE_ComposedMeshing:: side_global_vertices( size_t i_meshing ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: side_global_vertices" ) ;
   PEL_CHECK( static_cast<GE_Meshing*>( MESHINGS->at(i_meshing) )->valid_face() ) ;

   GE_Meshing* m = static_cast<GE_Meshing*>( MESHINGS->at(i_meshing) ) ;

   static size_t_vector result(0) ;
   result = m->face_vertices() ;
   
   // Local to global :
   for( size_t i=0 ; i<result.size() ; ++i )
   {
      result(i) = VERT_LOC_TO_GLOB( i_meshing, result(i) ) ;
   }
   
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_ComposedMeshing:: check_meshing( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ComposedMeshing:: check_meshing" ) ;

   bool check_consistency = GE_Mpolyhedron::check_consistency() ;
   if( !check_consistency ) GE_Mpolyhedron::set_check_consistency() ;

   // Build meshing :
   PEL_Vector* mesh_table = PEL_Vector::create( 0, 0 ) ;
   PEL_Vector* vect_pts = PEL_Vector::create( 0, 0 ) ;
   for( size_t i=0 ; i<MESHINGS->index_limit() ; ++i )
   {
      GE_Meshing* m = static_cast<GE_Meshing*>( MESHINGS->at(i) ) ;
      for( m->start_cell_iterator() ; m->valid_cell() ; m->go_next_cell() )
      {
         size_t_vector const& connect = m->cell_vertices() ;
         vect_pts->re_initialize(0) ;
         for( size_t k=0 ; k<connect.size() ; ++k )
         {
            GE_Point* pt = const_cast<GE_Point*>(
               VERTICES->point( VERT_LOC_TO_GLOB( i, connect(k) ) ) ) ;
            vect_pts->append( pt ) ;
         }
         GE_Mpolyhedron* cell =
            GE_Mpolyhedron::create( mesh_table,
                                    m->cell_polyhedron_name(),
                                    vect_pts ) ;
         mesh_table->append( cell ) ;
      }
   }
   vect_pts->destroy() ; vect_pts = 0 ;

   // Check that a vertices is a vertex of the cell or is not contains by
   // this later :
   for( size_t i=0 ; i<VERTICES->nb_points() ; ++i )
   {
      GE_Point const* pt = VERTICES->point(i) ;
      for( size_t im=0 ; im<mesh_table->index_limit() ; ++im )
      {
         GE_Mpolyhedron const* mesh =
              static_cast<GE_Mpolyhedron const*>( mesh_table->at(im) ) ;
         if( mesh->contains( pt ) )
         {
            bool a_vert = false ;
            for( size_t j=0 ; !a_vert && j<mesh->nb_vertices() ; ++j )
            {
               a_vert = ( pt==mesh->vertex(j) ) ;
            }
            if( !a_vert )
            {
               mesh_table->destroy() ; mesh_table = 0 ;
               PEL_Error::object()->raise_plain(
                  "GE_ComposedMeshing error : Overlapping meshings" ) ;
            }
         }
      }
   }
   mesh_table->destroy() ; mesh_table = 0 ;

   if( !check_consistency ) GE_Mpolyhedron::unset_check_consistency() ;
}
