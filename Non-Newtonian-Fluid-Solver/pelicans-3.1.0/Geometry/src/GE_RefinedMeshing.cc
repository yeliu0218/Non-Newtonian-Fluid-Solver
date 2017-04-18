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

#include <GE_RefinedMeshing.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>
#include <GE_SetOfPoints.hh>

#include <iostream>
#include <sstream>
#include <string>

using std::cout ;
using std::endl ;

GE_RefinedMeshing const* 
GE_RefinedMeshing:: PROTOTYPE = new GE_RefinedMeshing() ;

struct GE_RefinedMeshing_ERROR
{
   static void n0( GE_ReferencePolyhedronRefiner const* r0,
                   GE_ReferencePolyhedronRefiner const* r1,
                   GE_Color const* col ) ;
   static void n1( GE_Color const* col ) ;
} ;

//-------------------------------------------------------------------------
GE_RefinedMeshing:: GE_RefinedMeshing( void )
//-------------------------------------------------------------------------
   : GE_Meshing( "GE_RefinedMeshing" )
   , MESH_VERTICES( 0 )
   , CELL_SIDES( 0 )
   , global_sides(0,0)
   , global_vertices(0,0)
   , nb_meshes( 0 )
   , vertex_coord( 0 )
   , REFs( 0 )
   , REF_COLs( 0 )
{
}

//-------------------------------------------------------------------------
GE_RefinedMeshing*
GE_RefinedMeshing:: create_replica( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp,
                                    size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   GE_RefinedMeshing* result = 
                      new GE_RefinedMeshing( a_owner, exp, dim_space ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;   
}

//-------------------------------------------------------------------------
GE_RefinedMeshing:: GE_RefinedMeshing( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp,
                                       size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , initial( 0 )
   , MESH_VERTICES( 0 )
   , CELL_SIDES( 0 )
   , CELL_IT( false )
   , SIDE_IT( false )
   , sop( GE_SetOfPoints::create( this, dim_space ) )
   , global_sides(0,0)
   , global_vertices(0,0)
   , initial_side_poly(0)
   , initial_side_color(0)
   , SIDE_COLORS(0)
   , nb_meshes( dim_space+1 )
   , idx_vertex( (size_t)~0 )
   , idx_mesh( (size_t)~0 )
   , vertex_coord( dim_space )
   , verbose( false )
   , REFs( 0 )
   , REF_COLs( 0 )
{
   read_mesh_polyhedron( exp, dim_space, side_poly_name, mesh_poly_name ) ;

   side_ref = GE_Mpolyhedron::reference_polyhedron( side_poly_name ) ;
   mesh_ref = GE_Mpolyhedron::reference_polyhedron( mesh_poly_name ) ;
   
   PEL_ModuleExplorer const* sexp = exp->create_subexplorer( 0, "GE_Meshing" ) ;
   initial = GE_Meshing::create( this, sexp, dim_space ) ;
   sexp->destroy() ; sexp = 0 ;
   
   initial_side_poly  = PEL_Vector::create( this, initial->nb_faces() ) ;
   initial_side_color = PEL_Vector::create( this, initial->nb_faces() ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, 
                                 "list_of_GE_ReferencePolyhedronRefiner" ) ;
   se->start_module_iterator() ;
   for( ; se->is_valid_module() ; se->go_next_module() )
   {
      PEL_ModuleExplorer const* sse = se->create_subexplorer( se ) ;
      GE_ReferencePolyhedronRefiner const* refi = 
                   GE_ReferencePolyhedronRefiner::make( this, sse ) ;
      GE_Color const* col = 0 ;
      if( sse->has_entry( "color" ) )
      {
         col = GE_Color::object( sse->string_data( "color" ) ) ;
      }
      append_refiner( col, refi ) ;
   }
   se->destroy() ;
   
   build_vertices() ;
   build_sides() ;
   build_cells() ;
}

//------------------------------------------------------------------------
GE_RefinedMeshing:: ~GE_RefinedMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_RefinedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
  
   return sop->nb_points() ;
}

//------------------------------------------------------------------------
size_t
GE_RefinedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: nb_cells" ) ;
   size_t result = nb_meshes( nb_space_dimensions() ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t
GE_RefinedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: nb_faces" ) ;
   size_t result = nb_meshes( nb_space_dimensions()-1 ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   idx_vertex = 0 ;
}

//------------------------------------------------------------------------
bool
GE_RefinedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return idx_vertex<sop->nb_points() ;   
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   idx_vertex++ ;   
}

//------------------------------------------------------------------------
doubleVector const& 
GE_RefinedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
   GE_Point const* vert = sop->point( idx_vertex ) ;
   for( size_t i=0 ; i<nb_space_dimensions() ; i++ )
      vertex_coord( i ) =  vert->coordinate( i ) ;
   
   doubleVector const& result = vertex_coord ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   CELL_IT = true ;
   SIDE_IT = false ;

   idx_mesh=0 ;
   initial->start_cell_iterator() ;
   local=0 ;
   
   the_mesh_polyhedron_name = mesh_poly_name ;
   connectivity = true ;
   ref_poly = mesh_ref ;

   MESH_VERTICES.re_initialize( mesh_ref->nb_vertices() ) ;
   CELL_SIDES.re_initialize( mesh_ref->nb_faces() ) ;

   find_next_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   local++ ;
   last_side++ ;
   
   find_next_cell() ;
}

//------------------------------------------------------------------------
bool
GE_RefinedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: valid_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = CELL_IT && initial->valid_cell() ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
std::string const&
GE_RefinedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = the_mesh_polyhedron_name ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_RefinedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = MESH_VERTICES ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_RefinedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = CELL_SIDES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   CELL_IT = false ;
   SIDE_IT = true ;

   idx_mesh=0 ;
   initial->start_cell_iterator() ;
   local=0 ;
   
   the_mesh_polyhedron_name = side_poly_name ;
   connectivity = false ;
   last_side = 0 ;
   ref_poly = side_ref ;
   
   MESH_VERTICES.re_initialize( side_ref->nb_vertices() ) ;

   find_next_side() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   local++ ;
   last_side++ ;
   
   find_next_side() ;
}

//------------------------------------------------------------------------
bool
GE_RefinedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: valid_face" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = SIDE_IT && initial->valid_cell() ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
std::string const&
GE_RefinedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = the_mesh_polyhedron_name ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_RefinedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = MESH_VERTICES ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   initial->print( os, indent_width+3 ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: display( std::ostream& os, size_t indent_width )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: display" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const s = std::string( indent_width, ' ' ) ;
   os << s << "Initial mesh : " << endl ;
   initial->display( os, indent_width+2 ) ;
   os << endl ;
   os << s << "Refined mesh : " << endl ;
   GE_Meshing::display( os, indent_width+2  ) ;
   os << endl ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_RefinedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = sop->color( idx_vertex ) ;

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_RefinedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = initial->cell_color() ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_RefinedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = 
                static_cast<GE_Color const*>( SIDE_COLORS->at( last_side ) ) ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: build_vertices( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: build_vertices" ) ;
   
   bool first = true ;
   
//    bool check_consistency = GE_Mpolyhedron::check_consistency() ;
//    if( check_consistency ) GE_Mpolyhedron::unset_check_consistency() ;
   
   for( initial->start_vertex_iterator() ;
        initial->valid_vertex() ;
        initial->go_next_vertex() )
   {
      GE_Point* pt = GE_Point::create( sop, initial->vertex_coordinates() ) ;
      sop->append( pt, initial->vertex_color() ) ;
   }
   size_t dim = initial->nb_space_dimensions() ;
   
   size_t iSide=0 ;
   initial->start_face_iterator() ;
   for( ; initial->valid_face() ; initial->go_next_face() )
   {
      GE_Mpolyhedron const* poly =
         GE_Mpolyhedron::create( initial->face_polyhedron_name(),
                                 sop,  initial->face_vertices() ) ;
      initial_side_poly->set_at( iSide,
                                 const_cast<GE_Mpolyhedron*>(poly) ) ;
      initial_side_color->set_at( iSide,
                                  const_cast<GE_Color*>(initial->face_color()) ) ;
      iSide++ ;
   }
   
   GE_Point* real_pt = GE_Point::create( this, dim ) ;

   size_t iMesh=0 ;
   initial->start_cell_iterator() ;
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      GE_ReferencePolyhedronRefiner const* refi = refiner( initial ) ;
      if( first )
      {
         global_vertices.re_initialize( initial->nb_cells(),
                                        refi->nb_vertices() ) ;
         first = false ;
      }
      else
      {
         PEL_ASSERT( global_vertices.index_bound(1)==refi->nb_vertices() ) ;
      }
      GE_Mpolyhedron const* poly =
         GE_Mpolyhedron::create( initial->cell_polyhedron_name(),
                                 sop, initial->cell_vertices() ) ;
      for( size_t iVert=0 ; iVert<refi->nb_vertices() ; ++iVert )
      {
         GE_Point const* ref_pt = refi->vertex( iVert ) ;
         PEL_CHECK( dynamic_cast<GE_Point const*>(ref_pt)!=0 ) ;
         poly->apply_mapping( ref_pt, real_pt ) ;
         if( !sop->has( real_pt ) )
         {
            GE_Color const* col = local_color( real_pt ) ;
            sop->append( real_pt, col ) ;
            if( verbose )
            {
               PEL::out() << "Building vertex "
                          << sop->index( real_pt ) << " at " ;
               real_pt->print( PEL::out(), 1 ) ;
               PEL::out() << endl ;
            }
         }
         global_vertices( iMesh, iVert ) = sop->index( real_pt ) ;            
      }
      iMesh++ ;
   }
   
//    if( check_consistency ) GE_Mpolyhedron::set_check_consistency() ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: build_sides( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: build_sides" ) ;
   
   size_t dim = initial->nb_space_dimensions() ;
   GE_Point* bary = GE_Point::create( 0, dim ) ;
   PEL_BalancedBinaryTree* sides_tree = PEL_BalancedBinaryTree:: create( 0 ) ;
   
   bool first = true ;
   
   size_t iMesh=0 ;
   for( initial->start_cell_iterator() ;
        initial->valid_cell() ;
        initial->go_next_cell() )
   {
      GE_ReferencePolyhedronRefiner const* refi = refiner( initial ) ;
      if( refi->subface_reference_polyhedron()->nb_vertices() != side_ref->nb_vertices() )
                            raise_invalid_face_polyhedron( side_poly_name ) ;
      if( first )
      {
         size_t nb_initial_mesh = initial->nb_cells() ;
         
         global_sides.re_initialize( nb_initial_mesh,
                                     refi->nb_subfaces() ) ;
         SIDE_COLORS = PEL_Vector::create( this,
                                           refi->nb_subfaces()*nb_initial_mesh ) ;
         first = false ;
      }
      else
      {
         PEL_ASSERT( global_sides.index_bound(1)==refi->nb_subfaces() ) ;
      }
      for( size_t iSide=0 ; iSide<refi->nb_subfaces() ; iSide++ )
      {
         size_t nb_vert = refi->subface_reference_polyhedron()->nb_vertices() ;
         size_t_vector vertices(nb_vert) ;
         for( size_t iVert=0 ; iVert<nb_vert ; iVert++ )
         {
            size_t global_vert = 
                   global_vertices( iMesh, refi->subface_vertex(iSide,iVert) ) ;
            vertices(iVert) =  global_vert ;
            GE_Point const* vert = sop->point( global_vert ) ;
            for( size_t j=0 ; j< dim ; j++ )
            {
               double old = 0. ;
               if( iVert!=0 )
                  old = bary->coordinate(j) ;
               bary->set_coordinate( j, old+vert->coordinate(j)/nb_vert ) ;
            }
         }
         PEL_IndexSet* side_vertices_set =
            PEL_IndexSet::create( 0, vertices, sides_tree->count() ) ;
         bool const found = sides_tree->has( side_vertices_set ) ;
         
         size_t ind = PEL::bad_index() ;
         if( found )
         {
            PEL_IndexSet const* l =
               static_cast<PEL_IndexSet const*>(
                                    sides_tree->item( side_vertices_set ) ) ;
            ind = l->id() ;
            side_vertices_set->destroy() ; side_vertices_set = 0 ;
         }
         else
         {
            side_vertices_set->set_owner( sides_tree ) ;
            ind = side_vertices_set->id() ;
            sides_tree->extend( side_vertices_set ) ;
      
            GE_Color const* s_col = local_color( bary ) ;
            SIDE_COLORS->set_at( ind, const_cast<GE_Color*>( s_col ) ) ;
            if( verbose )
            {
               PEL::out() << "Building side " << ind << " at " ;
               bary->print( PEL::out(), 1 ) ;
               PEL::out() << " " << s_col->name() << endl ;
            }
            SIDE_COLORS->set_at(
               ind, const_cast<GE_Color*>( local_color( bary ) ) ) ;
         }
         global_sides( iMesh, iSide ) = ind ;
      }
      iMesh++;
   }
   nb_meshes( dim-1 ) = sides_tree->count() ;
   
   sides_tree->destroy() ; sides_tree=0 ;
   bary->destroy() ; bary = 0 ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: build_cells( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: build_cells" ) ;
   
   size_t dim = initial->nb_space_dimensions() ;

   size_t nbM = 0 ;
   for( initial->start_cell_iterator() ;
        initial->valid_cell() ;
        initial->go_next_cell() )
   {
      GE_ReferencePolyhedronRefiner const* refi = refiner( initial ) ;
      if( refi->subcell_reference_polyhedron()->nb_vertices() != mesh_ref->nb_vertices() )
      {
         raise_invalid_cell_polyhedron( mesh_poly_name ) ;
      }
      nbM+= refi->nb_subcells() ;
   }
   nb_meshes( dim ) = nbM ;
   
}

//------------------------------------------------------------------------
GE_Color const*
GE_RefinedMeshing:: local_color( GE_Point const* pt ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( pt!=0 ) ;
   
   GE_Color const* result = initial->cell_color() ;
   size_t_vector const& initialConnectivity = initial->cell_faces() ;
   for( size_t j=0 ; j<initialConnectivity.size() ; j++ )
   {
      GE_Mpolyhedron const* poly = static_cast<GE_Mpolyhedron const*>(
         initial_side_poly->at( initialConnectivity( j ) ) ) ;
      PEL_CHECK( dynamic_cast<GE_Mpolyhedron const*>( poly ) !=0 ) ;
      if( poly->contains( pt ) )
      {
         GE_Color const* col = static_cast<GE_Color const*>(
            initial_side_color->at( initialConnectivity( j ) ) ) ;
         PEL_CHECK( dynamic_cast<GE_Color const*>( col ) != 0 ) ;
         result = col ;
      }
   }
   PEL_CHECK( result!=0 ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: find_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: find_next_cell" ) ;
   PEL_CHECK( valid_cell() ) ;
   
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      GE_ReferencePolyhedronRefiner const* refi = refiner( initial ) ;
      PEL_ASSERT( refi->subcell_reference_polyhedron()->nb_vertices()==mesh_ref->nb_vertices() ) ;
      if( local < refi->nb_subcells() )
      {
         PEL_ASSERT( refi->subcell_reference_polyhedron()->nb_faces()== mesh_ref->nb_faces() ) ;
         for( size_t iv=0 ; iv<refi->subcell_reference_polyhedron()->nb_vertices() ; iv++ )
         {
            MESH_VERTICES(iv) = global_vertices( idx_mesh, refi->subcell_vertex( local, iv ) ) ;
         }
         for( size_t is=0 ; is<refi->subcell_reference_polyhedron()->nb_faces() ; ++is )
         {
            CELL_SIDES( is ) =
                  global_sides( idx_mesh, refi->subcell_face( local, is ) ) ;
            
         }
         break ;
      }
      idx_mesh ++ ;
      local = 0 ;
   }   
}

//------------------------------------------------------------------------
void
GE_RefinedMeshing:: find_next_side( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: find_next_side" ) ;
   PEL_CHECK( valid_face() ) ;
   
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      GE_ReferencePolyhedronRefiner const* refi = refiner( initial ) ;
      bool found = false ;
         
      for( ; local < refi->nb_subfaces() ; local++ )
      {
         if( global_sides( idx_mesh, local )==(int)last_side )
         {
            last_side= global_sides( idx_mesh, local ) ;
            for( size_t i=0 ; i<refi->subface_reference_polyhedron()->nb_vertices() ; i++ )
            {
               MESH_VERTICES(i) =
                     global_vertices( idx_mesh, refi->subface_vertex( local, i ) ) ;
            }
            found = true ;
            break ;
         }
            
      }
      if( found ) break ;
      idx_mesh ++ ;
      local = 0 ;
   }   
}

//----------------------------------------------------------------------------
void
GE_RefinedMeshing:: append_refiner( GE_Color const* col,
                                    GE_ReferencePolyhedronRefiner const* refi )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: append_refiner" ) ;
   
   size_t idx_poly = refi->reference_polyhedron()->id_number() ;
   if( col == 0 )
   {
      if( REFs == 0 )
      {
         REFs = PEL_Vector::create( this, idx_poly + 1 ) ;  
      }
      else if( idx_poly >= REFs->index_limit() )
      {
         REFs->resize( idx_poly + 1 ) ;
      }
      if( REFs->at( idx_poly ) != 0 )
      {
         GE_ReferencePolyhedronRefiner* refi0 = 
             static_cast<GE_ReferencePolyhedronRefiner*>( REFs->at( idx_poly ) ) ;
         GE_RefinedMeshing_ERROR::n0( refi0, refi, 0 ) ;
      }
      REFs->set_at( idx_poly, const_cast< GE_ReferencePolyhedronRefiner* >( refi ) ) ;
   }
   else
   {
      size_t idx_col = col->identifier() ;
      PEL_Vector* refs = 0 ;
      if( REF_COLs == 0 )
      {
         REF_COLs = PEL_Vector::create( this, idx_col + 1 ) ;
         refs = PEL_Vector::create( REF_COLs, idx_poly+1 ) ;
         REF_COLs->set_at( idx_col, refs ) ;
      }
      else if( idx_col >= REF_COLs->index_limit() )
      {
         REF_COLs->resize( idx_col + 1 ) ;
         refs = PEL_Vector::create( REF_COLs, idx_poly+1 ) ;
         REF_COLs->set_at( idx_col, refs ) ;
      }
      else
      {
         PEL_Object* oo = REF_COLs->at( idx_col ) ;
         if( oo == 0 )
         {
            refs = PEL_Vector::create( REF_COLs, idx_poly+1 ) ;
            REF_COLs->set_at( idx_col, refs ) ;
         }
         else
         {
            refs = static_cast< PEL_Vector* >( oo ) ;
         }
      }
      PEL_ASSERT( refs != 0 ) ;
      if( refs->at( idx_poly ) != 0 )
      {
         GE_ReferencePolyhedronRefiner* refi0 = 
             static_cast<GE_ReferencePolyhedronRefiner*>( refs->at( idx_poly ) ) ;
         GE_RefinedMeshing_ERROR::n0( refi0, refi, col ) ;
      }
      refs->set_at( idx_poly, 
                    const_cast< GE_ReferencePolyhedronRefiner* >( refi ) ) ;
   }
}

//----------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner const*
GE_RefinedMeshing:: refiner( GE_Meshing const* mm ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedMeshing:: refiner" ) ;
   
   GE_ReferencePolyhedronRefiner const* result = 0 ;
   size_t idx_poly = mm->cell_reference_polyhedron()->id_number() ;
   if( REF_COLs != 0 )
   {
      size_t idx_col = mm->cell_color()->identifier() ;
      if( idx_col < REF_COLs->index_limit() )
      {
         PEL_Object* oo = REF_COLs->at( idx_col ) ;
         if( oo != 0 )
         {
            PEL_Vector* refs = static_cast< PEL_Vector* >( oo ) ;
            if( idx_poly <refs->index_limit() )
            {
               oo = refs->at( idx_poly ) ;
               if( oo != 0 )
               {
                  PEL_ASSERT( result == 0 ) ;
                  result = static_cast< GE_ReferencePolyhedronRefiner* >( oo ) ;
               }
            }
         } 
      }
   }
   if( (result == 0) && (REFs != 0) )
   {
      if( idx_poly < REFs->index_limit() )
      {
         PEL_Object* oo = REFs->at( idx_poly ) ;
         if( oo != 0 )
         {
            result = static_cast< GE_ReferencePolyhedronRefiner* >( oo ) ;
         }
      }
   }
   if( result == 0 ) GE_RefinedMeshing_ERROR::n1( mm->cell_color() ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void
GE_RefinedMeshing_ERROR:: n0( GE_ReferencePolyhedronRefiner const* r0,
                              GE_ReferencePolyhedronRefiner const* r1,
                              GE_Color const* col )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** GE_RefinedMeshing error:" << std::endl ;
   msg << "    multiple possible refinements defined" << endl ;
   if( col != 0 )
   {
      msg << "    for cells of color \"" << col->name() << "\"" << endl ;
   }
   msg << "      - " << r0->type_name() << endl ;
   msg << "      - " << r1->type_name() ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
GE_RefinedMeshing_ERROR:: n1(  GE_Color const* col )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** GE_RefinedMeshing error:" << std::endl ;
   msg << "    incomplete definition of refinements" << endl ;
   msg << "    (for instance, no refinement has been found for" << endl ;
   msg << "     cells of color \"" << col->name() << "\")" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
