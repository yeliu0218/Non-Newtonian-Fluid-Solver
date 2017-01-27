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

#include <GE_ExtrudedMeshing.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <iostream>
using std::string ;
using std::endl ;

GE_ExtrudedMeshing const* 
GE_ExtrudedMeshing::prototype = new GE_ExtrudedMeshing() ;

//-------------------------------------------------------------------------
GE_ExtrudedMeshing:: GE_ExtrudedMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_ExtrudedMeshing" )
   , initial(0 )
   , local_mesh_vertices( 0 )
   , local_mesh_connectivity( 0 )
   , sop( 0 )
   , global_sides(0,0,0)
   , global_vertices(0,0,0)
   , side_connectivity(0,0)
   , SIDE_COLORS(0)
   , nb_meshes( 0 )
   , idx_vertex( PEL::bad_index() )
   , idx_mesh( PEL::bad_index() )
   , vertex_coord( 0 )
   , DEFAULT_ELEVATION( 0 )
   , SIDE_TREE( 0 )
{
}

//-------------------------------------------------------------------------
GE_ExtrudedMeshing*
GE_ExtrudedMeshing:: create_replica( PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp,
                                        size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;
   
   if( dim_space != 3 ) raise_invalid_nb_space_dimensions( dim_space, "3" ) ;
   
   GE_ExtrudedMeshing* result = 
                          new GE_ExtrudedMeshing( a_owner, exp, dim_space ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;   
}

//-------------------------------------------------------------------------
GE_ExtrudedMeshing:: GE_ExtrudedMeshing( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp,
                                         size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , initial(0 )
   , local_mesh_vertices( 0 )
   , local_mesh_connectivity( 0 )
   , sop( GE_SetOfPoints::create( this, dim_space ) )
   , global_sides(0,0,0)
   , global_vertices(0,0,0)
   , side_connectivity(0,0)
   , SIDE_COLORS( PEL_Vector::create( this, 0 ) )
   , nb_meshes( dim_space+1 )
   , idx_vertex( PEL::bad_index() )
   , idx_mesh( PEL::bad_index() )
   , vertex_coord( dim_space )
   , IC_COLORS( PEL_Vector::create( this, 0 ) )
   , ELEVATIONS( PEL_Vector::create( this, 0 ) )
   , DEFAULT_ELEVATION( 0 )
   , max_elevation_size(0)
   , SIDE_TREE( PEL_BalancedBinaryTree:: create( this ) )
{
   PEL_CHECK( dim_space == 3 ) ;

   read_mesh_polyhedron( exp, dim_space, side_poly_name, mesh_poly_name ) ;

   side_ref = GE_Mpolyhedron::reference_polyhedron( side_poly_name ) ;
   mesh_ref = GE_Mpolyhedron::reference_polyhedron( mesh_poly_name ) ;

   // Initial meshing:
   {
      PEL_ModuleExplorer* sexp =
         exp->create_subexplorer( 0, "GE_Meshing" ) ;
      initial = GE_Meshing::create( this, sexp, 2 ) ;
      sexp->destroy() ; sexp = 0 ;
   }
   
   // Elevation:
   {
      PEL_ModuleExplorer* sexp =
            exp->create_subexplorer( 0, "vertices_coordinate_2" ) ;
      sexp->start_entry_iterator() ;
      for( ; sexp->is_valid_entry() ; sexp->go_next_entry() )
      {
         std::string const& str = sexp->keyword() ;
         PEL_Data const* e = sexp->data( sexp, 0 ) ;
         if( e->data_type() != PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
                                 sexp, str, PEL_Data::DoubleVector ) ;
         }
         if( !e->value_can_be_evaluated( 0 ) )
         {
            PEL_Error::object()->raise_not_evaluable(
               sexp, str, e->undefined_variables( 0 ) ) ;
         }
         doubleVector const& elev = e->to_double_vector() ;
         if( elev.size() == 1 )
         {
            PEL_Error::object()->raise_data_error(
               sexp, str, "At least two elements are expected" ) ;
         }
         max_elevation_size = PEL::max( max_elevation_size, elev.size() ) ;
         if( str == "default" )
         {
            DEFAULT_ELEVATION = elev ;
         }
         else
         {
            GE_Color const* col = cell_color_in_initial_meshing( str ) ;
            IC_COLORS->append( const_cast<GE_Color*>( col ) ) ;
            ELEVATIONS->append( PEL_DoubleVector::create( ELEVATIONS, elev ) ) ;
         }
      }
      sexp->destroy() ; sexp = 0 ;
   }

   GE_Color::extend( "behind" ) ;
   GE_Color::extend( "front" ) ;
   behind_color = GE_Color::object( "behind" ) ;
   front_color = GE_Color::object( "front" ) ;
   
   build_vertices() ;
   build_sides() ;
   build_meshes() ;
}

//------------------------------------------------------------------------
GE_ExtrudedMeshing:: ~GE_ExtrudedMeshing( void )
//------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
size_t
GE_ExtrudedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: nb_vertices" ) ;
   
   return sop->nb_points() ;
}

//------------------------------------------------------------------------
size_t
GE_ExtrudedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: nb_cells" ) ;
   size_t result = nb_meshes( nb_space_dimensions() ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t
GE_ExtrudedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: nb_faces" ) ;
   size_t result = nb_meshes( nb_space_dimensions()-1 ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: start_vertex_iterator" ) ;
   idx_vertex = 0 ;
}

//------------------------------------------------------------------------
bool
GE_ExtrudedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: valid_vertex" ) ;
   return idx_vertex<sop->nb_points() ;   
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   idx_vertex++ ;   
}

//------------------------------------------------------------------------
doubleVector const& 
GE_ExtrudedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: vertex_coordinates" ) ;
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
GE_ExtrudedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   idx_mesh=0 ;
   initial->start_cell_iterator() ;
   local=0 ;
   
   mesh_iterator_dim = nb_space_dimensions() ;

   the_mesh_polyhedron_name = mesh_poly_name ;
   connectivity = true ;
   ref_poly = mesh_ref ; 

   find_next_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_ExtrudedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: valid_cell" ) ;

   bool result = connectivity && initial->valid_cell() ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   local++ ;
   last_side++ ;
   
   find_next_cell() ;
}

//------------------------------------------------------------------------
std::string const&
GE_ExtrudedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = the_mesh_polyhedron_name ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_ExtrudedMeshing:: cell_reference_polyhedron( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: cell_reference_polyhedron" ) ;
   PEL_CHECK_PRE( cell_reference_polyhedron_PRE() ) ;

   GE_ReferencePolyhedron const* result = ref_poly ;

   PEL_CHECK_POST( cell_reference_polyhedron_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExtrudedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = local_mesh_vertices ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExtrudedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = local_mesh_connectivity ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   idx_mesh=0 ;
   initial->start_cell_iterator() ;
   local=0 ;
   
   the_mesh_polyhedron_name = side_poly_name ;
   connectivity = false ;
   last_side = 0 ;

   find_next_side() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_ExtrudedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: valid_face" ) ;

   //??? connectivity pour savoir si boucle sur sides
   bool result = (!connectivity) && initial->valid_cell() ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   local++ ;
   last_side++ ;
   
   find_next_side() ;
}

//------------------------------------------------------------------------
std::string const&
GE_ExtrudedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = the_mesh_polyhedron_name ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_ExtrudedMeshing:: face_reference_polyhedron( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: face_reference_polyhedron" ) ;
   PEL_CHECK_PRE( face_reference_polyhedron_PRE() ) ;

   GE_ReferencePolyhedron const* result = ref_poly ;

   PEL_CHECK_POST( face_reference_polyhedron_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_ExtrudedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = local_mesh_vertices ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   initial->print( os, indent_width+3 ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: display( std::ostream& os, size_t indent_width )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: display" ) ;

   string space( indent_width, ' ' ) ;
   
   os << space << "Initial mesh : " << endl ;
   initial->display( os, indent_width+2 ) ;
   os << endl ;
   os << space << "Extruded mesh : " << endl ;
   GE_Meshing::display( os, indent_width+2  ) ;
   os << endl ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExtrudedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;

   GE_Color const* result = sop->color( idx_vertex ) ;

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExtrudedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = m_color ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExtrudedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = m_color ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: build_vertices( void )
//------------------------------------------------------------------------
{
   static std::string mesg = "GE_ExtrudedMeshing : inconsistent \"vertices_coordinate_2\" module" ;
   bool first = true ;
   size_t_vector nb_sweeping_steps( initial->nb_vertices() ) ;

   GE_SetOfPoints* ini_sop = GE_SetOfPoints::create( 0, 2 ) ;
   initial->start_vertex_iterator() ;
   for( ; initial->valid_vertex() ; initial->go_next_vertex() )
   {
      GE_Point* ini_pt = GE_Point::create( ini_sop, 
                                           initial->vertex_coordinates() ) ;
      ini_sop->append( ini_pt, initial->vertex_color() ) ;
   }
   
   GE_Point* pt = GE_Point::create( this, nb_space_dimensions() ) ;
   
   size_t iMesh=0 ;
   initial->start_cell_iterator() ;
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      doubleVector const& elev = elevation( initial->cell_color() ) ;
      size_t const nb_vert = initial->cell_reference_polyhedron()->nb_vertices() ;
      if( first )
      {
         global_vertices.re_initialize( nb_vert,
                                        max_elevation_size,
                                        initial->nb_cells() ) ;
         first = false ;
      }
      else
      {
         if( global_vertices.index_bound(0)<nb_vert )
         {
            global_vertices.raise_first_index_bound( nb_vert ) ;
         }
      }
      size_t_vector const& vert = initial->cell_vertices() ;
      
      for( size_t j=0 ; j<nb_vert ; j++ )
      {
         GE_Point const* ini = ini_sop->point( vert(j) ) ;
         pt->set_coordinate( 0, ini->coordinate( 0 ) ) ;
         pt->set_coordinate( 1, ini->coordinate( 1 ) ) ;
         for( size_t i=0 ; i<elev.size() ; i++ )
         {
            pt->set_coordinate( 2, elev(i) ) ;
            
            GE_Color const* extra_color = 0 ;
            if( i==0 )
            {
              extra_color  = behind_color ;
            }
            else if( i==elev.size()-1 )
            {
              extra_color  = front_color ;
            }
            GE_Color const* col = merge_color( extra_color,
                                               ini_sop->color( vert(j) ) ) ;
            if( i >= nb_sweeping_steps( vert(j) ) )
            {
	       if( sop->has( pt ) ) PEL_Error::object()->raise_plain( mesg ) ;
               sop->append( pt, col ) ;
               nb_sweeping_steps( vert(j) )++ ;
	    }
	    else
	    {
	       if( !sop->has( pt ) ) PEL_Error::object()->raise_plain( mesg ) ;
	    }
            global_vertices( j,i,iMesh ) = sop->index( pt ) ;
	 }
      }            
      iMesh++ ;
   }
   ini_sop->destroy() ; ini_sop=0 ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_ExtrudedMeshing:: merge_color( GE_Color const* first,
                                    GE_Color const* second )
//------------------------------------------------------------------------
{
   PEL_CHECK( second!=0 ) ;
   GE_Color const* result = second ;
   if( first!=0 )
   {
      if( second!=GE_Color::null_color() )
      {
         std::string new_name = first->name() + "_"  ;
         new_name += second->name() ;
         GE_Color::extend( new_name ) ;       
         result = GE_Color::object( new_name ) ;
      }
      else
      {
         result = first ;
      }
   }
   return result ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: recover_initial_sides( intArray2D& old_side_vertices,
                                              PEL_Vector * initial_face_color ) 
//------------------------------------------------------------------------
{
   size_t iSide=0 ;
   initial->start_face_iterator() ;
   for( ; initial->valid_face() ; initial->go_next_face() )
   {
      size_t_vector const& vert = initial->face_vertices() ;
      PEL_ASSERT( vert.size()==2 ) ;
      for( size_t i=0 ; i<vert.size() ; i++ )
         old_side_vertices( iSide, i ) = vert(i) ;
      initial_face_color->set_at( iSide,
                             const_cast<GE_Color*>( initial->face_color() ) ) ;
      iSide++;
   }
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: build_vertical_side( size_t_vector& loc,
                                            GE_Color const*& s_color,
                                            size_t old_side,
                                            size_t iz,
                                            size_t iMesh,
                                            size_t_vector const& vert,
                                            intArray2D const& old_side_vertices,
                                            PEL_Vector const* initial_face_color ) const
//------------------------------------------------------------------------
{
   s_color = static_cast<GE_Color const*>(
      initial_face_color->at( old_side ) ) ;
   size_t nb = old_side_vertices.index_bound(1) ;
   PEL_CHECK( nb==2 ) ;
   loc.re_initialize( 2*nb ) ;
   for( size_t i=0 ; i<nb ; i++ )
   {
      size_t k = old_side_vertices( old_side, i ) ;
      size_t v_loc = 0 ;
      for( ; v_loc<vert.size() ; v_loc++ )
      {
         if( vert(v_loc)==k ) break ;
      }
      PEL_ASSERT( v_loc<vert.size() ) ;
      loc(i) = global_vertices( v_loc, iz, iMesh ) ;
      loc(2*nb-1-i) = global_vertices( v_loc, iz+1, iMesh ) ;
   }
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: build_horizontal_side( size_t_vector& loc,
                                              GE_Color const*& s_color,
                                              bool up_side,
                                              size_t idx_z,
                                              size_t iMesh,
                                              doubleVector const& elev,
                                              size_t_vector const& vert ) const           
//------------------------------------------------------------------------
{
   size_t nbz = elev.size() ;
   PEL_CHECK( idx_z<nbz-1 ) ;
   bool first_mesh = idx_z==0 ;
   bool last_mesh = idx_z==nbz-2 ;
   size_t nbv = vert.size() ;
   
   loc.re_initialize( nbv ) ;
   if( !up_side && first_mesh )
   {
      s_color = behind_color ;
   }
   else if( up_side && last_mesh )
   {
      s_color = front_color ;
   }
   size_t iz = ( up_side ? idx_z+1 : idx_z ) ;
   
   for( size_t v_loc=0 ; v_loc<nbv ; v_loc++ )
   {
      loc( v_loc ) = global_vertices( v_loc, iz, iMesh ) ;
   }
}

//------------------------------------------------------------------------
size_t
GE_ExtrudedMeshing:: extend_side( size_t_vector const& loc,
                                  GE_Color const* s_color )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: extend_side" ) ;
   
   size_t result = PEL::bad_index() ;

   PEL_IndexSet* indx = PEL_IndexSet::create( 0,
                                              loc,
                                              SIDE_TREE->count() ) ;

   // Search if side has not already been registred
   bool const found = SIDE_TREE->has( indx ) ;
   if( found )
   {
      PEL_IndexSet const* l =
         static_cast<PEL_IndexSet const*>( SIDE_TREE->item( indx ) ) ;
      result = l->id() ;
      indx->destroy() ; indx = 0 ;
   }
   else
   {
      indx->set_owner( SIDE_TREE ) ;
      result = indx->id() ;
      SIDE_TREE->extend( indx ) ;
  
      size_t const bucket_size = 1024 ;
      size_t nb = SIDE_COLORS->index_limit() ;
      if( result >= nb )
      {
         size_t new_size = result+bucket_size ;
         SIDE_COLORS->resize( new_size ) ;
         side_connectivity.raise_first_index_bound( new_size ) ;
      }  
      SIDE_COLORS->set_at( result, const_cast<GE_Color*>( s_color ) ) ;
      for( size_t i=0 ; i<loc.size() ; i++ )
      {
         side_connectivity( result, i ) = loc(i) ;
      }
   }
   return( result ) ;      
}
//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: build_sides( void )
//------------------------------------------------------------------------
{
   intArray2D old_side_vertices( initial->nb_faces() , 2 ) ;
   PEL_Vector* initial_face_color =
                               PEL_Vector::create( 0, initial->nb_faces() ) ;
   recover_initial_sides( old_side_vertices, initial_face_color ) ;
   
   size_t iMesh=0 ;
   bool first = true ;
   
   initial->start_cell_iterator() ;
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      size_t_vector const& vert = initial->cell_vertices() ;
      size_t_vector const& conn = initial->cell_faces() ;
      doubleVector const& elev = elevation( initial->cell_color() ) ;
      
      size_t const nb_vert = vert.size() ;
      size_t const nbs = conn.size() ;
      PEL_ASSERT( nb_vert==nbs ) ;
         
      if( first )
      {
         global_sides.re_initialize( elev.size(),
                                     nbs+2,
                                     initial->nb_cells() ) ;
         side_connectivity.re_initialize( 0,
                                          PEL::max( nb_vert, (size_t)4 ) ) ;
         first = false ;
      }
      else
      {
         PEL_ASSERT( global_sides.index_bound(1)==nbs+2 ) ;
         if( elev.size() > global_sides.index_bound(0) )
         {
            global_sides.raise_first_index_bound( elev.size() ) ;
         }
      }
      size_t_vector loc(0) ;
      for( size_t iz=0 ; iz<elev.size()-1 ; iz++ )
      {
         for( size_t iSide=0 ; iSide<nbs+2 ; iSide++ )
         {
            GE_Color const* s_color = GE_Color::null_color() ;
            if( iSide<nbs )
            {
               size_t old_side = conn(iSide) ;
               build_vertical_side( loc, s_color, old_side, iz, iMesh, vert,
                                    old_side_vertices, initial_face_color ) ;
               
            }
            else
            {
               bool up_side = iSide==nbs+1 ;
               build_horizontal_side( loc, s_color, up_side, iz,
                                      iMesh, elev, vert ) ;
            }
            global_sides( iz, iSide, iMesh ) = extend_side( loc,
                                                            s_color ) ;
         }
      }
      iMesh++;
   }
   nb_meshes( nb_space_dimensions()-1 ) = SIDE_TREE->count() ;
   
   initial_face_color->destroy() ; initial_face_color=0 ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: build_meshes( void )
//------------------------------------------------------------------------
{
   size_t nbM = 0 ;
   initial->start_cell_iterator() ;
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      if( initial->cell_reference_polyhedron()->nb_vertices() != 4 )
      {
         PEL_Error::object()->raise_plain(
            "GE_ExtrudedMeshing :\n"
            "   a quadrangular 2D meshing is expected" ) ;
      }
      doubleVector const& elev = elevation( initial->cell_color() ) ;
      nbM+= elev.size()-1 ;
   }
   nb_meshes( nb_space_dimensions() ) = nbM ;
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: find_next_cell( void )
//------------------------------------------------------------------------
{
   for(  ; initial->valid_cell() ; initial->go_next_cell() )
   {
      doubleVector const& elev = elevation( initial->cell_color() ) ;

      local_mesh_vertices.re_initialize( cell_reference_polyhedron()->nb_vertices() ) ;
      local_mesh_connectivity.re_initialize( cell_reference_polyhedron()->nb_faces() ) ;

      if( local < elev.size()-1 )
      {
         size_t nb_vert = initial->cell_vertices().size() ;
         for( size_t i=0 ; i<nb_vert ; i++ )
         {
            local_mesh_vertices(i) = global_vertices( i, local, idx_mesh ) ;
            local_mesh_vertices(i+nb_vert) = global_vertices( i, local+1, idx_mesh ) ;
         }
         size_t nbs = initial->cell_faces().size() ;
         for( size_t i=0 ; i<nbs+2 ; i++ )
         {
            local_mesh_connectivity(i) = global_sides( local, i, idx_mesh ) ;
         }
         m_color = initial->cell_color() ;
         GE_Color const* extra_color = 0 ;
         if( local==0 )
         {
            extra_color  = behind_color ;
         }
         else if( local==elev.size()-2 )
         {
            extra_color = front_color ;
         }
         m_color = merge_color( extra_color, initial->cell_color() ) ;
         break ;
      }  
      idx_mesh ++ ;
      local = 0 ;
   }   
}

//------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: find_next_side( void )
//------------------------------------------------------------------------
{
   for(  ; initial->valid_cell() ; initial->go_next_cell() )
   {
      doubleVector const& elev = elevation( initial->cell_color() ) ;

      size_t nb_sides_ini = initial->cell_faces().size() ;
      size_t nbs = nb_sides_ini+2 ;
         
      bool found = false ;
         
      for( ; local < (elev.size()-1)*nbs ; local++ )
      {
         size_t iz = local/nbs ;
         size_t loc_side = local % nbs ;
         if( iz<elev.size()-1 )
         {
            if( global_sides( iz, loc_side, idx_mesh )==(int)last_side )
            {
               if( loc_side<nb_sides_ini )
               {
                  ref_poly = side_ref ;
               }
               else
               {
                  ref_poly = initial->cell_reference_polyhedron() ;
               }
               size_t nb = ref_poly->nb_vertices() ; 

   local_mesh_vertices.re_initialize( face_reference_polyhedron()->nb_vertices() ) ;

               for( size_t i=0 ; i<nb ; i++ )
               {
                  local_mesh_vertices(i) = side_connectivity( last_side, i ) ;
               }
               found = true ;
               m_color = static_cast<GE_Color const*>( SIDE_COLORS->at( last_side ) ) ;
               PEL_CHECK( dynamic_cast<GE_Color const*>( m_color )!=0 ) ;
            }
         }
         if( found ) break ;
      }
      if( found ) break ;
      idx_mesh ++ ;
      local = 0 ;
   }   
}

//------------------------------------------------------------------------
doubleVector const& 
GE_ExtrudedMeshing:: elevation( GE_Color const* cell ) const
//------------------------------------------------------------------------
{
   doubleVector const* result =0 ;

   size_t i = 0 ;
   for( ; i<IC_COLORS->index_limit() ; ++i )
   {
      GE_Color const* col = static_cast<GE_Color*>( IC_COLORS->at(i) ) ;
      if( col->is_matching( cell ) )
      {
         PEL_DoubleVector const* lst = 
                        static_cast<PEL_DoubleVector*>( ELEVATIONS->at(i) ) ;
         result = &lst->to_double_vector() ;
      }
   }
   if( result != 0 )
   {
      for( ++i ; i<IC_COLORS->index_limit() ; ++i )
      {
         GE_Color const* col = static_cast<GE_Color*>( IC_COLORS->at(i) ) ;
         if( col->is_matching( cell ) )
         {
	    PEL_Error::object()->raise_plain( "GE_ExtrudedMeshing : conflict " ) ;
         } 
      }
   }

   if( result==0 && DEFAULT_ELEVATION.size() > 1 )
   {
      result =  &DEFAULT_ELEVATION ;
   }

   if( result == 0 )
   {
      std::string mesg = "GE_ExtrudedMeshing : \n" ;
      mesg += "   \"vertices_coordinate_2\" module does not account\n" ;
      mesg += "   for cells of color \"" + cell->name() + "\"\n";
      PEL_Error::object()->raise_plain( mesg ) ;
   }
   
   return *result ;
}

//???????? quasiment dupliqué avec GE_RefinedMeshing
//-------------------------------------------------------------------------
GE_Color const*
GE_ExtrudedMeshing:: cell_color_in_initial_meshing( 
                                               std::string const& str ) const
//-------------------------------------------------------------------------
{
   GE_Color const* result = GE_Color::object( str ) ;

   bool ok = false ;

   initial->start_cell_iterator() ;
   for( ; initial->valid_cell() ; initial->go_next_cell() )
   {
      GE_Color const* col = initial->cell_color() ;
      if( col->is_matching( result ) )
      {
         ok = true ;
         break ;
      }
   }

   if( !ok )
   {
      std::string mesg = "GE_ExtrudedMeshing : \n" ;
      mesg += "   the meshing to be extruded does not have any cell\n" ;
      mesg += "   whose colors matches the color of name \"" + str +"\"";
      PEL_Error::object()->raise_plain( mesg ) ;
   }

   return( result ) ;
}

//-------------------------------------------------------------------------
void
GE_ExtrudedMeshing:: check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ExtrudedMeshing:: check_mesh_polyhedron" ) ;
   PEL_CHECK( check_mesh_polyhedron_PRE(
                            dim_space, face_poly_name, cell_poly_name ) ) ;

   // Check the cell reference polyhedron :
   GE_ReferencePolyhedron const* cell_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   if( cell_poly_ref->nb_vertices() != 8 )
   {
      raise_invalid_cell_polyhedron( cell_poly_name, 
                                     "   only hexahedral cells are allowed" ) ;
   }
   
   // Check the side reference polyhedron :
   GE_ReferencePolyhedron const* side_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( face_poly_name ) ;
   if( side_poly_ref->nb_vertices() != 4 )
   {
      raise_invalid_face_polyhedron( side_poly_ref->name(), 
                                     "   only quadrilateral faces are allowed" ) ;
   }
}
