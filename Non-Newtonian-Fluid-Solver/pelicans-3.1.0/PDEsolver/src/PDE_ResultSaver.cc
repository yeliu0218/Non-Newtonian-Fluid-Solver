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

#include <PDE_ResultSaver.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_DataOnMeshingWriter.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_IntVector.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_SimplePolygon2D.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ;
using std::endl ;
using std::string ;

struct PDE_ResultSaver_ERROR
{
   static void n0( std::string const& name ) ;
   static void n1( std::string const& name ) ;
   static void n2( void ) ;
   static void n3( void ) ;
} ;

//----------------------------------------------------------------------
PDE_ResultSaver:: PDE_ResultSaver( PDE_DomainAndFields* dom,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( dom )
   , MY_DOM( dom )
   , MY_EXP( exp->create_clone( this ) )
   , OPENED_CYCLE( false )
   , I_CYCLE( PEL::bad_index() )
   , VERTICES( dom->set_of_vertices() )
   , FIELDS( dom->set_of_discrete_fields() )
   , WRITERS( PEL_List::create( this ) )
   , MOD( 0 )
   , MOD_VARS( 0 )
   , MOD_FF( 0 )
   , NB_VARS( 0 )
   , NB_FIELDS( 0 )
   , EPS_DBL( exp->has_entry( "dbl_epsilon" ) ?
                 exp->double_data( "dbl_epsilon" ) : 1.E-4 )
   , MIN_DBL( exp->has_entry( "dbl_minimum" ) ?
                 exp->double_data( "dbl_minimum" ) : 1.E-8 )
   , MESHING_NAME( "" )
   , NB_ACTIVE_VERTS( 0 )
   , NB_ACTIVE_CELLS( 0 )
   , NB_ACTIVE_FACES( 0 )
   , ACTIVE_IDX_OF_VERT( 0 )
   , ACTIVE_IDX_OF_CELL( 0 )
   , PERIODIC( false )
   , DETECT_ACTIVE_ITEMS( ( dom->adapter_CHARMS() != 0 ) ||
                          ( dom->adapter_HN() != 0 ) ) 
   , cFE( dom->create_LocalFEcell( this ) )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , FIELD_LOCATION( InvalidLocation )
{
   cFE->include_color( GE_Color::halo_color() ) ;
   sFE->include_color( GE_Color::halo_color() ) ;
   bFE->include_color( GE_Color::halo_color() ) ;

   stringVector const& formats = exp->stringVector_data( "writers" ) ;
   for( size_t i=0 ; i<formats.size() ; ++i )
   {
      for( size_t j=i+1 ; j<formats.size() ; ++j )
      {
         if( formats(i) == formats(j) )
         {
            PEL_Error::object()->raise_bad_data_value(
               exp, "writers",
               "The writer \""+formats(i)+"\" is defined twice." ) ;
         }
      }
      PEL_DataOnMeshingWriter* ww = 
	      PEL_DataOnMeshingWriter::make( this, formats(i), exp ) ;
      WRITERS->append( ww ) ;
   }
   
   if( EPS_DBL <= 0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "dbl_epsilon", "A positive value is expected" ) ;
   }
   if( MIN_DBL <= 0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "dbl_minimum", "A positive value is expected" ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( owner() == dom ) ;
   PEL_CHECK_POST( !has_an_opened_cycle() ) ;
   PEL_CHECK_POST( cycle_number() == PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
PDE_ResultSaver:: ~PDE_ResultSaver( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_DomainAndFields const*
PDE_ResultSaver:: attached_domain( void ) const
//----------------------------------------------------------------------
{
   return( MY_DOM ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: start_cycle( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: start_cycle" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_SAVEOLD( size_t, cycle_number, cycle_number() ) ;
   
   if( has_an_opened_cycle() )
   {
      terminate_cycle() ;
   }
   OPENED_CYCLE = true ;

   if( I_CYCLE == PEL::bad_index() )
   {
      I_CYCLE = 1 ;
   }
   else
   {
      ++I_CYCLE ;
   }

   NB_VARS = 0 ;
   NB_FIELDS = 0 ;

   std::ostringstream mod_name ;
   mod_name << "cycle_" << I_CYCLE ;
   MOD = PEL_Module::create( this, mod_name.str() ) ;

   MOD->add_entry( "cycle_number", PEL_Int::create( MOD, I_CYCLE ) ) ;

   MOD_VARS = PEL_Module::create( MOD, "variables" ) ;
   MOD->add_module( MOD_VARS ) ;

   MOD_FF = PEL_Module::create( MOD, "fields" ) ;
   MOD->add_module( MOD_FF ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( has_an_opened_cycle() ) ;
   PEL_CHECK_POST( ( cycle_number() == 1 ) ||
                   ( cycle_number() == OLD( cycle_number ) + 1 ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: terminate_cycle( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: terminate_cycle" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_INV( invariant() ) ;

   MOD->add_entry( "nb_variables", PEL_Int::create( MOD, NB_VARS ) ) ;
   MOD->add_entry( "nb_fields", PEL_Int::create( MOD, NB_FIELDS ) ) ;

   PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( 0, MOD ) ;
   PEL_ListIterator* it_writers = PEL_ListIterator::create( 0, WRITERS ) ;
   for( it_writers->start() ; it_writers->is_valid() ; it_writers->go_next() )
   {
      PEL_DataOnMeshingWriter* ww =
                 static_cast<PEL_DataOnMeshingWriter*>( it_writers->item() ) ;
      ww->write_cycle( exp ) ;
   }
   it_writers->destroy() ; it_writers = 0 ;
   exp->destroy() ; exp = 0 ;

   destroy_possession( MOD ) ;
   MOD = 0 ; MOD_VARS = 0 ; MOD_FF = 0 ;
   OPENED_CYCLE = false ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !has_an_opened_cycle() ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultSaver:: has_an_opened_cycle( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: has_an_opened_cycle" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( OPENED_CYCLE ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultSaver:: cycle_number( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: cycle_number" ) ;

   size_t result = I_CYCLE ;

   PEL_CHECK_POST( ( result == PEL::bad_index() ) || ( result >= 1 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultSaver:: has_variable( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: has_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( MOD_VARS->has_entry( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_variable( double const& value,
                                 std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( !has_variable( name ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   ++NB_VARS ;
   if( MOD_VARS->has_entry( name ) )
   {
      PDE_ResultSaver_ERROR::n0( name ) ;
   }
   MOD_VARS->add_entry( name, PEL_Double::create( MOD_VARS, value ) ) ;

   PEL_CHECK_POST( has_variable( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_variable( doubleVector const& value,
                                 std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( !has_variable( name ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   ++NB_VARS ;
   if( MOD_VARS->has_entry( name ) )
   {
      PDE_ResultSaver_ERROR::n0( name ) ;
   }
   MOD_VARS->add_entry( name, PEL_DoubleVector::create( MOD_VARS, value ) ) ;

   PEL_CHECK_POST( has_variable( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_variable( doubleArray2D const& value,
                                 std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( !has_variable( name ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   ++NB_VARS ;
   if( MOD_VARS->has_entry( name ) )
   {
      PDE_ResultSaver_ERROR::n0( name ) ;
   }
   MOD_VARS->add_entry( name, PEL_DoubleArray2D::create( MOD_VARS, value ) ) ;

   PEL_CHECK_POST( has_variable( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_variable( int const& value,
                                 std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( !has_variable( name ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++NB_VARS ;
   if( MOD_VARS->has_entry( name ) )
   {
      PDE_ResultSaver_ERROR::n0( name ) ;
   }
   MOD_VARS->add_entry( name, PEL_Int::create( MOD_VARS, value ) ) ;

   PEL_CHECK_POST( has_variable( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_variable( intVector const& value,
                                 std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_variable" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( !has_variable( name ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++NB_VARS ;
   if( MOD_VARS->has_entry( name ) )
   {
      PDE_ResultSaver_ERROR::n0( name ) ;
   }
   MOD_VARS->add_entry( name, PEL_IntVector::create( MOD_VARS, value ) ) ;
   
   PEL_CHECK_POST( has_variable( name ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_time( double time )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_time" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   
   save_variable( time, "TIME" ) ;
   MOD->add_entry( "Time", PEL_Double::create( MOD, time ) ) ;
   
   PEL_CHECK_POST( has_variable( "TIME" ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_grid( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_grid" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   infer_active_grid_items() ;
   
   // *** allocation of local arrays representing the grid
   size_t const ndims = cFE->nb_space_dimensions() ;
   doubleArray2D vertices( ndims, NB_ACTIVE_VERTS ) ;
   intVector cell_nb_vertices( NB_ACTIVE_CELLS ) ;
   intVector cell_nb_faces(    NB_ACTIVE_CELLS ) ;
   intVector face_nb_vertices( NB_ACTIVE_FACES ) ;
   intArray2D cell2vertex( 0, NB_ACTIVE_CELLS ) ;
   intArray2D cell2face(   0, NB_ACTIVE_CELLS ) ;
   intArray2D face2vertex( 0, NB_ACTIVE_FACES ) ;
   intArray2D face2cell(   2, NB_ACTIVE_FACES,
                           PEL_DataOnMeshingWriter::index_for_trash() ) ;
   intVector vertex_color( NB_ACTIVE_VERTS ) ;
   intVector cell_color(   NB_ACTIVE_CELLS ) ;
   intVector face_color(   NB_ACTIVE_FACES ) ;

   // *** vertices *** arrays: vertices, vertex_color
   for( size_t iv=0 ; iv<VERTICES->nb_points() ; ++iv )
   {
      GE_Point const* vv = VERTICES->point( iv ) ;
      size_t ii = ACTIVE_IDX_OF_VERT( iv ) ;
      if( ii != PEL::bad_index() )
      {
         PEL_ASSERT( ii < NB_ACTIVE_VERTS ) ;
         for( size_t d=0 ; d<ndims ; ++d ) 
         {
            vertices( d, ii ) = vv->coordinate( d ) ;
         }
         vertex_color( ii ) = VERTICES->color( iv )->identifier() ;
      }
   }

   // *** cells *** arrays: cell_nb_vertices, cell2vertex, cell_color
   size_t icell = 0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next(), ++icell ) 
   {
      set_mesh_vertices( cFE->polyhedron(), icell, 
                         cell_nb_vertices, cell2vertex ) ;
      cell_color( icell ) = cFE->color()->identifier() ;
   }

   // *** faces *** arrays: cell_nb_faces, cell2face, 
   //                       face_nb_vertices, face2vertex, face2cell, face_color
   PDE_LocalFEcell const* cfe_0 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* cfe_1 = sFE->adjacent_localFEcell( 1 ) ;
   size_t iface = 0 ;
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next(), ++iface ) 
   {
      GE_Mpolyhedron const* poly = ( sFE->is_periodic() ? 
                                   sFE->polyhedron( 0 ) : sFE->polyhedron() ) ;
      set_mesh_vertices( poly, iface, face_nb_vertices, face2vertex ) ;

      GE_Color const* col = ( sFE->is_periodic() ? 
                              sFE->color( 0 ) : sFE->color() ) ;
      face_color( iface ) = col->identifier() ;

      icell = ACTIVE_IDX_OF_CELL( cfe_0->mesh_id() ) ;
      set_face_cell_link( iface, 0, icell, 
                          cell_nb_faces, cell2face, face2cell ) ;

      if( !sFE->is_periodic() )
      {
         icell = ACTIVE_IDX_OF_CELL( cfe_1->mesh_id() ) ;
         set_face_cell_link( iface, 1, icell,
                             cell_nb_faces, cell2face, face2cell ) ;
      }
      else
      {
         PEL_ASSERT( poly->nb_vertices() == 
                     sFE->polyhedron( 1 )->nb_vertices() ) ;
         ++iface ;
         
         poly = sFE->polyhedron( 1 ) ;
         set_mesh_vertices( poly, iface, face_nb_vertices, face2vertex ) ;
                              
         face_color( iface ) = sFE->color( 1 )->identifier() ;

         icell = ACTIVE_IDX_OF_CELL( cfe_1->mesh_id() ) ;
         set_face_cell_link( iface, 0, icell,
                             cell_nb_faces, cell2face, face2cell ) ;
      }
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next(), ++iface ) 
   {
      GE_Mpolyhedron const* poly = bFE->polyhedron() ;
      set_mesh_vertices( poly, iface, face_nb_vertices, face2vertex ) ;

      face_color( iface ) = bFE->color()->identifier() ;

      icell = ACTIVE_IDX_OF_CELL( bFE->adjacent_cell_id() ) ;
      set_face_cell_link( iface, 0, icell,
                          cell_nb_faces, cell2face, face2cell ) ;
   }

   save_grid( ndims, vertices, 
              cell_nb_vertices, cell2vertex,
              cell_nb_faces, cell2face,
              face_nb_vertices, face2vertex,
              face2cell,
              vertex_color, cell_color, face_color ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( grid_is_saved() ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultSaver:: grid_is_saved( void ) const
//----------------------------------------------------------------------
{
   return( !MESHING_NAME.empty() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultSaver:: nb_saved_cells( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: nb_saved_cells" ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   return( NB_ACTIVE_CELLS ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultSaver:: nb_saved_faces( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: nb_saved_faces" ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   return( NB_ACTIVE_FACES ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultSaver:: nb_saved_vertices( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: nb_saved_vertices" ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   return( NB_ACTIVE_VERTS ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultSaver:: vertex_index( GE_Point const* vertex ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: vertex_index" ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   PEL_CHECK_PRE( vertex != 0 ) ;
   PEL_CHECK_PRE( vertex->nb_coordinates() ==
                            attached_domain()->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( attached_domain()->set_of_vertices()->has( vertex ) ) ;

   size_t result = VERTICES->index( vertex ) ;
   if( ACTIVE_IDX_OF_VERT.size() != 0 ) result = ACTIVE_IDX_OF_VERT( result ) ;
   
   PEL_CHECK_POST( result == PEL::bad_index() || result < nb_saved_vertices() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void 
PDE_ResultSaver:: update_active_vertices( GE_SetOfPoints const* vertices,
                                          GE_Mpolyhedron const* poly,
                                          size_t_vector& active_idx,
                                          size_t& nb_active )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: update_active_vertices" ) ;
   PEL_CHECK_PRE( vertices != 0 ) ;
   PEL_CHECK_PRE( poly != 0 ) ;
   PEL_CHECK_PRE( active_idx.size() == vertices->nb_points() ) ;
   
   for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
   {
      GE_Point const* pt = poly->vertex( i ) ;
      size_t iv = vertices->index( pt ) ;
      if( active_idx( iv ) == PEL::bad_index()  )
      {
         active_idx( iv ) = nb_active++ ;
      }
   }   
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: set_face_cell_link( size_t iface, size_t il, size_t icell,
                                      intVector& cell_nb_faces, 
                                      intArray2D& cell2face,
                                      intArray2D& face2cell ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: set_face_cell_link" ) ;
   PEL_CHECK_PRE( icell < cell_nb_faces.size() ) ;
   PEL_CHECK_PRE( ( il == 0 ) || ( il == 1 ) ) ;
   PEL_CHECK_PRE( iface < face2cell.index_bound( 1 ) ) ;
   PEL_CHECK_PRE( face2cell.index_bound( 0 ) == 2 ) ;
   
   PEL_ASSERT( face2cell( il, iface ) == 
               PEL_DataOnMeshingWriter::index_for_trash() ) ;
   face2cell( il, iface ) = icell ;
   
   cell_nb_faces( icell )++ ;
   
   if( (size_t)cell_nb_faces( icell ) > cell2face.index_bound( 0 ) )
      cell2face.raise_first_index_bound( cell_nb_faces( icell ) ) ;
   
   cell2face( cell_nb_faces( icell ) - 1, icell ) = iface ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_grid( 
          size_t ndims, doubleArray2D const& vertices,
          intVector const& cell_nb_vertices, intArray2D const& cell2vertex,
          intVector const& cell_nb_faces,    intArray2D const& cell2face,
          intVector const& face_nb_vertices, intArray2D const& face2vertex,
          intArray2D const& face2cell,
          intVector const& vertex_color,
          intVector const& cell_color,
          intVector const& face_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_grid" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( ndims == attached_domain()->nb_space_dimensions() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertices.index_bound( 0 ) == (size_t)ndims ) ;
   PEL_CHECK_PRE( vertices.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 1 ) == cell_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_vertices.size() ; ++ic ),
         ( cell_nb_vertices( ic ) > 0 ) &&
         ( (size_t)cell_nb_vertices( ic ) <= cell2vertex.index_bound( 0 ) ))) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_vertices.size() ; ++ic ),
         FORALL( ( size_t iv=0 ; iv<(size_t)cell_nb_vertices( ic ) ; ++iv ),
            cell2vertex( iv, ic ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( cell2face.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( cell2face.index_bound( 1 ) == cell2vertex.index_bound( 1 ) );
   PEL_CHECK_PRE( cell2face.index_bound( 1 ) == cell_nb_faces.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_faces.size() ; ++ic ),
         ( cell_nb_faces( ic ) > 0 ) &&
         ( (size_t)cell_nb_faces( ic ) <= cell2face.index_bound( 0 ) ) ) )  ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_faces.size() ; ++ic ),
         FORALL( ( size_t f=0 ; f<(size_t)cell_nb_faces( ic ) ; ++f ),
            cell2face( f, ic ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 1 ) == face_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         ( face_nb_vertices( f ) > 0 ) && 
         ( (size_t)face_nb_vertices( f ) <= face2vertex.index_bound( 0 ) ) )) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         FORALL( ( size_t iv=0 ; iv<(size_t)face_nb_vertices( f ) ; ++iv ),
            face2vertex( iv, f ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( face2cell.index_bound( 0 ) == 2 ) ;
   PEL_CHECK_PRE( face2cell.index_bound( 1 ) == face_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         face2cell( 0, f ) >= 0  ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         ( face2cell( 1, f ) == PEL_DataOnMeshingWriter::index_for_trash() ) ||
         ( face2cell( 1, f ) >= 0 ) ) ) ;
   PEL_CHECK_PRE( vertex_color.size() == vertices.index_bound( 1 ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t iv=0 ; iv<vertex_color.size() ; ++iv ),
         vertex_color( iv ) >= 0 ) ) ;
   PEL_CHECK_PRE( cell_color.size() == cell_nb_vertices.size( ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_color.size() ; ++ic ),
         cell_color( ic ) >= 0 ) ) ;
   PEL_CHECK_PRE( face_color.size() == face_nb_vertices.size( ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_color.size() ; ++f ),
         face_color( f ) >= 0 ) ) ;

   static int i_meshing = -1 ;

   std::ostringstream name ;
   name << "meshing_" << ++i_meshing ;
   MESHING_NAME = name.str() ;
   
   PEL_Module* mm = PEL_DataOnMeshingWriter::create_meshing_module( 
              MOD, 
              MESHING_NAME, 
              (int)ndims, vertices,
              cell_nb_vertices, cell2vertex,
              cell_nb_faces, cell2face,
              face_nb_vertices, face2vertex, face2cell,
              GE_Color::color_table(), GE_Color::halo_color()->name(),
              GE_Color::color_table_connectivity(),
              vertex_color, cell_color, face_color ) ;

   MOD->add_module( mm ) ;

   NB_ACTIVE_VERTS = vertices.index_bound( 1 ) ;
   NB_ACTIVE_CELLS = cell2vertex.index_bound( 1 ) ;
   NB_ACTIVE_FACES = face2vertex.index_bound( 1 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( grid_is_saved() ) ;
   PEL_CHECK_POST( nb_saved_vertices() == vertices.index_bound( 1 ) ) ;
   PEL_CHECK_POST( nb_saved_cells() == cell2vertex.index_bound( 1 ) ) ;
   PEL_CHECK_POST( nb_saved_faces() == face2vertex.index_bound( 1 ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_fields( size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_fields(level)" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   doubleArray2D values( 0, 0 ) ;
   doubleVector default_values( 0 ) ; 

   MY_EXP->start_module_iterator() ;
   for( ; MY_EXP->is_valid_module() ; MY_EXP->go_next_module() )
   {
      PEL_ModuleExplorer* exp = MY_EXP->create_subexplorer( 0 ) ;
      std::string const& save_name = exp->string_data( "entry_name" ) ;
      std::string const& location = exp->string_data( "where_to_save" ) ;
      exp->test_data_in( "where_to_save", "at_vertices," "at_cell_centers" ) ;
      if( exp->has_entry( "field" ) )
      {
         std::string const& name = exp->string_data( "field" ) ;
         PDE_DiscreteField const* field = FIELDS->item( name ) ;
         PDE_ResultSaver::SavingLocation loc = 
               ( location == "at_vertices" ? PDE_ResultSaver::AtVertices :
                 PDE_ResultSaver::AtCellCenters ) ;
         prepare_for_field_saving( loc, save_name, field->nb_components(),
                                   values, default_values ) ;
         build_reconstruction( field, level, loc, values, default_values ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain( "Bad syntax in module "+exp->name() ) ;
      }
      exp->destroy() ;

      save_field( values, default_values ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: prepare_for_field_saving( SavingLocation location,
                                            std::string const& save_name,
                                            size_t nb_components,
                                            doubleArray2D& values,
                                            doubleVector& default_values ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: prepare_for_field_saving" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   PEL_CHECK_PRE( !field_saving_in_progress() ) ;
   PEL_CHECK_PRE( !save_name.empty() ) ;
   PEL_CHECK_PRE( nb_components != 0 ) ;
   PEL_CHECK_PRE( (location == AtVertices) || (location == AtCellCenters) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   FIELD_LOCATION = location ;
   FIELD_NAME = save_name ;

   if( location == AtVertices )
   {
      values.re_initialize( nb_components, NB_ACTIVE_VERTS, 
                            undefined_value() ) ;
   }
   else if( location == AtCellCenters )
   {
      values.re_initialize( nb_components, NB_ACTIVE_CELLS, 
                            undefined_value() ) ;
   }
   default_values.re_initialize( nb_components, undefined_value() ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( field_saving_in_progress() ) ;
   PEL_CHECK_POST( saving_location() == location ) ; 
   PEL_CHECK_POST( values.index_bound( 0 ) == nb_components ) ;
   PEL_CHECK_POST(
      IMPLIES( location == AtCellCenters,
               values.index_bound( 1 ) == nb_saved_cells() ) ) ;
   PEL_CHECK_POST(
      IMPLIES( location == AtVertices,
               values.index_bound( 1 ) == nb_saved_vertices() ) ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t ic=0 ; ic<values.index_bound( 0 ) ; ++ic ),
         FORALL( ( size_t im=0 ; im<values.index_bound( 1 ) ; ++im ),
            values( ic, im ) == PDE_ResultSaver::undefined_value() ) ) ) ;
   PEL_CHECK_POST( default_values.size() == nb_components ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t ic=0 ; ic<default_values.size() ; ++ic ),
            default_values( ic ) == PDE_ResultSaver::undefined_value() ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultSaver:: field_saving_in_progress( void ) const
//----------------------------------------------------------------------
{
   return( FIELD_LOCATION != InvalidLocation ) ;
}

//----------------------------------------------------------------------
PDE_ResultSaver::SavingLocation
PDE_ResultSaver:: saving_location( void ) const
//----------------------------------------------------------------------
{
   return( FIELD_LOCATION ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: save_field( doubleArray2D& values,
                              doubleVector const& default_values )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: save_field" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( grid_is_saved() ) ;
   PEL_CHECK_PRE( field_saving_in_progress() ) ;
   PEL_CHECK_PRE( values.index_bound( 0 ) == default_values.size() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( FIELD_LOCATION == AtCellCenters )
   {
      for( size_t ic=0 ; ic<values.index_bound( 0 ) ; ++ic )
         for( size_t j=0 ; j<values.index_bound( 1 ) ; ++j )
            PEL_ASSERT( values( ic, j ) != undefined_value() ) ;
   }
   else
   {
      PEL_ASSERT( FIELD_LOCATION == AtVertices ) ;
      for( size_t ic=0 ; ic<values.index_bound( 0 ) ; ++ic )
      {
         PEL_ASSERT( default_values( ic ) != PEL::bad_double() ) ;
         for( size_t j=0 ; j<values.index_bound( 1 ) ; ++j )
         {
            if( values( ic, j ) == undefined_value() )
            {
               PEL_ASSERT( PERIODIC ) ;
               values( ic, j ) = default_values( ic ) ;
               // default_values so that the min and max values given by the
               // postprocessing tool remain correct
            }
         }
      }      
   }
   
   ++NB_FIELDS ;
   if( MOD_FF->has_module( FIELD_NAME ) )
   {
      PDE_ResultSaver_ERROR::n1( FIELD_NAME ) ;
   }
   
   PEL_Module* fm = PEL_DataOnMeshingWriter::create_field_module(
                                 MOD_FF,  FIELD_NAME,
                                 MESHING_NAME,
                                 location2string(), values ) ;
   MOD_FF->add_module( fm ) ;

   FIELD_LOCATION = InvalidLocation ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !field_saving_in_progress() ) ;
}

//----------------------------------------------------------------------
double
PDE_ResultSaver:: undefined_value( void )
//----------------------------------------------------------------------
{
   return( PEL_DataOnMeshingWriter::undefined_value() ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: check_value_consistency_at_vertex(
                                     std::string const& banner,
                                     GE_Point const* vertex,
                                     double old_val, double new_val,
                                     double dbl_eps, double dbl_min )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: check_value_consistency_at_vertex" ) ;
   PEL_CHECK_PRE( vertex != 0 ) ;
   PEL_CHECK_PRE( dbl_eps > 0. ) ;
   PEL_CHECK_PRE( dbl_min > 0. ) ;
   
   if( old_val != undefined_value() &&
       !PEL::double_equality( old_val, new_val, dbl_eps, dbl_min ) )
   {
      std::ostringstream mesg ;
      if( !banner.empty() )  mesg << banner << std::endl ;
      mesg << std::setprecision( 7 )
           << std::setiosflags( std::ios::scientific ) ;
      mesg << "      the value at vertex " ;
      vertex->print( mesg, 0 ) ; mesg << endl ;
      mesg << "      is different when computed from two cells." << endl ;
      mesg << "         first  value : " << new_val << endl ;
      mesg << "         second value : " << old_val << endl ;
      mesg << "   If it is not an error, do adjust the values of" << endl ;
      mesg << "   the keywords \"dbl_minimum\" and \"dbl_epsilon\"." << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: build_reconstruction( PDE_DiscreteField const* field,
                                        size_t level,
                                        SavingLocation location,
                                        doubleArray2D& values,
                                        doubleVector& default_values)
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: build_reconstruction" ) ;
   PEL_CHECK( field != 0 ) ;
   PEL_CHECK( level < field->storage_depth() ) ;
   PEL_CHECK( has_an_opened_cycle() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( ACTIVE_IDX_OF_VERT.size() != VERTICES->nb_points() )
      PDE_ResultSaver_ERROR::n2() ;
   
   string mesg = "*** PDE_ResultSaver : error saving \""+field->name()+"\"" ;   
   size_t nb_cmps = field->nb_components() ;

   cFE->require_field_calculation( field, PDE_LocalFE::N ) ;
   size_t im = 0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next(), ++im )
   {
      if( im >= NB_ACTIVE_CELLS ) PDE_ResultSaver_ERROR::n2() ;
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;

      if( location == AtCellCenters )
      {
         cFE->set_calculation_point( poly->center() ) ;
         for( size_t ic=0 ; ic<nb_cmps ; ++ic )
         {
            if( cFE->nb_local_nodes( field ) != 0 )
            {
               values( ic, im ) = cFE->value_at_pt( field, level, ic ) ;
            }
            else
            {
               values( ic, im ) = 0.0 ;
            }
         }
      }
      else
      {
         PEL_ASSERT( location == AtVertices ) ;
         for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
         {
            GE_Point const* pt = poly->vertex( i ) ;
            PEL_CHECK( VERTICES->has( pt ) ) ;
            size_t const iv = ACTIVE_IDX_OF_VERT( VERTICES->index( pt ) ) ;
            if( iv >= NB_ACTIVE_VERTS) PDE_ResultSaver_ERROR::n2() ;
            cFE->set_calculation_point( pt ) ;
            for( size_t ic=0 ; ic<nb_cmps ; ic++ )
            {
               double val = 0.0 ;
               if( cFE->nb_local_nodes( field ) != 0 )
               {
                  val = cFE->value_at_pt( field, level, ic ) ;
                  if( val == PEL::bad_double() )
                  {
                     PEL::out() << "=> \"" << field->name() << "\"" << endl ;
                     PEL::out() << "------------------------------" << endl ;
                     cFE->print_current_mesh( PEL::out(), 3 ) ;
                     PEL::out() << "------------------------------" << endl ;
                     pt->print( PEL::out(), 3 ) ;
                  }
                  PEL_ASSERT( val != PEL::bad_double() ) ;
               }

               // the value at a vertex obtained from different meshes
               // should be the same
               check_value_consistency_at_vertex( mesg, pt, values(ic,iv), val, 
                                                  EPS_DBL, MIN_DBL ) ;
               
               if( ( default_values( ic ) == PEL::bad_double() ) ||
                   ( val > default_values( ic ) ) )
               {
                  default_values( ic ) = val ;
               }
               
               values( ic, iv ) = val ;
            }
         }
      }
   }
   
   PEL_CHECK( values.index_bound(0) == field->nb_components() ) ;
   PEL_CHECK( IMPLIES( location == AtCellCenters,
                       values.index_bound(1) == nb_saved_cells() ) ) ;
   PEL_CHECK( IMPLIES( location == AtVertices,
                       values.index_bound(1) == nb_saved_vertices() ) ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: set_mesh_vertices( GE_Mpolyhedron const* poly,
                                     size_t imesh,
                                     intVector& mesh_nb_vertices,
                                     intArray2D& mesh2vertex ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: set_mesh_vertices" ) ;
   
   size_t nvert = poly->nb_vertices() ;
   if( nvert > mesh2vertex.index_bound( 0 ) ) 
      mesh2vertex.raise_first_index_bound( nvert ) ;
   mesh_nb_vertices( imesh ) = nvert ;

   for( size_t i=0 ; i<nvert ; ++i ) 
   {
      GE_Point const* pt = poly->vertex( i ) ;
      size_t ii = ACTIVE_IDX_OF_VERT( VERTICES->index( pt ) ) ;
      mesh2vertex( i, imesh ) = ii ;
   }
}

//----------------------------------------------------------------------
void
PDE_ResultSaver:: infer_active_grid_items( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: infer_active_grid_items" ) ;
   
   NB_ACTIVE_VERTS = 0 ;
   NB_ACTIVE_CELLS = 0 ;
   NB_ACTIVE_FACES = 0 ;
   
   ACTIVE_IDX_OF_VERT.re_initialize( VERTICES->nb_points() ) ;
   ACTIVE_IDX_OF_VERT.set( PEL::bad_index() ) ;
   
   ACTIVE_IDX_OF_CELL.re_initialize( cFE->nb_meshes() ) ;
   ACTIVE_IDX_OF_VERT.set( PEL::bad_index() ) ;
   
   if( DETECT_ACTIVE_ITEMS )
   {
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         size_t id = cFE->mesh_id() ;
         if( id > ACTIVE_IDX_OF_CELL.size() ) ACTIVE_IDX_OF_CELL.resize( id+1) ;
         ACTIVE_IDX_OF_CELL( id ) = NB_ACTIVE_CELLS ;
         ++NB_ACTIVE_CELLS ;
         update_active_vertices( VERTICES, cFE->polyhedron(),
                                 ACTIVE_IDX_OF_VERT, NB_ACTIVE_VERTS ) ;
      }
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      { 
         ++NB_ACTIVE_FACES ;
         if( sFE->is_periodic() )
         {
            PERIODIC = true ;
            update_active_vertices( VERTICES, sFE->polyhedron( 0 ),
                                    ACTIVE_IDX_OF_VERT, NB_ACTIVE_VERTS ) ;
            update_active_vertices( VERTICES, sFE->polyhedron( 1 ),
                                    ACTIVE_IDX_OF_VERT, NB_ACTIVE_VERTS ) ;
            ++NB_ACTIVE_FACES ;
         }
      }
      for( bFE->start(); bFE->is_valid() ; bFE->go_next() )
      { 
         ++NB_ACTIVE_FACES ;
      }
   }
   else
   {
      NB_ACTIVE_VERTS = VERTICES->nb_points() ;
      for( size_t iv=0 ; iv<NB_ACTIVE_VERTS ; ++iv )
      {
         ACTIVE_IDX_OF_VERT( iv ) = iv ;
      }
      NB_ACTIVE_CELLS = cFE->nb_meshes() ;
      for( size_t ic=0 ; ic<NB_ACTIVE_CELLS ; ++ic )
      {
         ACTIVE_IDX_OF_CELL( ic ) = ic ;
      }
      NB_ACTIVE_FACES = bFE->nb_meshes() ;
      //??? for each periodic side, there are two faces
      //??? PDE_CursorFEside may return the number of periodic sides
      for( sFE->start(); sFE->is_valid() ; sFE->go_next() )
      { 
         ++NB_ACTIVE_FACES ;
         if( sFE->is_periodic() ) ++NB_ACTIVE_FACES ;
      }
   }
}

//----------------------------------------------------------------------
std::string const&
PDE_ResultSaver:: location2string( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultSaver:: location2string" ) ;
   
   static std::string at_vertices = "at_vertices" ;
   static std::string at_cell_centers = "at_cell_centers" ;
   
   if( FIELD_LOCATION == AtVertices )
   {
      return( at_vertices ) ;
   }
   else
   {
      PEL_ASSERT( FIELD_LOCATION == AtCellCenters ) ;
      return( at_cell_centers ) ;
   }
}

//internal--------------------------------------------------------------
void 
PDE_ResultSaver_ERROR:: n0( std::string const& name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** PDE_ResultSaver error:\n"
      "    several variables are saved with the same name \""+name+"\"" ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ResultSaver_ERROR:: n1( std::string const& name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** PDE_ResultSaver error:\n"
      "    several discrete fields are saved with the same name \""+name+"\"" ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ResultSaver_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PDE_ResultSaver error:" << endl ;
   msg << "    the grid has changed since the last call to save_grid()" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ResultSaver_ERROR:: n3( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PDE_ResultSaver error:" << endl ;
   msg << "    halo color is expected for \"averaged_at_vertices\" saving" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
