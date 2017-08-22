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

#include <PDE_CursorFEside.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <GE_Color.hh>
#include <GE_Customized_QR.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_CustomizedQR_provider.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Transform.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_GridFE.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_FaceFE.hh>
#include <PDE_HelperFIC.hh>

#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ostringstream ;

struct PDE_CursorFEside_ERROR
{
   static void n0( GE_Mpolyhedron const* poly ) ;
} ;

//----------------------------------------------------------------------
PDE_CursorFEside:: PDE_CursorFEside( PEL_Object* a_owner,
                                     PDE_GridFE const* a_grid,
                                     PDE_LocalFEcell* cfe_0,
                                     PDE_LocalFEcell* cfe_1 )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , GRID( a_grid )
   , GEO_STATE( a_grid->geometric_state() )
   , EXCLUDE_COLORS( PEL_Vector::create( this, 0 ) )
   , SIDES( a_grid->sides() )
   , IT( PEL_VectorIterator::create( this, a_grid->sides() ) )
   , SIDE( 0 )
   , OK( false )
   , CELL_FE_0( cfe_0 )
   , CELL_FE_1( cfe_1 )
   , NORMAL0( GE_Vector::create( this, a_grid->nb_space_dimensions() ) )
   , NORMAL1( GE_Vector::create( this, a_grid->nb_space_dimensions() ) )
   , DUMMY_VEC( GE_Vector::create( this, a_grid->nb_space_dimensions() ) )
   , QRP0( GE_CustomizedQR_provider::create( this ) )
   , QRP1( GE_CustomizedQR_provider::create( this ) )
   , DF_ROW( 0 )
   , DF_COL( 0 )
   , ROW_NODES( 0 )
   , COL_NODES( 0 )
   , ROW_LOC0( 0 )
   , ROW_LOC1( 0 )
   , COL_LOC0( 0 )
   , COL_LOC1( 0 )
   , SAME_ROW_AND_COL( false )
   , DIST( 0, 2 )
   , HELPER_FIC( 0 )
{
   PEL_LABEL( "PDE_CursorFEside:: PDE_CursorFEside" ) ;

   CELL_FE_0->set_owner( this ) ;
   CELL_FE_1->set_owner( this ) ;

   GRID->set_of_vertices()->attach_observer( this ) ;
   

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_CursorFEside:: ~PDE_CursorFEside( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: ~PDE_CursorFEside" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GRID->set_of_vertices()->detach_observer( this ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: update( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: update" ) ;

   DIST.re_initialize( 0, 2 ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( GRID->nb_space_dimensions() ) ;
}

//----------------------------------------------------------------------
PDE_LocalFEcell const*
PDE_CursorFEside:: adjacent_localFEcell( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: adjacent_localFEcell" ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;

   PDE_LocalFEcell const* result = ( i_adj==0 ? CELL_FE_0 : CELL_FE_1 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: require_field_calculation( PDE_DiscreteField const* ff,
                                              int order )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: require_field_calculation" ) ;
   PEL_CHECK_PRE( ff!=0 ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( order == PDE_LocalFE::node  ||
                  order == PDE_LocalFE::N  || 
                  order == PDE_LocalFE::dN || 
                  order == PDE_LocalFE::d2N ) ;

   CELL_FE_0->require_field_calculation( ff, order ) ;
   CELL_FE_1->require_field_calculation( ff, order ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( field_is_handled( ff ) ) ;
   PEL_CHECK_POST( field_calculation_is_handled( ff, order ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         adjacent_localFEcell( i_adj )->field_is_handled( ff ) &&
         adjacent_localFEcell( i_adj )->field_calculation_is_handled(
                                                           ff, order ) ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: field_calculation_is_handled( PDE_DiscreteField const* ff,
                                                 int order ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: field_calculation_is_handled" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( order == PDE_LocalFE::node  ||
                  order == PDE_LocalFE::N  || 
                  order == PDE_LocalFE::dN || 
                  order == PDE_LocalFE::d2N ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = CELL_FE_0->field_calculation_is_handled( ff, order ) ;
   PEL_CHECK( result == CELL_FE_1->field_calculation_is_handled( ff, order ) ) ;

   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         EQUIVALENT(
            result,
            adjacent_localFEcell( i_adj )->field_calculation_is_handled(
                                                         ff, order ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: field_is_handled( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: field_is_handled" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = CELL_FE_0->field_is_handled( ff ) ;
   PEL_CHECK( result == CELL_FE_1->field_is_handled( ff ) ) ;

   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         EQUIVALENT(
            result,
            adjacent_localFEcell( i_adj )->field_is_handled( ff ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: nb_handled_fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: PDE_CursorFEside:: nb_handled_fields" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t result = CELL_FE_0->nb_handled_fields() ;
   PEL_CHECK( result == CELL_FE_1->nb_handled_fields() ) ;

   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         result==adjacent_localFEcell( i_adj )->nb_handled_fields() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_CursorFEside:: handled_field( size_t iF ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: handled_field" ) ;
   PEL_CHECK_PRE( iF < nb_handled_fields() ) ;
   
   PDE_DiscreteField const* result = CELL_FE_0->handled_field( iF ) ;
   PEL_CHECK( result == CELL_FE_1->handled_field( iF ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         result==adjacent_localFEcell( i_adj )->handled_field( iF ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: nb_meshes( void ) const
//----------------------------------------------------------------------
{
   return( GRID->nb_sides() ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: exclude_color( GE_Color const* a_color ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: exclude_color" ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   PEL_CHECK_PRE( !is_excluded( a_color ) ) ;
   
   EXCLUDE_COLORS->append( const_cast<GE_Color*>( a_color ) ) ;
   
   PEL_CHECK_POST( is_excluded( a_color ) ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: include_color( GE_Color const* a_color ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: include_color" ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   PEL_CHECK_PRE( is_excluded( a_color ) ) ;
   
   EXCLUDE_COLORS->remove_at( EXCLUDE_COLORS->index_of( a_color ) ) ;
   
   PEL_CHECK_POST( !is_excluded( a_color ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: is_excluded( GE_Color const* a_color ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: is_excluded" ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   
   bool result = EXCLUDE_COLORS->has( a_color ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: go_i_th( size_t an_id_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: go_i_th" ) ;
   PEL_CHECK_PRE( an_id_mesh==PEL::bad_index() || an_id_mesh<nb_meshes() ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, SIDES ) ;
   }
   GEO_STATE = GRID->geometric_state() ;

   if( an_id_mesh != PEL::bad_index() )
   {
      IT->go_i_th( an_id_mesh ) ;
      OK = IT->is_valid() && !rejected_side() ;
   }
   else
   {
      OK = false ;
   }
   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( is_valid(), mesh_id() == an_id_mesh ) ) ;
   PEL_CHECK_POST( !valid_IP() ) ;
   PEL_CHECK_POST( calculation_point() == 0 ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside::start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, SIDES ) ;
   }
   GEO_STATE = GRID->geometric_state() ;

   IT->start() ;
   while( IT->is_valid() && rejected_side() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !valid_IP() ) ;
   PEL_CHECK_POST( calculation_point() == 0 ) ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: is_valid" ) ;
   
   bool result = OK && ( GEO_STATE == GRID->geometric_state() ) ;
   PEL_CHECK( OK==adjacent_localFEcell(0)->is_valid() ) ;
   PEL_CHECK( OK==adjacent_localFEcell(1)->is_valid() ) ;

   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         EQUIVALENT( result,
                     adjacent_localFEcell( i_adj )->is_valid() ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside::go_next" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   if( IT->container_has_been_modified() || 
       GEO_STATE != GRID->geometric_state() )
   {
      ostringstream mesg ;
      mesg << "PDE_CursorFEside : " << endl ;
      mesg << "   attempt to call \"go_next\" after" << endl ;
      mesg << "   a modification of the meshing" ;
      PEL_Error::object()->raise_internal( mesg.str() ) ;
   }

   IT->go_next() ;
   while( IT->is_valid() && rejected_side() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !valid_IP() ) ;
   PEL_CHECK_POST( calculation_point() == 0 ) ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: is_periodic( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: is_periodic" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = SIDE->is_periodic() ;

   return( result ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_CursorFEside:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: polyhedron" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = SIDE->polyhedron() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->dimension() == nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_CursorFEside:: polyhedron( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: polyhedron" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = SIDE->polyhedron() ;
   if( i_adj == 1 )
   {
      result = SIDE->periodic_neighbour()->polyhedron() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->dimension() == nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Transform const*
PDE_CursorFEside:: periodic_transform( size_t i_adj_source ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: periodic_transform" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj_source==0 || i_adj_source==1 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Transform const* result = SIDE->periodic_transform() ;
   if( i_adj_source == 1 )
   {
      result = SIDE->periodic_neighbour()->periodic_transform() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->inverse() == 
                   periodic_transform( 1 - i_adj_source ) ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: refinement_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( SIDE->refinement_level() ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_CursorFEside:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: color" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = SIDE->color() ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_CursorFEside:: color( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: color" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = 
                   ( i_adj==0 ? 
                     SIDE->color() : SIDE->periodic_neighbour()->color() ) ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: mesh_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: mesh_id" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = SIDE->side_id() ;

   PEL_CHECK_POST( result < nb_meshes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector const*
PDE_CursorFEside:: normal( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: normal" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Vector const* result = NORMAL0 ;

   PEL_CHECK_POST( result !=0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_components() == nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector const*
PDE_CursorFEside:: normal( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: normal" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Vector const* result = ( i_adj == 0 ? NORMAL0 : NORMAL1 ) ;

   PEL_CHECK_POST( result !=0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_components() == nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_CursorFEside:: distance_to_adjacent_finite_volume_center(
                                                    size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: distance_to_adjacent_finite_volume_center" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( DIST.index_bound( 0 ) < mesh_id()+1 )
   {
      size_t const new_size = mesh_id()+1 ;
      size_t const old_size = DIST.index_bound( 0 ) ;
      DIST.raise_first_index_bound( new_size ) ;
      for( size_t i=old_size ; i<new_size ; ++i )
      {
         DIST( i, 0 ) = PEL::bad_double() ;
         DIST( i, 1 ) = PEL::bad_double() ;
      }
   }
   
   double result = DIST( mesh_id(), i_adj ) ;
   if( result == PEL::bad_double() )
   {
      compute_current_distances() ;
      result = DIST( mesh_id(), i_adj ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_HelperFIC*
PDE_CursorFEside:: helper_FIC( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: helper_FIC" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   if( HELPER_FIC == 0 )
   {
      HELPER_FIC = new PDE_HelperFIC( const_cast<PDE_CursorFEside*>(this) ) ;
   }
   PDE_HelperFIC* result = HELPER_FIC ;
   
   result->set_side( SIDE,
                     adjacent_localFEcell( 0 )->mesh_id(),
                     adjacent_localFEcell( 1 )->mesh_id() ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->current_side_id() == mesh_id() ) ;
   PEL_CHECK_POST( !result->ready_for_interpolation() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: set_calculation_point( GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: set_calculation_point" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( polyhedron()->contains( pt ) ) ;

   CELL_FE_0->set_calculation_point( pt ) ;
   CELL_FE_1->set_calculation_point( pt ) ;
   
   PEL_CHECK_POST( calculation_point() != 0 ) ;
   PEL_CHECK_POST( calculation_point() != pt ) ;
   PEL_CHECK_POST( calculation_point()->nb_coordinates() ==
                                                  pt->nb_coordinates() ) ;
   PEL_CHECK_POST( 
      FORALL(
         ( size_t i=0 ; i<pt->nb_coordinates() ; ++i ),
         calculation_point()->coordinate( i ) == pt->coordinate( i ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         adjacent_localFEcell( i_adj )->calculation_point() != 0 ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         FORALL(
            ( size_t i=0 ; i<pt->nb_coordinates() ; ++i ),
            adjacent_localFEcell( i_adj )->calculation_point()->coordinate( i )
               == pt->coordinate( i ) ) ) ) ;     
}

//----------------------------------------------------------------------
GE_Point const*
PDE_CursorFEside:: calculation_point( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: calculation_point" ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;

   GE_Point const* result = CELL_FE_0->calculation_point() ;
   
   PEL_CHECK_POST( IMPLIES( result!=0, result->is_under_ownership_of( this ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, polyhedron()->contains( result ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0,
                            result->nb_coordinates()==nb_space_dimensions() ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         EQUIVALENT(
            result!=0,
            adjacent_localFEcell( i_adj )->calculation_point()!=0 ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: set_calculation_points( GE_Point const* pt_0,
                                           GE_Point const* pt_1 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: set_calculation_points" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( pt_0 != 0 ) ;
   PEL_CHECK_PRE( pt_0->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( polyhedron(0)->contains( pt_0 ) ) ;
   PEL_CHECK_PRE( pt_1 != 0 ) ;
   PEL_CHECK_PRE( pt_1->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( polyhedron(1)->contains( pt_1 ) ) ;

   CELL_FE_0->set_calculation_point( pt_0 ) ;
   CELL_FE_1->set_calculation_point( pt_1 ) ;

   PEL_CHECK_POST( calculation_point(0) != 0 ) ;
   PEL_CHECK_POST( calculation_point(0) != pt_0 ) ;
   PEL_CHECK_POST( calculation_point(0)->nb_coordinates() ==
                                                  pt_0->nb_coordinates() ) ;
   PEL_CHECK_POST( 
      FORALL(
         ( size_t i=0 ; i<pt_0->nb_coordinates() ; ++i ),
         calculation_point(0)->coordinate( i ) == pt_0->coordinate( i ) ) ) ;
   PEL_CHECK_POST( calculation_point(1) != 0 ) ;
   PEL_CHECK_POST( calculation_point(1) != pt_1 ) ;
   PEL_CHECK_POST( calculation_point(1)->nb_coordinates() ==
                                                  pt_1->nb_coordinates() ) ;
   PEL_CHECK_POST( 
      FORALL(
         ( size_t i=0 ; i<pt_1->nb_coordinates() ; ++i ),
         calculation_point(1)->coordinate( i ) == pt_1->coordinate( i ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         adjacent_localFEcell( i_adj )->calculation_point() != 0 ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<pt_0->nb_coordinates() ; ++i ),
         adjacent_localFEcell( 0 )->calculation_point()->coordinate( i )
               == pt_0->coordinate( i ) ) ) ;  
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<pt_1->nb_coordinates() ; ++i ),
         adjacent_localFEcell( 1 )->calculation_point()->coordinate( i )
               == pt_1->coordinate( i ) ) ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_CursorFEside:: calculation_point( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: calculation_point" ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj == 0 || i_adj == 1 ) ;

   GE_Point const* result = ( i_adj == 0 ? CELL_FE_0->calculation_point()
                                         : CELL_FE_1->calculation_point() ) ;
   
   PEL_CHECK_POST( IMPLIES( result!=0, result->is_under_ownership_of( this ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, polyhedron(i_adj)->contains( result ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0,
                            result->nb_coordinates()==nb_space_dimensions() ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<2 ; ++i ),
         EQUIVALENT(
            result!=0,
            adjacent_localFEcell( i )->calculation_point()!=0 ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: set_row_and_col_fields(
                                   PDE_DiscreteField const* row_field, 
                                   PDE_DiscreteField const* col_field )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: set_row_and_col_fields" ) ;
   PEL_CHECK_PRE( row_field !=0 && field_is_handled( row_field ) ) ;
   PEL_CHECK_PRE( col_field !=0 && field_is_handled( col_field ) ) ;  
   PEL_CHECK_PRE( is_valid() && !valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CELL_FE_0->set_row_and_col_fields( row_field, col_field ) ;
   CELL_FE_1->set_row_and_col_fields( row_field, col_field ) ;
   
   DF_ROW = row_field ;
   DF_COL = col_field ;

   // Row nodes and local connectivity :
   determine_local( ROW_NODES,
                    CELL_FE_0->row_field_node_connectivity(),
                    CELL_FE_1->row_field_node_connectivity(),
                    ROW_LOC0, ROW_LOC1 ) ;
   
   // Col nodes and local connectivity :
   if( col_field==row_field )
   {
      SAME_ROW_AND_COL = true ;      
   }
   else
   {
      SAME_ROW_AND_COL = false ;
      determine_local( COL_NODES,
                       CELL_FE_0->col_field_node_connectivity(),
                       CELL_FE_1->col_field_node_connectivity(),
                       COL_LOC0, COL_LOC1 ) ;
   }

   PEL_CHECK_POST( field( PDE_LocalFE::row ) == row_field ) ;
   PEL_CHECK_POST( field( PDE_LocalFE::col ) == col_field ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         adjacent_localFEcell( i_adj )->field( PDE_LocalFE::row ) == row_field &&
         adjacent_localFEcell( i_adj )->field( PDE_LocalFE::col ) == col_field ) ) ;
}


//----------------------------------------------------------------------
void
PDE_CursorFEside:: determine_local( size_t_vector& n,
                                    size_t_vector const& n0,
                                    size_t_vector const& n1,
                                    size_t_vector& loc0,
                                    size_t_vector& loc1 ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: determine_local" ) ;
   
   size_t nb = n0.size() ;
   for( size_t i=0 ; i<n1.size() ; i++ )
      if( ! n0.has( n1(i) ) )
          nb++ ;
   
   PEL_CHECK( n.size() == loc0.size() ) ;
   PEL_CHECK( n.size() == loc1.size() ) ;

   if( n.size() != nb ) 
   {
      n.re_initialize( nb ) ;
      loc0.re_initialize( nb ) ;
      loc1.re_initialize( nb ) ;
   }
   
   size_t idx = 0 ;
   
   for( size_t i=0 ; i<n0.size() ; ++i )
   {
      n(idx) = n0(i) ;
      loc0(idx) = i ;
      loc1(idx) = PEL::bad_index() ;
      idx++ ;
   }
   for( size_t i=0 ; i<n1.size() ; ++i )
   {
      size_t index = n0.index_of( n1(i) ) ;
      if( index<n0.size() )
      {
         loc1(index) = i ;
      }
      else
      {
         n(idx) = n1(i) ;
         loc0(idx) = PEL::bad_index() ;
         loc1(idx) = i ;
         idx++ ;
      }
   }
   PEL_CHECK( idx==nb ) ;
}


//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_CursorFEside:: field( PDE_LocalFE::field_id sf ) const
//----------------------------------------------------------------------
{
   return( sf==PDE_LocalFE::row ? DF_ROW : DF_COL ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_CursorFEside:: row_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: row_field_node_connectivity" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t_vector const& result = ROW_NODES ;
   
   PEL_CHECK_POST( result.size()==nb_basis_functions( PDE_LocalFE::row ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_CursorFEside:: row_field_node_side_2_cell_index( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: row_field_node_side_2_cell_index" ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t_vector const& result = ( i_adj==0 ? ROW_LOC0 : ROW_LOC1 ) ;
   
   PEL_CHECK_POST( result.size()==nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result.size() ; ++i ),
         result(i)==PEL::bad_index() ||
         adjacent_localFEcell(i_adj)->row_field_node_connectivity()(result(i))
                                     ==row_field_node_connectivity()(i) ) ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_CursorFEside:: col_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: col_field_node_connectivity" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t_vector const& result = ( SAME_ROW_AND_COL ? ROW_NODES : COL_NODES ) ;
   
   PEL_CHECK_POST( result.size()==nb_basis_functions( PDE_LocalFE::col ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_CursorFEside:: col_field_node_side_2_cell_index( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: col_field_node_side_2_cell_index" ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t_vector const& result = ( SAME_ROW_AND_COL ?
                                   ( i_adj==0 ? ROW_LOC0 : ROW_LOC1 ) :
                                   ( i_adj==0 ? COL_LOC0 : COL_LOC1 ) ) ;
   
   PEL_CHECK_POST( result.size()==nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result.size() ; ++i ),
         result(i)==PEL::bad_index() ||
         adjacent_localFEcell(i_adj)->col_field_node_connectivity()(result(i))
                                     ==col_field_node_connectivity()(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CursorFEside:: nb_basis_functions( PDE_LocalFE::field_id sf ) const
//----------------------------------------------------------------------
{
   return( SAME_ROW_AND_COL || sf==PDE_LocalFE::row ? ROW_NODES.size() : COL_NODES.size() ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: start_IP_iterator( GE_QRprovider const* qrp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: start_IP_iterator" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( qrp != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   customize_providers( qrp ) ;
   CELL_FE_0->start_IP_iterator( QRP0 ) ;
   CELL_FE_1->start_IP_iterator( QRP1 ) ;

   PEL_CHECK_POST( valid_IP() ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: go_next_IP( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: go_next_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   CELL_FE_0->go_next_IP() ;
   CELL_FE_1->go_next_IP() ;
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: valid_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: valid_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = CELL_FE_0->valid_IP() ;
   PEL_CHECK( EQUIVALENT( result, CELL_FE_1->valid_IP() ) ) ;

   PEL_CHECK_POST(
             EQUIVALENT( result, adjacent_localFEcell( 0 )->valid_IP() ) ) ;
   PEL_CHECK_POST(
             EQUIVALENT( result, adjacent_localFEcell( 1 )->valid_IP() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_CursorFEside:: weight_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: weight_of_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_PRE( valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double result = CELL_FE_0->weight_of_IP() ;
   
   PEL_CHECK_POST(
      FORALL(
         ( size_t i_adj=0 ; i_adj<2 ; ++i_adj ),
         FORMAL(
            result == adjacent_localFEcell( i_adj )->weight_of_IP() ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_CursorFEside:: weight_of_IP( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: weight_of_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_PRE( valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double result = ( i_adj == 0 ? CELL_FE_0->weight_of_IP() : 
                                  CELL_FE_1->weight_of_IP() ) ;
   
   PEL_CHECK_POST(
         FORMAL(
            result == adjacent_localFEcell( i_adj )->weight_of_IP() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_CursorFEside:: coordinates_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: coordinates_of_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_PRE( valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Point const* result = CELL_FE_0->coordinates_of_IP() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( FORMAL( 
       result->is_equal( adjacent_localFEcell( 0 )->coordinates_of_IP() ) ) ) ;
   PEL_CHECK_POST( FORMAL( 
       result->is_equal( adjacent_localFEcell( 1 )->coordinates_of_IP() ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_CursorFEside:: coordinates_of_IP( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: coordinates_of_IP" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_PRE( valid_IP() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Point const* result = ( i_adj == 0 ? CELL_FE_0->coordinates_of_IP() : 
                                           CELL_FE_1->coordinates_of_IP() ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( FORMAL( 
       result->is_equal( adjacent_localFEcell( 0 )->coordinates_of_IP() ) ) ) ;
   PEL_CHECK_POST( FORMAL( 
       result->is_equal( adjacent_localFEcell( 1 )->coordinates_of_IP() ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: print_current_mesh" ) ;
   PEL_CHECK_PRE( os ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   SIDE->print( os, indent_width ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: compute_current_distances( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: compute_current_distances" ) ;
   PEL_CHECK( DIST.index_bound(0) > mesh_id() ) ;

   size_t ii = mesh_id() ;
   GE_Mpolyhedron const* s0_poly = SIDE->polyhedron() ;
   GE_Mpolyhedron const* s1_poly = ( SIDE->is_periodic() ? 
                         SIDE->periodic_neighbour()->polyhedron() : s0_poly ) ;
   GE_Mpolyhedron const* c0_poly = CELL_FE_0->polyhedron() ;
   GE_Mpolyhedron const* c1_poly = CELL_FE_1->polyhedron() ;

   GE_Point const* pt0 = c0_poly->finite_volume_center() ;
   if( pt0 == 0 ) PDE_CursorFEside_ERROR::n0( c0_poly ) ;
   
   DUMMY_VEC->re_initialize( s0_poly->vertex( 0 ), pt0 ) ;
   DIST( ii, 0 ) = DUMMY_VEC->dot_product( NORMAL0 ) ;

   GE_Point const* pt1 = c1_poly->finite_volume_center() ;
   if( pt1 == 0 ) PDE_CursorFEside_ERROR::n0( c1_poly ) ;
   
   DUMMY_VEC->re_initialize( s1_poly->vertex( 0 ), pt1 ) ;
   if( SIDE->is_periodic() )
   {
      DIST( ii, 1 ) = - DUMMY_VEC->dot_product( NORMAL1 ) ;
   }
   else
   {
      DIST( ii, 1 ) = - DUMMY_VEC->dot_product( NORMAL0 ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: do_internals_on_current_mesh( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: do_internals_on_current_mesh" ) ;

   if( OK )
   {
      SIDE = the_side( IT ) ;

      CELL_FE_0->go_i_th( SIDE->adjacent_cell( 0 )->id_number() ) ;
      if( !SIDE->is_periodic() )
      {
         CELL_FE_1->go_i_th( SIDE->adjacent_cell( 1 )->id_number() ) ;
      }
      else
      {
         PDE_FaceFE const* ss = SIDE->periodic_neighbour() ;
         CELL_FE_1->go_i_th( ss->adjacent_cell( 0 )->id_number() ) ;         
      }

      NORMAL0->set( SIDE->polyhedron()->unit_normal() ) ;
      DUMMY_VEC->re_initialize( SIDE->polyhedron()->center(),
                                CELL_FE_0->polyhedron()->center() ) ;
      if( NORMAL0->dot_product( DUMMY_VEC )<0. ) 
      {
         NORMAL0->scale( -1. ) ;
      }
      
      if( SIDE->is_periodic() )
      {
         PDE_FaceFE const* ss = SIDE->periodic_neighbour() ;
         NORMAL1->set( ss->polyhedron()->unit_normal() ) ;
         DUMMY_VEC->re_initialize( CELL_FE_1->polyhedron()->center(),
                                   ss->polyhedron()->center() ) ;
         if( NORMAL1->dot_product( DUMMY_VEC )<0. ) 
         {
            NORMAL1->scale( -1. ) ;
         }
      }
   }
   else
   {
      SIDE = 0 ;
      CELL_FE_0->go_i_th( PEL::bad_index() ) ;
      CELL_FE_1->go_i_th( PEL::bad_index() ) ;
   }
}

//----------------------------------------------------------------------
bool
PDE_CursorFEside:: rejected_side( void ) const
//----------------------------------------------------------------------
{
   bool result = is_excluded( the_side( IT )->color() )
              || !the_side( IT )->is_active() ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_FaceFE const*
PDE_CursorFEside:: the_side( PEL_VectorIterator const* it )
//----------------------------------------------------------------------
{
   return( static_cast<PDE_FaceFE*>( it->item() ) ) ;
}

//----------------------------------------------------------------------
void
PDE_CursorFEside:: customize_providers( GE_QRprovider const* qrp ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CursorFEside:: customize_providers" ) ;

   GE_Mpolyhedron const* c0 = CELL_FE_0->polyhedron() ;
   GE_Customized_QR* qr0 =
      QRP0->quadrature_rule_to_be_customized( c0->reference_polyhedron() )  ;

   GE_Mpolyhedron const* c1 = CELL_FE_1->polyhedron() ;
   GE_Customized_QR* qr1 =
      QRP1->quadrature_rule_to_be_customized( c1->reference_polyhedron() )  ;

   GE_Mpolyhedron const* s0 = SIDE->polyhedron() ;
   GE_Mpolyhedron const* s1 = ( SIDE->is_periodic() ? 
                              SIDE->periodic_neighbour()->polyhedron() : s0 ) ;
   
   GE_ReferencePolyhedron const* rpoly = s0->reference_polyhedron() ;
   PEL_ASSERT( rpoly == s1->reference_polyhedron() ) ;
   GE_QuadratureRule const* qrs = qrp->quadrature_rule( rpoly ) ;

   double const de0 = s0->measure()/rpoly->measure() ;
   double const de1 = s1->measure()/rpoly->measure() ;
   double const dc0 = c0->measure()/c0->reference_polyhedron()->measure() ;
   double const dc1 = c1->measure()/c1->reference_polyhedron()->measure() ;
   double const w0 = de0/dc0 ;
   double const w1 = de1/dc1 ;
   
   GE_Point* pt_real = GE_Point::create( 0, nb_space_dimensions() ) ;
   GE_Point* c_pt_ref  = GE_Point::create( 0, nb_space_dimensions() ) ;
   for( size_t i=0 ; i<qrs->nb_points() ; ++i )
   {
      GE_Point const* s_pt_ref = qrs->point( i ) ;
      s0->apply_mapping( s_pt_ref, pt_real ) ;
      c0->apply_inverse_mapping( pt_real, c_pt_ref ) ;
      qr0->insert_point( c_pt_ref, w0*qrs->weight(i) ) ;
      if( s1 != s0 )
      {
         s1->apply_mapping( s_pt_ref, pt_real ) ;
      }
      c1->apply_inverse_mapping( pt_real, c_pt_ref ) ;
      qr1->insert_point( c_pt_ref, w1*qrs->weight(i) ) ;
   }
   qr0->finalize() ;
   qr1->finalize() ;
   
   pt_real->destroy() ; pt_real = 0 ;     //??? a chaque fois
   c_pt_ref->destroy() ; c_pt_ref = 0 ;   //??? a chaque fois
}

//------------------------------------------------------------------------
bool
PDE_CursorFEside:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
PDE_CursorFEside_ERROR:: n0( GE_Mpolyhedron const* poly )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_CursorFEside:" << endl << endl ;
   msg << "    the \"finite volume center\" has not been defined" << endl ;
   msg << "    for the following polyhedron:" << endl << endl ;
   poly->print( msg, 4 ) ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
