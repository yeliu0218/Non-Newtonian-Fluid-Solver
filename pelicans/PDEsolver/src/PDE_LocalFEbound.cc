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

#include <PDE_LocalFEbound.hh>

#include <GE_Customized_QR.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_SetOfPoints.hh>
#include <GE_SimplePolygon2D.hh>
#include <GE_Vector.hh>

#include <GE_Matrix.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_GridFE.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <iostream>
#include <sstream>
#include <string>

using std::cout ; using std::endl ;
using std::ostringstream ;

struct PDE_LocalFEbound_ERROR
{
   static void n0( GE_Mpolyhedron const* poly ) ;
} ;

//----------------------------------------------------------------------
PDE_LocalFEbound:: PDE_LocalFEbound( PEL_Object* a_owner,
                                     PDE_GridFE const* a_grid )
//----------------------------------------------------------------------
   : PDE_LocalFEmulti( a_owner, a_grid->nb_space_dimensions() )
   , BOUNDS( a_grid->bounds() )
   , IT( PEL_VectorIterator::create( this, a_grid->bounds() ) )
   , BOUND( 0 )
   , OK( false )
   , GRID( a_grid )
   , GEO_STATE( a_grid->geometric_state() )
   , PT( GE_Point::create( this, a_grid->nb_space_dimensions() ) )
   , BD_PT_REF( GE_Point::create( this, a_grid->nb_space_dimensions()-1 ) )
   , CELL_PT_REF( GE_Point::create( this, a_grid->nb_space_dimensions() ) )
   , CELL_tr_JAC( GE_Matrix::create( this,
                                     a_grid->nb_space_dimensions(),
                                     a_grid->nb_space_dimensions() ) )
   , HESSIAN( new doubleArray3D( a_grid->nb_space_dimensions(), 
                                      a_grid->nb_space_dimensions(), 
                                      a_grid->nb_space_dimensions() ) )
   , DUMMY_VEC( GE_Vector::create( this, a_grid->nb_space_dimensions() ) )
   , DIST( 0 )
     
{
   GRID->set_of_vertices()->attach_observer( this ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_LocalFEbound:: ~PDE_LocalFEbound( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: ~PDE_LocalFEbound" ) ;

   GRID->set_of_vertices()->detach_observer( this ) ;
   delete HESSIAN ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEbound:: update( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: update" ) ;

   DIST.re_initialize( 0 ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEbound:: nb_meshes( void ) const
//----------------------------------------------------------------------
{
   return( GRID->nb_bounds() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEbound:: go_i_th( size_t an_id_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: go_i_th" ) ;
   PEL_CHECK_PRE( go_i_th_PRE( an_id_mesh ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, BOUNDS ) ;
   }
   GEO_STATE = GRID->geometric_state() ;
   if( an_id_mesh != PEL::bad_index() )
   {
      IT->go_i_th( an_id_mesh ) ;
      OK = IT->is_valid() && !rejected_bound() ;
   }
   else
   {
      OK = false ;
   }

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( go_i_th_POST( an_id_mesh ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEbound:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound::start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, BOUNDS ) ;
   }
   GEO_STATE = GRID->geometric_state() ;

   IT->start() ;
   while( IT->is_valid() && rejected_bound() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( start_POST() ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEbound:: is_valid( void ) const
//----------------------------------------------------------------------
{
   bool result = OK && ( GEO_STATE == GRID->geometric_state() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEbound:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound::go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( IT->container_has_been_modified() || 
       GEO_STATE != GRID->geometric_state() )
   {
      ostringstream mesg ;
      mesg << "PDE_LocalFEbound : " << endl ;
      mesg << "   attempt to call \"go_next\" after" << endl ;
      mesg << "   a modification of the meshing" ;
      PEL_Error::object()->raise_internal( mesg.str() ) ;
   }

   IT->go_next() ;
   while( IT->is_valid() && rejected_bound() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( go_next_POST() ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_LocalFEbound:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: polyhedron" ) ;
   PEL_CHECK_PRE( polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = BOUND->polyhedron() ;

   PEL_CHECK_POST( polyhedron_POST( result ) ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEbound:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: refinement_level" ) ;
   PEL_CHECK_PRE( refinement_level_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( BOUND->refinement_level() ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_LocalFEbound:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: color" ) ;
   PEL_CHECK_PRE( color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = BOUND->color() ;

   PEL_CHECK_POST( color_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEbound:: mesh_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: mesh_id" ) ;
   PEL_CHECK_PRE( mesh_id_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = BOUND->id_number() ;

   PEL_CHECK_POST( mesh_id_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector const*
PDE_LocalFEbound:: outward_normal( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: outward_normal" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Vector const* result = BOUND->unit_outward_normal() ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( PEL::toler( result->norm()-1.0 ) ) ;
   PEL_CHECK_POST( result->nb_components()==nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_LocalFEbound:: adjacent_cell_polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: adjacent_cell_polyhedron" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = BOUND->adjacent_cell()->polyhedron() ;

   PEL_CHECK_POST( polyhedron_POST( result ) ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEbound:: adjacent_cell_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: adjacent_cell_id" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = BOUND->adjacent_cell()->id_number() ;

   return( result ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_LocalFEbound:: adjacent_cell_color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: adjacent_cell_color" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = BOUND->adjacent_cell()->color() ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEbound:: distance_to_adjacent_finite_volume_center( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: distance_to_adjacent_finite_volume_center" );
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( DIST.size() < mesh_id()+1 )
   {
      size_t const new_size = mesh_id()+1 ;
      size_t const old_size = DIST.size() ;
      DIST.resize( new_size ) ;
      for( size_t i=old_size ; i<new_size ; ++i )
      {
         DIST( i ) = PEL::bad_double() ;
      }
   }
   
   double result = DIST( mesh_id() ) ;
   if( result == PEL::bad_double() )
   {
      compute_current_distance() ;
      result =  DIST( mesh_id() ) ;
   }

   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_LocalFEbound:: set_calculation_point( GE_Point const* pt )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: set_calculation_point" ) ;
   PEL_CHECK_PRE( set_calculation_point_PRE( pt ) ) ;

   set_CP( pt ) ;

   GE_Mpolyhedron const* bd_poly = BOUND->polyhedron() ;

   PDE_CellFE const* cell = BOUND->adjacent_cell() ;
   GE_Mpolyhedron const* cell_poly = cell->polyhedron() ;

   PDE_BoundFE const* adj_bd = 0 ;
   GE_Mpolyhedron const* adj_bd_poly = 0 ;
   PDE_CellFE const* adj_bd_cell = 0 ;
   GE_Mpolyhedron const* adj_bd_cell_poly = 0 ;
   if( BOUND->has_adjacent_bound() )
   {
      PEL_ASSERT( nb_nested_meshes( BOUND ) <= 1 ) ; //?????
      PEL_ASSERT( nb_nested_meshes( cell ) <= 1 ) ;  //?????
      adj_bd = BOUND->adjacent_bound() ;
      adj_bd_poly = adj_bd->polyhedron() ;
      adj_bd_cell = adj_bd->adjacent_cell() ;
      adj_bd_cell_poly = adj_bd_cell->polyhedron() ;
      PEL_ASSERT( adj_bd_poly->reference_polyhedron() == 
                  bd_poly->reference_polyhedron() ) ;
   }

   size_t nb_bound_levels = nb_nested_meshes( BOUND ) ;
   size_t nb_cell_levels = nb_nested_meshes( cell ) ;
   if( nb_bound_levels!=0 && nb_cell_levels!=0 ) 
      PEL_ASSERT( nb_bound_levels == nb_cell_levels ) ;
   size_t nb_levels = PEL::max( nb_bound_levels, nb_cell_levels ) ;

   bd_poly->apply_inverse_mapping( pt, BD_PT_REF ) ;

   connect_CP( BOUND, 0, BD_PT_REF, 0, 0 ) ;

   cell_poly->apply_inverse_mapping( pt, CELL_PT_REF ) ;
   cell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
   CELL_tr_JAC->compute_determinant() ;
      
   doubleArray3D const* hess = 0 ;
   if( hessian_of_mapping_required() )
   {
      bool nonzero = false ;
      cell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero ) ;
      if( nonzero ) hess = HESSIAN ;
   }

   connect_CP( cell, 0, CELL_PT_REF, CELL_tr_JAC, hess ) ;

   if( adj_bd != 0 )
   {
      // le BD_PT_REF est le meme sur les deux frontières adjacentes
      connect_CP( adj_bd, 0, BD_PT_REF, 0, 0 ) ;
      
      adj_bd_poly->apply_mapping( BD_PT_REF, PT ) ;
      adj_bd_cell_poly->apply_inverse_mapping( PT, CELL_PT_REF ) ;
      adj_bd_cell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
      CELL_tr_JAC->compute_determinant() ;
     
      hess = 0 ;
      if( hessian_of_mapping_required() )
      {
         bool nonzero = false ;
         adj_bd_cell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero ) ;
         if( nonzero ) hess = HESSIAN ;
      }

      connect_CP( adj_bd_cell, 0, CELL_PT_REF, CELL_tr_JAC, hess ) ;
   }

   PDE_BoundFE const* bound = BOUND ;
   for( size_t cl=1 ; cl<nb_levels ; ++cl )
   {
      PDE_BoundFE const* pbound = bound->parent() ;
      GE_Mpolyhedron const* pbdpoly = pbound->polyhedron() ;
      pbdpoly->apply_inverse_mapping( pt, BD_PT_REF ) ;

      connect_CP( BOUND, cl, BD_PT_REF, 0, 0 ) ;

      PDE_CellFE const* pcell = pbound->adjacent_cell() ;
      GE_Mpolyhedron const* pcell_poly = pcell->polyhedron() ;
      pcell_poly->apply_inverse_mapping( pt, CELL_PT_REF ) ;
      pcell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
      CELL_tr_JAC->compute_determinant() ;

      hess = 0 ;
      if( hessian_of_mapping_required() )
      {
         bool nonzero = false ;
         pcell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero ) ;
         if( nonzero ) hess = HESSIAN ;
      }

      connect_CP( cell, cl, CELL_PT_REF, CELL_tr_JAC, hess ) ;

      bound = pbound ;         
   }

   PEL_CHECK_POST( set_calculation_point_POST( pt ) ) ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEbound:: print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: print_current_mesh" ) ;
   PEL_CHECK_PRE( print_current_mesh_PRE( os, indent_width ) ) ;

   BOUND->print( os, indent_width ) ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEbound:: compute_current_distance( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: compute_current_distance" ) ;
   PEL_CHECK( DIST.size() > mesh_id() ) ;

   size_t ii = mesh_id() ;
   GE_Mpolyhedron const* pbound = polyhedron() ;
   GE_Mpolyhedron const* pcell = BOUND->adjacent_cell()->polyhedron() ;

   GE_Point const* pt = pcell->finite_volume_center() ;
   if( pt == 0 ) PDE_LocalFEbound_ERROR::n0( pcell ) ;
   DUMMY_VEC->re_initialize( pbound->vertex( 0 ), pt ) ;
   DIST( ii ) = DUMMY_VEC->dot_product( outward_normal() ) ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEbound:: do_internals_on_current_mesh( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: do_internals_on_current_mesh" ) ;

   if( OK )
   {
      BOUND = the_bound( IT ) ;

      clear_local_fields() ;

      PDE_CellFE const* adjacent_cell = BOUND->adjacent_cell() ;
      for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
      {
         PDE_DiscreteField const* ff = handled_field( iF ) ;
        
         if( adjacent_cell->has_discretization( ff ) )
         {
            append_local_field( iF, adjacent_cell ) ;
         }
         else if( BOUND->has_discretization( ff ) )
         {
            append_local_field( iF, BOUND ) ;
         }
      }

      if( BOUND->has_adjacent_bound() )
      {
         PDE_BoundFE const* bd = BOUND->adjacent_bound() ;
         PDE_CellFE const* cell = bd->adjacent_cell() ;
         for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
         {
            PDE_DiscreteField const* ff = handled_field( iF ) ;

            if( cell->has_discretization( ff ) )
            {
               append_local_field( iF, cell ) ;
            }
            else if( bd->has_discretization( ff ) )
            {
               append_local_field( iF, bd ) ;
            }
         }
      }

      terminate_local_fields() ;
   }
   else
   {
      BOUND = 0 ;
   }

   reset_calculation_point() ;

   reset_row_and_col_fields() ;

   reset_IP_iterator() ;
}

//---------------------------------------------------------------------------
bool
PDE_LocalFEbound:: is_in_mesh( PDE_DiscreteField const* ff,
                               size_t i,
                               PDE_BasisFunction const* bf ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: is_in_mesh" ) ;

   bool result =  bf->is_located_on_bound( BOUND ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_QuadratureRule const*
PDE_LocalFEbound:: quadrature_rule( GE_QRprovider const* qrp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: quadrature_rule" ) ;
   PEL_CHECK( quadrature_rule_PRE( qrp ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* result = PDE_LocalFEmulti::quadrature_rule( qrp ) ;

   PEL_CHECK( quadrature_rule_POST( result, qrp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEbound:: compute_itg_pts( GE_QuadratureRule const* qr )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEbound:: compute_itg_pts" ) ;
   PEL_CHECK( is_valid() ) ;
   PEL_CHECK( !valid_IP() ) ;

   set_nb_IPs( qr->nb_points() ) ;

   GE_Mpolyhedron const* bd_poly = BOUND->polyhedron() ;

   PDE_CellFE const* cell = BOUND->adjacent_cell() ;
   GE_Mpolyhedron const* cell_poly = cell->polyhedron() ;

   PDE_BoundFE const* adj_bd = 0 ;
   GE_Mpolyhedron const* adj_bd_poly = 0 ;
   PDE_CellFE const* adj_bd_cell = 0 ;
   GE_Mpolyhedron const* adj_bd_cell_poly = 0 ;
   if( BOUND->has_adjacent_bound() )
   {
      PEL_ASSERT( nb_nested_meshes( BOUND ) <= 1 ) ; //?????
      PEL_ASSERT( nb_nested_meshes( cell ) <= 1 ) ;  //?????
      adj_bd = BOUND->adjacent_bound() ;
      adj_bd_poly = adj_bd->polyhedron() ;
      adj_bd_cell = adj_bd->adjacent_cell() ;
      adj_bd_cell_poly = adj_bd_cell->polyhedron() ;
   }

   size_t nb_bound_levels = nb_nested_meshes( BOUND ) ;
   size_t nb_cell_levels = nb_nested_meshes( cell ) ;
   if( nb_bound_levels!=0 && nb_cell_levels!=0 ) 
      PEL_ASSERT( nb_bound_levels == nb_cell_levels ) ;
   size_t nb_levels = PEL::max( nb_bound_levels, nb_cell_levels ) ;

   double coef = bd_poly->measure() / qr->sum_of_weights() ;

   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      GE_Point const* bd_pt_ref = qr->point( i ) ;
      bd_poly->apply_mapping( bd_pt_ref, PT ) ;
      double weight = qr->weight(i) * coef ;

      append_IP( PT, weight ) ;
      connect_IP( BOUND, 0, bd_pt_ref, 0, 0 ) ;
      
      cell_poly->apply_inverse_mapping( PT, CELL_PT_REF ) ;
      cell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
      CELL_tr_JAC->compute_determinant() ;
      
      doubleArray3D const* hess = 0 ;
      if( hessian_of_mapping_required() )
      {
         bool nonzero = false ;
         cell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero ) ;
         if( nonzero ) hess = HESSIAN ;
      }

      connect_IP( cell, 0, CELL_PT_REF, CELL_tr_JAC, hess ) ;
 
      if( adj_bd != 0 )
      {
         // le bd_pt_ref est le meme sur les deux frontières adjacentes
         connect_IP( adj_bd, 0, bd_pt_ref, 0, 0 ) ;
      
         adj_bd_poly->apply_mapping( bd_pt_ref, PT ) ;
         adj_bd_cell_poly->apply_inverse_mapping( PT, CELL_PT_REF ) ;
         adj_bd_cell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
         CELL_tr_JAC->compute_determinant() ;
      
         hess = 0 ;
         if( hessian_of_mapping_required() )
         {
            bool nonzero = false ;
            adj_bd_cell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero ) ;
            if( nonzero ) hess = HESSIAN ;
         }

         connect_IP( adj_bd_cell, 0, CELL_PT_REF, CELL_tr_JAC, hess ) ;
      }

      PDE_BoundFE const* bound = BOUND ;
      for( size_t cl=1 ; cl<nb_levels ; ++cl )
      {
         PDE_BoundFE const* pbound = bound->parent() ;
         GE_Mpolyhedron const* pbdpoly = pbound->polyhedron() ;
         pbdpoly->apply_inverse_mapping( PT, BD_PT_REF ) ;

         connect_IP( BOUND, cl, BD_PT_REF, 0, 0 ) ;

         PDE_CellFE const* pcell = pbound->adjacent_cell() ;
         GE_Mpolyhedron const* pcell_poly = pcell->polyhedron() ;
         pcell_poly->apply_inverse_mapping( PT, CELL_PT_REF ) ;
         pcell_poly->build_tr_mapping_derivative( CELL_PT_REF, CELL_tr_JAC ) ;
         CELL_tr_JAC->compute_determinant() ;
         hess = 0 ;
         if( hessian_of_mapping_required() )
         {
            bool nonzero = false ;
            pcell_poly->build_mapping_hessian( CELL_PT_REF, HESSIAN, nonzero );
            if( nonzero ) hess = HESSIAN ;
         }

         connect_IP( cell, cl, CELL_PT_REF, CELL_tr_JAC, hess ) ;

         //????? si il y a un maillage adjacent ???

         bound = pbound ;         
      }
   }
}

//---------------------------------------------------------------------------
bool
PDE_LocalFEbound:: rejected_bound( void ) const
//---------------------------------------------------------------------------
{
   bool result =  is_excluded( the_bound( IT )->color() )
               || is_excluded( the_bound( IT )->adjacent_cell()->color() )
               || !the_bound( IT )->is_active() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_BoundFE const*
PDE_LocalFEbound:: the_bound( PEL_VectorIterator const* it )
//---------------------------------------------------------------------------
{
   return( static_cast<PDE_BoundFE*>( it->item() ) ) ;
}

//internal--------------------------------------------------------------
void 
PDE_LocalFEbound_ERROR:: n0( GE_Mpolyhedron const* poly )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_LocalFEbound:" << endl << endl ;
   msg << "    the \"finite volume center\" has not been defined" << endl ;
   msg << "    for the following polyhedron:" << endl << endl ;
   poly->print( msg, 4 ) ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
