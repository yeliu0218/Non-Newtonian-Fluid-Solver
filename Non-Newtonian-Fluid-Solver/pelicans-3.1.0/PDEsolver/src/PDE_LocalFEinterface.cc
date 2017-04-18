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

#include <PDE_LocalFEinterface.hh>

#include <GE_QuadratureRule.hh>
#include <GE_Color.hh>
#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <LA_DenseMatrix.hh>

#include <PDE_CellFE.hh>
#include <PDE_BoundFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_GridFE.hh>
#include <PDE_FaceFE.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ostringstream ;

//----------------------------------------------------------------------
PDE_LocalFEinterface:: PDE_LocalFEinterface( PEL_Object* a_owner,
                                             PDE_GridFE const* grid,
                                             GE_Color const* inward_domain,
                                             GE_Color const* outward_domain )
//----------------------------------------------------------------------
   : PDE_LocalFEmulti( a_owner, grid->nb_space_dimensions() )
   , IT( 0 )
   , SIDE( 0 )
   , OK( false )
   , NB_SIDES( PEL::bad_index() )
   , IN_COLOR( inward_domain )
   , OUT_COLOR( outward_domain )
   , INTERFACE( PEL_Vector::create( this, 0 ) )
   , CELL_tr_JAC( GE_Matrix::create( this,
                                     grid->nb_space_dimensions(),
                                     grid->nb_space_dimensions() ) )
   , HESSIAN( new doubleArray3D( grid->nb_space_dimensions(), 
                                 grid->nb_space_dimensions(), 
                                 grid->nb_space_dimensions() ) )
   , CELL_PT_REF( GE_Point::create( this, grid->nb_space_dimensions() ) )
   , PT( GE_Point::create( this, grid->nb_space_dimensions() ) )
   , NORMAL( GE_Vector::create( this, grid->nb_space_dimensions() ) )
{
   PEL_CHECK( grid!=0 && inward_domain!=0 && outward_domain!=0 ) ;
   PEL_CHECK( !inward_domain->is_matching( outward_domain ) ) ;

   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, grid->cells() ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_CellFE* cell = static_cast<PDE_CellFE*>( it->item() ) ;
      if( inward_domain->is_matching( cell->color() ) )
      {
         PEL_VectorIterator* itS = 
                     PEL_VectorIterator::create( 0, cell->faces() ) ;
         for( itS->start() ; itS->is_valid() ; itS->go_next() )
         {
            PDE_FaceFE* side = static_cast<PDE_FaceFE*>( itS->item() ) ;
            if( side->nb_adjacent_cells()==2 )
            {
               PDE_CellFE const* other_cell = side->adjacent_cell_other_than( cell ) ;
               if( outward_domain->is_matching( other_cell->color() ) )
               {
                  INTERFACE->append( side ) ;
               }
            }
         }
         itS->destroy() ; itS=0 ;
      }
   }
   it->destroy() ; it=0 ;

   NB_SIDES = INTERFACE->count() ;

   IT = PEL_VectorIterator::create( this, INTERFACE ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_LocalFEinterface:: ~PDE_LocalFEinterface( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: ~PDE_LocalFEinterface" ) ;
   PEL_CHECK_INV( invariant() ) ;
   delete HESSIAN ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEinterface:: nb_meshes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: nb_meshes" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( NB_SIDES ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEinterface:: go_i_th( size_t an_id_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: go_i_th" ) ;
   PEL_CHECK_PRE( go_i_th_PRE( an_id_mesh ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( an_id_mesh != PEL::bad_index() )
   {
      IT->go_i_th( an_id_mesh ) ;
      OK =  IT->is_valid() && !rejected_side() ;
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
PDE_LocalFEinterface:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface::start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   OK = false ;

   IT->start() ;
   while( IT->is_valid() && rejected_side() )
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
PDE_LocalFEinterface:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( OK ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEinterface:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface::go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   IT->go_next() ;
   while( IT->is_valid() && rejected_side() )
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
PDE_LocalFEinterface:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: polyhedron" ) ;
   PEL_CHECK_PRE( polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = SIDE->polyhedron() ;

   PEL_CHECK_POST( polyhedron_POST( result ) ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEinterface:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: refinement_level" ) ;
   PEL_CHECK_PRE( refinement_level_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( SIDE->refinement_level() ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_LocalFEinterface:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: color" ) ;
   PEL_CHECK_PRE( color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = SIDE->color() ;

   PEL_CHECK_POST( color_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEinterface:: mesh_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: mesh_id" ) ;
   PEL_CHECK_PRE( mesh_id_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = IT->index_of_item() ;

   PEL_CHECK_POST( mesh_id_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector const*
PDE_LocalFEinterface:: outward_normal( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface::outward_normal" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Vector const* result = NORMAL ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( PEL::toler( result->norm()-1.0 ) ) ;
   PEL_CHECK_POST( result->nb_components()==nb_space_dimensions() ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_LocalFEinterface:: set_calculation_point( GE_Point const* pt )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: set_calculation_point" ) ;
   PEL_CHECK_PRE( set_calculation_point_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   set_CP( pt ) ;

   GE_Mpolyhedron const* cell_poly = CELL->polyhedron() ;

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

   connect_CP( CELL, 0, CELL_PT_REF, CELL_tr_JAC, HESSIAN ) ;

   PEL_CHECK_POST( set_calculation_point_POST( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEinterface:: print_current_mesh( std::ostream& os, 
                                           size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: print_current_mesh" ) ;
   PEL_CHECK_PRE( print_current_mesh_PRE( os, indent_width ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   SIDE->print( os, indent_width ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEinterface:: do_internals_on_current_mesh( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: do_internals_on_current_mesh" ) ;

   if( OK )
   {
      clear_local_fields() ;
      
      SIDE = the_side( IT ) ;
      CELL = SIDE->adjacent_cell( 0 ) ;
      if( !IN_COLOR->is_matching( CELL->color() ) )
      {
         PEL_CHECK( OUT_COLOR->is_matching( CELL->color() ) ) ;
         CELL = SIDE->adjacent_cell( 1 ) ;
         PEL_CHECK( IN_COLOR->is_matching( CELL->color() ) ) ;
      }
      for( size_t iF=0 ; iF<nb_handled_fields() ; iF++ )
      {
         PDE_DiscreteField const* ff = handled_field( iF ) ;
         if( !CELL->has_discretization( ff ) ) 
            raise_no_discretization_error( ff ) ;
         append_local_field( iF, CELL ) ;
      }

      NORMAL->set( SIDE->polyhedron()->unit_normal() ) ;
      GE_Vector* uu = GE_Vector::create( 0,
                                         CELL->polyhedron()->center(),
                                         SIDE->polyhedron()->vertex(0) ) ;
      if( NORMAL->dot_product(uu)>0.0 )
      {
         NORMAL->scale( -1.0 ) ;
      }
      uu->destroy() ;      

      terminate_local_fields() ;      
   }

   reset_calculation_point() ;

   reset_row_and_col_fields() ;

   reset_IP_iterator() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEinterface:: compute_itg_pts( GE_QuadratureRule const* qr )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEinterface:: compute_itg_pts" ) ;
   PEL_CHECK_INV( invariant() ) ;

   set_nb_IPs( qr->nb_points() ) ;

   GE_Mpolyhedron const* side_poly = SIDE->polyhedron() ;
   GE_Mpolyhedron const* cell_poly = CELL->polyhedron() ;

   double coef = side_poly->measure() / qr->sum_of_weights() ;

   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      GE_Point const* side_pt_ref = qr->point( i ) ;
      side_poly->apply_mapping( side_pt_ref, PT ) ;
      double weight = qr->weight( i ) * coef ;

      append_IP( PT, weight ) ;

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

      connect_IP( CELL, 0, CELL_PT_REF, CELL_tr_JAC, HESSIAN ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFEinterface:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( PDE_LocalFE::invariant() ) ;
   PEL_ASSERT( NB_SIDES == INTERFACE->count() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFEinterface:: rejected_side( void ) const
//------------------------------------------------------------------------
{
   bool result = is_excluded( the_side( IT )->color() ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
PDE_FaceFE const*
PDE_LocalFEinterface:: the_side( PEL_VectorIterator const* it )
//------------------------------------------------------------------------
{
   return( static_cast<PDE_FaceFE*>( it->item() ) ) ;
}
