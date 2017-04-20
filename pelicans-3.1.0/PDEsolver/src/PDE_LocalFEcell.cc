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

#include <PDE_LocalFEcell.hh>

#include <GE_Customized_QR.hh>
#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_BoundFE.hh>
#include <PDE_CFootFinder.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_FacesOfCellFE.hh>
#include <PDE_GridFE.hh>
#include <PDE_MeshFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_FaceFE.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <size_t_vector.hh>

#include <string>
#include <sstream>
#include <iostream>

using std::cout ; using std::endl ;
using std::ostringstream ;

//----------------------------------------------------------------------
PDE_LocalFEcell:: PDE_LocalFEcell( PEL_Object* a_owner,
                                   PDE_GridFE const* a_grid )
//----------------------------------------------------------------------
   : PDE_LocalFEsingle( a_owner, a_grid->nb_space_dimensions() )
   , CELLS( a_grid->cells() )
   , IT( PEL_VectorIterator::create( this, a_grid->cells() ) )
   , CELL( 0 )
   , OK( false )
   , GRID( a_grid )
   , GEO_STATE( a_grid->geometric_state() )
   , FFINDER( 0 )
   , tr_JAC( GE_Matrix::create( this,
                                a_grid->nb_space_dimensions(),
                                a_grid->nb_space_dimensions() ) )
   , HESSIAN( new doubleArray3D( a_grid->nb_space_dimensions(), 
                                 a_grid->nb_space_dimensions(), 
                                 a_grid->nb_space_dimensions() ) )
   , PT( GE_Point::create( this, a_grid->nb_space_dimensions() ) )
   , PT_REF( GE_Point::create( this, a_grid->nb_space_dimensions() ) )
   , ADJ_CELL( 0 )
   , ADJ_SIDE( 0 )
   , ADJ_BOUND( 0 )
{
   PEL_LABEL( "PDE_LocalFEcell:: PDE_LocalFEcell" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_LocalFEcell:: ~PDE_LocalFEcell( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: ~PDE_LocalFEcell" ) ;
   delete HESSIAN ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEcell:: nb_meshes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: nb_meshes" ) ;

   size_t result = GRID->nb_cells() ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: go_i_th( size_t an_id_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: go_i_th" ) ;
   PEL_CHECK_PRE( go_i_th_PRE( an_id_mesh ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, CELLS ) ;
   }
   GEO_STATE = GRID->geometric_state() ;

   if( an_id_mesh != PEL::bad_index() )
   {
      IT->go_i_th( an_id_mesh ) ;
      OK = IT->is_valid() && !rejected_cell() ;
   }
   else
   {
      OK = false ;
   }

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( go_i_th_POST( an_id_mesh ) ) ;
   PEL_CHECK_POST( IMPLIES( 
                   is_valid() && has_foot_finder(),
                   ( foot_finder()->polyhedron_of_head_cell() == polyhedron() ) &&
                    !foot_finder()->search_has_been_performed() ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() )
   {
      destroy_possession( IT ) ;
      IT = PEL_VectorIterator::create( this, CELLS ) ;
   }
   GEO_STATE = GRID->geometric_state() ;
   
   IT->start() ;
   while( IT->is_valid() && rejected_cell() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( start_POST() ) ;
   PEL_CHECK_POST( IMPLIES( 
                   is_valid() && has_foot_finder(),
                   ( foot_finder()->polyhedron_of_head_cell() == polyhedron() ) &&
                    !foot_finder()->search_has_been_performed() ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEcell:: is_valid( void ) const
//----------------------------------------------------------------------
{
   bool result = OK && ( GEO_STATE == GRID->geometric_state() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( IT->container_has_been_modified() || 
       GEO_STATE != GRID->geometric_state() )
   {
      ostringstream mesg ;
      mesg << "PDE_LocalFEcell : " << endl ;
      mesg << "   attempt to call \"go_next\" after" << endl ;
      mesg << "   a modification of the meshing" ;
      PEL_Error::object()->raise_internal( mesg.str() ) ;
   }

   IT->go_next() ;
   while( IT->is_valid() && rejected_cell() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( go_next_POST() ) ;
   PEL_CHECK_POST( IMPLIES( 
                   is_valid() && has_foot_finder(),
                   ( foot_finder()->polyhedron_of_head_cell() == polyhedron() ) &&
                    !foot_finder()->search_has_been_performed() ) ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_LocalFEcell:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: polyhedron" ) ;
   PEL_CHECK_PRE( polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = CELL->polyhedron() ;

   PEL_CHECK_POST( polyhedron_POST( result ) ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEcell:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: refinement_level" ) ;
   PEL_CHECK_PRE( refinement_level_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( CELL->refinement_level() ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_LocalFEcell:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: color" ) ;
   PEL_CHECK_PRE( color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = CELL->color() ;

   PEL_CHECK_POST( color_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEcell:: mesh_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: mesh_id" ) ;
   PEL_CHECK_PRE( mesh_id_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = CELL->id_number() ;

   PEL_CHECK_POST( mesh_id_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEcell:: adjacent_cell_ids( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: adjacent_cell_ids" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t_vector& result = ADJ_CELL ;
   
   result.resize( 0 ) ;
   PDE_FacesOfCellFE* its = CELL->create_faces_iterator( 0 ) ;
   for( its->start() ; its->is_valid() ; its->go_next() )
   {
      PDE_FaceFE const* s = its->item() ;
      if( s->is_active() && ! s->has_adjacent_bound() ) 
      {
         PDE_CellFE const* other_cell = 0 ;
         if( s->is_periodic() )
         {
            other_cell = s->periodic_neighbour()->adjacent_cell( 0 ) ;
         }
         else
         {
            other_cell = s->adjacent_cell_other_than( CELL ) ;
         }
         if( other_cell->is_active() )
         {
            result.append( other_cell->id_number() ) ;
         }
      }
   }
   its->destroy() ;
     
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEcell:: adjacent_side_ids( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: adjacent_side_ids" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t_vector& result = ADJ_SIDE ;
   
   result.resize( 0 ) ;
   PDE_FacesOfCellFE* its = CELL->create_faces_iterator( 0 ) ;
   for( its->start() ; its->is_valid() ; its->go_next() )
   {
      PDE_FaceFE const* s = its->item() ;
      if( s->is_active() && ! s->has_adjacent_bound() )
      {
         size_t ii = s->side_id() ;
         if( ii == PEL::bad_index() )
         {
            PEL_CHECK( s->is_periodic() ) ;
            ii = s->periodic_neighbour()->side_id() ;
            PEL_CHECK( ii != PEL::bad_index() ) ;
         }
         result.append( ii ) ;
      }
   }
   its->destroy() ;
     
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEcell:: adjacent_bound_ids( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: adjacent_bound_ids" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t_vector& result = ADJ_BOUND ;
   
   result.resize( 0 ) ;
   PDE_FacesOfCellFE* its = CELL->create_faces_iterator( 0 ) ;
   for( its->start() ; its->is_valid() ; its->go_next() )
   {
      PDE_FaceFE const* s = its->item() ;
      if( s->is_active() && s->has_adjacent_bound() )
      {
         PDE_BoundFE const* b = s->adjacent_bound() ;
         result.append( b->id_number() ) ;
      }
   }
   its->destroy() ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_LocalFEcell:: reference_element( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: reference_element" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( field_is_handled( ff ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const ee = CELL->index_of_reference_element( ff ) ;
   PDE_ReferenceElement const* result = CELL->reference_element( ee ) ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFEcell:: set_calculation_point( GE_Point const* pt )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: set_calculation_point" ) ;
   PEL_CHECK_PRE( set_calculation_point_PRE( pt ) ) ;

   set_CP( pt ) ;

   if( jacobian_of_mapping_required() )
   {
      GE_Mpolyhedron const* poly = CELL->polyhedron() ;
      poly->apply_inverse_mapping( pt, PT_REF ) ;
      poly->build_tr_mapping_derivative( PT_REF, tr_JAC ) ;
      tr_JAC->compute_determinant() ;

      doubleArray3D const* hess = 0 ;
      if( hessian_of_mapping_required() )
      {
         bool nonzero = false ;
         poly->build_mapping_hessian( PT_REF, HESSIAN, nonzero ) ;
         if( nonzero ) hess = HESSIAN ;
      }

      connect_CP( 0, PT_REF, tr_JAC, hess ) ;

      size_t nb_cell_levels = nb_nested_meshes() ;
      PDE_MeshFE const* pcell = CELL->parent() ;
      for( size_t cl=1 ; cl<nb_cell_levels ; ++cl )
      {
         GE_Mpolyhedron const* ppoly = pcell->polyhedron() ;
         ppoly->apply_inverse_mapping( pt, PT_REF ) ;
         ppoly->build_tr_mapping_derivative( PT_REF, tr_JAC ) ;
         tr_JAC->compute_determinant() ;
         hess = 0 ;
         if( hessian_of_mapping_required() )
         {
            bool nonzero = false ;
            poly->build_mapping_hessian( PT_REF, HESSIAN, nonzero ) ;
            if( nonzero ) hess = HESSIAN ;
         }
         connect_CP( cl, PT_REF, tr_JAC, hess ) ;
         pcell = pcell->parent() ;
      }
   }

   PEL_CHECK_POST( set_calculation_point_POST( pt ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: set_foot_finder( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: set_foot_finder" ) ;
   PEL_CHECK_PRE( exp!=0 &&
                  exp->has_entry( "concrete_name" ) ) ;
   PEL_CHECK_PRE( !has_foot_finder() ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   
   std::string const& concrete_name = exp->string_data( "concrete_name" ) ;
   FFINDER = PDE_CFootFinder::create( this, concrete_name, exp, GRID ) ;
   
   PEL_CHECK_POST( has_foot_finder() ) ;
   PEL_CHECK_POST( foot_finder()->owner()==this ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEcell:: has_foot_finder( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: has_foot_finder" ) ;
   
   return( FFINDER!=0 ) ;
}

//----------------------------------------------------------------------
PDE_CFootFinder*
PDE_LocalFEcell:: foot_finder( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: foot_finder" ) ;
   PEL_CHECK_PRE( has_foot_finder() ) ;
   
   PDE_CFootFinder* result = FFINDER ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: synchronize_foot_finder( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: synchronize_foot_finder" ) ;
   PEL_CHECK_PRE( is_valid() && has_foot_finder() ) ;
   
   FFINDER->synchronize( this, CELL ) ;
   
   PEL_CHECK_POST( foot_finder()->polyhedron_of_head_cell() == polyhedron() ) ;
   PEL_CHECK_POST( !foot_finder()->search_has_been_performed() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: print_current_mesh( std::ostream& os, 
                                      size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: print_current_mesh" ) ;
   PEL_CHECK_PRE( print_current_mesh_PRE( os, indent_width ) ) ;

   std::string space( indent_width, ' ' ) ;

   CELL->print( os, indent_width ) ;
   if( nb_nested_meshes() > 1 )
   {
      os << space 
         << nb_nested_meshes() << " superimposed nested cells" << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEcell:: do_internals_on_current_mesh( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: do_internals_on_current_mesh" ) ;

   if( OK )
   {
      CELL = the_cell( IT ) ;
      set_single_mesh( CELL ) ;
      if( FFINDER != 0 ) synchronize_foot_finder() ;
   }
   else
   {
      CELL = 0 ;
   }

   reset_calculation_point() ;

   reset_row_and_col_fields() ;

   reset_IP_iterator() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEcell:: is_in_mesh( PDE_BasisFunction const* bf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: is_in_mesh" ) ;

   bool result = bf->is_located_in_cell( CELL ) ;
   return( result ) ;

}

//----------------------------------------------------------------------
GE_QuadratureRule const*
PDE_LocalFEcell:: quadrature_rule( GE_QRprovider const* qrp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: quadrature_rule" ) ;
   PEL_CHECK( quadrature_rule_PRE( qrp ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* result = PDE_LocalFEsingle::quadrature_rule( qrp ) ;

   PEL_CHECK( quadrature_rule_POST( result, qrp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFEcell:: compute_itg_pts( GE_QuadratureRule const* qr )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: compute_itg_pts" ) ;
   PEL_CHECK( is_valid() ) ;
   PEL_CHECK( !valid_IP() ) ;
   
   set_nb_IPs( qr->nb_points() ) ;

   GE_Mpolyhedron const* poly = CELL->polyhedron() ;
   size_t nb_cell_levels = nb_nested_meshes() ;

   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      GE_Point const* pt_ref = qr->point( i ) ;
      poly->apply_mapping( pt_ref, PT ) ;

      poly->build_tr_mapping_derivative( pt_ref, tr_JAC ) ;
      tr_JAC->compute_determinant() ;
      double weight = qr->weight( i ) * PEL::abs( tr_JAC->determinant() ) ;

      doubleArray3D const* hess = 0 ;
      if( hessian_of_mapping_required() )
      {
         bool nonzero = false ;
         poly->build_mapping_hessian( pt_ref, HESSIAN, nonzero ) ;
         if( nonzero ) hess = HESSIAN ;
      }

//????? calculer toujours PT et le stocker... cf PDE_LocalFEbound, multi
      append_IP( weight ) ;

      if( jacobian_of_mapping_required() )
      {
         connect_IP( 0, pt_ref, tr_JAC, hess ) ;

         PDE_MeshFE const* pcell = CELL->parent() ;
         for( size_t cl=1 ; cl<nb_cell_levels ; ++cl )
         {
            GE_Mpolyhedron const* ppoly = pcell->polyhedron() ;
            ppoly->apply_inverse_mapping( PT, PT_REF ) ;
            ppoly->build_tr_mapping_derivative( PT_REF, tr_JAC ) ;
            tr_JAC->compute_determinant() ;
            hess = 0 ;
            if( hessian_of_mapping_required() )
            {
               bool nonzero = false ;
               poly->build_mapping_hessian( pt_ref, HESSIAN, nonzero ) ;
               if( nonzero ) hess = HESSIAN ;
            }
            connect_IP( cl, PT_REF, tr_JAC, hess ) ;
            pcell = pcell->parent() ;
         }
      }
   }
}

//---------------------------------------------------------------------------
bool
PDE_LocalFEcell:: rejected_cell( void ) const
//---------------------------------------------------------------------------
{
   bool result =  is_excluded( the_cell( IT )->color() )
               || !the_cell( IT )->is_active() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_CellFE*
PDE_LocalFEcell:: the_cell( PEL_VectorIterator const* it )
//---------------------------------------------------------------------------
{
   return( static_cast<PDE_CellFE*>( it->item() ) ) ;
}

