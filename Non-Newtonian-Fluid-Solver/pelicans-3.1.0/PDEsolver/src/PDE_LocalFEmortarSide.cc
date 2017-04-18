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

#include <PDE_LocalFEmortarSide.hh>

#include <GE_QuadratureRule.hh>
#include <GE_Color.hh>
#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <LA_DenseMatrix.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_InterfaceBuilder.hh>
#include <PDE_MortarSideFE.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <iostream>
using std::endl ;


//----------------------------------------------------------------------
PDE_LocalFEmortarSide:: PDE_LocalFEmortarSide( PEL_Object* a_owner,
					       size_t nb_sp_dims,
					       PEL_Vector const* mortar_sides )
//----------------------------------------------------------------------
   : PDE_LocalFEmulti( a_owner, nb_sp_dims )
   , NB_SIDES( mortar_sides->count() )
   , IT( PEL_VectorIterator::create( this, mortar_sides ) )
   , SIDE( 0 )
   , OK( false )     
   , PT( GE_Point::create( this, nb_sp_dims ) )
   , BD_PT_REF( GE_Point::create( this, nb_sp_dims-1 ) )
   , CELL_PT_REF( GE_Point::create( this, nb_sp_dims ) )
   , CELL_tr_JAC( GE_Matrix::create( this, nb_sp_dims, nb_sp_dims ) )
   , CELL_HESSIAN( new doubleArray3D( nb_sp_dims, nb_sp_dims, nb_sp_dims  ) )
{
}

//----------------------------------------------------------------------
PDE_LocalFEmortarSide:: ~PDE_LocalFEmortarSide( void )
//----------------------------------------------------------------------
{
   delete CELL_HESSIAN ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmortarSide:: nb_meshes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: nb_meshes" ) ;

   return( NB_SIDES ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: go_i_th( size_t an_id_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: go_i_th" ) ;
   PEL_CHECK_PRE( go_i_th_PRE( an_id_mesh ) ) ;
   PEL_CHECK_INV( invariant() ) ;

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
   PEL_CHECK_POST( go_i_th_POST( an_id_mesh ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide::start" ) ;

   IT->start() ;
   while( IT->is_valid() && rejected_side() )
   {
      IT->go_next() ;
   }
   OK = IT->is_valid() ;

   do_internals_on_current_mesh() ;

   PEL_CHECK_POST( start_POST() ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmortarSide:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( OK ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide::go_next" ) ;
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
PDE_LocalFEmortarSide:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: polyhedron" ) ;
   PEL_CHECK_PRE( polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = SIDE->polyhedron() ;

   PEL_CHECK_POST( polyhedron_POST( result ) ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmortarSide:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: refinement_level" ) ;
   PEL_CHECK_PRE( refinement_level_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( SIDE->refinement_level() ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_LocalFEmortarSide:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: color" ) ;
   PEL_CHECK_PRE( color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = SIDE->color() ;

   PEL_CHECK_POST( color_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmortarSide:: mesh_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: mesh_id" ) ;
   PEL_CHECK_PRE( mesh_id_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = SIDE->id_number() ;

   PEL_CHECK_POST( mesh_id_POST( result ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: set_calculation_point( GE_Point const* pt )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: set_calculation_point" ) ;
   PEL_CHECK_PRE( set_calculation_point_PRE( pt ) ) ;

   set_CP( pt ) ;

   SIDE->polyhedron()->apply_inverse_mapping( pt, BD_PT_REF ) ;

   connect_CP( SIDE, 0, BD_PT_REF, 0, 0 ) ;
	            
   for( size_t d=0 ; d<2 ; ++d ) 
   { 
      connect_with_one_domain_meshes( false, pt, SIDE->domain_bounds( d ) ) ;
   }

   PEL_CHECK_POST( set_calculation_point_POST( pt ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: print_current_mesh( std::ostream& os, 
                                            size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: print_current_mesh" ) ;
   PEL_CHECK_PRE( print_current_mesh_PRE( os, indent_width ) ) ;

   SIDE->print( os, indent_width ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: do_internals_on_current_mesh( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: do_internals_on_current_mesh" ) ;

   if( OK )
   {
      SIDE = the_side( IT ) ;

      clear_local_fields() ;
      
      for( size_t domain_id=0 ; domain_id<2 ; ++domain_id )
      {
         append_one_domain_local_fields( SIDE->domain_bounds( domain_id ) ) ;
      }

      for( size_t iF=0 ; iF<nb_handled_fields() ; iF++ )
      {
         PDE_DiscreteField const* ff = handled_field( iF ) ;
         if( SIDE->has_discretization( ff ) ) 
         {
            append_local_field( iF, SIDE ) ;
         }
      }

      terminate_local_fields() ;
   }

   reset_calculation_point() ;

   reset_row_and_col_fields() ;

   reset_IP_iterator() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmortarSide:: is_in_mesh( PDE_DiscreteField const* ff,
                                    size_t i,
                                    PDE_BasisFunction const* bf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: is_in_mesh" ) ;

   bool result = SIDE->polyhedron()->contains( local_node_location( ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: append_one_domain_local_fields( 
                                            PEL_Vector const* domain_bounds ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: append_one_domain_local_fields" ) ;

   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, domain_bounds ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_BoundFE const* bd = static_cast<PDE_BoundFE*>( it->item() ) ;
      PDE_CellFE const* cell = bd->adjacent_cell() ;
   
      for( size_t iF=0 ; iF<nb_handled_fields() ; iF++ )
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
   it->destroy() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: compute_itg_pts( GE_QuadratureRule const* qr )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: compute_itg_pts" ) ;

   set_nb_IPs( qr->nb_points() ) ;

   GE_Mpolyhedron const* side_poly = SIDE->polyhedron() ;
   double coef = side_poly->measure()/ qr->sum_of_weights() ;

   for( size_t i=0 ; i<qr->nb_points() ; ++i )
   {
      GE_Point const* pt_ref = qr->point( i ) ;
      side_poly->apply_mapping( pt_ref, PT ) ;
      double weight = qr->weight(i) * coef ;

      append_IP( PT, weight ) ;
      connect_IP( SIDE, 0, pt_ref, 0, 0 ) ;
	            
      for( size_t d=0 ; d<2 ; ++d ) 
      { 
         connect_with_one_domain_meshes( true, PT, SIDE->domain_bounds( d ) ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEmortarSide:: connect_with_one_domain_meshes(
                                             bool is_IP,
                                             GE_Point const* pt,
                                             PEL_Vector const* domain_bounds )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmortarSide:: connect_with_one_domain_meshes" ) ;

   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, domain_bounds ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_BoundFE const* bd = static_cast<PDE_BoundFE*>( it->item() ) ;
      if( bd->polyhedron()->contains( pt ) )
      {
	 bd->polyhedron()->apply_inverse_mapping( pt, BD_PT_REF ) ;
         if( is_IP )
	    connect_IP( bd, 0, BD_PT_REF, 0, 0 ) ;
         else
	    connect_CP( bd, 0, BD_PT_REF, 0, 0 ) ;

	 PDE_CellFE const* cell = bd->adjacent_cell() ;
	 PEL_CHECK( cell->polyhedron()->contains( pt ) ) ;
	 cell->polyhedron()->apply_inverse_mapping( pt, CELL_PT_REF ) ;
	 cell->polyhedron()->build_tr_mapping_derivative( CELL_PT_REF, 
							  CELL_tr_JAC ) ;
         CELL_tr_JAC->compute_determinant() ;

         doubleArray3D const* hess = 0 ;
         if( hessian_of_mapping_required() )
         {
            bool nonzero = false ;
            cell->polyhedron()->build_mapping_hessian( 
                                        CELL_PT_REF, CELL_HESSIAN, nonzero ) ;
            if( nonzero ) hess = CELL_HESSIAN ;
         }

         if( is_IP )
	    connect_IP( cell, 0, CELL_PT_REF, CELL_tr_JAC, CELL_HESSIAN ) ;
         else
 	    connect_CP( cell, 0, CELL_PT_REF, CELL_tr_JAC, CELL_HESSIAN ) ;

//??? il se peut qu'un point d'integration soit au sommet d'un bound,
//??? auquel cas il appartienta deux bounds. On l'affecte a la
//??? premiere trouvee.
	 break ;
      }
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmortarSide:: rejected_side( void ) const
//----------------------------------------------------------------------
{
   bool result = is_excluded( the_side( IT )->color() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_MortarSideFE const*
PDE_LocalFEmortarSide:: the_side( PEL_VectorIterator const* it )
//----------------------------------------------------------------------
{
   return( static_cast<PDE_MortarSideFE*>( it->item() ) ) ;
}

