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

#include <PDE_HelperFIC.hh>

#include <PEL_Error.hh>
#include <PEL_VectorIterator.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_FaceFE.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;

struct PDE_HelperFIC_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( void ) ;
   static void n3( GE_Mpolyhedron const* poly ) ;
} ;

//----------------------------------------------------------------------------
PDE_HelperFIC:: PDE_HelperFIC( PEL_Object* a_owner )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ID_NUMBER( 2 )
   , NB_DIMS( PEL::bad_index() )
   , SIDE( 0 )
   , MESHES( 3, (PDE_MeshFE*)0 )
   , X_REF( PEL::bad_index() )
   , Y_REF( PEL::bad_index() ) 
{
}

//----------------------------------------------------------------------------
PDE_HelperFIC:: ~PDE_HelperFIC( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
size_t
PDE_HelperFIC:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------------
{
   return( NB_DIMS ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_HelperFIC:: current_side_id( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: current_side_id" ) ;
   
   size_t result = ( SIDE!=0 ? SIDE->side_id() : PEL::bad_index() ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_HelperFIC:: parent_cell_id( size_t i_adj ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: parent_cell_id" ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   
   PDE_CellFE* rcell = SIDE->adjacent_cell( i_adj ) ;
   
   // it is essential that i_adj represents the 
   // same adjacent cell as PDE_CursorFEside::adjacent_localFEcell( i_adj )
   PEL_ASSERT( rcell->id_number() == ID_NUMBER( i_adj ) ) ;
   
   bool call_test = ( rcell->parent() != 0 ) && 
                    ( rcell->refinement_level()==SIDE->refinement_level() ) ;
   if( !call_test ) PDE_HelperFIC_ERROR::n2() ;
   
   size_t result = rcell->parent()->id_number() ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_HelperFIC:: prepare_for_interpolation( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: prepare_for_interpolation" ) ;
   PEL_CHECK_PRE( !ready_for_interpolation() ) ;
   
   size_t level = SIDE->refinement_level() ;
   if( level == 0 ) PDE_HelperFIC_ERROR::n0() ;
   PEL_ASSERT( SIDE->nb_adjacent_cells() == 2 ) ;
   
   MESHES[ 0 ] = 0 ;
   PDE_CellFE* cell_fine = 0 ;
   PDE_CellFE* cc = SIDE->adjacent_cell( 0 ) ;
   if( cc->refinement_level() == level )
   {
      MESHES[ 0 ] = SIDE->adjacent_cell( 1 ) ;
      cell_fine = cc ;
   }
   else
   {
      MESHES[ 0 ] = cc ;
      cell_fine = SIDE->adjacent_cell( 1 ) ;
   }
   PEL_ASSERT( MESHES[ 0 ]->refinement_level() == level-1 ) ;
   PEL_ASSERT( cell_fine->refinement_level() == level ) ;
   
   MESHES[ 1 ] = cell_fine->parent() ;
   PEL_ASSERT( MESHES[ 1 ] != 0 ) ;
   
   PDE_CellFE* cell_2 = 0 ;
   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, 
                                                        cell_fine->faces() ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_FaceFE* ff = static_cast<PDE_FaceFE*>( it->item() ) ;
      if( ff->nb_adjacent_cells() != 2 ) PDE_HelperFIC_ERROR::n1() ;
      PDE_CellFE* ocell = ff->adjacent_cell_other_than( cell_fine ) ;
      if( !( ocell->refinement_level()==level &&
             ocell->parent() == MESHES[ 1 ] ) )
      {
         // ff n'est pas une face interne à MESHES[ 1 ]
         if( ff != SIDE )
         {
            PEL_ASSERT( cell_2 == 0 ) ;
            cell_2 = ocell ;
         }
         if( cell_2 == 0 )
         {
            PEL_ASSERT( ff == SIDE ) ;
         }
      }
   }
   if( cell_2->refinement_level() ==  level )
   {
      cell_2 = cell_2->parent() ;
   }
   else
   {
      if( cell_2->refinement_level() != level-1 )
      {
         PEL_Error::object()->raise_plain( "bad meshing for FIC" ) ;
      }
   }
   MESHES[ 2 ] = cell_2 ;
   
   it->destroy() ;
   
   PEL_CHECK_POST( ready_for_interpolation() ) ;
}

//---------------------------------------------------------------------------
bool
PDE_HelperFIC:: ready_for_interpolation( void ) const
//---------------------------------------------------------------------------
{
   return( MESHES[0] != 0 ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_HelperFIC:: node( size_t i_pt,
                      PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: node" ) ;
   PEL_CHECK_PRE( ready_for_interpolation() ) ;
   PEL_CHECK_PRE( IMPLIES( nb_space_dimensions()==2, i_pt<3 )  ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   
   PDE_MeshFE const* mesh = MESHES[ i_pt ] ;
   
   size_t ee = mesh->index_of_reference_element( ff ) ;
   PEL_ASSERT( mesh->nb_basis_functions( ee ) == 1 ) ;
   
   PDE_BasisFunction const* bf = mesh->basis_function( ee, 0 ) ;
   
   size_t result = bf->node_of_DOF( ff ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
PDE_HelperFIC:: point_of_node( size_t i_pt ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: point_of_node" ) ;
   PEL_CHECK_PRE( ready_for_interpolation() ) ;
   PEL_CHECK_PRE( IMPLIES( nb_space_dimensions()==2, i_pt<3 )  ) ;
   
   GE_Mpolyhedron const* poly = MESHES[ i_pt ]->polyhedron() ;
   GE_Point const* result = poly->finite_volume_center() ;
   if( result == 0 ) PDE_HelperFIC_ERROR::n3( poly ) ; 
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_HelperFIC:: set_calculation_point( GE_Point const* pt )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: set_calculation_point" ) ;
   PEL_CHECK_PRE( ready_for_interpolation() ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates() == nb_space_dimensions() ) ;
   
   GE_Point const* V0 = point_of_node( 0 ) ;
   GE_Point const* V1 = point_of_node( 1 ) ;
   GE_Point const* V2 = point_of_node( 2 ) ;
   
   double V0M_x = pt->coordinate( 0 ) - V0->coordinate( 0 ) ;
   double V0M_y = pt->coordinate( 1 ) - V0->coordinate( 1 ) ;
   
   double V0V1_x = V1->coordinate( 0 ) - V0->coordinate( 0 ) ;
   double V0V1_y = V1->coordinate( 1 ) - V0->coordinate( 1 ) ;
   
   double V0V2_x = V2->coordinate( 0 ) - V0->coordinate( 0 ) ;
   double V0V2_y = V2->coordinate( 1 ) - V0->coordinate( 1 ) ;

   double det = V0V1_x * V0V2_y - V0V1_y * V0V2_x ;
   
   X_REF = ( V0M_x  * V0V2_y - V0M_y  * V0V2_x ) / det ;
   Y_REF = ( V0V1_x * V0M_y  - V0V1_y * V0M_x  ) / det ;
}

//----------------------------------------------------------------------------
double
PDE_HelperFIC:: interpolated_value( doubleVector const& value_at_nodes ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: interpolated_value" ) ;
   PEL_CHECK_PRE( ready_for_interpolation() ) ;
   PEL_CHECK_PRE( IMPLIES( nb_space_dimensions() == 2,
                           value_at_nodes.size() == 3 ) ) ;
   
   double result = value_at_nodes( 0 ) * ( 1.0 - X_REF - Y_REF ) +
                   value_at_nodes( 1 ) * ( X_REF ) +
                   value_at_nodes( 2 ) * ( Y_REF ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_HelperFIC:: set_side( PDE_FaceFE const* a_side,
                          size_t id_cell_0,
                          size_t id_cell_1 )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_HelperFIC:: set_side" ) ;
   PEL_CHECK( a_side != 0 ) ;
   PEL_CHECK( a_side->is_active() ) ;
   
   ID_NUMBER( 0 ) = id_cell_0 ;
   ID_NUMBER( 1 ) = id_cell_1 ;
   
   SIDE = a_side ;
   NB_DIMS = a_side->polyhedron()->nb_space_dimensions() ;
   
   if( NB_DIMS != 2 )
   {
      PEL_Error::object()->raise_plain( 
            "*** PDE_HelperFIC is only available for 2D geometries" ) ;
   }

   MESHES[0] = 0 ;
   
   X_REF = PEL::bad_double() ;
   Y_REF = PEL::bad_double() ;
   
   PEL_CHECK( !ready_for_interpolation() ) ;
}

//internal--------------------------------------------------------------
void
PDE_HelperFIC_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** PDE_HelperFIC:" << std::endl << std::endl ;
   msg << "    the interpolation facilities are only" << std::endl ;
   msg << "    meaningful from a refined level" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_HelperFIC_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** PDE_HelperFIC:" << std::endl << std::endl ;
   msg << "    the facilities handing refinements next to" << std::endl ;
   msg << "    a boundary remain to be implemented" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_HelperFIC_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** PDE_HelperFIC:" << std::endl << std::endl ;
   msg << "    invalid call to \"parent_cell_id\"" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_HelperFIC_ERROR:: n3( GE_Mpolyhedron const* poly )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** PDE_HelperFIC:" << std::endl << std::endl ;
   msg << "    the \"finite volume center\" has not been defined" << endl ;
   msg << "    for the following polyhedron:" << endl << endl ;
   poly->print( msg, 4 ) ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

