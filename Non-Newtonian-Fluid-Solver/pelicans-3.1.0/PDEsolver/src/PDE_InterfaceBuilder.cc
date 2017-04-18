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

#include <PDE_InterfaceBuilder.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>

#include <GE_Meshing.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_BasisFunctionMortarSide.hh>
#include <PDE_BoundFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DiscOnMeshFE.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_GridFE.hh>
#include <PDE_MortarSideFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>
#include <sstream>

//----------------------------------------------------------------------------
PDE_InterfaceBuilder*
PDE_InterfaceBuilder:: create( PEL_Object* a_owner,
                               PEL_ModuleExplorer* exp,
                               PDE_DomainBuilder const* builder_0,
                               PDE_DomainBuilder const* builder_1 )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: create" ) ;

   PDE_InterfaceBuilder* result = 
               new PDE_InterfaceBuilder( a_owner, exp, builder_0, builder_1 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_InterfaceBuilder:: PDE_InterfaceBuilder( PEL_Object* a_owner,
                                             PEL_ModuleExplorer* exp,
                                             PDE_DomainBuilder const* builder_0,
                                             PDE_DomainBuilder const* builder_1 )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DIM( exp->int_data( "nb_space_dimensions" ) )
   , EXP( exp )
   , FIELDS( PDE_SetOfDiscreteFields::create( this ) )
   , SIDES( 0 )
   , D0_NAME( builder_0->name() )
   , D1_NAME( builder_1->name() )
   , DISCS( PEL_Vector::create( this, GE_ReferencePolyhedron::nb_objects() ) ) 
{
   VERTICES = GE_SetOfPoints::create( this, DIM ) ;

   build_disc_ref_poly( exp ) ;
   
   build_meshing() ;

   if( EXP->has_module( "fields" ) ) 
   {      
      PEL_ModuleExplorer* ee = EXP->create_subexplorer( 0, "fields" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ;
         build_one_field( se ) ;
         se->destroy() ;
      }
      ee->destroy() ;
   }

   connect_with_domain( 0, builder_0->finite_element_grid()->bounds() ) ;

   connect_with_domain( 1, builder_1->finite_element_grid()->bounds() ) ;
}

//----------------------------------------------------------------------------
PDE_InterfaceBuilder:: ~PDE_InterfaceBuilder( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
size_t
PDE_InterfaceBuilder:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------------
{
   return( DIM ) ;
}

//----------------------------------------------------------------------------
PDE_SetOfDiscreteFields*
PDE_InterfaceBuilder:: set_of_discrete_fields( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: set_of_discrete_fields" ) ;

   PDE_SetOfDiscreteFields* result = FIELDS ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner()==result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_SetOfPoints*
PDE_InterfaceBuilder:: set_of_vertices( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: set_of_vertices" ) ;

   GE_SetOfPoints* result = VERTICES ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_Vector const*
PDE_InterfaceBuilder:: mortar_sides( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: mortar_sides" ) ;

   PEL_Vector const* result = SIDES ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; ++i ),
                     dynamic_cast<PDE_MortarSideFE*>( result->at(i) ) != 0 )) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; ++i ),
                     result->at(i)->owner() == set_of_vertices() )) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_InterfaceBuilder:: build_disc_ref_poly( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: build_disc_ref_poly" ) ;
   
   if( EXP->has_module( "fields" ) ) 
   {
      PEL_ModuleExplorer* ee = EXP->create_subexplorer( 0, "fields" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ;
         
         PDE_DiscreteField* ff = 
             PDE_DiscreteField::create( 0,
                                        se->string_data( "name" ),
                                        se->int_data( "nb_components" ),
                                        se->int_data( "storage_depth" ) ) ;
         
         PDE_ReferenceElement const* elm = 
            PDE_ReferenceElement::object( se->string_data( "element_name" ) ) ;
         PEL_Vector* elms = PEL_Vector::create( 0, 1 ) ;
         elms->set_at( 0, const_cast< PDE_ReferenceElement* >( elm ) ) ;
         
         FIELDS->append( ff, false, elms ) ;

         GE_ReferencePolyhedron const* rpoly = elm->reference_polyhedron() ;
         size_t i = rpoly->id_number() ;
         
         PEL_Object* oo = DISCS->at( i ) ;
         if( oo == 0 )
         {
            DISCS->set_at( i, PDE_DiscOnMeshFE::create( DISCS, ff, elm ) ) ;
         }
         else
         {
            PDE_DiscOnMeshFE* disc = 
                                   static_cast<PDE_DiscOnMeshFE*>( oo ) ;
            disc->add_discretization( ff, elm ) ;
         }
         se->destroy() ; se = 0 ;
      }
      ee->destroy() ;
   }
}

//----------------------------------------------------------------------------
void
PDE_InterfaceBuilder:: build_meshing( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: build_meshing" ) ;
   PEL_CHECK( VERTICES->nb_points() == 0 ) ;

   bool const c = GE_Mpolyhedron::check_consistency() ;
   if( EXP->has_entry( "check_meshing_consistency" ) )
   {
      bool const check_consistency =
                              EXP->bool_data( "check_meshing_consistency" ) ;
      if( check_consistency )
      {
         GE_Mpolyhedron::set_check_consistency() ;
      }
      else
      {
         GE_Mpolyhedron::unset_check_consistency() ;
      }
   }
   
   PEL_ModuleExplorer const* ee = EXP->create_subexplorer( 0, "GE_Meshing" ) ;

   GE_Meshing* T = GE_Meshing::create( 0, ee, DIM ) ;
//????????????? tester T

   GE_Point* pt = GE_Point::create( 0, DIM ) ;
   for( T->start_vertex_iterator() ; T->valid_vertex() ; T->go_next_vertex() )
   {
      pt->set_coordinates( T->vertex_coordinates() ) ;
      VERTICES->append( pt, T->vertex_color() ) ;
   }
   pt->destroy() ; pt=0 ;
   ee->destroy() ; ee=0 ;

   SIDES = PEL_Vector::create( this, T->nb_cells() ) ;
   size_t i_side = 0 ;
   for( T->start_cell_iterator() ; T->valid_cell() ; T->go_next_cell() )
   {
      GE_Mpolyhedron* poly = GE_Mpolyhedron::create( 
                                             T->cell_polyhedron_name(),
                                             VERTICES,
                                             T->cell_vertices() ) ;
      PEL_Object* oo = DISCS->at( poly->reference_polyhedron()->id_number() ) ;
      PDE_DiscOnMeshFE* disc = static_cast< PDE_DiscOnMeshFE* >( oo ) ;
      PDE_MortarSideFE* side = PDE_MortarSideFE::create( VERTICES, 
                                                         i_side,
                                                         poly,
                                                         T->cell_color(), 
                                                         0,
                                                         disc ) ;
      SIDES->set_at( i_side, side ) ;
      i_side++ ;
   }
   PEL_CHECK( SIDES->count() == T->nb_cells() ) ;
   T->destroy() ;

   if( c )
   {
      GE_Mpolyhedron::set_check_consistency() ;
   }
   else
   {
      GE_Mpolyhedron::unset_check_consistency() ;
   }
}

//----------------------------------------------------------------------------
void
PDE_InterfaceBuilder:: build_one_field( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceBuilder:: build_one_field" ) ;
   
   PDE_ReferenceElement const* elm = 
      PDE_ReferenceElement::object( exp->string_data( "element_name" ) ) ;

   GE_SetOfPoints* nodes = GE_SetOfPoints::create( 0, DIM ) ;
   size_t_array2D loc2glob( SIDES->count(), elm->nb_nodes() ) ;
   PDE_MortarSideFE* side = 0 ;
   GE_Point* pt = GE_Point::create( 0, DIM ) ;
   for( size_t im=0 ; im<SIDES->count() ; ++im )
   {
      side = static_cast<PDE_MortarSideFE*>( SIDES->at( im ) ) ;
      for( size_t i=0 ; i<elm->nb_nodes() ; ++i )
      {
         GE_Point const* pt_ref = elm->node_location( i ) ;
         side->polyhedron()->apply_mapping( pt_ref, pt ) ;
         size_t idx ;
         if( !nodes->has( pt ) )
         { 
            idx = nodes->nb_points() ;
            nodes->append( pt ) ;
         }
         else
         {
            idx = nodes->index( pt ) ;
         }
         PEL_CHECK( idx == nodes->index( pt ) ) ;
         loc2glob( im, i ) = idx ;
      }
   }
   pt->destroy() ; pt=0 ;
   
   PDE_DiscreteField* ff = FIELDS->item( exp->string_data( "name" ) ) ;
   ff->set_nb_nodes( nodes->nb_points() ) ;

   std::vector< PDE_BasisFunctionMortarSide* > 
             built_bfs( nodes->nb_points(), (PDE_BasisFunctionMortarSide*)0 ) ;
   for( size_t im=0 ; im<SIDES->count() ; ++im )
   {
      side = static_cast<PDE_MortarSideFE*>( SIDES->at( im ) ) ;

      size_t ee = side->index_of_element( elm ) ;
      for( size_t ln=0 ; ln<elm->nb_nodes() ; ++ln )
      {
         size_t gn = loc2glob( im, ln ) ;
         PDE_BasisFunctionMortarSide* bf = side->basis_function( ee, ln ) ;
         if( bf == 0 )
         {
            bf = built_bfs[ gn ] ;
            if( bf == 0 )
            {
               bf = PDE_BasisFunctionMortarSide::create( VERTICES ) ;
               bf->set_active() ;
               built_bfs[ gn ] = bf ;
               bf->append_one_field( ff, gn ) ;
            }
            bf->extend_pieces( side, ee, ln ) ;
            side->set_basis_function( ee, ln, bf ) ;
         }
         else if( built_bfs[ gn ] == 0 )
         {
            bf->append_one_field( ff, gn ) ;
            built_bfs[ gn ] = bf ;
            PEL_ASSERT( bf->node_of_DOF( ff ) == gn ) ;
         }
      }
   }

   nodes->destroy() ; nodes=0 ;
}

//----------------------------------------------------------------------------
void
PDE_InterfaceBuilder:: connect_with_domain( size_t domain_id,
                                            PEL_Vector const* domain_bounds )
//----------------------------------------------------------------------------
{
   PEL_ASSERT( DIM == 2 ) ;

   double alpha, beta ;
   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, domain_bounds ) ;

   PEL_VectorIterator* itm = PEL_VectorIterator::create( 0, SIDES ) ;
   size_t nb_connections = 0 ;
   for( itm->start() ; itm->is_valid() ; itm->go_next() )
   {
      PDE_MortarSideFE* side = static_cast<PDE_MortarSideFE*>( itm->item() ) ;
      GE_Mpolyhedron const* side_poly = side->polyhedron() ;
      GE_Point const* A0 = side_poly->vertex( 0 ) ;
      GE_Point const* A1 = side_poly->vertex( 1 ) ;
      double length = A0->distance( A1 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         PDE_BoundFE* bd = static_cast<PDE_BoundFE*>( it->item() ) ;
         GE_Mpolyhedron const* bd_poly = bd->polyhedron() ;
         GE_Point const* B0 = bd_poly->vertex( 0 ) ;
         GE_Point const* B1 = bd_poly->vertex( 1 ) ;
         size_t ov = intersection_nature( A0, A1, B0, B1, alpha, beta ) ;
         bool connection = ( ov==4 || ov==3 ) ;
         if( ov==4 )
         {
            GE_Point const* pt1 = ( alpha == 0.0 ) ? A0 : A1  ;
            GE_Point const* pt2 = ( beta  == 0.0 ) ? B0 : B1  ;
            connection = ( (pt1->distance(pt2)/length) > 0.001 ) ; 
         }
         if( connection )
         {
            nb_connections++ ;
            side->append_domain_bound( domain_id, bd ) ;
         }
      }
      if( nb_connections==0 )
      {
         std::ostringstream message ;
         message << "invalid interface between domains \"" << D0_NAME
                 << "\" and " << D1_NAME << "\"" << std::endl ;
         message << "side of vertices : " << std::endl ;
         A0->print( message, 3 ) ; message << std::endl ;
         A1->print( message, 3 ) ; message << std::endl ;
         message << "is not connected to both domains" << std::endl ;
         PEL_Error::object()->raise_plain( message.str() ) ;
      }
   }
 
   itm->destroy() ;
   it->destroy() ;
}

// result = 0 : pas d'intersection
// result = 1 : pas d'intersection et paralleles 
// result = 2 : pas d'intersection et colineaires
// result = 3 : l'un est inclu dans l'autre
// result = 4 : colineaires et d'intersection non vide
//              alpha=0 : V0 de polyA est dans polyB
//              alpha=1 : V1 de polyA est dans polyB
//              beta=0  : V0 de polyB est dans polyA
//              beta=1  : V1 de polyB est dans polyA
// result = 5 : non colineaires et intersection reduite a un point
// intersection point is given by (1-`alpha')*`GE_PolygonSide::endpoint'(0)
// +`alpha'*`GE_PolygonSide::endpoint'(1) and 
// ( 1 - `beta' )*`pt0' + `beta'*`pt1'
//----------------------------------------------------------------------------
size_t
PDE_InterfaceBuilder:: intersection_nature( GE_Point const* A0, GE_Point const* A1,
                                       GE_Point const* B0, GE_Point const* B1,
                                       double& alpha, double& beta )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PolygonSide2D:: intersected_by" ) ;
   PEL_CHECK_PRE( ( A0 != 0 ) && ( A0->nb_coordinates() == 2 ) ) ;
   PEL_CHECK_PRE( ( A1 != 0 ) && ( A1->nb_coordinates() == 2 ) ) ;
   PEL_CHECK_PRE( ( B0 != 0 ) && ( B0->nb_coordinates() == 2 ) ) ;
   PEL_CHECK_PRE( ( B1 != 0 ) && ( B1->nb_coordinates() == 2 ) ) ;

   double epsilon = 1.E-08 ; //PEL::epsMachine ;

   double aX = A0->coordinate( 0 ) ;
   double aZ = A0->coordinate( 1 ) ;
   double bX = A1->coordinate( 0 ) ;
   double bZ = A1->coordinate( 1 ) ;
   double cX = B0->coordinate( 0 ) ;
   double cZ = B0->coordinate( 1 ) ;
   double dX = B1->coordinate( 0 ) ;
   double dZ = B1->coordinate( 1 ) ;

   double num_AB = (dX-cX)*(aZ-cZ)-(dZ-cZ)*(aX-cX) ;
   double num_CD = (bX-aX)*(aZ-cZ)-(bZ-aZ)*(aX-cX) ;
   double denomi = (bX-aX)*(dZ-cZ)-(bZ-aZ)*(dX-cX) ;

   size_t result = 0 ;

   if( PEL::abs( denomi )>epsilon ) // Denomi non nul ... the sides are not parallel
   {
      alpha = num_AB/denomi ;
      beta = num_CD/denomi ;
      if( alpha>-epsilon && alpha-1.<epsilon && 
          beta>-epsilon && beta-1.<epsilon )
         result = 5 ;
      else result = 0 ;
   }
   else // Sides are parallel
   {
      if( PEL::abs( num_AB )>epsilon && PEL::abs( num_CD )>epsilon )
         result = 1 ; // Sides are strictly parallel
      else // Sides are colinear
      {
         if( PEL::abs( aX - bX )>epsilon ) // [AB] is not vertical projection on the horizontal axis
         {
            double abXmin = PEL::min( aX, bX ), abXmax = PEL::max( aX, bX ) ;
            double cdXmin = PEL::min( cX, dX ), cdXmax = PEL::max( cX, dX ) ;
            if( cdXmax-abXmin<-epsilon || abXmax-cdXmin<-epsilon )
               result = 2 ; // Sides do not overlap
            else if( (abXmin-cdXmin<epsilon && cdXmax-abXmax<epsilon) || 
                     (cdXmin-abXmin<epsilon && abXmax-cdXmax<epsilon) )
               result = 3 ; // One side lies enterely in the other one
            else
            {
               result = 4 ; // each side share one enpoint of the other one
               if( cdXmin-aX<epsilon && aX-cdXmax<epsilon ) alpha = 0. ; // A is in CD
               else if ( cdXmin-bX<epsilon && bX-cdXmax<epsilon ) alpha = 1. ; // B is in CD
               if( abXmin-cX<epsilon && cX-abXmax<epsilon ) beta = 0. ; // C is in AB
               else if ( abXmin-dX<epsilon && dX-abXmax<epsilon ) beta = 1. ; // D is in AB
            }
         }
         else // [AB] is not horizontal projection on the vertical axis
         {
            double abZmin = PEL::min( aZ, bZ ), abZmax = PEL::max( aZ, bZ ) ;
            double cdZmin = PEL::min( cZ, dZ ), cdZmax = PEL::max( cZ, dZ ) ;
            if( cdZmax-abZmin<-epsilon || abZmax-cdZmin<-epsilon )
               result = 2 ; // Sides do not overlap
            else if( (abZmin-cdZmin<epsilon && cdZmax-abZmax<epsilon) || 
                     (cdZmin-abZmin<epsilon && abZmax-cdZmax<epsilon) )
               result = 3 ;// One side lies enterely in the other one
            else
            {
               result = 4 ; // each side share one enpoint of the other one
               if( cdZmin-aZ<epsilon && aZ-cdZmax<epsilon ) alpha = 0. ; // A is in CD
               else if ( cdZmin-bZ<epsilon && bZ-cdZmax<epsilon ) alpha = 1. ; // B is in CD
               if( abZmin-cZ<epsilon && cZ-abZmax<epsilon ) beta = 0. ; // C is in AB
               else if ( abZmin-dZ<epsilon && dZ-abZmax<epsilon ) beta = 1. ; // D is in AB
            }
         }
      }
   }

   PEL_CHECK_POST( IMPLIES( result==4, (alpha==0.0 || alpha==1.0) ) ) ;
   PEL_CHECK_POST( IMPLIES( result==4, (beta==0.0 || beta==1.0) ) ) ;
   return( result ) ;
}


