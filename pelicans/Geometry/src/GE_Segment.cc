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

#include <GE_Segment.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>

#include <iostream>

GE_Segment const* GE_Segment::MODEL = new GE_Segment() ;

//-----------------------------------------------------------------------------
GE_Segment:: GE_Segment( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Segment" )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Segment:: GE_Segment" ) ;
//   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Segment* 
GE_Segment:: create_replica( PEL_Object* a_owner, 
                             PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Segment* result = new GE_Segment( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Segment:: GE_Segment( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
   , UNIT_NORMAL( 0 )
   , MEASURE( 0. )
{
   PEL_LABEL( "GE_Segment:: GE_Segment" ) ;
   
   if( nb_space_dimensions()==2 )
   {
      UNIT_NORMAL = GE_Vector::create( this, 2 ) ;
   }
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Segment:: ~GE_Segment( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: ~GE_Segment" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
bool
GE_Segment:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   bool ok = ( nb_space_dimensions()>=1 ) ;
   return( ok ) ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_Segment:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static std::string const result = "GE_Segment" ;
   
   PEL_CHECK_POST( result == "GE_Segment" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Segment:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   double result = MEASURE ;
   PEL_CHECK_POST( measure_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_Segment:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = center() ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Segment:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: contains" ) ;
   PEL_CHECK_PRE( contains_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   
   double dot  = 0. ;
   double snorm2 = 0. ;
   for( size_t d = 0; d<nb_space_dimensions(); d++ )
   {
      double const v0 = vertex(0)->coordinate(d) ;
      double const v0v1 = vertex(1)->coordinate(d)-v0 ;
      double const v0pt = pt->coordinate(d)-v0 ;
      dot += v0pt*v0v1 ;
      snorm2 += v0pt*v0pt ;
   }
   
   double const alpha = dot/MEASURE/MEASURE ;
   double const epsilon_ref = reference_polyhedron()->epsilon() ;

   bool is_in = ( PEL::toler(alpha+epsilon_ref) ||
                  PEL::toler(alpha-1.-epsilon_ref) ||
                  ( alpha+epsilon_ref>0. && alpha-epsilon_ref<1.) ) ;
   if( is_in && nb_space_dimensions()>1 )
   {
      double const epsilon_real = epsilon_ref*MEASURE ;
      double const normP_Pproject =
           PEL::sqrt( PEL::abs( snorm2-dot*dot/MEASURE/MEASURE ) ) ;
      is_in = PEL::toler( normP_Pproject, epsilon_real ) ;
   }
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_Segment:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   GE_Vector const* result = 0 ;

   result = UNIT_NORMAL ;

   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Segment:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;
   
   static GE_ReferencePolyhedron const* result = GE_ReferenceSegment::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Segment:: apply_mapping( GE_Point const* pt_ref,
                            GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double x = pt_ref->coordinate( 0 ) ;
   
   double umx = 1.0-x ;
   
   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      pt->set_coordinate( d, umx * vertex( 0 )->coordinate( d ) +
                                  x * vertex( 1 )->coordinate( d ) ) ;
   }

   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Segment:: apply_inverse_mapping( GE_Point const* pt,
                                    GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   double x = 0. ;
   for( size_t d = 0; d<nb_space_dimensions(); d++ )
   {
      double const v0 = vertex( 0 )->coordinate( d )  ;
      x += ( pt->coordinate(d)-v0 )*( vertex(1)->coordinate(d)-v0 ) ;
   }
   pt_ref->set_coordinate( 0, x/MEASURE/MEASURE ) ;

   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Segment:: build_mapping_derivative( GE_Point const* pt_ref,
                                       GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;
   
   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      jac->set_item( d, 0, vertex( 1 )->coordinate( d ) - 
                           vertex( 0 )->coordinate( d ) ) ;
   }
}

//----------------------------------------------------------------------------
void
GE_Segment:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                          GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;
   
   size_t dim = nb_space_dimensions() ;
   for( size_t d=0 ; d<dim ; ++d )
   {
      tjac->set_item( 0, d, vertex( 1 )->coordinate( d ) - 
                            vertex( 0 )->coordinate( d ) ) ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Segment:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;

   MEASURE = computed_measure() ;
   if( nb_space_dimensions()==2 )
   {
      GE_Point const* pt0 = vertex( 0 ) ;
      GE_Point const* pt1 = vertex( 1 ) ;
      UNIT_NORMAL->set_component( 0, 
                    ( pt0->coordinate( 1 ) - pt1->coordinate( 1 ) )/MEASURE ) ;
      UNIT_NORMAL->set_component( 1, 
                    ( pt1->coordinate( 0 ) - pt0->coordinate( 0 ) )/MEASURE ) ;
   }
}

//-----------------------------------------------------------------------------
double
GE_Segment:: computed_measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Segment:: computed_measure" ) ;
   PEL_CHECK( !is_prototype() ) ;
   
   return( vertex(0)->distance( vertex(1) ) ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Segment:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result != 0 && result == center() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Segment:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const 
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceSegment::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Segment:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT( PEL::equal( MEASURE, computed_measure() ) ) ;
      PEL_ASSERT( IMPLIES( check_consistency(), is_consistent() ) ) ;
      if( nb_space_dimensions()==2 )
      {
         GE_Vector* V0V1 =
            GE_Vector::create( 0, vertex( 1 ), vertex( 0 ) ) ;
         PEL_ASSERT( UNIT_NORMAL!=0 &&
                     UNIT_NORMAL->nb_components()==2 ) ;
         PEL_ASSERT( PEL::toler( UNIT_NORMAL->dot_product( V0V1 ) ) ) ;
         V0V1->destroy() ;
      }
      else
      {
         PEL_ASSERT( UNIT_NORMAL==0 ) ;
      }
   }
   return( true ) ;
}
