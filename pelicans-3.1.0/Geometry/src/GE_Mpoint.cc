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

#include <GE_Mpoint.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferencePoint.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>

#include <iostream>

GE_Mpoint const* GE_Mpoint::MODEL = new GE_Mpoint() ;
GE_Vector const* GE_Mpoint::UNIT_NORMAL_MODEL = 0 ;

//-----------------------------------------------------------------------------
GE_Mpoint:: GE_Mpoint( void )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( "GE_Mpoint" )
{
   PEL_LABEL( "GE_Mpoint:: GE_Mpoint" ) ;
   GE_Vector* n = GE_Vector::create( this, 1 ) ;
   n->set_component( 0, 1. ) ;
   UNIT_NORMAL_MODEL = n ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Mpoint* 
GE_Mpoint:: create_replica( PEL_Object* a_owner, 
                            PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, vertices ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpoint* result = new GE_Mpoint( a_owner, vertices ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( create_replica_POST( a_owner, vertices, result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Mpoint:: GE_Mpoint( PEL_Object* a_owner, PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : GE_Mpolyhedron( a_owner, vertices )
{
   PEL_LABEL( "GE_Mpoint:: GE_Mpoint" ) ;
   update() ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Mpoint:: ~GE_Mpoint( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: ~GE_Mpoint" ) ;
   if( is_prototype() )
   {
      MODEL = 0 ;
      UNIT_NORMAL_MODEL = 0 ;
   }
}

//-----------------------------------------------------------------------------
std::string const&
GE_Mpoint:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static std::string const result = "GE_Mpoint" ;
   
   PEL_CHECK_POST( result == "GE_Mpoint" ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpoint:: is_consistent( std::ostream& os, bool verbose ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: is_consistent" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   bool const ok = true ;
   return( ok ) ;
}

//-----------------------------------------------------------------------------
double
GE_Mpoint:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   return( 1. ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Mpoint:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = 0 ;
   
   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpoint:: contains( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: contains" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   bool const is_in = PEL::toler( pt->distance( vertex( 0 ) ) ) ;
   return( is_in ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_Mpoint:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   GE_Vector const* result = 0 ;
   if( nb_space_dimensions()==1 )
   {
      PEL_CHECK( UNIT_NORMAL_MODEL!=0 ) ;
      result = UNIT_NORMAL_MODEL ;
   }
   else
   {
      result = GE_Mpolyhedron::unit_normal() ;
   }
   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Mpoint:: reference_polyhedron( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( reference_polyhedron_PRE() ) ;
   
   static GE_ReferencePolyhedron const* result = GE_ReferencePoint::object() ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpoint:: apply_mapping( GE_Point const* pt_ref,
                           GE_Point* pt  ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: apply_mapping" ) ;
   PEL_CHECK_PRE( apply_mapping_PRE( pt_ref, pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   pt->set( vertex( 0 ) ) ;
   PEL_CHECK_POST( apply_mapping_POST( pt_ref, pt ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpoint:: apply_inverse_mapping( GE_Point const* pt,
                                   GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: apply_inverse_mapping" ) ;
   PEL_CHECK_PRE( apply_inverse_mapping_PRE( pt, pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;
   PEL_CHECK_POST( apply_inverse_mapping_POST( pt, pt_ref ) ) ;
}

//----------------------------------------------------------------------------
void
GE_Mpoint:: build_mapping_derivative( GE_Point const* pt_ref,
                                      GE_Matrix* jac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: build_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_mapping_derivative_PRE( pt_ref, jac ) ) ;
   
   jac->nullify() ;
}

//----------------------------------------------------------------------------
void
GE_Mpoint:: build_tr_mapping_derivative( GE_Point const* pt_ref,
                                         GE_Matrix* tjac ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: build_tr_mapping_derivative" ) ;
   PEL_CHECK_PRE( build_tr_mapping_derivative_PRE( pt_ref, tjac ) ) ;
   
   tjac->nullify() ;
}

//-----------------------------------------------------------------------------
void
GE_Mpoint:: update_internal( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpoint:: update_internal" ) ;
   PEL_CHECK( update_internal_PRE() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpoint:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   // more than the parent postcondition
   PEL_ASSERT( result == 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpoint:: reference_polyhedron_POST( GE_ReferencePolyhedron const* result ) const 
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferencePoint::object() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpoint:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::invariant() ) ;
   if( !is_prototype() && !is_updating() )
   {
      PEL_ASSERT( UNIT_NORMAL_MODEL!=0 &&
                  UNIT_NORMAL_MODEL->nb_components()==1 &&
                  UNIT_NORMAL_MODEL->component(0)==1. ) ;
   }
   return( true ) ;
}
