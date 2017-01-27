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

#include <GE_SegmentPolyhedron2D_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_Triangle.hh>
#include <GE_Vector.hh>

#include <PEL_Vector.hh>


GE_SegmentPolyhedron_INT const*
GE_SegmentPolyhedron2D_INT::PROTOTYPE = new GE_SegmentPolyhedron2D_INT() ;

//----------------------------------------------------------------------
void
GE_SegmentPolyhedron2D_INT:: check_intersection(
                                      GE_Point const* S0,
                                      GE_Point const* S1,
                                      GE_Mpolyhedron const* M )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: check_intersection" ) ;
   PEL_CHECK_PRE( check_intersection_PRE( S0, S1, M ) ) ;

   reset( S0, S1, M ) ;
   ONE_INTER = false ;

   size_t const nb_verts = M->nb_vertices() ;
   if( nb_verts == 3 )
   {
      check_intersection( S0, S1,
                          M->vertex( 0 ), M->vertex( 1 ), M->vertex( 2 ) ) ;
   }
   else
   {
      for( size_t i=0 ; !ONE_INTER && i<nb_verts ; ++i )
      {
         size_t ii = ( i==nb_verts-1 ? 0 : i+1 ) ;
         check_intersection( S0, S1,
                             M->vertex( i ), M->center(), M->vertex( ii ) ) ;
      }
   }

   declare_intersection_checked() ;

   PEL_CHECK_POST( check_intersection_POST( S0, S1, M ) ) ;
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron2D_INT:: one_single_intersection( void ) const
//----------------------------------------------------------------------
{
   return( ONE_INTER ) ;
}

//----------------------------------------------------------------------
void 
GE_SegmentPolyhedron2D_INT:: intersection_point( GE_Point* pt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: intersection_point" ) ;
   PEL_CHECK_PRE( intersection_point_PRE( pt ) ) ;

   pt->copy( PI ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron_INT*
GE_SegmentPolyhedron2D_INT:: create_replica(
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* a_mod_exp ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, a_mod_exp ) ) ;

   GE_SegmentPolyhedron2D_INT* result =
                         new GE_SegmentPolyhedron2D_INT( a_owner, a_mod_exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, a_mod_exp ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron2D_INT:: GE_SegmentPolyhedron2D_INT(
                     PEL_Object* a_owner, PEL_ModuleExplorer const* a_mod_exp )
//-----------------------------------------------------------------------------
   : GE_SegmentPolyhedron_INT( a_owner, 3 )
   , ONE_INTER( false )
   , PI( GE_Point::create( this, (size_t) 3 ) )
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: GE_SegmentPolyhedron2D_INT" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron2D_INT:: GE_SegmentPolyhedron2D_INT( void )
//-----------------------------------------------------------------------------
   : GE_SegmentPolyhedron_INT( "GE_SegmentPolyhedron2D_INT", 3 )
   , ONE_INTER( false )
   , PI( 0 )
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: GE_SegmentPolyhedron2D_INT" ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron2D_INT:: ~GE_SegmentPolyhedron2D_INT( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: ~GE_SegmentPolyhedron2D_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
GE_SegmentPolyhedron2D_INT:: check_intersection(
                                      GE_Point const* S0,
                                      GE_Point const* S1,
                                      GE_Point const* V0,
                                      GE_Point const* V1,
                                      GE_Point const* V2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron2D_INT:: check_intersection" ) ;
   PEL_CHECK( S0 != 0 && S0->nb_coordinates() == 3 ) ;
   PEL_CHECK( S1 != 0 && S1->nb_coordinates() == 3 ) ;
   PEL_CHECK( V0 != 0 && V0->nb_coordinates() == 3 ) ;
   PEL_CHECK( V1 != 0 && V1->nb_coordinates() == 3 ) ;
   PEL_CHECK( V2 != 0 && V2->nb_coordinates() == 3 ) ;
   PEL_CHECK( ONE_INTER == false ) ;

   double const epsilon = 1.e-8 ;
   
   ONE_INTER = false ;

   static GE_ReferenceTriangle const* poly_ref =
                                     GE_ReferenceTriangle::object() ;
   static GE_Vector* S0S1 =
                GE_Vector::create( PEL_Root::object(), (size_t) 3 ) ;
   static GE_Vector* S0V0 =
                GE_Vector::create( PEL_Root::object(), (size_t) 3 ) ;
   static GE_Vector* V0V1 =
                GE_Vector::create( PEL_Root::object(), (size_t) 3 ) ;
   static GE_Vector* V0V2 =
                GE_Vector::create( PEL_Root::object(), (size_t) 3 ) ;
   static GE_Vector* n =
                GE_Vector::create( PEL_Root::object(), (size_t) 3 ) ;
   static GE_Point* pt_ref =
                 GE_Point::create( PEL_Root::object(), (size_t) 2 ) ;
   
   S0S1->re_initialize( S1, S0 ) ;
   S0V0->re_initialize( V0, S0 ) ;
   V0V1->re_initialize( V1, V0 ) ;
   V0V2->re_initialize( V2, V0 ) ;
   n->set_as_cross_product( V0V1, V0V2 ) ;

   double const bary_seg = n->dot_product( S0V0 )/n->dot_product( S0S1 ) ;
   if( bary_seg>-epsilon && bary_seg<1.+epsilon )
   {
      PI->set_as_barycenter( bary_seg, S0, S1 ) ;
      GE_Triangle::apply_inverse_mapping( V0, V1, V2, n->norm()/2., PI,
                                          pt_ref ) ;
      ONE_INTER = poly_ref->contains( pt_ref ) ;
   }
}
