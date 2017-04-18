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

#include <GE_SegmentPolyhedron1D_INT.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_SegmentSegment_INT.hh>
#include <GE_SegmentSegment2_INT.hh>

GE_SegmentPolyhedron_INT const*
GE_SegmentPolyhedron1D_INT::PROTOTYPE = new GE_SegmentPolyhedron1D_INT() ;

//----------------------------------------------------------------------
void 
GE_SegmentPolyhedron1D_INT:: check_intersection( GE_Point const* S0,
                                                 GE_Point const* S1,
                                                 GE_Mpolyhedron const* M )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: check_intersection" ) ;
   PEL_CHECK_PRE( check_intersection_PRE( S0, S1, M ) ) ;

   reset( S0, S1, M ) ;   
   ONE_INTER = false ;

   INTERSECTOR->compute_intersection( S0, S1, 
                                      M->vertex( 0 ), M->vertex( 1 ) ) ;
   GE_SegmentSegment_INT::IntersectionType IT
                                   = INTERSECTOR->intersection_type() ;
   
   if( IT==GE_SegmentSegment_INT::one_intersection )
   {
      ONE_INTER = true ;
      ALPHA = INTERSECTOR->alpha() ;
   }
   declare_intersection_checked() ;

   PEL_CHECK_INV( invariant() ) ;  
   PEL_CHECK_POST( check_intersection_POST( S0, S1, M ) ) ; 
}

//----------------------------------------------------------------------
void 
GE_SegmentPolyhedron1D_INT:: intersection_point( GE_Point* pt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: intersection_point" ) ;
   PEL_CHECK_PRE( intersection_point_PRE( pt ) ) ;

   pt->set_as_barycenter( 
             ALPHA, segment_first_vertex(), segment_second_vertex() ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentPolyhedron1D_INT:: one_single_intersection( void ) const
//----------------------------------------------------------------------
{
   return( ONE_INTER ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron_INT*
GE_SegmentPolyhedron1D_INT:: create_replica(
               PEL_Object* a_owner, PEL_ModuleExplorer const* a_mod_exp ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, a_mod_exp ) ) ;

   GE_SegmentPolyhedron1D_INT* result =
                       new GE_SegmentPolyhedron1D_INT( a_owner, a_mod_exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, a_mod_exp ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron1D_INT:: GE_SegmentPolyhedron1D_INT(
                     PEL_Object* a_owner, PEL_ModuleExplorer const* a_mod_exp )
//-----------------------------------------------------------------------------
   : GE_SegmentPolyhedron_INT( a_owner, 2 )
   , INTERSECTOR( 0 )
   , ONE_INTER( false )
   , ALPHA( PEL::bad_double() ) 
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: GE_SegmentPolyhedron1D_INT" ) ;

   if( a_mod_exp != 0 )
   {
      GE_PointPoint_INT* pt_pt_int = 0 ;
      if( a_mod_exp->has_module( "GE_PointPoint_INT" ) )
      {
         PEL_ModuleExplorer* se =
                  a_mod_exp->create_subexplorer( 0, "GE_PointPoint_INT" ) ;
         pt_pt_int =
            GE_PointPoint_INT::create( this,
                                       se->string_data( "concrete_name" ),
                                       se ) ;
         se->destroy() ; se = 0 ;
      }

      GE_PointSegment_INT* pt_seg_int = 0 ;
      if( a_mod_exp->has_module( "GE_PointSegment_INT" ) )
      {
         PEL_ModuleExplorer* se =
                  a_mod_exp->create_subexplorer( 0, "GE_PointSegment_INT" ) ;
         pt_seg_int =
            GE_PointSegment_INT::create( this,
                                         se->string_data( "concrete_name" ),
                                         se, pt_pt_int ) ;
         se->destroy() ; se = 0 ;
      }
      
      PEL_ModuleExplorer* se =
                  a_mod_exp->create_subexplorer( 0, "GE_SegmentSegment_INT" ) ;
      INTERSECTOR =
         GE_SegmentSegment_INT::create( this,
                                        se->string_data( "concrete_name" ),
                                        se, pt_pt_int, pt_seg_int ) ;
      se->destroy() ; se = 0 ;
   }
   else
   {
      INTERSECTOR = GE_SegmentSegment2_INT::create( this,
						    1.E-12, 1.E-12, 1.E-12 ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron1D_INT:: GE_SegmentPolyhedron1D_INT( void )
//-----------------------------------------------------------------------------
   : GE_SegmentPolyhedron_INT( "GE_SegmentPolyhedron1D_INT", 2 )
   , INTERSECTOR( 0 )
   , ONE_INTER( false )
   , ALPHA( PEL::bad_double() ) 
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: GE_SegmentPolyhedron1D_INT" ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_SegmentPolyhedron1D_INT:: ~GE_SegmentPolyhedron1D_INT( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron1D_INT:: ~GE_SegmentPolyhedron1D_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}
