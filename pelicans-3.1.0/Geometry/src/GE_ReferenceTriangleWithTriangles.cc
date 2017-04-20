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

#include <GE_ReferenceTriangleWithTriangles.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceTriangle.hh>

#include <iostream>

GE_ReferenceTriangleWithTriangles const*
GE_ReferenceTriangleWithTriangles:: PROTOTYPE = 
                                    new GE_ReferenceTriangleWithTriangles() ;

//---------------------------------------------------------------------------
GE_ReferenceTriangleWithTriangles:: GE_ReferenceTriangleWithTriangles( void )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( "GE_ReferenceTriangleWithTriangles" ) 
{
}

//---------------------------------------------------------------------------
GE_ReferenceTriangleWithTriangles const*
GE_ReferenceTriangleWithTriangles:: object_2( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangleWithTriangles:: object_2" ) ;

   static GE_ReferenceTriangleWithTriangles* result = 0 ;
   if( result == 0 )
   {
      result = new GE_ReferenceTriangleWithTriangles( PEL_Root::object(), 2 ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   PEL_CHECK_POST( result->nb_subintervals_per_edge() == 2 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferenceTriangleWithTriangles const*
GE_ReferenceTriangleWithTriangles:: create_replica( 
                                           PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceTriangleWithTriangles:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   size_t nbs = exp->int_data( "nb_subintervals_per_edge" ) ;

   GE_ReferenceTriangleWithTriangles* result = 
               new GE_ReferenceTriangleWithTriangles( a_owner, nbs ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferenceTriangleWithTriangles:: GE_ReferenceTriangleWithTriangles( 
                                              PEL_Object* a_owner,
                                              size_t nb_sub_per_edge )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( a_owner, 
                                    GE_ReferenceTriangle::object(),
                                    GE_ReferenceTriangle::object(),
                                    GE_ReferenceSegment::object() )
   , NBS( nb_sub_per_edge )
{
   PEL_LABEL( 
   "GE_ReferenceTriangleWithTriangles:: GE_ReferenceTriangleWithTriangles" ) ;
                                                    
   size_t n_vertices = (NBS+1)*(NBS+2)/2 ;
   size_t n_subfaces = 3*NBS*(NBS+1)/2 ;
   size_t n_subcells = NBS*NBS ;
   size_t n_subfaces_per_face = NBS ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;
   
   size_t_array2D global_vertex( NBS+1, NBS+1 ) ;
   size_t i_vert = 0 ;
   for( size_t ix=0 ; ix!=NBS+1 ; ++ix )
   {
      for( size_t iy=0 ; ix+iy!=NBS+1 ; ++iy )
      {
         double xx = ((double) ix)/((double) NBS ) ;
         double yy = ((double) iy)/((double) NBS ) ;
         set_vertex( i_vert, GE_Point::create( this, xx, yy ) ) ;

         global_vertex( ix, iy ) = i_vert ;
         i_vert++ ;
      }
   }
   PEL_ASSERT( i_vert == (NBS+1)*(NBS+2)/2 ) ;
   
   size_t_vector vertex_indices( 3 ) ;
   
   size_t i_cell = 0 ;
   for( size_t ix=0 ; ix!=NBS ; ++ix )
   {
      for( size_t iy=0 ; ix+iy!=NBS ; ++iy )
      {
         vertex_indices( 0 ) = global_vertex( ix,   iy   ) ;
         vertex_indices( 1 ) = global_vertex( ix+1, iy   ) ;
         vertex_indices( 2 ) = global_vertex( ix,   iy+1 ) ;
         
         set_subcell( i_cell, vertex_indices ) ;
         i_cell++ ;

         if( ix+iy != NBS-1 )
         {
            vertex_indices( 0 ) = global_vertex( ix+1, iy+1 ) ;
            vertex_indices( 1 ) = global_vertex( ix,   iy+1 ) ;
            vertex_indices( 2 ) = global_vertex( ix+1, iy   ) ;

            set_subcell( i_cell, vertex_indices ) ;
            i_cell++ ;
         }
      }
   }
   PEL_ASSERT( i_cell == NBS*NBS ) ;
   
   PEL_CHECK( nb_subintervals_per_edge() == nb_sub_per_edge ) ;
}

//---------------------------------------------------------------------------
GE_ReferenceTriangleWithTriangles:: ~GE_ReferenceTriangleWithTriangles( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
size_t
GE_ReferenceTriangleWithTriangles:: nb_subintervals_per_edge( void ) const
//---------------------------------------------------------------------------
{
   return( NBS ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceTriangleWithTriangles:: compute_location_in_subcell( 
                                        size_t ic, 
                                        GE_Point const* pt_cell,
                                        GE_Point* pt_subcell ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( 
   "GE_ReferenceTriangleWithTriangles:: compute_location_in_subcell") ;
   PEL_CHECK_PRE( compute_location_in_subcell_PRE( ic, pt_cell, pt_subcell ) );
   
   PEL_ASSERT( NBS == 2  ) ;

   if( ic == 0 )
   {
      pt_subcell->set_coordinate( 0, 2.0*pt_cell->coordinate( 0 ) ) ;
      pt_subcell->set_coordinate( 1, 2.0*pt_cell->coordinate( 1 ) ) ;
   }
   else if( ic == 3 )
   {
      pt_subcell->set_coordinate( 0, 2.0*pt_cell->coordinate( 0 )-1.0 ) ;
      pt_subcell->set_coordinate( 1, 2.0*pt_cell->coordinate( 1 ) ) ;
   }
   else if( ic == 1 )
   {
      pt_subcell->set_coordinate( 0, 1.0-2.0*pt_cell->coordinate( 0 ) ) ;
      pt_subcell->set_coordinate( 1, 1.0-2.0*pt_cell->coordinate( 1 ) ) ;
   }
   else if( ic == 2 )
   {
      pt_subcell->set_coordinate( 0, 2.0*pt_cell->coordinate( 0 ) ) ;
      pt_subcell->set_coordinate( 1, 2.0*pt_cell->coordinate( 1 )-1.0 ) ;
   }
}
