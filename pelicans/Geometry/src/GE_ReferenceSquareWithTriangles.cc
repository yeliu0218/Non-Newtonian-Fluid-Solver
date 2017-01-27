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

#include <GE_ReferenceSquareWithTriangles.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <stringVector.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_ReferenceTriangle.hh>

#include <iostream>

GE_ReferenceSquareWithTriangles const*
GE_ReferenceSquareWithTriangles:: PROTOTYPE = 
                                  new GE_ReferenceSquareWithTriangles() ;

//---------------------------------------------------------------------------
GE_ReferenceSquareWithTriangles:: GE_ReferenceSquareWithTriangles( void )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( "GE_ReferenceSquareWithTriangles" ) 
{
}

//---------------------------------------------------------------------------
GE_ReferenceSquareWithTriangles const*
GE_ReferenceSquareWithTriangles:: create_replica( 
                                         PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) const 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquareWithTriangles:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   GE_ReferenceSquareWithTriangles const* result = 
               new GE_ReferenceSquareWithTriangles( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferenceSquareWithTriangles:: GE_ReferenceSquareWithTriangles( 
                                              PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( a_owner, 
                                    GE_ReferenceSquare::object(),
                                    GE_ReferenceTriangle::object(),
                                    GE_ReferenceSegment::object() ) 
{
   PEL_LABEL( 
   "GE_ReferenceSquareWithTriangles:: GE_ReferenceSquareWithTriangles" ) ;
   
   std::string const& strat = exp->string_data( "strategy" ) ;
   exp->test_data_in( "strategy","/,\\,X,x,26_acute_triangles" ) ;
   if( strat == "/" )
   {
      initialize_right() ;
   }
   else if( strat == "\\" )
   {
      initialize_left() ;
   }
   else if( strat == "X" || strat == "x" )
   {
      initialize_cross() ;
   }
   else if( strat == "26_acute_triangles" )
   {
      initialize_26_acute_triangles() ;
   }
   else
   {
      PEL_Error::object()->raise_data_error( exp, "strategy",
                                             "invalid strategy: " + strat ) ;
   }
}

//---------------------------------------------------------------------------
GE_ReferenceSquareWithTriangles:: ~GE_ReferenceSquareWithTriangles( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_ReferenceSquareWithTriangles:: compute_location_in_subcell( 
                                        size_t ic, 
                                        GE_Point const* pt_cell,
                                        GE_Point* pt_subcell ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL("GE_ReferenceSquareWithTriangles:: compute_location_in_subcell") ;
   PEL_CHECK_PRE( compute_location_in_subcell_PRE( ic, pt_cell, pt_subcell ) );

   PEL_Error::object()->raise_not_implemented( this, 
                                  "compute_location_in_subcell" ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceSquareWithTriangles:: initialize_right( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquareWithTriangles:: initialize_right" ) ;

   size_t n_vertices = 4 ;
   size_t n_subfaces = 5 ;
   size_t n_subcells = 2 ;
   size_t n_subfaces_per_face = 1 ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;

   set_vertex( 0, GE_Point::create( this, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 1.0, 1.0 ) ) ;
   set_vertex( 3, GE_Point::create( this, 0.0, 1.0 ) ) ;

   size_t_vector vertex_indices( 3 ) ;

   // sub-triangle 0 
   vertex_indices( 0 ) = 0 ;
   vertex_indices( 1 ) = 1 ;
   vertex_indices( 2 ) = 2 ;
   set_subcell( 0, vertex_indices ) ;

   // sub-triangle 1
   vertex_indices( 0 ) = 2 ;
   vertex_indices( 1 ) = 3 ;
   vertex_indices( 2 ) = 0 ;
   set_subcell( 1, vertex_indices ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceSquareWithTriangles:: initialize_left( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquareWithTriangles:: initialize_left" ) ;

   size_t n_vertices = 4 ;
   size_t n_subfaces = 5 ;
   size_t n_subcells = 2 ;
   size_t n_subfaces_per_face = 1 ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;

   set_vertex( 0, GE_Point::create( this, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 1.0, 1.0 ) ) ;
   set_vertex( 3, GE_Point::create( this, 0.0, 1.0 ) ) ;

   size_t_vector vertex_indices( 3 ) ;

   // sub-triangle 0 
   vertex_indices( 0 ) = 0 ;
   vertex_indices( 1 ) = 1 ;
   vertex_indices( 2 ) = 3 ;
   set_subcell( 0, vertex_indices ) ;

   // sub-triangle 1
   vertex_indices( 0 ) = 1 ;
   vertex_indices( 1 ) = 2 ;
   vertex_indices( 2 ) = 3 ;
   set_subcell( 1, vertex_indices ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceSquareWithTriangles:: initialize_cross( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceSquareWithTriangles:: initialize_cross" ) ;

   size_t n_vertices = 5 ;
   size_t n_subfaces = 8 ;
   size_t n_subcells = 4 ;
   size_t n_subfaces_per_face = 1 ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;

   set_vertex( 0, GE_Point::create( this, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 1.0, 1.0 ) ) ;
   set_vertex( 3, GE_Point::create( this, 0.0, 1.0 ) ) ;
   set_vertex( 4, GE_Point::create( this, 0.5, 0.5 ) ) ;

   size_t_vector vertex_indices( 3 ) ;

   // sub-triangle 0 
   vertex_indices( 0 ) = 4 ;
   vertex_indices( 1 ) = 0 ;
   vertex_indices( 2 ) = 1 ;
   set_subcell( 0, vertex_indices ) ;

   // sub-triangle 1
   vertex_indices( 0 ) = 4 ;
   vertex_indices( 1 ) = 1 ;
   vertex_indices( 2 ) = 2 ;
   set_subcell( 1, vertex_indices ) ;

   // sub-triangle 2
   vertex_indices( 0 ) = 4 ;
   vertex_indices( 1 ) = 2 ;
   vertex_indices( 2 ) = 3 ;
   set_subcell( 2, vertex_indices ) ;

   // sub-triangle 3
   vertex_indices( 0 ) = 4 ;
   vertex_indices( 1 ) = 3 ;
   vertex_indices( 2 ) = 0 ;
   set_subcell( 3, vertex_indices ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceSquareWithTriangles:: initialize_26_acute_triangles( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( 
   "GE_ReferenceSquareWithTriangles:: initialize_26_acute_triangles" ) ;

   size_t n_vertices = 20 ;
   size_t n_subfaces = 45 ;
   size_t n_subcells = 26 ;
   size_t n_subfaces_per_face = 3 ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;

   set_vertex( 0,  GE_Point::create( this, 0.    , 0.    ) ) ;
   set_vertex( 1,  GE_Point::create( this, 1./3. , 0.0   ) ) ;
   set_vertex( 2,  GE_Point::create( this, 1./5. , 1./4. ) ) ;
   set_vertex( 3,  GE_Point::create( this, 0.    , 1./3. ) ) ;
   set_vertex( 4,  GE_Point::create( this, 1./4. , 1./2. ) ) ;
   set_vertex( 5,  GE_Point::create( this, 0.    , 2./3. ) ) ;
   set_vertex( 6,  GE_Point::create( this, 1./5. , 3./4. ) ) ;
   set_vertex( 7,  GE_Point::create( this, 0.    , 1.    ) ) ;
   set_vertex( 8,  GE_Point::create( this, 1./3. , 1.    ) ) ;
   set_vertex( 9,  GE_Point::create( this, 1./2. , 1./3. ) ) ;
   set_vertex( 10, GE_Point::create( this, 1./2. , 2./3. ) ) ;
   set_vertex( 11, GE_Point::create( this, 2./3. , 0.    ) ) ;
   set_vertex( 12, GE_Point::create( this, 3./4. , 1./2. ) ) ;
   set_vertex( 13, GE_Point::create( this, 2./3. , 1.    ) ) ;
   set_vertex( 14, GE_Point::create( this, 4./5. , 1./4. ) ) ;
   set_vertex( 15, GE_Point::create( this, 4./5. , 3./4. ) ) ;
   set_vertex( 16, GE_Point::create( this, 1.    , 0.    ) ) ;
   set_vertex( 17, GE_Point::create( this, 1.    , 1./3. ) ) ;
   set_vertex( 18, GE_Point::create( this, 1.    , 2./3. ) ) ;
   set_vertex( 19, GE_Point::create( this, 1.    , 1.    ) ) ;

   size_t_vector vert_indices( 3 ) ;

   // sub-triangle 0
   vert_indices( 0 ) = 0 ; vert_indices( 1 ) = 1 ; vert_indices( 2 ) = 2  ; 
   set_subcell( 0, vert_indices ) ;

   // sub-triangle 1
   vert_indices( 0 ) = 0 ; vert_indices( 1 ) = 2 ; vert_indices( 2 ) = 3  ; 
   set_subcell( 1, vert_indices ) ;

   // sub-triangle 2
   vert_indices( 0 ) = 3 ; vert_indices( 1 ) = 2 ; vert_indices( 2 ) = 4  ; 
   set_subcell( 2, vert_indices ) ;

   // sub-triangle 3
   vert_indices( 0 ) = 3 ; vert_indices( 1 ) = 4 ; vert_indices( 2 ) = 5  ; 
   set_subcell( 3, vert_indices ) ;

   // sub-triangle 4
   vert_indices( 0 ) = 5 ; vert_indices( 1 ) = 4  ; vert_indices( 2 ) = 6  ; 
   set_subcell( 4, vert_indices ) ;

   // sub-triangle 5
   vert_indices( 0 ) = 5 ; vert_indices( 1 ) = 6  ; vert_indices( 2 ) = 7  ; 
   set_subcell( 5, vert_indices ) ;

   // sub-triangle 6
   vert_indices( 0 ) = 7 ; vert_indices( 1 ) = 6  ; vert_indices( 2 ) = 8  ; 
   set_subcell( 6, vert_indices ) ;

   // sub-triangle 7
   vert_indices( 0 ) = 1 ; vert_indices( 1 ) = 9  ; vert_indices( 2 ) = 2  ; 
   set_subcell( 7, vert_indices ) ;

   // sub-triangle 8
   vert_indices( 0 ) = 2 ; vert_indices( 1 ) = 9  ; vert_indices( 2 ) = 4  ; 
   set_subcell( 8, vert_indices ) ;

   // sub-triangle 9
   vert_indices( 0 ) = 4 ; vert_indices( 1 ) = 9  ; vert_indices( 2 ) = 10 ; 
   set_subcell( 9, vert_indices ) ;

   // sub-triangle 10
   vert_indices( 0 ) = 4 ; vert_indices( 1 ) = 10 ; vert_indices( 2 ) = 6  ; 
   set_subcell( 10, vert_indices ) ;

   // sub-triangle 11
   vert_indices( 0 ) = 6 ; vert_indices( 1 ) = 10 ; vert_indices( 2 ) = 8  ; 
   set_subcell( 11, vert_indices ) ;

   // sub-triangle 12
   vert_indices( 0 ) = 1 ; vert_indices( 1 ) = 11 ; vert_indices( 2 ) = 9  ; 
   set_subcell( 12, vert_indices ) ;

   // sub-triangle 13
   vert_indices( 0 ) = 9 ; vert_indices( 1 ) = 12 ; vert_indices( 2 ) = 10 ; 
   set_subcell( 13, vert_indices ) ;

   // sub-triangle 14
   vert_indices( 0 ) = 10 ; vert_indices( 1 ) = 13 ; vert_indices( 2 ) = 8  ; 
   set_subcell( 14, vert_indices ) ;

   // sub-triangle 15
   vert_indices( 0 ) = 11 ; vert_indices( 1 ) = 14 ; vert_indices( 2 ) = 9  ; 
   set_subcell( 15, vert_indices ) ;

   // sub-triangle 16
   vert_indices( 0 ) = 9  ; vert_indices( 1 ) = 14 ; vert_indices( 2 ) = 12 ; 
   set_subcell( 16, vert_indices ) ;

   // sub-triangle 17
   vert_indices( 0 ) = 10 ; vert_indices( 1 ) = 12 ; vert_indices( 2 ) = 15 ; 
   set_subcell( 17, vert_indices ) ;

   // sub-triangle 18
   vert_indices( 0 ) = 10 ; vert_indices( 1 ) = 15 ; vert_indices( 2 ) = 13 ; 
   set_subcell( 18, vert_indices ) ;

   // sub-triangle 18
   vert_indices( 0 ) = 11 ; vert_indices( 1 ) = 16 ; vert_indices( 2 ) = 14 ; 
   set_subcell( 19, vert_indices ) ;

   // sub-triangle 19
   vert_indices( 0 ) = 14 ; vert_indices( 1 ) = 17 ; vert_indices( 2 ) = 12 ; 
   set_subcell( 20, vert_indices ) ;

   // sub-triangle 20
   vert_indices( 0 ) = 12 ; vert_indices( 1 ) = 18 ; vert_indices( 2 ) = 15 ; 
   set_subcell( 21, vert_indices ) ;

   // sub-triangle 21
   vert_indices( 0 ) = 15 ; vert_indices( 1 ) = 19 ; vert_indices( 2 ) = 13 ; 
   set_subcell( 22, vert_indices ) ;

   // sub-triangle 22
   vert_indices( 0 ) = 14 ; vert_indices( 1 ) = 16 ; vert_indices( 2 ) = 17 ; 
   set_subcell( 23, vert_indices ) ;

   // sub-triangle 23
   vert_indices( 0 ) = 12 ; vert_indices( 1 ) = 17 ; vert_indices( 2 ) = 18 ; 
   set_subcell( 24, vert_indices ) ;

   // sub-triangle 24
   vert_indices( 0 ) = 15 ; vert_indices( 1 ) = 18 ; vert_indices( 2 ) = 19 ; 
   set_subcell( 25, vert_indices ) ;
}
