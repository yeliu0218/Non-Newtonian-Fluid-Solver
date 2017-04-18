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

#include <GE_ReferenceCubeWithTetrahedra.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_ReferenceCube.hh>
#include <GE_ReferenceTetrahedron.hh>
#include <GE_ReferenceTriangle.hh>

#include <iostream>

GE_ReferenceCubeWithTetrahedra const*
GE_ReferenceCubeWithTetrahedra:: PROTOTYPE = 
                                  new GE_ReferenceCubeWithTetrahedra() ;

//---------------------------------------------------------------------------
GE_ReferenceCubeWithTetrahedra:: GE_ReferenceCubeWithTetrahedra( void )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( "GE_ReferenceCubeWithTetrahedra" ) 
{
}

//---------------------------------------------------------------------------
GE_ReferenceCubeWithTetrahedra const*
GE_ReferenceCubeWithTetrahedra:: create_replica( 
                                         PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) const 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceCubeWithTetrahedra:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   GE_ReferenceCubeWithTetrahedra const* result = 
               new GE_ReferenceCubeWithTetrahedra( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
GE_ReferenceCubeWithTetrahedra:: GE_ReferenceCubeWithTetrahedra( 
                                              PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : GE_ReferencePolyhedronRefiner( a_owner, 
                                    GE_ReferenceCube::object(),
                                    GE_ReferenceTetrahedron::object(),
                                    GE_ReferenceTriangle::object() ) 
{
   PEL_LABEL( 
   "GE_ReferenceCubeWithTetrahedra:: GE_ReferenceCubeWithTetrahedra" ) ;
   
   initialize() ;
}

//---------------------------------------------------------------------------
GE_ReferenceCubeWithTetrahedra:: ~GE_ReferenceCubeWithTetrahedra( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_ReferenceCubeWithTetrahedra:: compute_location_in_subcell( 
                                        size_t ic, 
                                        GE_Point const* pt_cell,
                                        GE_Point* pt_subcell ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL("GE_ReferenceCubeWithTetrahedra:: compute_location_in_subcell") ;
   PEL_CHECK_PRE( compute_location_in_subcell_PRE( ic, pt_cell, pt_subcell ) );

   PEL_Error::object()->raise_not_implemented( this, 
                                  "compute_location_in_subcell" ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferenceCubeWithTetrahedra:: initialize( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferenceCubeWithTetrahedra:: initialize_right" ) ;

   size_t n_vertices = 8 ;
   size_t n_subfaces = 18 ;
   size_t n_subcells = 6 ;
   size_t n_subfaces_per_face = 1 ;

   set_dimensions( n_vertices, n_subfaces, n_subcells, n_subfaces_per_face ) ;

   set_vertex( 0, GE_Point::create( this, 0.0, 0.0, 0.0 ) ) ;
   set_vertex( 1, GE_Point::create( this, 1.0, 0.0, 0.0 ) ) ;
   set_vertex( 2, GE_Point::create( this, 0.0, 1.0, 0.0 ) ) ;
   set_vertex( 3, GE_Point::create( this, 1.0, 1.0, 0.0 ) ) ;
   set_vertex( 4, GE_Point::create( this, 0.0, 0.0, 1.0 ) ) ;
   set_vertex( 5, GE_Point::create( this, 1.0, 0.0, 1.0 ) ) ;
   set_vertex( 6, GE_Point::create( this, 0.0, 1.0, 1.0 ) ) ;
   set_vertex( 7, GE_Point::create( this, 1.0, 1.0, 1.0 ) ) ;

   size_t_vector vertex_indices( 4 ) ;

   // sub-tetrahedron 0 
   vertex_indices( 0 ) = 0 ; 
   vertex_indices( 1 ) = 1 ; 
   vertex_indices( 2 ) = 4 ; 
   vertex_indices( 3 ) = 6 ;
   set_subcell( 0, vertex_indices ) ;
   
   // sub-tetrahedron 1 
   vertex_indices( 0 ) = 1 ; 
   vertex_indices( 1 ) = 4 ; 
   vertex_indices( 2 ) = 5 ; 
   vertex_indices( 3 ) = 6 ;
   set_subcell( 1, vertex_indices ) ;
   
   // sub-tetrahedron 2 
   vertex_indices( 0 ) = 1 ; 
   vertex_indices( 1 ) = 5 ; 
   vertex_indices( 2 ) = 6 ; 
   vertex_indices( 3 ) = 7 ;
   set_subcell( 2, vertex_indices ) ;
   
   // sub-tetrahedron 3 
   vertex_indices( 0 ) = 1 ; 
   vertex_indices( 1 ) = 3 ; 
   vertex_indices( 2 ) = 6 ; 
   vertex_indices( 3 ) = 7 ;
   set_subcell( 3, vertex_indices ) ;
   
   // sub-tetrahedron 4 
   vertex_indices( 0 ) = 0 ; 
   vertex_indices( 1 ) = 1 ; 
   vertex_indices( 2 ) = 2 ; 
   vertex_indices( 3 ) = 6 ;
   set_subcell( 4, vertex_indices ) ;
   
   // sub-tetrahedron 5 
   vertex_indices( 0 ) = 1 ; 
   vertex_indices( 1 ) = 2 ; 
   vertex_indices( 2 ) = 3 ; 
   vertex_indices( 3 ) = 6 ;
   set_subcell( 5, vertex_indices ) ;
}

