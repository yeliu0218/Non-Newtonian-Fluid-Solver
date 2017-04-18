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

#include <GE_ReferencePolyhedron.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <iostream>

size_t GE_ReferencePolyhedron::NB_INSTANCES = 0 ;
   
//-----------------------------------------------------------------------------
GE_ReferencePolyhedron:: GE_ReferencePolyhedron( std::string const& a_name,
                                                 size_t const a_nb_vertices,
                                                 size_t const a_nb_faces,
                                                 size_t const a_dimension,
                                                 double const a_measure )
//-----------------------------------------------------------------------------
   : PEL_Object( PEL_Root::object() )
   , NAME( a_name )
   , ID( NB_INSTANCES ) 
   , POLY_DIM( a_dimension )
   , NB_VERTS( a_nb_vertices )
   , VERTS( PEL_Vector::create( this, a_nb_vertices ) )
   , NB_FACES( a_nb_faces )
   , NB_FVERTS( a_nb_faces )
   , IDX_FVERTS( 0, a_nb_faces )
   , F_NORMALS( PEL_Vector::create( this, a_nb_faces ) )
   , CENTER( 0 )
   , MEASURE( a_measure )
{
   ++NB_INSTANCES ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron:: ~GE_ReferencePolyhedron( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_ReferencePolyhedron:: name( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( NAME ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: id_number( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: id_number" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t result = ID ;
   
   PEL_CHECK_POST( result < nb_objects() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: nb_objects( void )
//-----------------------------------------------------------------------------
{
   return( NB_INSTANCES ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: dimension( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: dimension" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( POLY_DIM ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: nb_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( NB_VERTS ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_ReferencePolyhedron:: vertex( size_t i ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: vertex" ) ;
   PEL_CHECK_PRE( i<nb_vertices() ) ;
   PEL_CHECK_INV( invariant() ) ;
   GE_Point const* result =
                   static_cast<GE_Point const*>( VERTS->at(i) ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->nb_coordinates()==dimension() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_ReferencePolyhedron:: measure( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: measure" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( MEASURE ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: nb_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: nb_faces" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( NB_FACES ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: nb_face_vertices( size_t i_face ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: nb_face_vertices" ) ;
   PEL_CHECK_PRE( i_face < nb_faces() ) ;

   size_t result = NB_FVERTS( i_face ) ;

   PEL_CHECK_POST( result < nb_vertices() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_ReferencePolyhedron:: face_vertex( size_t i_face, size_t i_vert ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: face_vertex" ) ;
   PEL_CHECK_PRE( i_face < nb_faces() ) ;
   PEL_CHECK_PRE( i_vert < nb_face_vertices( i_face ) ) ;

   size_t result = IDX_FVERTS( i_vert, i_face ) ;

   PEL_CHECK_POST( result < nb_vertices() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_ReferencePolyhedron:: face_outward_normal( size_t i_face ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: face_outward_normal" ) ;
   PEL_CHECK_PRE( i_face < nb_faces() ) ;

   GE_Vector const* result =
                    static_cast<GE_Vector const*>( F_NORMALS->at( i_face ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_components() == dimension() ) ;
   PEL_CHECK_POST( FORMAL( result->norm() == 1. ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_ReferencePolyhedron:: center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: center" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( CENTER==0 )
   {
      CENTER = GE_Point::create( PEL_Root::object(), dimension() ) ;
      for( size_t i=0 ; i<dimension() ; ++i )
      {
         double coord = 0. ;
         for( size_t j=0 ; j<nb_vertices() ; ++j )
         {
            coord += vertex(j)->coordinate(i) ;
         }
         CENTER->set_coordinate( i, coord/nb_vertices() ) ;
      }
      
   }
   GE_Point const* result = CENTER ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==PEL_Root::object() ) ;
   PEL_CHECK_POST( result->nb_coordinates()==dimension() ) ;
   PEL_CHECK_POST( contains( result ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_ReferencePolyhedron:: epsilon( void )
//-----------------------------------------------------------------------------
{
   return( 1.E-5 ) ;
}

//-----------------------------------------------------------------------------
void
GE_ReferencePolyhedron:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string const space( indent_width, ' ' ) ;
   os << space << name() << " : " << std::endl ;
   for( size_t i=0 ; i<nb_vertices() ; i++ )
   {
      vertex( i )->print( os, indent_width+3 ) ;
      os << std::endl ;
   }
}

//-----------------------------------------------------------------------------
void
GE_ReferencePolyhedron:: set_vertex( size_t i, GE_Point* pt )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: set_vertex" ) ;
   PEL_CHECK( i < nb_vertices() ) ;
   PEL_CHECK( pt != 0 ) ;
   PEL_CHECK( pt->owner() == this ) ;
   PEL_CHECK( pt->nb_coordinates() == dimension() ) ;
   PEL_CHECK( 
      FORALL( ( size_t j=0 ; j<dimension() ; ++j ),
              pt->coordinate( j ) == 0.0 ||
              pt->coordinate( j ) == 1.0 ) ) ;

   VERTS->set_at( i, pt ) ;
}

//-----------------------------------------------------------------------------
void
GE_ReferencePolyhedron:: set_face_normal( size_t i_face, GE_Vector* n )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: set_face_normal" ) ;
   PEL_CHECK( i_face < nb_faces() ) ;
   PEL_CHECK( n != 0 ) ;
   PEL_CHECK( n->owner() == this ) ;
   PEL_CHECK( n->nb_components() == dimension() ) ;

   F_NORMALS->set_at( i_face, n ) ;
}

//-----------------------------------------------------------------------------
void 
GE_ReferencePolyhedron:: append_face_vertex( size_t i_face, size_t vertex_id )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron:: append_face_vertex" ) ;
   PEL_CHECK( i_face < nb_faces() ) ;
   PEL_CHECK( vertex_id < nb_vertices() ) ;

   size_t idx = NB_FVERTS( i_face ) ;
   ++NB_FVERTS( i_face ) ;
   if( idx >= IDX_FVERTS.index_bound( 0 ) )
      IDX_FVERTS.raise_first_index_bound( idx+1 ) ;
   IDX_FVERTS( idx, i_face ) = vertex_id ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: project_PRE( GE_Point const* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref!=0 &&
               pt_ref->nb_coordinates()==dimension() &&
               !contains( pt_ref ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: project_POST( GE_Point const* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( contains( pt_ref ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: contains_PRE( GE_Point const* pt_ref,
                                       double tol ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates()==dimension() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: face_contains_PRE( size_t i_face, 
                                            GE_Point const* pt_ref,
                                            double tol ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( i_face < nb_faces() ) ;
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates()==dimension() ) ;
   PEL_ASSERT( tol > 0. ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: build_neighbor_PRE( GE_Point const* pt_ref,
                                             size_t ic,
                                             GE_Point const* neighbor ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref!=0 &&
               pt_ref->nb_coordinates()==dimension() &&
               contains( pt_ref ) ) ;
   PEL_ASSERT( ic<dimension() ) ;
   PEL_ASSERT( neighbor!=0 &&
               neighbor->nb_coordinates()==dimension() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_ReferencePolyhedron:: build_neighbor_POST( GE_Point const* pt_ref,
                                              size_t ic,
                                              GE_Point const* neighbor,
                                              double delta ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( contains( neighbor ) ) ;
   PEL_ASSERT( PEL::equal( neighbor->coordinate(ic),
                           pt_ref->coordinate(ic)+delta ) ) ;
   return( true ) ;
}
