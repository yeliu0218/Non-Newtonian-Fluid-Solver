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

#include <GE_Mpolyhedron.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::cout ; using std::endl ;
using std::ostringstream ;

bool GE_Mpolyhedron::CHECK_CONSISTENCY = true ;
bool GE_Mpolyhedron::UPDATING = false ;

struct GE_Mpolyhedron_ERROR
{
   static void n0( GE_Mpolyhedron const* poly_1,
                   GE_Mpolyhedron const* poly_2 ) ;
} ;

//-----------------------------------------------------------------------------
GE_Mpolyhedron*
GE_Mpolyhedron:: create( std::string const& a_name,
                         GE_SetOfPoints* a_set_of_vertices,
                         size_t_vector const& a_vertices_index_table )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_set_of_vertices!=0 ) ;
   PEL_CHECK_PRE( a_vertices_index_table.size()==
                         reference_polyhedron(a_name)->nb_vertices() ) ;
   PEL_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<a_vertices_index_table.size() ; ++i ),
              a_vertices_index_table(i)<a_set_of_vertices->nb_points() ) ) ;
   PEL_CHECK_PRE( a_set_of_vertices->dimension()>=
                        reference_polyhedron(a_name)->dimension() ) ;

   GE_Mpolyhedron const* prototype =
      static_cast<GE_Mpolyhedron const*>( plugins_map()->item( a_name ) ) ;

   PEL_Vector* vert = PEL_Vector::create( 0, a_vertices_index_table.size() ) ;
   for( size_t i=0; i<a_vertices_index_table.size(); i++ )
   {
      vert->set_at( i, const_cast<GE_Point*>(
                    a_set_of_vertices->point( a_vertices_index_table(i) ) ) ) ;
   }

   GE_Mpolyhedron* result = prototype->create_replica( a_set_of_vertices, vert ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_set_of_vertices ) ;
   PEL_CHECK_POST( result->name()==a_name ) ;
   PEL_CHECK_POST( result->reference_polyhedron()==
                                            reference_polyhedron(a_name) ) ;
   PEL_CHECK_POST( result->dimension()==
                               reference_polyhedron(a_name)->dimension() ) ;
   PEL_CHECK_POST( result->nb_space_dimensions()==
                                           a_set_of_vertices->dimension() ) ;
   PEL_CHECK_POST( result->nb_vertices()==a_vertices_index_table.size() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<result->nb_vertices() ; ++i ),
              result->vertex(i)==
                   a_set_of_vertices->point( a_vertices_index_table(i) ) ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Mpolyhedron*
GE_Mpolyhedron:: create( PEL_Object* a_owner,
                         std::string const& a_name,
                         PEL_Vector const* vertices )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( vertices!=0 ) ;
   PEL_CHECK_PRE( vertices->count()==
                                reference_polyhedron(a_name)->nb_vertices() ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=0 ; i<vertices->count() ; i++ ),
	            dynamic_cast<GE_Point const*>( vertices->at( i ) )!=0 ) ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=1 ; i<vertices->count() ; i++ ),
     static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates()==
     static_cast<GE_Point const*>( vertices->at( i ) )->nb_coordinates() ) ) ;

   GE_Mpolyhedron const* prototype =
      static_cast<GE_Mpolyhedron const*>( plugins_map()->item( a_name ) ) ;

   GE_Mpolyhedron* result = prototype->create_replica( a_owner,
                                              vertices->create_clone( 0 ) ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->name()==a_name ) ;
   PEL_CHECK_POST( result->reference_polyhedron()==
                                            reference_polyhedron(a_name) ) ;
   PEL_CHECK_POST( result->dimension()==
                               reference_polyhedron(a_name)->dimension() ) ;
   PEL_CHECK_POST( result->nb_vertices()==vertices->count() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<result->nb_vertices() ; i++ ),
    result->vertex(i)== static_cast<GE_Point const*>( vertices->at( i ) ) ) ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Mpolyhedron:: GE_Mpolyhedron( std::string const& a_name )
//-----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , VERTS( 0 )
   , CENTER( 0 )
{
   PEL_LABEL( "GE_Mpolyhedron:: GE_Mpolyhedron" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK( is_prototype() ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Mpolyhedron:: GE_Mpolyhedron( PEL_Object* a_owner,
                                 PEL_Vector* vertices )
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , VERTS( vertices )
   , CENTER( 0 )
{
   PEL_LABEL( "GE_Mpolyhedron:: GE_Mpolyhedron" ) ;
   PEL_CHECK_PRE( !is_prototype() ) ;
   PEL_CHECK_PRE( vertices != 0 ) ;
   PEL_CHECK_PRE( vertices->owner() == 0 ) ;

   vertices->set_owner( this ) ;

   PEL_CHECK_POST( vertices->owner() == this ) ;
}

//-----------------------------------------------------------------------------
GE_Mpolyhedron:: ~GE_Mpolyhedron( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: ~GE_Mpolyhedron" ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: update( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: update" ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   UPDATING = true ;
   update_internal() ;
   UPDATING = false ;

   if( CENTER!=0 )
   {
      apply_mapping( reference_polyhedron()->center(), CENTER ) ;
   }

   if( CHECK_CONSISTENCY )
   {
      if( !is_consistent( std::cout, false ) )
      {
         std::ostringstream message ;
         message << std::endl ;
         message << "The polyhedron : " << std::endl ;
         UPDATING = true ; // Bypass INVARIANT for printing...
         print( message, 5 ) ;
         UPDATING = false ;
         message << "is not valid." << std::endl << std::endl ;
         message << "The following checks have been performed :" << std::endl ;
         is_consistent( message, true ) ;
         PEL_Error::object()->raise_plain( message.str() ) ;
      }
   }
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: check_consistency( void )
//-----------------------------------------------------------------------------
{
   return( CHECK_CONSISTENCY ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: unset_check_consistency( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: unset_check_consistency" ) ;
   CHECK_CONSISTENCY = false ;
   PEL_CHECK_POST( !check_consistency() ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: set_check_consistency( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: set_check_consistency" ) ;
   CHECK_CONSISTENCY = true ;
   PEL_CHECK_POST( check_consistency() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: is_equal( PEL_Object const* other ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( three_way_comparison( other ) == 0 ) ;

   PEL_CHECK_INV( invariant() ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
int
GE_Mpolyhedron:: three_way_comparison( PEL_Object const* other ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* poly =
                          static_cast<GE_Mpolyhedron const*>( other ) ;

   size_t const nbVert = nb_vertices() ;

   int result = nbVert - poly->nb_vertices() ;
   if( result == 0 )
   {
      for( size_t i=0 ; i<nbVert && (result == 0) ; i++ )
      {
         result = vertex(i)->three_way_comparison( poly->vertex( i ) ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_Mpolyhedron:: hash_code( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: hash_code" ) ;
   PEL_Error::object()->raise_not_implemented( this, "hash_code" ) ;
   return( 0 ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: reorder_vertices_according_to( GE_Mpolyhedron const* other )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: reorder_vertices_according_to" ) ;

   size_t nbvs = nb_vertices() ;
   bool ok = ( other->nb_vertices() == nbvs ) ;

   if( !ok ) GE_Mpolyhedron_ERROR::n0( this, other ) ;

   for( size_t i=0 ; i<nbvs ; ++i )
   {
      GE_Point const* other_v = other->vertex( i ) ;
      ok = false ;
      for( size_t k=i ; k<nbvs ; ++k )
      {
         GE_Point* vv = static_cast<GE_Point*>( VERTS->at( k ) ) ;
         if( vv->distance( other_v ) < 1.e-8 )
         {
            if( k != i )
            {
               VERTS->set_at( k, VERTS->at( i ) ) ;
               VERTS->set_at( i, vv ) ;
            }
            ok = true ;
            break ;
         }
      }
      if( !ok ) GE_Mpolyhedron_ERROR::n0( this, other ) ;
   }

   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<nb_vertices() ; i++ ),
              vertex(i)->distance( other->vertex(i) ) < 1.e-8 ) ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_Mpolyhedron:: nb_space_dimensions( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: nb_space_dimensions" ) ;
   PEL_CHECK( !is_prototype() ) ;

   PEL_CHECK( dynamic_cast<GE_Point const*>( VERTS->at( 0 ) ) != 0 ) ;
   size_t result = static_cast<GE_Point*>( VERTS->at( 0 ) )->nb_coordinates() ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_Mpolyhedron:: dimension( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: dimension" ) ;
   return( reference_polyhedron()->dimension() ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_Mpolyhedron:: nb_vertices( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: nb_vertices" ) ;
   return( reference_polyhedron()->nb_vertices() ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_Mpolyhedron:: nb_faces( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: nb_faces" ) ;
   return( reference_polyhedron()->nb_faces() ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Mpolyhedron:: vertex( size_t i ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: vertex" ) ;
   PEL_CHECK_PRE( i<nb_vertices() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() ) ;

   PEL_CHECK( dynamic_cast<GE_Point const*>( VERTS->at( i ) ) != 0 ) ;
   GE_Point const* result = static_cast<GE_Point const*>( VERTS->at(i) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->nb_coordinates() == nb_space_dimensions() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_Mpolyhedron:: inter_vertices_maximum_distance( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: inter_vertices_maximum_distance" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   size_t const nbVertices = nb_vertices() ;

   double maxLength = 0.0 ;
   for( size_t iV1=0 ; iV1<nbVertices ; iV1++ )
   {
      for( size_t iV2=iV1+1 ; iV2<nbVertices ; iV2++ )
      {
         double d = vertex(iV1)->distance( vertex(iV2) ) ;
         if( d>maxLength )
         {
            maxLength = d ;
         }
      }
   }

   return( maxLength ) ;
}

//-----------------------------------------------------------------------------
double
GE_Mpolyhedron:: inter_vertices_maximum_distance( size_t dir ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: inter_vertices_maximum_distance(dir)" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( dir<nb_space_dimensions() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   size_t const nbVertices = nb_vertices() ;
   double xMin =  PEL::max_double() ;
   double xMax = -PEL::max_double() ;

   for( size_t i=0 ; i<nbVertices ; ++i )
   {
      double const x = vertex(i)->coordinate(dir) ;
      xMin = PEL::min( xMin, x ) ;
      xMax = PEL::max( xMax, x ) ;
   }

   return( xMax-xMin ) ;
}

//-----------------------------------------------------------------------------
double
GE_Mpolyhedron:: equivalent_ball_diameter( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: equivalent_ball_diameter" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   size_t const nbDims = nb_space_dimensions() ;

   double h = 0.0 ;
   if ( nbDims==2 )
   {
      h = 2.0*PEL::sqrt(measure()/PEL::pi()) ;
   }
   else if ( nbDims==3)
   {
      h = 2.0*PEL::pow( 3.0*measure()/4.0/PEL::pi(), 1.0/3.0 ) ;
   }
   else if ( nbDims==1 )
   {
      h = measure()/PEL::pi() ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "Invalid number of space dimensions" ) ;
   }

   // The subsequent formula could be used. It works for any dimension
   // but it is less readable
   //
   // h = pow( pow(2,nbDims-1)*nbDims*measure()/Pi/max(1,nbDims-1),
   //          1.0/nbDims ) );

   return( h ) ;
}

//-----------------------------------------------------------------------------
double
GE_Mpolyhedron:: reference_distance( GE_Point const* pt1,
                                     GE_Point const* pt2 ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: reference_distance" ) ;
   PEL_CHECK_PRE( pt1!=0 && pt1->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_CHECK_PRE( pt2!=0 && pt2->nb_coordinates()==nb_space_dimensions() ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_space_dimensions() ; i++ )
   {
      double const w = inter_vertices_maximum_distance(i) ;
      if( w>1.E-8 )
      {
         double v = (pt1->coordinate(i)-pt2->coordinate(i))/w ;
         result += v*v ;
      }
   }
   result = PEL::sqrt( result ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Mpolyhedron:: center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: center" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !is_prototype() && !is_updating() ) ;

   // lazy initialization of center
   // fake this approach : creation of a pointer-to-non-const that
   // points to the same object as this does [Meyers p. 89]
   if( CENTER == 0 )
   {
      GE_Mpolyhedron* fakeThis = const_cast<GE_Mpolyhedron*>(this) ;
      CENTER = GE_Point::create( fakeThis, nb_space_dimensions() ) ;
      apply_mapping( reference_polyhedron()->center(), CENTER ) ;
   }
   GE_Point const* result = CENTER ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_Mpolyhedron:: finite_volume_center( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: finite_volume_center" ) ;
   PEL_CHECK_PRE( finite_volume_center_PRE() ) ;

   GE_Point const* result = 0 ;

   PEL_CHECK_POST( finite_volume_center_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Vector const*
GE_Mpolyhedron:: unit_normal( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: unit_normal" ) ;
   PEL_CHECK_PRE( unit_normal_PRE() ) ;

   PEL_Error::object()->raise_not_implemented( this, "unit_normal" ) ;

   GE_Vector const* result = 0 ;

   PEL_CHECK_POST( unit_normal_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Mpolyhedron:: reference_polyhedron( std::string const& a_name )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: reference_polyhedron" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   GE_Mpolyhedron const* prototype =
      static_cast<GE_Mpolyhedron const*>( plugins_map()->item( a_name ) ) ;
   PEL_CHECK( prototype->is_prototype() ) ;

   GE_ReferencePolyhedron const* result = prototype->reference_polyhedron() ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}


//??? this function would be useless if build_mapping_hessian was coded
//??? in the derived classes (and it would be faster as it is for
//??? build_mappinf and build_(tr_)mapping_derivative
//----------------------------------------------------------------------------
void
GE_Mpolyhedron:: geobfs_2_mapping_hessian( doubleArray3D const& d2Ngeom,
                                           doubleArray3D* hessian ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: geobfs_2_mapping_hessian" ) ;
   PEL_CHECK_PRE( d2Ngeom.index_bound( 0 ) == nb_vertices() ) ;
   PEL_CHECK_PRE( d2Ngeom.index_bound( 1 ) == dimension() ) ;
   PEL_CHECK_PRE( d2Ngeom.index_bound( 2 ) == dimension() ) ;
   PEL_CHECK_PRE( hessian != 0 ) ;
   PEL_CHECK_PRE( hessian->index_bound( 0 ) == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( hessian->index_bound( 1 ) == dimension() ) ;
   PEL_CHECK_PRE( hessian->index_bound( 2 ) == dimension() ) ;

   size_t const nbv = nb_vertices() ;
   size_t const refe_dim = dimension() ;
   size_t const real_dim = nb_space_dimensions() ;

   hessian->set( 0.0 ) ;
   for( size_t iv = 0 ; iv<nbv ; ++iv )
   {
      GE_Point const* v = vertex( iv ) ;
      for( size_t i1=0 ; i1<real_dim ; ++i1 )
      {
         double const x_vert = v->coordinate( i1 ) ;
         for( size_t d1=0 ; d1<refe_dim ; ++d1 )
         {
            for( size_t d2=0 ; d2<refe_dim ; ++d2 )
            {
               (*hessian)( i1, d1, d2 ) += x_vert * d2Ngeom( iv, d1, d2 ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------------
void
GE_Mpolyhedron:: build_mapping_hessian( GE_Point const* pt_ref,
                                        doubleArray3D* hessian,
                                        bool& nonzero_hessian ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: build_mapping_hessian" ) ;
   PEL_CHECK_PRE( build_mapping_hessian_PRE( pt_ref, hessian, nonzero_hessian ) ) ;

   nonzero_hessian = false ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string const space( indent_width, ' ' ) ;
   os << space << name() << " : " ;
   if( !is_prototype() )
   {
      os << std::endl ;
      for( size_t i=0 ; i<nb_vertices() ; i++ )
      {
         os << space << "   V" << i << " : " ;
         vertex( i )->print( os, 0 ) ;
         os << std::endl ;
      }
   }
   else
   {
      os << "prototype" << std::endl ;
   }
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: is_updating( void )
//-----------------------------------------------------------------------------
{
   return( UPDATING ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: is_prototype( void ) const
//-----------------------------------------------------------------------------
{
   return( VERTS == 0 ) ;
}

//-----------------------------------------------------------------------------
void
GE_Mpolyhedron:: display_check( std::ostream& os,
                                std::string const& displayed_name,
                                bool success ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Mpolyhedron:: display_check" ) ;

   os << std::setw(40) << displayed_name << " :" ;
   if( success )
   {
      os << "  OK" ;
   }
   else
   {
      os << " FAIL" ;
   }
   os << endl ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: finite_volume_center_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( dimension() == nb_space_dimensions() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: finite_volume_center_POST( GE_Point const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result == 0 ||
               result->nb_coordinates() == nb_space_dimensions() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: contains_PRE( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt!=0 && pt->nb_coordinates()==nb_space_dimensions() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: measure_POST( double result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result>=0. ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: unit_normal_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( dimension()==nb_space_dimensions()-1 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: unit_normal_POST( GE_Vector const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->nb_components()==nb_space_dimensions() ) ;
   PEL_ASSERT( PEL::equal( result->norm(), 1. ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: apply_mapping_PRE( GE_Point const* pt_ref,
                                    GE_Point* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( pt != 0 ) ;
   PEL_ASSERT( pt->nb_coordinates() == nb_space_dimensions() ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: apply_mapping_POST( GE_Point const* pt_ref,
                                     GE_Point* pt ) const
//-----------------------------------------------------------------------------
{
   //??????????????????????????????????
   // SUPPRIMER cette postcondition et celle de inverse_mapping
   // le contains utilise lui-meme apply_mapping etc
//   PEL_ASSERT( IMPLIES( reference_polyhedron()->contains( pt_ref ),
//                        contains( pt ) ) ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: apply_inverse_mapping_PRE( GE_Point const* pt,
                                            GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt != 0 ) ;
   PEL_ASSERT( pt->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: apply_inverse_mapping_POST( GE_Point const* pt,
                                             GE_Point* pt_ref ) const
//-----------------------------------------------------------------------------
{
//   PEL_ASSERT( IMPLIES( contains( pt ),
//                        reference_polyhedron()->contains( pt_ref ) ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: build_mapping_derivative_PRE( GE_Point const* pt_ref,
                                               GE_Matrix* tjac ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref!= 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( tjac!=0 ) ;
   PEL_ASSERT( tjac->nb_rows()==nb_space_dimensions() ) ;
   PEL_ASSERT( tjac->nb_cols()==dimension() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: build_tr_mapping_derivative_PRE( GE_Point const* pt_ref,
                                                  GE_Matrix* tjac ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref!= 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( tjac!=0 ) ;
   PEL_ASSERT( tjac->nb_rows()==dimension() ) ;
   PEL_ASSERT( tjac->nb_cols()==nb_space_dimensions() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: build_mapping_hessian_PRE( GE_Point const* pt_ref,
                                            doubleArray3D* hessian,
                                            bool& nonzero_hessian ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( pt_ref!= 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
//   PEL_ASSERT( reference_polyhedron()->contains( pt_ref ) ) ;
   PEL_ASSERT( hessian != 0 ) ;
   PEL_ASSERT( hessian->index_bound( 0 ) == nb_space_dimensions() ) ;
   PEL_ASSERT( hessian->index_bound( 1 ) == dimension() ) ;
   PEL_ASSERT( hessian->index_bound( 2 ) == dimension() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: update_internal_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( !is_prototype() && is_updating() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: reference_polyhedron_PRE( void ) const
//-----------------------------------------------------------------------------
{
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: reference_polyhedron_POST(
                                   GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: create_replica_PRE( PEL_Object* a_owner,
                                     PEL_Vector* vertices ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_prototype() ) ;
   PEL_ASSERT( vertices != 0 ) ;
   PEL_ASSERT( vertices->owner() == 0 ) ;
   PEL_ASSERT( vertices->count() == nb_vertices() ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=0 ; i<vertices->count() ; i++ ),
         dynamic_cast<GE_Point const*>( vertices->at( i ) )!=0 ) ) ;
   PEL_ASSERT(
      static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates()==1 ||
      static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates()==2 ||
      static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates()==3 ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=1 ; i<vertices->count() ; i++ ),
         static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates()==
         static_cast<GE_Point const*>( vertices->at( i ) )->nb_coordinates() ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: create_replica_POST( PEL_Object* a_owner,
                                      PEL_Vector* vertices,
                                      GE_Mpolyhedron* result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_prototype() ) ;
   PEL_ASSERT( result->reference_polyhedron() == reference_polyhedron() ) ;
   PEL_ASSERT( result->nb_vertices() == nb_vertices() ) ;
   PEL_ASSERT( result->nb_faces() == nb_faces() ) ;
   PEL_ASSERT( result->name() == name() ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=0 ; i<result->nb_vertices() ; i++ ),
         result->vertex(i) == vertices->at(i) ) ) ;
   PEL_ASSERT(
      result->nb_space_dimensions() ==
      static_cast<GE_Point const*>( vertices->at( 0 ) )->nb_coordinates() ) ;
   PEL_ASSERT( IMPLIES( check_consistency(), result->is_consistent() ) ) ;
   PEL_ASSERT( vertices->owner() == result ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Mpolyhedron:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;

   if( !is_updating() )
   {
      if( !is_prototype() )
      {
         PEL_ASSERT( VERTS!=0 ) ;
         PEL_ASSERT( FORALL( ( size_t i=0;i<VERTS->count(); i++ ),
 		       dynamic_cast<GE_Point const*>(VERTS->at(i))!=0 ) ) ;

         PEL_ASSERT( dimension()<=nb_space_dimensions() ) ;

         PEL_ASSERT( IMPLIES( CENTER!=0,
                              CENTER->nb_coordinates()==
                                                    nb_space_dimensions() &&
                              CENTER->owner()==this ) ) ;
      }
      else
      {
         PEL_ASSERT( VERTS==0 ) ;
         PEL_ASSERT( CENTER==0 ) ;
      }
   }
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_Mpolyhedron:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
             PEL_ObjectRegister::create( PEL_Root::object(),
                                         "GE_Mpolyhedron descendant" ) ;
   return( result ) ;
}

//internal-----------------------------------------------------------------
void
GE_Mpolyhedron_ERROR:: n0( GE_Mpolyhedron const* poly_1,
                           GE_Mpolyhedron const* poly_2 )
//internal-----------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "GE_Mpolyhedron : " << endl ;
   mesg << "   the following two polyhedra do not have" << endl ;
   mesg << "   coincident vertices" << endl ;
   poly_1->print( mesg, 3 ) ;
   poly_2->print( mesg, 3 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

