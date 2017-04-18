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

#include <GE_SimplePolygon2D.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_PointIterator.hh>
#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_Polygon2D.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SegmentSegment_INT.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_Vector.hh>

#include <size_t_vector.hh>
#include <size_t_array2D.hh>

#include <iostream>
#include <sstream>



//------------------------------------------------------------------------------
double const AREA_NOT_COMPUTED = -PEL::max_double() ;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
GE_SimplePolygon2D:: update( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: update" ) ;
   SIGNED_AREA = AREA_NOT_COMPUTED ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D*
GE_SimplePolygon2D:: create_clone( PEL_Object* a_owner ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: create_clone" ) ;
   GE_SimplePolygon2D* result = new GE_SimplePolygon2D( a_owner, this ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==dimension() ) ;
   PEL_CHECK_POST( result->size()==size() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0; i<result->size(); i++ ),
              result->vertex(i)==vertex(i) ) ) ;
   PEL_CHECK_POST( result->is_set_point_point_intersector()==
                                           is_set_point_point_intersector() ) ;
   PEL_CHECK_POST( IMPLIES( is_set_point_point_intersector(),
                            result->point_point_intersector()==
                                           point_point_intersector() ) ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D*
GE_SimplePolygon2D:: create( PEL_Object* a_owner,
                             PEL_Sequence const* a_vertex_table ) 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: create" ) ;
   PEL_CHECK_PRE( a_vertex_table!=0 ) ;
   PEL_CHECK_PRE( a_vertex_table->index_limit()==a_vertex_table->count() ) ;
   PEL_CHECK_PRE( a_vertex_table->count()>2 ) ;
   PEL_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_vertex_table->index_limit() ; ++i ),
         dynamic_cast<GE_Point const*>( a_vertex_table->at(i) )!=0 &&
         static_cast<GE_Point const*>( a_vertex_table->at(i) )->nb_coordinates()==2 ) ) ;
   GE_SimplePolygon2D* result = new GE_SimplePolygon2D( a_owner,
                                                        a_vertex_table ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==2 ) ;
   PEL_CHECK_POST( result->size()==a_vertex_table->count() ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0; i<result->size(); i++ ),
         result->vertex(i)==static_cast<GE_Point const*>( a_vertex_table->at(i) ) ) ) ;
   PEL_CHECK_POST( !result->is_set_point_point_intersector() ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D*
GE_SimplePolygon2D:: create( PEL_Object* a_owner,
                             doubleVector const& coordinate_table ) 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: create" ) ;
   PEL_ASSERT( coordinate_table.size()%2==0 ) ;
   PEL_CHECK_PRE( coordinate_table.size()>4 ) ;
   GE_SimplePolygon2D* result = new GE_SimplePolygon2D( a_owner,
                                                        coordinate_table ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==2 ) ;
   PEL_CHECK_POST( result->size()==coordinate_table.size()/2 ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0; i<result->size() ; i++ ),
         result->vertex(i)->owner()==result &&
         FORALL(
            ( size_t j=0; j<2 ; j++ ),
            result->vertex(i)->coordinate(j)==coordinate_table(2*i+j) ) ) ) ;
   PEL_CHECK_POST( !result->is_set_point_point_intersector() ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D:: ~GE_SimplePolygon2D( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: ~GE_SimplePolygon2D" ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: is_consistent( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: is_consistent" ) ;
   return( GE_Polygon::is_consistent() &&
           IMPLIES( is_set_segment_segment_intersector(), is_simple() ) ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: is_simple( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: is_simple" ) ;
   PEL_CHECK( is_set_segment_segment_intersector() ) ;
   
   bool result = true ;
   if( size()>=4 )
   {
      GE_SegmentSegment_INT* seg_seg_algo = segment_segment_intersector() ;
      for( size_t i=0 ; result && i<size()-1 ; ++i )
      {
         GE_Point const* P1 = vertex(i) ;
         GE_Point const* P2 = vertex(i+1) ;
         size_t j = i+2 ;
         while( result && j<size()-1 )
         {
            GE_Point const* Q1 = vertex(j) ;
            GE_Point const* Q2 = next_vertex(j) ;
            result = !seg_seg_algo->has_intersection( P1, P2, Q1, Q2 ) ;
            j++ ;
         }
      }
   }
   return( result ) ;
}



//------------------------------------------------------------------------------
int
GE_SimplePolygon2D:: orientation( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: orientation" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( SIGNED_AREA==AREA_NOT_COMPUTED )
   {
      compute_area() ;
   }
   int const result = SIGNED_AREA<0 ? -1 : 1 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result==1 || result==-1 ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
double
GE_SimplePolygon2D:: area( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: area" ) ;
   PEL_CHECK_PRE( area_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( SIGNED_AREA==AREA_NOT_COMPUTED )
   {
      compute_area() ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( PEL::abs( SIGNED_AREA ) ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: has_in_interior( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: has_in_interior" ) ;
   PEL_CHECK_PRE( has_in_interior_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = has_on_boundary( pt ) ;

   if( !result )
   {
      double x = pt->coordinate( 0 ), z = pt->coordinate( 1 ) ;
      GE_Point const* P0 = vertex( size()-1 ) ;
      double x0 = P0->coordinate( 0 ) ;
      double z0 = P0->coordinate( 1 ) ;
      bool zFlag0 = ( z0>=z ) ;
      GE_PointIterator* pt_it = create_vertex_iterator( 0 ) ;
      for( pt_it->start() ; pt_it->is_valid() ; pt_it->go_next() )
      {
         GE_Point const* P1 = pt_it->item() ;
         double x1 = P1->coordinate( 0 ) ;
         double z1 = P1->coordinate( 1 ) ;
         bool zFlag1 = ( z1>=z ) ;
         if( zFlag0!=zFlag1 )
         {
            if( ( (z1-z)*(x0-x1)>=(x1-x)*(z0-z1) ) == zFlag1 )
	    { 
               result = !result ;
	    }
         }
         zFlag0 = zFlag1 ;
         x0 = x1 ;
         z0 = z1 ;
      }
      pt_it->destroy() ; pt_it = 0 ;
   }
  
   PEL_CHECK_INV( invariant() ) ;

   return( result ) ;

// Ray crossing algorithm
}



//------------------------------------------------------------------------------
void
GE_SimplePolygon2D:: suppress_flat_angles( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: suppress_flat_angles" ) ;
   PEL_CHECK_PRE( is_set_point_segment_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_PointSegment_INT* inter = point_segment_intersector() ;
   size_t i=0 ;
   while( i<size() )
   {
      GE_Point const* prec = previous_vertex( i ) ;
      GE_Point const* curr = vertex( i ) ;
      GE_Point const* next = next_vertex( i ) ;
      if( inter->point_in_segment( curr, prec, next ) )
      {
         remove_vertex_at( i ) ;
      }
      else
      {
         i++ ;
      }
   }
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
void
GE_SimplePolygon2D:: suppress_flat_angles( double sine_of_flat_angle )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: suppress_flat_angles(sinus)" ) ;
   PEL_CHECK_INV( invariant() ) ;

   static GE_Vector* ps = GE_Vector::create( PEL_Root::object(), 2 ) ;
   static GE_Vector* cs = GE_Vector::create( PEL_Root::object(), 2 ) ;
   PEL_CHECK( ps!=0 && ps->nb_components()==2 ) ;
   PEL_CHECK( cs!=0 && cs->nb_components()==2 ) ;
   double const val_min = sine_of_flat_angle*sine_of_flat_angle ;
   size_t i=0 ;
   while( i<size() )
   {
      ps->re_initialize( previous_vertex( i ), vertex( i ) ) ;
      cs->re_initialize( vertex( i ), next_vertex( i ) ) ;
      double c = ps->cosine( cs ) ;      
      if( PEL::abs( 1.-c*c ) < val_min )
      {
         remove_vertex_at( i ) ;
      }
      else
      {
         i++ ;
      }
   }
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
PEL_Vector*
GE_SimplePolygon2D:: create_triangulation( PEL_Object* seq_owner ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: create_triangulation" ) ;
   PEL_CHECK_PRE( size()>=3 ) ;
   PEL_CHECK_PRE( is_set_segment_segment_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Vector* result = PEL_Vector::create( seq_owner, 0 ) ;

   PEL_List* polygon = PEL_List::create(0) ;
   for( size_t i=0 ; i<size() ; ++i )
   {
      polygon->append( const_cast<GE_Point*>( vertex(i) ) ) ;
   }
   
   while( polygon->count()>=3 )
   {
      if( polygon->count()>3 )
      {
         set_ear_first_vertex( polygon ) ;
      }
      PEL_Vector* triangle_pts = PEL_Vector::create( result, 3 ) ;
      triangle_pts->set_at(
         0, static_cast<GE_Point*>( polygon->at( polygon->count()-1 ) ) ) ;
      triangle_pts->set_at(
         1, static_cast<GE_Point*>( polygon->at(0) ) ) ;
      triangle_pts->set_at(
         2, static_cast<GE_Point*>( polygon->at(1) ) ) ;
      polygon->remove_at(0) ;
      result->append( triangle_pts ) ;
   }
   polygon->destroy() ; polygon = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==seq_owner ) ;
   PEL_CHECK_POST( result->index_limit()==size()-2 ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<size()-2 ; ++i ),
         dynamic_cast<PEL_Vector*>( result->at(i) )!=0 &&
         static_cast<PEL_Vector*>( result->at(i) )->index_limit()==3 &&
         static_cast<PEL_Vector*>( result->at(i) )->owner()==result &&
         FORALL(
            ( size_t j=0 ; j<3 ; j++ ),
            dynamic_cast<GE_Point*>( static_cast<PEL_Vector*>( result->at(i) )->at(j) )!=0 ) ) ) ;

   return( result ) ;

// Use the so-called "two-ears theorem" : except for triangles every 
// simple polygon has at least two non-overlapping ears
// and the resulting polygon is also simple. In that way each simple
// polygon has a triangulation made of size()-2 triangles. This algorithm in 
// o(n*n) in general case but in o(n) for weakly convex polygons.
}



//------------------------------------------------------------------------------
void
GE_SimplePolygon2D:: compute_area( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: compute_area" ) ;
   PEL_CHECK( SIGNED_AREA==AREA_NOT_COMPUTED ) ;

   GE_Point const* pt = vertex( size() - 1 ) ;
   double x_prec = pt->coordinate( 0 ), z_prec = pt->coordinate( 1 ) ;
   double x_curr, z_curr ;

   SIGNED_AREA  = 0.  ;
   GE_PointIterator* it = create_vertex_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      pt = it->item() ;
      x_curr = pt->coordinate( 0 ) ;
      z_curr = pt->coordinate( 1 ) ;
      SIGNED_AREA += z_curr*x_prec - x_curr*z_prec ;
      x_prec = x_curr ;
      z_prec = z_curr ;
   }
   it->destroy() ; it = 0 ;
   SIGNED_AREA *= 0.5 ;

   PEL_CHECK( SIGNED_AREA!=AREA_NOT_COMPUTED ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: invariant( void ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Polygon::invariant() ) ;
   PEL_ASSERT( dimension()==2 ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: area_PRE( void ) const
//------------------------------------------------------------------------------
{
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: has_in_interior_PRE( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( is_set_point_segment_intersector() ) ;
   PEL_ASSERT( pt!=0 && pt->nb_coordinates()==2 ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
void
GE_SimplePolygon2D:: set_ear_first_vertex( PEL_List* polygon ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: set_ear_first_vertex" ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK( polygon->count()==polygon->index_limit() ) ;
   PEL_CHECK( polygon->count()>=3 ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<polygon->count() ; ++i ),
                      dynamic_cast<GE_Point const*>( polygon->at(i) )!=0 ) ) ;

   size_t counter = 0 ;
   while( !is_first_vertex_an_ear_vertex( polygon ) )
   {
      PEL_Object* first_pt = polygon->at(0) ;
      polygon->remove_at(0) ;
      polygon->append( first_pt ) ;
      if( counter++>polygon->count() )
      {
         std::ostringstream message ;
         message << "Unable to triangulate the simple polygon : "
                 << std::endl ;
         print( message, 5 ) ;
         message << "No ear found for the triangulation resulting simple polygon : "
                 << std::endl ;
         GE_SimplePolygon2D* p = GE_SimplePolygon2D::create( 0, polygon ) ;
         p->print( message, 5 ) ;
         p->destroy() ; p = 0 ;
         PEL_Error::object()->raise_plain( message.str() ) ;
      }
   }
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: is_first_vertex_an_ear_vertex( PEL_List const* polygon ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: is_first_vertex_an_ear_vertex" ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK( polygon->count()==polygon->index_limit() ) ;
   PEL_CHECK( polygon->count()>=3 ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<polygon->count() ; ++i ),
                      dynamic_cast<GE_Point const*>( polygon->at(i) )!=0 ) ) ;

   bool result = is_first_vertex_a_convex_vertex( polygon ) ;
   size_t const nb_points = polygon->count() ;
   
   if( result && nb_points>3 )
   {
      GE_Point const* P =
         dynamic_cast<GE_Point const*>( polygon->at(0) ) ;
      GE_Point const* P1 =
         dynamic_cast<GE_Point const*>( polygon->at(1) ) ;
      GE_Point const* P2 =
         dynamic_cast<GE_Point const*>( polygon->at(2) ) ;
      GE_Point const* M1 =
         dynamic_cast<GE_Point const*>( polygon->at(nb_points-1) ) ;
      GE_Point const* M2 =
         dynamic_cast<GE_Point const*>( polygon->at(nb_points-2) ) ;
      GE_SegmentSegment_INT* seg_seg_algo = segment_segment_intersector() ;

      // Avoid creation of concave angles :
      if( result )
      {
         double det0 = -(P->coordinate(1)-M1->coordinate(1))*
                                 (P1->coordinate(0)-M1->coordinate(0))
                       +(P->coordinate(0)-M1->coordinate(0))*
                                 (P1->coordinate(1)-M1->coordinate(1)) ;
         double det1 = -(P->coordinate(1)-M1->coordinate(1))*
                                 (M2->coordinate(0)-M1->coordinate(0))
                       +(P->coordinate(0)-M1->coordinate(0))*
                                 (M2->coordinate(1)-M1->coordinate(1)) ;
         double det2 = -(P1->coordinate(1)-M1->coordinate(1))*
                                 (M2->coordinate(0)-M1->coordinate(0))
                       +(P1->coordinate(0)-M1->coordinate(0))*
                                 (M2->coordinate(1)-M1->coordinate(1)) ;
         result = !( det0*det1>.0 && det0*det2<0. ) ;
      }
      if( result )
      {
         double det0 = -(P->coordinate(1)-P1->coordinate(1))*
                                 (M1->coordinate(0)-P1->coordinate(0))
                       +(P->coordinate(0)-P1->coordinate(0))*
                                 (M1->coordinate(1)-P1->coordinate(1)) ;
         double det1 = -(P->coordinate(1)-P1->coordinate(1))*
                                 (P2->coordinate(0)-P1->coordinate(0))
                       +(P->coordinate(0)-P1->coordinate(0))*
                                 (P2->coordinate(1)-P1->coordinate(1)) ;
         double det2 = -(M1->coordinate(1)-P1->coordinate(1))*
                                 (P2->coordinate(0)-P1->coordinate(0))
                       +(M1->coordinate(0)-P1->coordinate(0))*
                                 (P2->coordinate(1)-P1->coordinate(1)) ;
         result = !( det0*det1>0. && det0*det2<0. ) ;
      }    

      // Find for edges cutting the new triangle :
      for( size_t j=2 ; result && j<polygon->count()-3 ; ++j )
      {
         GE_Point const* Q1 = dynamic_cast<GE_Point const*>( polygon->at(j) ) ;
         GE_Point const* Q2 = dynamic_cast<GE_Point const*>( polygon->at(j+1) ) ;
         result = !seg_seg_algo->has_intersection( P1, M1, Q1, Q2 ) ;
      }
   }
   return( result ) ;
   
   // A vertex is an ear vertex if it is a convex vertex and if no size of
   // the polygon cut the triangle built with this vertex and its two neighbours
}



//------------------------------------------------------------------------------
bool
GE_SimplePolygon2D:: is_first_vertex_a_convex_vertex(  PEL_List const* polygon ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SimplePolygon2D:: is_first_vertex_a_convex_vertex" ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK( polygon->count()==polygon->index_limit() ) ;
   PEL_CHECK( polygon->count()>=3 ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<polygon->count() ; ++i ),
                      dynamic_cast<GE_Point const*>( polygon->at(i) )!=0 ) ) ;
   
   GE_Point const* cur_vert =
                       static_cast<GE_Point const*>( polygon->at( 0 ) ) ;
   GE_Point const* pre_vert =
      static_cast<GE_Point const*>( polygon->at( polygon->count()-1 ) ) ;
   GE_Point const* nex_vert =
                       static_cast<GE_Point const*>( polygon->at( 1 ) ) ;
   double const x = cur_vert->coordinate( 0 ) ; 
   double const z = cur_vert->coordinate( 1 ) ;
   double const ori = (double) orientation() ;
   bool result = ( ( ( pre_vert->coordinate( 1 )-z )*
                     ( nex_vert->coordinate( 0 )-x )-
                     ( pre_vert->coordinate( 0 )-x )*
                     ( nex_vert->coordinate( 1 )-z ) ) * ori ) >0. ;

   return( result ) ;
   
   // A vertex is convex if a right ( left ) turn is made while going from
   // its previous vertex to its next vertex if the polygon is clockwise 
   // ( counter-clockwise ) ordered.
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D:: GE_SimplePolygon2D( PEL_Object* a_owner,
                                         PEL_Sequence const* vertex_table )
//------------------------------------------------------------------------------
      : GE_Polygon( a_owner, 2, vertex_table )
{
   PEL_LABEL( "GE_SimplePolygon2D:: GE_SimplePolygon2D" ) ;

   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D:: GE_SimplePolygon2D( PEL_Object* a_owner,
                                         doubleVector const& coordinate_table )
//------------------------------------------------------------------------------
   : GE_Polygon( a_owner, 2, coordinate_table )
{
   PEL_LABEL( "GE_SimplePolygon2D:: GE_SimplePolygon2D" ) ;

   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_SimplePolygon2D:: GE_SimplePolygon2D( PEL_Object* a_owner,
                                         GE_SimplePolygon2D const* other )
//------------------------------------------------------------------------------
   : GE_Polygon( a_owner, other )
{
   PEL_LABEL( "GE_SimplePolygon2D:: GE_SimplePolygon2D" ) ;

   update() ;

   PEL_CHECK_INV( invariant() ) ;
}
