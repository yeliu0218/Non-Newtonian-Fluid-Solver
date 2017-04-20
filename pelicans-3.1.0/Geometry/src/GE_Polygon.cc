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

#include <GE_Polygon.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_PointIterator.hh>
#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_SegmentSegment_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_List.hh>
#include <PEL_Sequence.hh>

#include <doubleVector.hh>

#include <iostream>
#include <string>


//------------------------------------------------------------------------------
size_t
GE_Polygon:: size( void ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: size" ) ;
  PEL_CHECK_INV( invariant() ) ;

  return( LIST_OF_POINTS->count() ) ;
}



//------------------------------------------------------------------------------
size_t
GE_Polygon:: dimension( void ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: dimension" ) ;
  PEL_CHECK_INV( invariant() ) ;
  return( DIM ) ;
}



//------------------------------------------------------------------------------
GE_PointIterator*
GE_Polygon:: create_vertex_iterator( PEL_Object* owner ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: create_vertex_iterator" ) ;
  PEL_CHECK_INV( invariant() ) ;
  
  GE_PointIterator* result = GE_PointIterator::create( owner, LIST_OF_POINTS ) ;

  PEL_CHECK_INV( invariant() ) ;
  PEL_CHECK_POST( result != 0 ) ;
  PEL_CHECK_POST( result->owner()==owner ) ;
  PEL_CHECK_POST( EQUIVALENT( size()!=0 , result->is_valid() ) ) ;

  return( result ) ;
}



//------------------------------------------------------------------------------
GE_Point const*
GE_Polygon:: vertex( size_t i ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: vertex" ) ;
  PEL_CHECK_PRE( i < size() ) ; 
  PEL_CHECK_INV( invariant() ) ;
  
  GE_Point const* result =
                      static_cast<GE_Point const*>( LIST_OF_POINTS->at( i ) ) ;

  PEL_CHECK_INV( invariant() ) ;
  PEL_CHECK_POST( result != 0 ) ;
  PEL_CHECK_POST( result->nb_coordinates()==dimension() ) ;

  return( result ) ;
}



//------------------------------------------------------------------------------
GE_Point const*
GE_Polygon:: next_vertex( size_t i ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: next_vertex" ) ;
  PEL_CHECK_PRE( i < size() ) ;
  PEL_CHECK_INV( invariant() ) ;
  
  GE_Point const* result = static_cast<GE_Point const*>(
                                LIST_OF_POINTS->at( next_vertex_index( i ) ) ) ;

  PEL_CHECK_INV( invariant() ) ;
  PEL_CHECK_POST( result != 0 ) ;
  PEL_CHECK_POST( result == vertex( next_vertex_index( i ) ) ) ;

  return( result ) ;
}



//------------------------------------------------------------------------------
GE_Point const*
GE_Polygon:: previous_vertex( size_t i ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: previous_vertex" ) ;
  PEL_CHECK_PRE( i < size() ) ; 
  PEL_CHECK_INV( invariant() ) ;
  
  GE_Point const* result = static_cast<GE_Point const*>(
                              LIST_OF_POINTS->at( previous_vertex_index(i) ) ) ;

  PEL_CHECK_INV( invariant() ) ;
  PEL_CHECK_POST( result != 0 ) ;
  PEL_CHECK_POST( result == vertex( previous_vertex_index( i ) ) ) ;

  return( result ) ;
}



//------------------------------------------------------------------------------
size_t
GE_Polygon:: next_vertex_index( size_t i ) const
//------------------------------------------------------------------------------
{
  PEL_LABEL( "GE_Polygon:: next_vertex_index" ) ;
  PEL_CHECK_PRE( i < size() ) ; 
  PEL_CHECK_INV( invariant() ) ;

  size_t result = i + 1 ;
  if( result == size() ) result = 0 ;

  PEL_CHECK_INV( invariant() ) ;
  PEL_CHECK_POST( result < size() ) ;
  PEL_CHECK_POST( vertex( result ) == next_vertex( i ) ) ;

  return( result ) ; 
}



//------------------------------------------------------------------------------
size_t
GE_Polygon:: previous_vertex_index( size_t i ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: previous_vertex_index" ) ;
   PEL_CHECK_PRE( i < size() ) ; 
   PEL_CHECK_INV( invariant() ) ;

   size_t result = size()-1 ;
   if( i > 0 ) result = i-1 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result < size() ) ;
   PEL_CHECK_POST( vertex( result ) == previous_vertex( i ) ) ;

   return( result ) ; 
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: is_equal( PEL_Object const* other ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: is_equal" ) ;
   PEL_CHECK_PRE( other!=0 ) ; 
   PEL_CHECK_PRE( dynamic_cast<GE_Polygon const*>( other )!=0 ) ; 
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Polygon const* polygon = static_cast<GE_Polygon const*>( other ) ;

   bool result = ( size() == polygon->size() ) ;
   result = result && ( dimension()==polygon->dimension() ) ;

   if( result )
   {

      // First search for a vertex of other result to first vertex of self.
      result = false ;
      size_t iS = 0 ;
      for( size_t i=0; !result && i<size() ; i++ && iS++ )
      {
         if( PT_PT_INTERSECTOR==0 )
         {
            result = polygon->vertex( i )->is_equal( vertex( 0 ) ) ;
         }
         else
         {
            result = PT_PT_INTERSECTOR->points_are_close(
                                                     polygon->vertex( i ),
                                                     vertex( 0 ) ) ;
         }
      }
     
      // Then test the points sequence.
      if( result )
      {
         for( size_t j=1; result && j<size(); j++ )
         {
            iS = polygon->next_vertex_index( iS ) ;
            if( PT_PT_INTERSECTOR==0 )
            {
               result = vertex( j )->is_equal( polygon->vertex( iS ) ) ;
            }
            else
            {
               result = PT_PT_INTERSECTOR->points_are_close(
                                       vertex( j ),polygon->vertex( iS ) ) ;
               
            }
         }
      }
   }

   PEL_CHECK_INV( invariant() ) ;

   return( result ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: insert_vertex_at( size_t i, GE_Point const* pt )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: insert_vertex_at" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( i <= size() ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0 && i!=size(),
               !point_point_intersector()->points_are_close(
                                               pt, previous_vertex( i ) ) ) ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0 && i!=size(),
               !point_point_intersector()->points_are_close(
                                               pt, vertex( i ) ) ) ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0 && i==size(),
               !point_point_intersector()->points_are_close(
                                               pt, previous_vertex( size()-1 ) ) ) ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0 && i==size(),
               !point_point_intersector()->points_are_close( pt, vertex( 0 ) ) ) ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;
   
   if( i==size() )
   {
      LIST_OF_POINTS->append( const_cast<GE_Point*>( pt ) ) ;
   }
   else
   {
      LIST_OF_POINTS->insert_at( i, const_cast<GE_Point*>( pt ) ) ;
   }
   update() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == OLD( size ) + 1 ) ;
   PEL_CHECK_POST( vertex(i)==pt ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: append_vertex( GE_Point const* pt )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: append_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE( IMPLIES( is_set_point_point_intersector() && size()>0,
                           !point_point_intersector()->points_are_close(
                                                    pt, vertex( size()-1 ) ) ) ) ;
   PEL_CHECK_PRE( IMPLIES( is_set_point_point_intersector() && size()>0,
                           !point_point_intersector()->points_are_close(
                                                    pt, vertex( 0 ) ) ) ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;

   LIST_OF_POINTS->append( const_cast<GE_Point*>( pt ) ) ;
   update() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == OLD( size ) + 1 ) ;
   PEL_CHECK_POST( vertex( OLD( size ) )==pt ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: remove_vertex_at( size_t i )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: remove_vertex_at" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( i < size() ) ; 
   PEL_SAVEOLD( size_t, size, size() ) ;

   LIST_OF_POINTS->remove_at( i ) ;
   update() ;
   
   PEL_CHECK_POST( size() == OLD( size ) - 1 ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: replace_vertex( size_t i, GE_Point const* pt )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: replace_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( i < size() ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0,
               !point_point_intersector()->points_are_close(
                                               pt, previous_vertex( i ) ) ) ) ;
   PEL_CHECK_PRE(
      IMPLIES( is_set_point_point_intersector() && size()>0,
               !point_point_intersector()->points_are_close(
                                               pt, next_vertex( i ) ) ) ) ;
   PEL_SAVEOLD( size_t, size, size() ) ;

   LIST_OF_POINTS->set_at( i, const_cast<GE_Point*>( pt ) ) ;
   update() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( size() == OLD( size ) ) ;
   PEL_CHECK_POST( vertex(i)==pt ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: remove_vertices( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: remove_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;

   while( size()!=0 )
   {
      LIST_OF_POINTS->remove_at( 0 ) ;
   }
   update() ;
   
   PEL_CHECK_POST( size() == 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: has_on_boundary( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: has_on_boundary" ) ;
   
   PEL_CHECK_PRE( is_set_point_segment_intersector() ) ;
   PEL_CHECK_PRE( pt!=0 && pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_PointSegment_INT* pt_seg_algo = point_segment_intersector() ;
   
   bool result = false ;
   GE_PointIterator* pt_it = create_vertex_iterator( 0 ) ;
   GE_Point const* P1 = vertex( size()-1 ) ;
   for( pt_it->start() ; pt_it->is_valid() && !result ; pt_it->go_next() )
   {
      GE_Point const* P2 = pt_it->item() ;
      result = pt_seg_algo->point_in_segment( pt, P1, P2 ) ;
      P1 = P2 ;
   }
   pt_it->destroy() ; pt_it = 0 ;
   PEL_CHECK_INV( invariant() ) ;

   return( result ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: project_on_boundary( GE_Point* pt ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: project_on_boundary" ) ;
   PEL_CHECK_PRE( project_on_boundary_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Point* pt_sol = GE_Point::create( 0, dimension() ) ;
   GE_Point* pt_pro = GE_Point::create( 0, dimension() ) ;
   GE_PointIterator* pt_it = create_vertex_iterator( 0 ) ;

   double dist_min = PEL::max_double() ;
   GE_Point const* P1 = vertex( size()-1 ) ;
   for( pt_it->start() ; pt_it->is_valid() ; pt_it->go_next() )
   {
      GE_Point const* P2 = pt_it->item() ;
      double prod_scal = 0. ;
      double v0v1 = 0. ;
      for( size_t j=0; j<dimension(); ++j )
      {
         double const v0 = P1->coordinate(j)  ;
         double const  x = P2->coordinate(j)-v0 ;
         v0v1 += x*x ;
         prod_scal += ( pt->coordinate(j)-v0 )*x ;
      }
      double const alpha = PEL::max( 0., PEL::min( 1., prod_scal/v0v1 ) ) ;
      pt_pro->set_as_barycenter( alpha, P1, P2 ) ;
      double dist = pt->distance( pt_pro ) ;
      if( dist<dist_min )
      {
         dist_min = dist ;
         pt_sol->copy( pt_pro ) ;
      }
      P1 = P2 ;
   }
   pt->copy( pt_sol ) ;
   pt_sol->destroy() ; pt_sol = 0 ;
   pt_pro->destroy() ; pt_pro = 0 ;
   pt_it->destroy() ; pt_it = 0 ;
   
   PEL_CHECK_POST( project_on_boundary_POST( pt ) ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: set_point_point_intersector(
                                        GE_PointPoint_INT* a_pt_pt_intersector )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: set_point_point_intersector" ) ;
   PEL_CHECK_PRE( a_pt_pt_intersector!=0 ) ;
   PEL_CHECK_PRE( !is_set_point_point_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PT_PT_INTERSECTOR = a_pt_pt_intersector ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_set_point_point_intersector() ) ;
   PEL_CHECK_POST( point_point_intersector()==a_pt_pt_intersector ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: is_set_point_point_intersector( void ) const
//------------------------------------------------------------------------------
{
   return( PT_PT_INTERSECTOR!=0 ) ;
}



//------------------------------------------------------------------------------
GE_PointPoint_INT*
GE_Polygon:: point_point_intersector( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: point_point_intersector" ) ;
   PEL_CHECK_PRE( is_set_point_point_intersector() ) ;
   GE_PointPoint_INT* result = PT_PT_INTERSECTOR ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: set_point_segment_intersector(
                                     GE_PointSegment_INT* a_pt_seg_intersector )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: set_point_segment_intersector" ) ;
   PEL_CHECK_PRE( a_pt_seg_intersector!=0 ) ;
   PEL_CHECK_PRE( !is_set_point_segment_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PT_SEG_INTERSECTOR = a_pt_seg_intersector ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_set_point_segment_intersector() ) ;
   PEL_CHECK_POST( point_segment_intersector()==a_pt_seg_intersector ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: is_set_point_segment_intersector( void ) const
//------------------------------------------------------------------------------
{
   return( PT_SEG_INTERSECTOR!=0 ) ;
}



//------------------------------------------------------------------------------
GE_PointSegment_INT*
GE_Polygon:: point_segment_intersector( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: point_segment_intersector" ) ;
   PEL_CHECK_PRE( is_set_point_segment_intersector() ) ;
   GE_PointSegment_INT* result = PT_SEG_INTERSECTOR ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: set_segment_segment_intersector(
                                  GE_SegmentSegment_INT* a_seg_seg_intersector )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: set_segment_segment_intersector" ) ;
   PEL_CHECK_PRE( a_seg_seg_intersector!=0 ) ;
   PEL_CHECK_PRE( !is_set_segment_segment_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;
   SEG_SEG_INTERSECTOR = a_seg_seg_intersector ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_set_segment_segment_intersector() ) ;
   PEL_CHECK_POST( segment_segment_intersector()==a_seg_seg_intersector ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: is_set_segment_segment_intersector( void ) const
//------------------------------------------------------------------------------
{
   return( SEG_SEG_INTERSECTOR!=0 ) ;
}



//------------------------------------------------------------------------------
GE_SegmentSegment_INT*
GE_Polygon:: segment_segment_intersector( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: segment_segment_intersector" ) ;
   PEL_CHECK_PRE( is_set_segment_segment_intersector() ) ;
   GE_SegmentSegment_INT* result = SEG_SEG_INTERSECTOR ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon:: print( std::ostream& os, size_t indent_width ) const 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string space( indent_width, ' ' ) ;
   os << space << "Polygon :" << std::endl ;
   for( size_t i=0 ; i<size() ; ++i )
   {
      vertex(i)->print( os, indent_width+3 ) ;
      os << std::endl ;
   }
}



//------------------------------------------------------------------------------
GE_Polygon:: GE_Polygon( PEL_Object* a_owner, size_t dim )
//------------------------------------------------------------------------------
  : PEL_Object( a_owner ),
    LIST_OF_POINTS( PEL_List::create( this ) ),
    DIM( dim ),
    PT_PT_INTERSECTOR( 0 ),
    PT_SEG_INTERSECTOR( 0 ),
    SEG_SEG_INTERSECTOR( 0 )
{
   PEL_LABEL( "GE_Polygon:: GE_Polygon" ) ;
   
   PEL_CHECK( dim==2 || dim==3 ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}




//------------------------------------------------------------------------------
GE_Polygon:: GE_Polygon( PEL_Object* a_owner,
                         size_t dim,
                         PEL_Sequence const* vertex_table )
//------------------------------------------------------------------------------
  : PEL_Object( a_owner ),
    LIST_OF_POINTS( PEL_List::create( this ) ),
    DIM( dim ),
    PT_PT_INTERSECTOR( 0 ),
    PT_SEG_INTERSECTOR( 0 ),
    SEG_SEG_INTERSECTOR( 0 )
{
   PEL_LABEL( "GE_Polygon:: GE_Polygon" ) ;

   PEL_CHECK( dim==2 || dim==3 ) ;
   PEL_CHECK( vertex_table!=0 ) ;
   PEL_CHECK( vertex_table->index_limit()==vertex_table->count() ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<vertex_table->index_limit() ; ++i ),
         dynamic_cast<GE_Point const*>( vertex_table->at(i) )!=0 &&
         static_cast<GE_Point const*>( vertex_table->at(i) )->nb_coordinates()==dim ) ) ;
   
   for( size_t i=0; i<vertex_table->index_limit(); i++ )
   {
      GE_Point* pt = static_cast<GE_Point*>( vertex_table->at( i ) ) ;
      LIST_OF_POINTS->append( pt ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon:: GE_Polygon( PEL_Object* a_owner, GE_Polygon const* other )
//------------------------------------------------------------------------------
  : PEL_Object( a_owner ),
    LIST_OF_POINTS( PEL_List::create( this ) ),
    DIM( other->DIM ),
    PT_PT_INTERSECTOR( other->PT_PT_INTERSECTOR ),
    PT_SEG_INTERSECTOR( other->PT_SEG_INTERSECTOR ),
    SEG_SEG_INTERSECTOR( other->SEG_SEG_INTERSECTOR )
{
   PEL_LABEL( "GE_Polygon:: GE_Polygon" ) ;

   PEL_CHECK( other!=0 ) ;
   for( size_t i=0; i<other->size(); i++ )
   {
      LIST_OF_POINTS->append( const_cast<GE_Point*>( other->vertex( i ) ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon:: GE_Polygon( PEL_Object* a_owner,
                         size_t dim,
                         doubleVector const& coordinate_table )
//------------------------------------------------------------------------------
  : PEL_Object( a_owner ),
    LIST_OF_POINTS( PEL_List::create( this ) ),
    DIM( dim ),
    PT_PT_INTERSECTOR( 0 ),
    PT_SEG_INTERSECTOR( 0 ),
    SEG_SEG_INTERSECTOR( 0 )
{
   PEL_LABEL( "GE_Polygon:: GE_Polygon" ) ;

   PEL_CHECK( dim==2 || dim==3 ) ;
   PEL_CHECK( coordinate_table.size()%dim==0 ) ;
   
   for( size_t i=0; i<coordinate_table.size(); i=i+DIM )
   {
      GE_Point* pt = GE_Point::create( this, DIM ) ;
      for( size_t j=0; j<DIM; j++ )
      {
         pt->set_coordinate( j, coordinate_table(i+j) ) ;
      }
      LIST_OF_POINTS->append( pt ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon:: ~GE_Polygon( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: ~GE_Polygon" ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: is_consistent( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon:: is_consistent" ) ;
   bool result = true ;
   if( PT_PT_INTERSECTOR!=0 && size()>=2 )
   {
      for( size_t i=0 ; result && i<size() ; ++i )
      {
         result = !PT_PT_INTERSECTOR->points_are_close(
                                                   vertex(i), next_vertex(i) ) ;
      }
   }
   return( result ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: area_PRE( void ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( is_set_point_point_intersector() ) ;
   PEL_ASSERT( is_set_point_segment_intersector() ) ;
   PEL_ASSERT( is_set_segment_segment_intersector() ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: has_in_interior_PRE( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( is_set_point_point_intersector() ) ;
   PEL_ASSERT( is_set_point_segment_intersector() ) ;
   PEL_ASSERT( is_set_segment_segment_intersector() ) ;
   PEL_ASSERT( pt!=0 && pt->nb_coordinates()==dimension() ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: project_on_boundary_PRE( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( pt!=0 && pt->nb_coordinates()==dimension() ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: project_on_boundary_POST( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_set_point_segment_intersector(),
                        has_on_boundary( pt ) ) ) ;
   return( true ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon:: invariant( void ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( LIST_OF_POINTS!=0 ) ;
   PEL_ASSERT( LIST_OF_POINTS->index_limit()==LIST_OF_POINTS->count() ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=0 ; i<LIST_OF_POINTS->index_limit() ; ++i ),
         dynamic_cast<GE_Point const*>( LIST_OF_POINTS->at(i) )!=0 ) ) ;
   PEL_ASSERT( DIM==2 || DIM==3 ) ;
   PEL_ASSERT( is_consistent() ) ;
   return( true ) ;
}
