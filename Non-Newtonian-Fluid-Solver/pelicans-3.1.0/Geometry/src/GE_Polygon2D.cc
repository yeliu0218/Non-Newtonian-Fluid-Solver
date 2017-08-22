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

#include <GE_Polygon2D.hh>

#include <GE_Point.hh>
#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_SimplePolygon2D.hh>
#include <GE_SegmentSegment_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <size_t_vector.hh>

#include <iostream>
#include <sstream>

//------------------------------------------------------------------------------
double const AREA_NOT_COMPUTED = -PEL::max_double() ;
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
void
GE_Polygon2D:: update( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: update" ) ;
   AREA = AREA_NOT_COMPUTED ;
   if( SIMPLE_POLYGONS!=0 )
   {
      for( size_t i=0 ; i<SIMPLE_POLYGONS->index_limit() ; ++i )
      {
         destroy_possession( SIMPLE_POLYGONS->at( i ) ) ;
      }
      destroy_possession( SIMPLE_POLYGONS ) ;
      SIMPLE_POLYGONS = 0 ;
   }
}



//------------------------------------------------------------------------------
GE_Polygon2D*
GE_Polygon2D:: create( PEL_Object* a_owner ) 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: create" ) ;
   GE_Polygon2D* result = new GE_Polygon2D( a_owner ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==2 ) ;
   PEL_CHECK_POST( result->size()==0 ) ;
   PEL_CHECK_POST( !result->is_set_point_point_intersector() ) ;
   PEL_CHECK_POST( !result->is_set_point_segment_intersector() ) ;
   PEL_CHECK_POST( !result->is_set_segment_segment_intersector() ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D*
GE_Polygon2D:: create( PEL_Object* a_owner,
                       PEL_Sequence const* vertex_table ) 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: create" ) ;
   PEL_CHECK_PRE( vertex_table!=0 ) ;
   PEL_CHECK_PRE( vertex_table->index_limit()==vertex_table->count() ) ;
   PEL_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<vertex_table->index_limit() ; ++i ),
         dynamic_cast<GE_Point const*>( vertex_table->at(i) )!=0 &&
         static_cast<GE_Point const*>( vertex_table->at(i) )->nb_coordinates()==2 ) ) ;
   GE_Polygon2D* result = new GE_Polygon2D( a_owner, vertex_table ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==2 ) ;
   PEL_CHECK_POST( result->size()==vertex_table->count() ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0; i<result->size(); i++ ),
         result->vertex(i)==static_cast<GE_Point const*>( vertex_table->at(i) ) ) ) ;
   PEL_CHECK_POST( !result->is_set_point_point_intersector() ) ;
   PEL_CHECK_POST( !result->is_set_point_segment_intersector() ) ;
   PEL_CHECK_POST( !result->is_set_segment_segment_intersector() ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D*
GE_Polygon2D:: create( PEL_Object* a_owner,
                       doubleVector const& coordinate_table ) 
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: create" ) ;
   PEL_CHECK_PRE( coordinate_table.size()%2==0 ) ;
   GE_Polygon2D* result = new GE_Polygon2D( a_owner, coordinate_table ) ;
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
   PEL_CHECK_POST( !result->is_set_point_segment_intersector() ) ;
   PEL_CHECK_POST( !result->is_set_segment_segment_intersector() ) ;
   return( result ) ;
}




//------------------------------------------------------------------------------
GE_Polygon2D*
GE_Polygon2D:: create_clone( PEL_Object* a_owner ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: create_clone" ) ;
   GE_Polygon2D* result = new GE_Polygon2D( a_owner, this ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->dimension()==2 ) ;
   PEL_CHECK_POST( result->size()==size() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0; i<result->size(); i++ ),
                            result->vertex(i)==vertex(i) ) ) ;
   PEL_CHECK_POST( result->is_equal( this ) ) ;
   PEL_CHECK_POST( result->is_set_point_point_intersector()==
                                           is_set_point_point_intersector() ) ;
   PEL_CHECK_POST( IMPLIES( is_set_point_point_intersector(),
                            result->point_point_intersector()==
                                           point_point_intersector() ) ) ;
   PEL_CHECK_POST( result->is_set_point_segment_intersector()==
                                           is_set_point_segment_intersector() ) ;
   PEL_CHECK_POST( IMPLIES( is_set_point_segment_intersector(),
                            result->point_segment_intersector()==
                                           point_segment_intersector() ) ) ;
   PEL_CHECK_POST( result->is_set_segment_segment_intersector()==
                                           is_set_segment_segment_intersector() ) ;
   PEL_CHECK_POST( IMPLIES( is_set_segment_segment_intersector(),
                            result->segment_segment_intersector()==
                                           segment_segment_intersector() ) ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
double
GE_Polygon2D:: area( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: area" ) ;
   PEL_CHECK_PRE( area_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( AREA==AREA_NOT_COMPUTED )
   {
      compute_area() ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( AREA ) ;
}



//------------------------------------------------------------------------------
bool
GE_Polygon2D:: has_in_interior( GE_Point const* pt ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: has_in_interior" ) ;
   PEL_CHECK_PRE( has_in_interior_PRE( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   bool result = false ;
   if( size()>=3 )
   {
      PEL_Vector const* poly_table = split_of_simple_polygons() ;
      for( size_t i=0 ; !result && i<poly_table->index_limit() ; ++i )
      {
         GE_SimplePolygon2D const* poly =
                static_cast<GE_SimplePolygon2D const*>( poly_table->at(i) ) ;
         result = poly->has_in_interior( pt ) ;
      }
   }
   return( result ) ;
}



//------------------------------------------------------------------------------
PEL_Vector const*
GE_Polygon2D:: split_of_simple_polygons( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: split_of_simple_polygons" ) ;
   PEL_CHECK( size()>=3 ) ;
   PEL_CHECK( is_set_point_point_intersector() ) ;
   PEL_CHECK( is_set_point_segment_intersector() ) ;
   PEL_CHECK( is_set_segment_segment_intersector() ) ;
   
   if( SIMPLE_POLYGONS==0 )
   {
      SIMPLE_POLYGONS =
         create_split_of_simple_polygons( const_cast<GE_Polygon2D*>( this ) ) ;
   }
   PEL_Vector const* result = SIMPLE_POLYGONS ;
   if( result==0 )
   {
      std::ostringstream message ;
      message << std::endl ;
      message << "  Unable to split in simple polygon the polygon : "
              << std::endl << std::endl ;
      print( message, 6 ) ;
      PEL_Error::object()->raise_plain( message.str() ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==this ) ;
   PEL_CHECK( FORALL( ( size_t i=0; i<result->count(); i++ ),
                      result->at( i )->owner()==this ) ) ;
   PEL_CHECK( FORALL( ( size_t i=0; i<result->count(); i++ ),
                      dynamic_cast<GE_SimplePolygon2D const*>( result->at( i ) )!=0 ) ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_point_point_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->point_point_intersector()==point_point_intersector() ) ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_point_segment_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->point_segment_intersector()==point_segment_intersector() ) ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_segment_segment_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->segment_segment_intersector()==segment_segment_intersector() ) ) ;
   return( result ) ;   
}



//------------------------------------------------------------------------------
PEL_Vector*
GE_Polygon2D:: create_split_of_simple_polygons( PEL_Object* a_owner ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: create_split_of_simple_polygons" ) ;
   PEL_CHECK_PRE( size()>=3 ) ;
   PEL_CHECK_PRE( is_set_point_point_intersector() ) ;
   PEL_CHECK_PRE( is_set_point_segment_intersector() ) ;
   PEL_CHECK_PRE( is_set_segment_segment_intersector() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Polygon2D* p = create_clone( 0 ) ;

   GE_PointPoint_INT* pt_pt_algo = point_point_intersector() ;
   GE_PointSegment_INT* pt_seg_algo = point_segment_intersector() ;
   GE_SegmentSegment_INT* seg_seg_algo = segment_segment_intersector() ;

   static PEL_List* ptL = PEL_List::create( PEL_Root::object() ) ;
   static GE_Point* newPt = GE_Point::create( PEL_Root::object(),
                                              (size_t) 2 ) ;
   PEL_CHECK( ptL!=0 && ptL->count()==0 ) ;
   PEL_CHECK( newPt!=0 && newPt->nb_coordinates()==2 ) ;
   
   GE_Point const* curPt = 0 ;
   GE_Point const* nexPt = 0 ;
   GE_Point const* iniPt = 0 ;

   PEL_Vector* result = PEL_Vector::create( a_owner, 0 ) ;

   size_t test = 0 ;
   size_t i = 0, iStart = 0 ;
   size_t_vector inS( p->size() ) ;
   inS.set( 1 ) ;
   while( iStart < p->size() )
   {
      i = iStart ;
      curPt = p->vertex( iStart ) ;
      iniPt = curPt ;
      inS( iStart ) = 0 ;
      bool cutting_a_polygon = true ;
      while( cutting_a_polygon )
      {
         if( ++test>2*p->size() )
         {
            if( ptL->count()>0 )
            {
               ptL->clear() ;
            }
            p->destroy() ; p = 0 ;
            for( i=0 ; i<result->index_limit() ; i++ )
            {
               PEL_Object* pol = result->at(i) ;
               if( a_owner==0 )
               {
                  pol->destroy() ;
               }
               else
               {
                  a_owner->destroy_possession( pol ) ;
               }
            }    
            if( a_owner==0 )
            {
               result->destroy() ;
            }
            else
            {
               a_owner->destroy_possession( result ) ;
            }
            result = 0 ;
            return( result ) ;
         }
         ptL->append( const_cast<GE_Point*>( curPt ) ) ;

         inS( i ) = 0 ;
         size_t iprev = i ;
         i = p->next_vertex_index( iprev ) ;
         nexPt = p->vertex( i ) ;
         double alpha = 1.1 ;  

	 bool intersect = false ;
	 bool inlined = false ;
	 size_t iSide = p->size() ;
	 for( size_t kInt=0; kInt<p->size() ; kInt++ )
	 {
            if( kInt!=iprev && kInt!=i )
	    {
	       size_t jInt = p->next_vertex_index( kInt ) ;
               seg_seg_algo->compute_intersection( curPt, nexPt,
                                                   p->vertex( kInt ),
                                                   p->vertex( jInt ) ) ;
               GE_SegmentSegment_INT::IntersectionType inter_type =
                                     seg_seg_algo->intersection_type() ;
               if( inter_type==GE_SegmentSegment_INT::one_intersection )
               {
                  double const al = seg_seg_algo->alpha() ;
                  newPt->set_as_barycenter( al, curPt, nexPt ) ;
                  if(  al<alpha &&
                       !pt_pt_algo->points_are_close( newPt, curPt ) )
                  {
                     intersect = true ;
                     alpha = al ;
                     iSide = kInt ;
                     inlined = false ;
                  }
               }
               else if( inter_type==GE_SegmentSegment_INT::colinear )
               {
                  if( 1.<alpha &&
                      pt_seg_algo->point_in_segment( nexPt,
                                                     p->vertex( kInt ),
                                                     p->vertex( jInt ) ) )
                  {
                     intersect = true ;
                     alpha = 1. ;
                     iSide = kInt ;
                     inlined = false ;
                  }
               }
               else if( !intersect &&
                        inter_type==GE_SegmentSegment_INT::colinear_one_in_the_other  )
	       {
		  inlined = true ;
	       }
	    }
	 }
            
	 if( intersect )
	 {
            newPt->set_as_barycenter( alpha, curPt, nexPt ) ;
            if( pt_pt_algo->points_are_close( newPt, curPt ) ) // alpha<epsilon
            {
               // ok !
               inS( iprev ) = 0 ;
            }
            else if( pt_pt_algo->points_are_close( newPt, nexPt ) ) // alpha>1.-epsilon
            {
               curPt = nexPt ;
               inS( i ) = 0 ;
            }
            else if( pt_pt_algo->points_are_close( newPt, p->vertex(iSide) ) ) // beta<epsilon
            {
               curPt = p->vertex(iSide) ;
               inS( iSide ) = 0 ;
            }
            else if( pt_pt_algo->points_are_close( newPt, p->next_vertex(iSide) ) ) // beta>1-epsilon
            {
               curPt = p->next_vertex(iSide) ;
               inS( p->next_vertex_index( iSide ) ) = 0 ;
            }
            else
            {
               curPt = newPt->create_clone( result ) ;
               p->insert_vertex_at( i, curPt ) ;
	       // un insert_at à la main dans inS ...
	       size_t length = inS.size() ;
	       inS.resize( length+1 ) ;
	       for( size_t l=inS.size()-1; l>i; l-- )
	       {
		  inS( l ) = inS( l-1 ) ;
	       }
               inS( i ) = 0 ;
            }
            i = iSide ;
            if( !pt_pt_algo->points_are_close( newPt, nexPt ) &&
                !pt_pt_algo->points_are_close( newPt, p->vertex(iSide) ) )
            {
               i = p->next_vertex_index( iSide ) ;
            }
         }
         else if( inlined )
	 {
            inS( i ) = 0 ;
            ptL->clear() ;
            break ;
	 }
         else
         {
	    curPt = nexPt ;
         }
         cutting_a_polygon = ( !pt_pt_algo->points_are_close( iniPt, curPt ) ) ;
      }

      if( ptL->count()>2 )
      {
         GE_SimplePolygon2D* pol = GE_SimplePolygon2D::create( a_owner, ptL ) ;
         pol->set_point_point_intersector( pt_pt_algo ) ;
         if( is_set_point_segment_intersector() )
         {
            pol->set_point_segment_intersector( point_segment_intersector() ) ;
         }
         if( is_set_segment_segment_intersector() )
         {
            pol->set_segment_segment_intersector( segment_segment_intersector() ) ;
         }
         result->append( pol ) ;
      }
      iStart = inS.index_of( 1 ) ; // Initialisation before a new start
      if( ptL->count()>0 )
      {
         ptL->clear() ;
      }
   }
   p->destroy() ; p = 0 ;
 
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->owner()==a_owner ) ) ;
   PEL_CHECK_POST(
      FORALL( ( i=0; result!=0 && i<result->count(); i++ ),
              result->at( i )->owner()==a_owner ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( i=0; result!=0 && i<result->count(); i++ ),
         dynamic_cast<GE_SimplePolygon2D const*>( result->at( i ) )!=0 ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( i=0 ; result!=0 && i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_point_point_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->point_point_intersector()==point_point_intersector() ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( i=0 ; result!=0 && i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_point_segment_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->point_segment_intersector()==point_segment_intersector() ) );
   PEL_CHECK_POST(
      FORALL(
         ( i=0 ; result!=0 && i<result->count() ; i++ ),
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->is_set_segment_segment_intersector() &&
         static_cast<GE_SimplePolygon2D const*>( result->at( i ) )->segment_segment_intersector()==segment_segment_intersector() ) ) ;
   return( result ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D:: GE_Polygon2D( PEL_Object* a_owner )
//------------------------------------------------------------------------------
      : GE_Polygon( a_owner, 2 ),
        SIMPLE_POLYGONS( 0 ),
        AREA( AREA_NOT_COMPUTED )
{
   PEL_LABEL( "GE_Polygon2D:: GE_Polygon2D" ) ;

   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D:: GE_Polygon2D( PEL_Object* a_owner, 
                             PEL_Sequence const* vertex_table )
//------------------------------------------------------------------------------
      : GE_Polygon( a_owner, 2, vertex_table ),
        SIMPLE_POLYGONS( 0 ),
        AREA( AREA_NOT_COMPUTED )
{
   PEL_LABEL( "GE_Polygon2D:: GE_Polygon2D" ) ;
   
   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D:: GE_Polygon2D( PEL_Object* a_owner, 
                             doubleVector const& coordinate_table )
//------------------------------------------------------------------------------
      : GE_Polygon( a_owner, 2, coordinate_table ),
        SIMPLE_POLYGONS( 0 ),
        AREA( AREA_NOT_COMPUTED )
{
   PEL_LABEL( "GE_Polygon2D:: GE_Polygon2D" ) ;

   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D:: GE_Polygon2D( PEL_Object* a_owner, 
			     GE_Polygon const* other )
//------------------------------------------------------------------------------
      : GE_Polygon( a_owner, other ),
        SIMPLE_POLYGONS( 0 ),
        AREA( AREA_NOT_COMPUTED )
{
   PEL_LABEL( "GE_Polygon2D:: GE_Polygon2D" ) ;
   
   update() ;

   PEL_CHECK_INV( invariant() ) ;
}



//------------------------------------------------------------------------------
GE_Polygon2D:: ~GE_Polygon2D( void )
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: ~GE_Polygon2D" ) ;
}



//------------------------------------------------------------------------------
void
GE_Polygon2D:: compute_area( void ) const
//------------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon2D:: compute_area" ) ;
   PEL_CHECK_INV( invariant() ) ;
   AREA = 0. ;
   if( size()>=3 )
   {
      PEL_Vector const* poly_table = split_of_simple_polygons() ;
      for( size_t i=0 ; i<poly_table->index_limit() ; ++i )
      {
         GE_SimplePolygon2D const* poly =
                static_cast<GE_SimplePolygon2D const*>( poly_table->at(i) ) ;
         AREA += poly->area() ;
      }
   }
}



//------------------------------------------------------------------------------
bool
GE_Polygon2D:: invariant( void ) const
//------------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Polygon::invariant() ) ;
   PEL_ASSERT( dimension()==2 ) ;

   return( true ) ;
}
