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

#include <GE_Polygon_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <PEL_assertions.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Polygon.hh>
#include <GE_Polygon2D.hh>
#include <GE_SimplePolygon2D.hh>
#include <GE_Vector.hh>

#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_SegmentSegment_INT.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <size_t_array2D.hh>

#include <iomanip>
#include <iostream>
#include <string>
 
using std::cout ;
using std::endl ;
using std::string ;

//--------------------------------------------------------------------------
GE_Polygon_TEST*
GE_Polygon_TEST::unique_instance = new GE_Polygon_TEST() ;
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
GE_Polygon_TEST:: GE_Polygon_TEST( void ) 
//--------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Polygon", "GE_Polygon_TEST" )
{
   PEL_LABEL( "GE_Polygon_TEST:: GE_Polygon_TEST" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//--------------------------------------------------------------------------
GE_Polygon_TEST:: ~GE_Polygon_TEST( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: ~GE_Polygon_TEST" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//--------------------------------------------------------------------------
void
GE_Polygon_TEST:: process_all_tests( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: process_all_tests" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( data_deck_explorer()->has_module( "GE_SimplePolygon2D" ) )
   {
      PEL_ModuleExplorer* exp =
         data_deck_explorer()->create_subexplorer( 0, "GE_SimplePolygon2D" ) ;
      exp->start_module_iterator() ;
      for( ; exp->is_valid_module() ; exp->go_next_module() )
      {
          PEL_ModuleExplorer const* texp = exp->create_subexplorer( 0 ) ;
          process_simple_polygon2D_tests( texp ) ;
	  texp->destroy() ; texp = 0 ;
      }
      exp->destroy() ;
   }

   if( data_deck_explorer()->has_module( "GE_Polygon2D" ) )
   {
      PEL_ModuleExplorer* exp =
         data_deck_explorer()->create_subexplorer( 0, "GE_Polygon2D" ) ;
      exp->start_module_iterator() ;
      for( ; exp->is_valid_module() ; exp->go_next_module() )
      {
          PEL_ModuleExplorer const* texp = exp->create_subexplorer( 0 ) ;
          process_polygon2D_tests( texp ) ;
	  texp->destroy() ; texp = 0 ;
      }
      exp->destroy() ;
   }

   PEL_CHECK_INV( invariant() ) ;
}



//--------------------------------------------------------------------------
void
GE_Polygon_TEST:: process_simple_polygon2D_tests(
                                             PEL_ModuleExplorer const* exp )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: process_simple_polygon2D_tests" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const n = "GE_SimplePolygon2D/"+exp->name() ;

   // Polygon building :
   doubleVector const& coordinates =
                          exp->doubleVector_data( "vertices_coordinates" ) ;
   GE_SimplePolygon2D* polygon =
                              GE_SimplePolygon2D::create( 0, coordinates ) ;
   set_geometrical_algorithms( exp, polygon ) ;

   // Chain of vertices :
   bool const ok_chain = test_vertices_chain( exp, polygon ) ;
   notify_one_test_result( n+"/vertices",  ok_chain ) ;

   // Area :
   if( exp->has_entry( "area" ) )
   {
      bool const ok_area = test_area( exp, polygon ) ;
      notify_one_test_result( n+"/area", ok_area ) ;
   }
   
   // Orientation :
   if( exp->has_entry( "orientation" ) )
   {
      bool const ok_orient = test_orientation( exp, polygon ) ;
      notify_one_test_result( n+"/orientation", ok_orient ) ;
   }

   // Triangulation :
   bool const ok_tria = test_triangulation( exp, polygon ) ;
   notify_one_test_result( n+"/triangulation", ok_tria ) ;

   // Inclusion :
   if( exp->has_entry( "points_for_inclusion_test" ) )
   {
      bool const ok_inclu = test_inclusion( exp, polygon ) ;
      notify_one_test_result( n+"/inclusion", ok_inclu ) ;      
   }

   // Bound inclusion :
   if( exp->has_entry( "points_for_inclusion_on_boundary_test" ) )
   {
      bool const ok_inclu_on_bd = test_inclusion_on_boundary( exp,
                                                              polygon ) ;
      notify_one_test_result( n+"/boundary_inclusion", ok_inclu_on_bd ) ;
   }

   // Simplification :
   bool ok_simpli = test_simplification( polygon ) ;
   notify_one_test_result( n+"/simplification", ok_simpli ) ;

   polygon->destroy() ; polygon = 0 ;

   PEL_CHECK_INV( invariant() ) ;
}



//--------------------------------------------------------------------------
void
GE_Polygon_TEST:: process_polygon2D_tests( PEL_ModuleExplorer const* exp )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: process_polygon2D_tests" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const n = "GE_Polygon2D/"+exp->name() ;

   // Polygon building :
   doubleVector const& coordinates =
                      exp->doubleVector_data( "vertices_coordinates" ) ;
   GE_Polygon2D* polygon = GE_Polygon2D::create( 0, coordinates ) ;
   set_geometrical_algorithms( exp, polygon ) ;

   // Chain of vertices :
   bool const ok_chain = test_vertices_chain( exp, polygon ) ;
   notify_one_test_result( n+"/vertices",  ok_chain ) ;

   // Splitting in simple polygon test :
   if( exp->has_entry( "split_of_simple_polygons" ) )
   {
      bool const ok_split = test_split_of_simple_polygons( exp, polygon ) ;
      notify_one_test_result( n+"/splitting_in_simple_polygons", ok_split ) ;
   }
   else
   {
      PEL_Vector* r = polygon->create_split_of_simple_polygons( polygon ) ;
      if( r==0 )
      {
         out() << std::endl ;
         out() << "  Unable to split in simple polygon the polygon : "
                 << std::endl << std::endl ;
         polygon->print( out(), 6 ) ;
         out() << std::endl ;
         notify_one_test_result( n+"/splitting_in_simple_polygons", false ) ;
      }
      else
      {
         notify_one_test_result( n+"/splitting_in_simple_polygons", true ) ;
      }
   }

   // Area :
   if( exp->has_entry( "area" ) )
   {
      bool const ok_area = test_area( exp, polygon ) ;
      notify_one_test_result( n+"/area", ok_area ) ;
   }

   polygon->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
}



//--------------------------------------------------------------------------
void
GE_Polygon_TEST:: set_geometrical_algorithms( PEL_ModuleExplorer const* exp,
                                              GE_Polygon* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: set_geometrical_algorithms" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   // Point-Point :
   GE_PointPoint_INT* pt_pt_algo = 0 ;
   if( exp->has_module( "GE_PointPoint_INT" ) )
   {
      PEL_ModuleExplorer const* pt_pt_exp =
                       exp->create_subexplorer( 0, "GE_PointPoint_INT" ) ;
      std::string const& pt_pt_name =
                               pt_pt_exp->string_data( "concrete_name" ) ;
      
      pt_pt_algo = GE_PointPoint_INT::create( polygon,
                                              pt_pt_name,
                                              pt_pt_exp ) ;
      pt_pt_exp->destroy() ; pt_pt_exp = 0 ;
      polygon->set_point_point_intersector( pt_pt_algo ) ;
   }

   // Point-Segment :
   GE_PointSegment_INT* pt_seg_algo = 0 ;
   if( exp->has_module( "GE_PointSegment_INT" ) )
   {
      PEL_ModuleExplorer const* pt_seg_exp =
                       exp->create_subexplorer( 0, "GE_PointSegment_INT" ) ;
      std::string const& pt_seg_name =
                               pt_seg_exp->string_data( "concrete_name" ) ;
      
      pt_seg_algo = GE_PointSegment_INT::create( polygon,
                                                 pt_seg_name,
                                                 pt_seg_exp,
                                                 pt_pt_algo ) ;
      pt_seg_exp->destroy() ; pt_seg_exp = 0 ;
      polygon->set_point_segment_intersector( pt_seg_algo ) ;
   }

   // Segment-Segment :
   GE_SegmentSegment_INT* seg_seg_algo = 0 ;
   if( exp->has_module( "GE_SegmentSegment_INT" ) )
   {
      PEL_ModuleExplorer const* seg_seg_exp =
                    exp->create_subexplorer( 0, "GE_SegmentSegment_INT" ) ;
      std::string const& seg_seg_name =
                              seg_seg_exp->string_data( "concrete_name" ) ;
      
      seg_seg_algo = GE_SegmentSegment_INT::create( polygon,
                                                    seg_seg_name,
                                                    seg_seg_exp,
                                                    pt_pt_algo,
                                                    pt_seg_algo ) ;
      seg_seg_exp->destroy() ; seg_seg_exp = 0 ;
      polygon->set_segment_segment_intersector( seg_seg_algo ) ;
   }
   
      
}

//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_vertices_chain( PEL_ModuleExplorer const* exp,
                                       GE_Polygon const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_vertices_chain" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool ok = true ;

   doubleVector const& coordinates =
                           exp->doubleVector_data( "vertices_coordinates" ) ;
   
   size_t dim = polygon->dimension() ;
   PEL_List* list = PEL_List::create( 0 ) ;
   PEL_Vector* vector = PEL_Vector::create( 0, polygon->size() ) ;

   for( size_t i=0; i<polygon->size(); i++ )
   {
      ok = ok && ( dim==polygon->vertex( i )->nb_coordinates() ) ;

      GE_Point* vertex_i = const_cast<GE_Point*>( polygon->vertex( i ) ) ;

      for( size_t d=0; d<dim; d++ )
      {
         ok = ok && PEL::equal( vertex_i->coordinate( d ),
                                coordinates( dim*i+d ) ) ;
      }
      list->append( vertex_i ) ;
      vector->set_at( i, vertex_i ) ;

      ok = ok && 
          ( polygon->next_vertex_index( polygon->previous_vertex_index( i ) )==i ) ;
      ok = ok && ( polygon->next_vertex( i )
                        ==polygon->vertex( polygon->next_vertex_index( i ) ) ) ;

      ok = ok && ( polygon->previous_vertex( i )
                        ==polygon->vertex( polygon->previous_vertex_index( i ) ) ) ;
   }

   GE_Polygon2D* poly_vector = GE_Polygon2D::create( 0, vector ) ; 
   ok = ok && poly_vector->is_equal( polygon ) ;
   GE_Polygon2D* poly_list = GE_Polygon2D::create( 0, list ) ; 
   ok = ok && poly_list->is_equal( polygon ) ;

   list->destroy() ;
   vector->destroy() ;
   poly_list->destroy() ;
   poly_vector->destroy() ;

   PEL_CHECK_INV( invariant() ) ;
   
   return( ok ) ;
}



//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_area( PEL_ModuleExplorer const* exp,
                             GE_Polygon const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_area" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   double given_area = exp->double_data( "area" ) ;
   double computed_area = polygon->area() ;
   bool ok = PEL::toler( given_area-computed_area, 1.E-8 ) ;

   if( !ok )
   {
      std::ios_base::fmtflags original_flags = out().flags() ;
      out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      out() << std::setprecision( 9 ) ;
      out() << "   area computed: " << computed_area << std::endl ;
      out() << "   area expected: " << given_area << std::endl ;
      out().flags( original_flags ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   
   return( ok ) ;
}




//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_orientation( PEL_ModuleExplorer const* exp,
                                    GE_SimplePolygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_orientation_of_sp2D" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   int const given_orientation = exp->int_data( "orientation" ) ;
   int const computed_orientation = polygon->orientation() ;
   bool ok = given_orientation==computed_orientation ;

   GE_SimplePolygon2D* inverse_poly = polygon->create_clone( 0 ) ;
   inverse_poly->remove_vertices() ;
   for( size_t i=0 ; i<polygon->size() ; ++i )
   {
      inverse_poly->insert_vertex_at( 0, polygon->vertex(i) ) ;
   }

   ok = ok && PEL::equal( polygon->orientation()*inverse_poly->orientation(), -1 ) ;

   inverse_poly->destroy() ; inverse_poly = 0 ;

   PEL_CHECK_INV( invariant() ) ;
   
   return( ok ) ;
}


//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_triangulation(
                             PEL_ModuleExplorer const* exp,
                             GE_SimplePolygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_triangulation" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_SimplePolygon2D* p = polygon->create_clone(0) ;
   p->suppress_flat_angles() ;
   PEL_Sequence const* connectivity = p->create_triangulation( 0 ) ;
   bool ok = ( connectivity->index_limit()==p->size()-2 ) ;
   p->destroy() ; p = 0 ;
   
   double computed_area = 0. ;
   for( size_t i=0 ; ok && i<connectivity->index_limit() ; ++i )
   {
      PEL_Vector const* tr_pts =
              static_cast<PEL_Vector const*>( connectivity->at(i) ) ;
      ok = ok && ( tr_pts->index_limit()==3 ) ;
      if( ok )
      {
         GE_Mpolyhedron const* tr = GE_Mpolyhedron::create( 0,
                                                            "GE_Triangle",
                                                            tr_pts ) ;
         computed_area += tr->measure() ;
         tr->destroy() ; tr = 0 ;
      }
   }
   ok = ok && PEL::toler( computed_area-polygon->area(), 1.E-8 ) ;
   if( !ok )
   {
      out() << "  Polygon area = " << polygon->area() << std::endl ;
      out() << "  Triangulation area = " << computed_area << std::endl ;
      polygon->print( out(), 2 ) ;
      out() << "  Triangulation : " << std::endl ;
      for( size_t i=0 ; i<connectivity->index_limit() ; ++i )
      {
         PEL_Vector const* tr_pts =
                     static_cast<PEL_Vector const*>( connectivity->at(i) ) ;
         out() << "   " ;
         for( size_t j=0 ; j<tr_pts->index_limit() ; ++j )
         {
            GE_Point const* pt =
                             static_cast<GE_Point const*>( tr_pts->at(j) ) ;
            pt->print( out(), 2 ) ;
         }
         out() << std::endl ;
      }
   }

   connectivity->destroy() ; connectivity = 0 ;
   
   PEL_CHECK_INV( invariant() ) ;

   return( ok ) ;
}


//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_inclusion( PEL_ModuleExplorer const* exp,
                                  GE_SimplePolygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_inclusion" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( exp->has_entry( "points_for_inclusion_test" ) ) ;
   PEL_CHECK( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Sequence const* connectivity = polygon->create_triangulation( 0 ) ;
   PEL_List* list_of_triangles = PEL_List::create( 0 ) ;
   for( size_t i=0 ; i<connectivity->index_limit() ; ++i )
   {
      PEL_Vector const* tr_pts =
              static_cast<PEL_Vector const*>( connectivity->at(i) ) ;
      GE_Mpolyhedron const* tr =
            GE_Mpolyhedron::create( list_of_triangles,
                                    "GE_Triangle",
                                    tr_pts ) ;
      list_of_triangles->append( const_cast<GE_Mpolyhedron*>( tr ) ) ;
   }   

   doubleVector const& coordinates = 
      exp->doubleVector_data( "points_for_inclusion_test" ) ;
   size_t points_nb = (size_t)( coordinates.size()/2. ) ;
   bool ok = true ;
   for( size_t i=0; i<points_nb; ++i )
   {
      GE_Point* pt = GE_Point::create( 0, coordinates( 2*i ),
                                          coordinates( 2*i+1 ) ) ;
      bool inside = polygon->has_in_interior( pt ) ;
      PEL_ListIterator* it = list_of_triangles->create_iterator( 0 ) ;
      bool really_inside = false ;
      for( it->start(); !really_inside && it->is_valid(); it->go_next() )
      {
         GE_Mpolyhedron const* tr =
                     static_cast<GE_Mpolyhedron const*>( it->item() ) ;
         really_inside = tr->contains( pt ) ;
      }
      ok = ok && ( inside==really_inside ) ;
      it->destroy() ; it = 0 ;
      pt->destroy() ; pt = 0 ;
   }

   connectivity->destroy() ; connectivity = 0 ;
   list_of_triangles->destroy() ; list_of_triangles = 0 ;

   PEL_CHECK_INV( invariant() ) ;

   return( ok ) ;
}


//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_inclusion_on_boundary( 
                                  PEL_ModuleExplorer const* exp,
                                  GE_SimplePolygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_triangulation_on_boundary" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( exp->has_entry( "points_for_inclusion_on_boundary_test" ) ) ;
   PEL_CHECK_PRE( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   doubleVector const& coordinates = 
      exp->doubleVector_data( "points_for_inclusion_on_boundary_test" ) ;
   size_t points_nb = (size_t)( coordinates.size()/2. ) ;

   GE_PointSegment_INT* PT_SEG_ALGO = polygon->point_segment_intersector() ;
   
   bool ok = true ;
   for( size_t i=0; i<points_nb; i++ )
   {
      GE_Point* pt = GE_Point::create( 0, coordinates( 2*i ),
                                          coordinates( 2*i+1 ) ) ;
      bool on_boundary = polygon->has_on_boundary( pt ) ;
      bool really_on_boundary =
         PT_SEG_ALGO->point_in_segment( pt,
                                        polygon->vertex( polygon->size()-1 ),
                                        polygon->vertex( 0 ) ) ;
      for( size_t j=0; !really_on_boundary && j<polygon->size()-1; j++ )
      {
         really_on_boundary =
            PT_SEG_ALGO->point_in_segment( pt,
                                           polygon->vertex( j ),
                                           polygon->next_vertex( j ) ) ;
      }
      ok = ok && ( on_boundary==really_on_boundary ) ;
      pt->destroy() ;
   }

   PEL_CHECK_INV( invariant() ) ;

   return( ok ) ;
}



//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_simplification(
                                   GE_SimplePolygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_simplification" ) ;
   PEL_CHECK_PRE( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Vector* ps = GE_Vector::create( 0, 2 ) ;
   GE_Vector* cs = GE_Vector::create( 0, 2 ) ;

   GE_SimplePolygon2D* simple_polygon = polygon->create_clone( 0 ) ;
   double sine_of_flat_angle = 1.E-6 ;
   simple_polygon->suppress_flat_angles( sine_of_flat_angle ) ;   
   bool ok = true ;
   double const val_min = sine_of_flat_angle*sine_of_flat_angle ;
   for( size_t i=0 ; ok && i<simple_polygon->size() ; ++i )
   {
      ps->re_initialize( simple_polygon->previous_vertex( i ),
                         simple_polygon->vertex( i ) ) ;
      cs->re_initialize( simple_polygon->vertex( i ),
                         simple_polygon->next_vertex( i ) ) ;
      double const c = ps->cosine( cs ) ;      
      ok = PEL::abs( 1.-c*c )>=val_min ;
   }
   ok = ok && ( PEL::toler( polygon->area()-simple_polygon->area(), 1.E-8 ) ) ;

   ps->destroy() ; ps = 0 ;
   cs->destroy() ; cs = 0 ;
   simple_polygon->destroy() ; simple_polygon = 0 ;

   PEL_CHECK_INV( invariant() ) ;

   return( ok ) ;
}



//--------------------------------------------------------------------------
bool
GE_Polygon_TEST:: test_split_of_simple_polygons(
                                         PEL_ModuleExplorer const* exp,
                                         GE_Polygon2D const* polygon ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Polygon_TEST:: test_split_of_simple_polygons" ) ;
   PEL_CHECK( exp!=0 ) ;
   PEL_CHECK( exp->has_entry( "split_of_simple_polygons" ) ) ;
   PEL_CHECK_PRE( polygon!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Vector* sp = polygon->create_split_of_simple_polygons( 0 ) ;
   bool ok =  ( sp!=0 );
   if( !ok )
   {
      out() << std::endl ;
      out() << "  Unable to split in simple polygon the polygon : "
                << std::endl << std::endl ;
      polygon->print( out(), 6 ) ;
      out() << std::endl ;
   }
   else
   {
      doubleArray2D const& sp_coords =
         exp->doubleArray2D_data( "split_of_simple_polygons" ) ;

      ok = ( sp->index_limit()==sp_coords.index_bound(0) ) ;

      for( size_t i=0 ; ok && i<sp->index_limit() ; ++i )
      {
         GE_SimplePolygon2D const* poly =
            static_cast<GE_SimplePolygon2D const*>( sp->at(i) ) ;
         ok = ( 2*poly->size()<=sp_coords.index_bound(1) ) ;
         for( size_t j=0 ; ok && j<poly->size() ; ++j )
         {
            GE_Point const* pt = poly->vertex(j) ;
            ok = PEL::toler( pt->coordinate(0)-sp_coords(i,2*j), 1.E-8 ) &&
               PEL::toler( pt->coordinate(1)-sp_coords(i,2*j+1), 1.E-8 ) ;
         }
         for( size_t j=2*poly->size() ; ok && j<sp_coords.index_bound(1) ; ++j )
         {
            ok = PEL::equal( sp_coords(i,j), 0. ) ;
         }
      }
      if( !ok )
      {
         out() << "Number of simple polygons : "
                   << sp->index_limit() << std::endl ;
         for( size_t i=0 ; i<sp->index_limit() ; ++i )
         {
            GE_SimplePolygon2D const* poly =
               static_cast<GE_SimplePolygon2D const*>( sp->at(i) ) ;
            out() << "Simple polygon #" << i << " : " << std::endl ;
            poly->print( out(), 5 ) ;
         }
      }

      sp->destroy_items_and_clear() ;
      sp->destroy() ; sp = 0 ;
   }
   
   PEL_CHECK_INV( invariant() ) ;

   return( ok ) ;
}
