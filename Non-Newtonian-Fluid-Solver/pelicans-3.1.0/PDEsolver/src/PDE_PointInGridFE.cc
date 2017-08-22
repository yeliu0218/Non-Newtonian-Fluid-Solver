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

#include <PDE_PointInGridFE.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <boolVector.hh>
#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SegmentPolyhedron_INT.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>

#include <sstream>

using std::ostringstream ;
using std::endl ;

struct PDE_PointInGridFE_ERROR
{
   static void n0( GE_Point const* pt,
                   PDE_LocalFEcell const* cFE,
                   PDE_CursorFEside* sFE,
                   PDE_LocalFEbound* bFE ) ; 
   static void n1( GE_Point const* pt ) ; 
} ;

//----------------------------------------------------------------------
PDE_PointInGridFE*
PDE_PointInGridFE:: create( PEL_Object* a_owner, 
                            PDE_DomainAndFields const* dom,
                            PEL_ModuleExplorer const* exp )  
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: create" ) ;
   PEL_CHECK_PRE( dom!= 0 ) ;
   PEL_CHECK_PRE( dom->nb_space_dimensions()==2 ||
                  dom->nb_space_dimensions()==3 ) ;
 
   PDE_PointInGridFE* result = new PDE_PointInGridFE( a_owner, dom, exp ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == dom->nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_PointInGridFE:: PDE_PointInGridFE( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DIM( dom->nb_space_dimensions() )
   , INTERSECTOR( 0 )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
{
   PEL_LABEL( "PDE_PointInGridFE:: PDE_PointInGridFE" ) ;
   PEL_CHECK_INV( invariant() ) ;

   // Intersector:
   {
      std::string name ;
      PEL_ModuleExplorer* se = 0 ;
      if( exp != 0 )
      {
         se = exp->create_subexplorer( 0, "GE_SegmentPolyhedron_INT" ) ;
         name = se->string_data( "concrete_name" ) ;
      }
      else
      {
         name = (  DIM == 2 ? "GE_SegmentPolyhedron1D_INT"
                            : "GE_SegmentPolyhedron2D_INT" ) ;
      }
      INTERSECTOR = GE_SegmentPolyhedron_INT::make( this, name, DIM, se ) ;
      
      if( se != 0 ) se->destroy() ;
   }

   // Cell face iterators:
   if( sFE->is_excluded( GE_Color::halo_color() ) )
   {
      sFE->include_color( GE_Color::halo_color() ) ;
   }
   if( bFE->is_excluded( GE_Color::halo_color() ) )
   {
      bFE->include_color( GE_Color::halo_color() ) ;
   }
}

//----------------------------------------------------------------------
PDE_PointInGridFE:: ~PDE_PointInGridFE( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: ~PDE_PointInGridFE" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_PointInGridFE:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( DIM ) ;
}

//----------------------------------------------------------------------
bool 
PDE_PointInGridFE:: is_in_grid( GE_Point const* pt,
                                PDE_LocalFEcell* cFE ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: is_in_grid" ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( cFE != 0 ) ;
   PEL_CHECK_PRE( cFE->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( cFE->is_valid() ) ;
   PEL_CHECK_PRE( ! cFE->is_excluded( GE_Color::halo_color() ) ) ;
   
   bool result = false ;
   boolVector visited_cells(0) ;
   size_t next_cell_id = PEL::bad_index() ;
   
   bool cont = true ;
   while( cont && !result )
   {
      result = ( cFE->polyhedron()->contains( pt ) ) ;
      if( !result )
      {
         size_t const id_cell = cFE->mesh_id() ;
         set_visited( visited_cells, id_cell ) ;
         
         // Search intersection with sides:
         next_cell_id = cell_neighbour_id( pt, cFE, visited_cells ) ;
         if( next_cell_id != PEL::bad_index() )
         {
            cont = true ;
            cFE->go_i_th( next_cell_id ) ;
         }
         else
         {
            // Search intersection with bounds:
            next_cell_id = cell_id_after_hole( pt, cFE, visited_cells ) ;
            if( next_cell_id != PEL::bad_index() )
            {
               cont = true ;
               cFE->go_i_th( next_cell_id ) ;
            }
            else
            {
               cont = false ;
            }
         }
      }
   }

//    if( !result )
//    {
//       search_in_all_cells( pt, cFE, result ) ;
//    }
   
   PEL_CHECK_POST( IMPLIES( result, cFE->polyhedron()->contains( pt ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_PointInGridFE:: cell_neighbour_id(
                                GE_Point const* pt,
                                PDE_LocalFEcell const* cFE,
                                boolVector const& visited_cells ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: cell_neighbour_id" ) ;
   PEL_CHECK( pt != 0 ) ;
   PEL_CHECK( cFE != 0 ) ;
   PEL_CHECK( cFE->is_valid() ) ;
   PEL_CHECK( is_visited( visited_cells, cFE->mesh_id() ) ) ;

   size_t result = PEL::bad_index() ;

   bool found = false ;
   size_t const id_cell = cFE->mesh_id() ;
   GE_Point const* cell_center = cFE->polyhedron()->center() ;
   size_t_vector const& sid = cFE->adjacent_side_ids() ;
   for( size_t i=0 ; !found && i<sid.size() ; ++i )
   {
      sFE->go_i_th( sid( i ) ) ;
      if( !sFE->is_periodic() )
      {
         INTERSECTOR->check_intersection( cell_center, pt,
                                          sFE->polyhedron() ) ;
         if( INTERSECTOR->one_single_intersection() )
         {
            size_t id = sFE->adjacent_localFEcell(0)->mesh_id() ;
            if( id==id_cell )
            {
               id = sFE->adjacent_localFEcell(1)->mesh_id() ;
            }
            if( !is_visited( visited_cells, id ) )
            {
               found = true ;
               result = id ;
            }
         }
      }
      else
      {
         size_t const id =
            ( id_cell == sFE->adjacent_localFEcell(0)->mesh_id() ? 0 : 1 ) ;
         INTERSECTOR->check_intersection( cell_center, pt,
                                          sFE->polyhedron(id) ) ;
         if( INTERSECTOR->one_single_intersection() )
         {
            PDE_PointInGridFE_ERROR::n1( pt ) ;
         }
      }
   }

   PEL_CHECK( IMPLIES( result != PEL::bad_index(),
                       !is_visited( visited_cells, result ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_PointInGridFE:: cell_id_after_hole(
                                GE_Point const* pt,
                                PDE_LocalFEcell const* cFE,
                                boolVector const& visited_cells ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: cell_id_after_hole" ) ;
   PEL_CHECK( pt != 0 ) ;
   PEL_CHECK( cFE != 0 ) ;
   PEL_CHECK( cFE->is_valid() ) ;
   PEL_CHECK( is_visited( visited_cells, cFE->mesh_id() ) ) ;

   size_t result = PEL::bad_index() ;
   
   bool found = false ;
   GE_Point const* cell_center = cFE->polyhedron()->center() ;
   size_t_vector const& bid = cFE->adjacent_bound_ids() ;
   for( size_t i=0 ; !found && i<bid.size() ; ++i )
   {
      bFE->go_i_th( bid( i ) ) ;
      INTERSECTOR->check_intersection( cell_center, pt,
                                       bFE->polyhedron() ) ;
      found = ( INTERSECTOR->one_single_intersection() ) ;
   }

   if( !found ) PDE_PointInGridFE_ERROR::n0( pt, cFE, sFE, bFE ) ;

   // Search next cell after the hole:
   GE_Point const* bd_pt = bFE->polyhedron()->center() ;
   found = false ;
   for( bFE->start() ; !found && bFE->is_valid() ; bFE->go_next() )
   {
      size_t const id = bFE->adjacent_cell_id() ;
      if( !is_visited( visited_cells, id ) )
      {
         INTERSECTOR->check_intersection( bd_pt, pt,
                                          bFE->polyhedron() ) ;
         if( INTERSECTOR->one_single_intersection() )
         {
            found = true ;
            result = id ;
         }
      }
   }
   
   PEL_CHECK( IMPLIES( result != PEL::bad_index(),
                       !is_visited( visited_cells, result ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void 
PDE_PointInGridFE:: search_in_all_cells( GE_Point const* pt,
                                         PDE_LocalFEcell* cFE,
                                         bool& found ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: search_in_all_cells" ) ;
   PEL_CHECK( pt != 0 ) ;
   PEL_CHECK( pt->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK( cFE != 0 ) ;
   PEL_CHECK( cFE->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK( !found ) ;

   std::ostringstream msg ;
   msg << "*** PDE_PointInGridFE: heuristical search failure" << std::endl ;
   msg << "    point: " ;
   pt->print( msg, 0 ) ;
   msg << std::endl ;
   msg << "    switch to complete search..." ;
   PEL_Error::object()->display_info( msg.str() ) ;

   cFE->start() ;
   while( !found && cFE->is_valid() )
   {
      found = ( cFE->polyhedron()->contains( pt ) ) ;
      if( !found ) cFE->go_next() ;
   }
   
   PEL_CHECK_POST( IMPLIES( found, cFE->polyhedron()->contains( pt ) ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_PointInGridFE:: is_visited( boolVector const& visited_cells,
                                size_t id_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: is_visited" ) ;
   return( id_cell<visited_cells.size() && visited_cells(id_cell) ) ;
}

//----------------------------------------------------------------------
void
PDE_PointInGridFE:: set_visited( boolVector& visited_cells,
                                 size_t id_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE:: set_visited" ) ;
   PEL_CHECK( !is_visited( visited_cells, id_cell ) ) ;
   if( id_cell>=visited_cells.size() ) visited_cells.resize( id_cell+1 ) ;
   visited_cells(id_cell) = true ;
   PEL_CHECK_POST( is_visited( visited_cells, id_cell ) ) ;
}

//internal--------------------------------------------------------------
void
PDE_PointInGridFE_ERROR:: n0( GE_Point const* pt,
                              PDE_LocalFEcell const* cFE,
                              PDE_CursorFEside* sFE,
                              PDE_LocalFEbound* bFE )
//internal--------------------------------------------------------------
{
   // No intersection: geometry failure
   ostringstream msg ;
   msg << "*** PDE_PointInGridFE: internal failure" << endl << endl ;
   msg << "The following contradiction occurs." << endl << endl ;
   msg << "1. The point " << endl ;
   msg << "      PT=" ; pt->print( msg, 0 ) ; msg << endl ;
   msg << "   is detected not to be in the CURRENT CELL." << endl ;
   msg << "2. The segment between PT and the center of the" << endl ;
   msg << "   CURRENT CELL does not intersect the boundary" << endl ;
   msg << "   of that CURRENT CELL." << endl << endl ;
   msg << "Characteristics of the CURRENT CELL:" << endl << endl ;
   cFE->print_current_mesh( msg, 3 ) ;
   msg << endl ;
   
   size_t_vector const& sid = cFE->adjacent_side_ids() ;
   if( sid.size() == 0 )
      msg << "   no side" << endl ;
   for( size_t i=0 ; i<sid.size() ; ++i )
   {
      msg << "   Side " << i << endl ;
      sFE->go_i_th( sid( i ) ) ;
      sFE->print_current_mesh( msg, 6 ) ;
   }
   
   size_t_vector const& bid = cFE->adjacent_bound_ids() ;
   if( sid.size() == 0 )
      msg << "   no bound" << endl ;
   for( size_t i=0 ; bid.size() ; ++i )
   {
      msg << "   Bound " << i << endl ;
      bFE->go_i_th( bid( i ) ) ;
      sFE->print_current_mesh( msg, 6 ) ;
   }
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_PointInGridFE_ERROR:: n1( GE_Point const* pt )
//internal--------------------------------------------------------------
{
   // No intersection: geometry failure
   ostringstream msg ;
   msg << "*** PDE_PointInGridFE error" << endl << endl ;
   msg << "    A periodic boundary side is crossed." << endl ;
   msg << "    (the case of periodic domain with hole s not implemented)" <<  endl ;;
   msg << "    Fails to find the point " ; pt->print( msg, 0 ) ;
   msg << " in the grid." ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

