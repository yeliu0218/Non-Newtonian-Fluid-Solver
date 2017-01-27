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

#include <PDE_CFootFinder.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_GridFE.hh>

#include <PDE_ForwardEulerFinder.hh>

#include <iostream>
#include <sstream>

using std::string ;
using std::endl ;
using std::ostringstream ;

//--------------------------------------------------------------------------
PDE_CFootFinder*
PDE_CFootFinder:: create( PDE_LocalFEcell* fem,
                          std::string const& concrete_name,
                          PEL_ModuleExplorer const* exp,
                          PDE_GridFE const* grid )
//--------------------------------------------------------------------------
{
   PEL_CHECK( fem != 0 ) ;
   PEL_CHECK( !concrete_name.empty() ) ;
   PEL_CHECK( exp != 0 ) ;
   PEL_CHECK( grid != 0 ) ;

   PDE_CFootFinder* result = 0 ;
   if( concrete_name=="PDE_ForwardEulerFinder" )
   {
      result = PDE_ForwardEulerFinder:: create( fem, exp, grid ) ;
   }
   else
   {
      std::string allowed_values ;
      allowed_values += "   \"PDE_ForwardEulerFinder\"\n" ;
      PEL_Error::object()->raise_bad_data_value( exp,
                                                 "concrete_name",
                                                 allowed_values ) ;
   }
   
   PEL_CHECK( result != 0 ) ;
   PEL_CHECK( result->owner() == fem ) ;
   PEL_CHECK( !result->search_has_been_performed() ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PDE_CFootFinder:: PDE_CFootFinder( PDE_LocalFEcell* fem,
                                   PEL_ModuleExplorer const* exp )
//--------------------------------------------------------------------------
   : PEL_Object( fem )
   , SEARCH_DONE( false )
   , HEAD_CELL( 0 )
   , FOOT_CELL( 0 )
   , FOOT_BOUND( 0 )
{
   PEL_LABEL( "PDE_CFootFinder:: PDE_CFootFinder" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//--------------------------------------------------------------------------
PDE_CFootFinder:: ~PDE_CFootFinder( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: ~PDE_CFootFinder" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//--------------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_CFootFinder:: polyhedron_of_head_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: polyhedron_of_head_cell" ) ;
   
   GE_Mpolyhedron const* result = ( HEAD_CELL == 0 ? 
                                       0 : HEAD_CELL->polyhedron() ) ;
   
   return( result ) ;
}

//--------------------------------------------------------------------------
size_t
PDE_CFootFinder:: mesh_id_of_head_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: mesh_id_of_head_cell" ) ;
   PEL_CHECK_PRE( polyhedron_of_head_cell() != 0 ) ;

   return( HEAD_CELL->id_number() ) ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder:: search_foot( PDE_DiscreteField const* aa,
                                              size_t level,
                                              double time_step,
                                              GE_Point const* head )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: search_foot" ) ;
   PEL_CHECK_PRE( polyhedron_of_head_cell() != 0 ) ;
   PEL_CHECK_PRE( aa != 0 ) ;
   PEL_CHECK_PRE( level < aa->storage_depth() ) ;
   PEL_CHECK_PRE( head != 0 ) ;
   PEL_CHECK_PRE( polyhedron_of_head_cell()->contains( head ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   FOOT_CELL = 0 ;
   FOOT_BOUND = 0 ;
   SEARCH_DONE = false ;   
   HEAD = head ;
   search_foot( aa, level, time_step, head, HEAD_CELL ) ;
   SEARCH_DONE = true ;

   PEL_CHECK_POST( search_has_been_performed() ) ;
   PEL_CHECK_POST( IMPLIES( foot_has_been_found(), foot()!=0 ) ) ;
   PEL_CHECK_POST( IMPLIES( foot_has_been_found(), foot_mesh()!=0 ) ) ;
   PEL_CHECK_POST( 
      IMPLIES( foot_has_been_found(), 
               foot_mesh()->polyhedron()->contains(foot()) )) ;
   PEL_CHECK_POST( 
      IMPLIES( foot_has_been_found(),
               foot_is_on_boundary() || foot_is_interior() ) ) ;
   PEL_CHECK_POST( IMPLIES( foot_is_interior(),
                            foot_cell()==foot_mesh() ) ) ;
   PEL_CHECK_POST( IMPLIES( foot_is_on_boundary(),
                            foot_bound()==foot_mesh() ) ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: search_has_been_performed( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: search_has_been_performed" ) ;

   return( SEARCH_DONE ) ;
}

//--------------------------------------------------------------------------
double
PDE_CFootFinder:: value_at_foot( PDE_DiscreteField const* f,
                                 size_t level,
                                 size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: value_at_foot" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( f != 0 ) ;
   PEL_CHECK_PRE( level < f->storage_depth() ) ;
   PEL_CHECK_PRE( ic < f->nb_components() ) ;
   PEL_CHECK_PRE( foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( foot_mesh()->value( f, level, foot_ref(), ic ) ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_has_been_found( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_has_been_found" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( FOOT_CELL!=0 || FOOT_BOUND!=0 ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_is_interior( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_is_interior" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( FOOT_CELL!=0 ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_is_on_boundary( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_is_on_boundary" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( FOOT_BOUND!=0 ) ;
}

//--------------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_CFootFinder:: polyhedron_of_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: polyhedron_of_foot_cell" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Mpolyhedron const* result = foot_mesh()->polyhedron() ;
   
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
GE_Color const*
PDE_CFootFinder:: color_of_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: color_of_foot_cell" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Color const* result = foot_mesh()->color() ;
   
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
size_t
PDE_CFootFinder:: mesh_id_of_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: mesh_id_of_foot_cell" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( foot_mesh()->id_number() ) ;
}

//--------------------------------------------------------------------------
double
PDE_CFootFinder:: CFL_number( size_t dim ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: CFL_number" ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( foot_has_been_found() ) ;
   PEL_CHECK_PRE( dim<foot()->nb_coordinates() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double d_mean =
          HEAD_CELL->polyhedron()->inter_vertices_maximum_distance( dim ) ;
   if( foot_is_interior() )
   {
      d_mean =
         0.5*( d_mean
              +foot_cell()->polyhedron()->inter_vertices_maximum_distance(dim) ) ;
   }
   return( PEL::abs( HEAD->coordinate( dim ) - foot()->coordinate( dim ) )
           /d_mean ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: is_ersatz_foot_investigated( void ) const
//--------------------------------------------------------------------------
{
   return( false ) ;
}

//--------------------------------------------------------------------------
double
PDE_CFootFinder:: value_at_ersatz_foot( PDE_DiscreteField const* f,
                                         size_t level,
                                         size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: value_at_ersatz_foot" ) ;
   PEL_CHECK_PRE( is_ersatz_foot_investigated() ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( f != 0 ) ;
   PEL_CHECK_PRE( level < f->storage_depth() ) ;
   PEL_CHECK_PRE( ic < f->nb_components() ) ;
   PEL_CHECK_PRE( !foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( ersatz_foot_cell()->value( f, level,
                                      ersatz_foot_ref(),
                                      ic ) ) ;
}

//--------------------------------------------------------------------------
GE_Mpolyhedron const*
PDE_CFootFinder:: polyhedron_of_ersatz_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: polyhedron_of_ersatz_foot_cell" ) ;
   PEL_CHECK_PRE( is_ersatz_foot_investigated() ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( !foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Mpolyhedron const* result = ersatz_foot_cell()->polyhedron() ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
GE_Color const*
PDE_CFootFinder:: color_of_ersatz_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: color_of_ersatz_foot_cell" ) ;
   PEL_CHECK_PRE( is_ersatz_foot_investigated() ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( !foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Color const* result = ersatz_foot_cell()->color() ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
size_t
PDE_CFootFinder:: mesh_id_of_ersatz_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: mesh_id_of_ersatz_foot_cell" ) ;
   PEL_CHECK_PRE( is_ersatz_foot_investigated() ) ;
   PEL_CHECK_PRE( search_has_been_performed() ) ;
   PEL_CHECK_PRE( !foot_has_been_found() ) ;
   PEL_CHECK_INV( invariant() ) ;
    
   return( ersatz_foot_cell()->id_number() ) ;
}

//--------------------------------------------------------------------------
GE_Point const*
PDE_CFootFinder:: ersatz_foot( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: ersatz_foot" ) ;
   PEL_CHECK_PRE( ersatz_foot_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Point const* result = 0 ;
   
   PEL_CHECK_POST( ersatz_foot_POST( result ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
GE_Point const*
PDE_CFootFinder:: ersatz_foot_ref( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: ersatz_foot_ref" ) ;
   PEL_CHECK_PRE( ersatz_foot_ref_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Point const* result = 0 ;
   
   PEL_CHECK_POST( ersatz_foot_ref_POST( result ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder:: set_found_cell( PDE_CellFE const* a_cell )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: set_found_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !search_has_been_performed() ) ;
   PEL_CHECK( a_cell!=0 ) ;
   PEL_CHECK( FOOT_BOUND == 0 ) ;
   FOOT_CELL = a_cell ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder:: set_found_bound( PDE_BoundFE const* a_bound )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: set_found_bound" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( !search_has_been_performed() ) ;
   PEL_CHECK( a_bound!=0 ) ;
   PEL_CHECK( FOOT_CELL == 0 ) ;
   FOOT_BOUND = a_bound ;
}

//--------------------------------------------------------------------------
PDE_BoundFE const*
PDE_CFootFinder:: foot_bound( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_bound" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( search_has_been_performed() ) ;
   PEL_CHECK( foot_is_on_boundary() ) ;
   PEL_CHECK( FOOT_CELL == 0 ) ;

   PDE_BoundFE const* result = FOOT_BOUND ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PDE_CellFE const*
PDE_CFootFinder:: foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( search_has_been_performed() ) ;
   PEL_CHECK( foot_is_interior() ) ;
   PEL_CHECK( FOOT_BOUND == 0 ) ;

   PDE_CellFE const* result = FOOT_CELL ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PDE_MeshFE const*
PDE_CFootFinder:: foot_mesh( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: foot_mesh" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( search_has_been_performed() ) ;
   PEL_CHECK( foot_has_been_found() ) ;

   PDE_MeshFE const* result = 0 ;
   ( FOOT_CELL != 0 ? result = FOOT_CELL : result = FOOT_BOUND ) ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PDE_CellFE const*
PDE_CFootFinder:: ersatz_foot_cell( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder:: ersatz_foot_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( ersatz_foot_cell_PRE() ) ;

   PDE_CellFE const* result = 0 ;

   PEL_CHECK_POST( ersatz_foot_cell_POST( result ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder:: raise_periodicity_error( void ) const
//--------------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PDE_CFootFinder error" << endl ;
   mesg << "    domains with periodicity are not handled" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//--------------------------------------------------------------------------
bool 
PDE_CFootFinder:: search_foot_PRE( PDE_DiscreteField const* aa,
                                   size_t level,
                                   double time_step,
                                   GE_Point const* head,
                                   PDE_CellFE const* head_cell ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( !search_has_been_performed() ) ;
   PEL_ASSERT( aa != 0 ) ;
   PEL_ASSERT( level < aa->storage_depth() ) ;
   PEL_ASSERT( head_cell != 0 ) ;
   PEL_ASSERT( head != 0 ) ;
   PEL_ASSERT( head_cell->polyhedron()->contains( head ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool PDE_CFootFinder:: search_foot_POST( PDE_DiscreteField const* aa,
                                         size_t level,
                                         double time_step,
                                         GE_Point const* head,
                                          PDE_CellFE const* head_cell ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( foot_has_been_found(), foot()!=0 ) ) ;
   PEL_ASSERT( IMPLIES( foot_has_been_found(), foot_mesh()!=0 ) ) ;
   PEL_ASSERT( 
      IMPLIES( foot_has_been_found(), 
               foot_mesh()->polyhedron()->contains(foot()) )) ;
   PEL_ASSERT( 
      IMPLIES( foot_has_been_found(),
               foot_is_on_boundary() || foot_is_interior() ) ) ;
   PEL_ASSERT( IMPLIES( foot_is_interior(), foot_cell()==foot_mesh() ) ) ;
   PEL_ASSERT( IMPLIES( foot_is_on_boundary(), foot_bound()==foot_mesh() ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_PRE( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( search_has_been_performed() ) ;
   PEL_ASSERT( foot_has_been_found() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_POST( GE_Point const* result ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( foot_mesh()->polyhedron()->contains( result ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_ref_PRE( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( search_has_been_performed() ) ;
   PEL_ASSERT( foot_has_been_found() ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: foot_ref_POST( GE_Point const* result ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( foot_mesh()->polyhedron()->reference_polyhedron()->contains( result ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_PRE( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_ersatz_foot_investigated() ) ;
   PEL_ASSERT( search_has_been_performed() ) ;
   PEL_ASSERT( !foot_has_been_found() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_POST( GE_Point const* result ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( ersatz_foot_cell()->polyhedron()->contains( result ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_ref_PRE( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_ersatz_foot_investigated() ) ;
   PEL_ASSERT( search_has_been_performed() ) ;
   PEL_ASSERT( !foot_has_been_found() ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_ref_POST( GE_Point const* result ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( ersatz_foot_cell()->polyhedron()->reference_polyhedron()->contains( result ) ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_cell_PRE( void ) const 
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_ersatz_foot_investigated() ) ;
   PEL_ASSERT( search_has_been_performed() ) ;
   PEL_ASSERT( !foot_has_been_found() ) ;

   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: ersatz_foot_cell_POST( PDE_CellFE const* result ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_CFootFinder:: invariant( void ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder:: synchronize( PDE_LocalFEcell const* fem,
                               PDE_CellFE const* initial_cell )
//--------------------------------------------------------------------------
{
   PEL_CHECK( initial_cell != 0 ) ;
   
   HEAD_CELL = initial_cell ;
   FOOT_CELL = 0 ;
   FOOT_BOUND = 0 ;
   SEARCH_DONE = false ;
   
   PEL_CHECK( !search_has_been_performed() ) ;
}
