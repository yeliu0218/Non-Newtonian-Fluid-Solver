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

#include <GE_SegmentPolyhedron_INT.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <iostream>

//----------------------------------------------------------------------
GE_SegmentPolyhedron_INT*
GE_SegmentPolyhedron_INT:: make( PEL_Object* a_owner,
                                 std::string const& a_name,
                                 size_t nb_space_dim,
                                 PEL_ModuleExplorer const* a_mod_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: make" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( nb_space_dim == 2 || nb_space_dim == 3 ) ;

   GE_SegmentPolyhedron_INT const* proto =
      static_cast<GE_SegmentPolyhedron_INT const*>(
                                       plugins_map()->item( a_name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
   if( proto->nb_space_dimensions() != nb_space_dim )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_SegmentPolyhedron_INT error:\n"
         "    object of type \""+a_name+"\"\n"
         "    is not valid for this number of space dimensions" ) ;
   }
   
   GE_SegmentPolyhedron_INT* result =
                          proto->create_replica( a_owner, a_mod_exp ) ;
   PEL_CHECK( !result->is_a_prototype() ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dim ) ;
   PEL_CHECK_POST( !result->intersection_checked() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
GE_SegmentPolyhedron_INT:: GE_SegmentPolyhedron_INT(
                        std::string const& a_name, size_t nb_space_dim )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , PROTO( true )
   , DIM( nb_space_dim )
   , INTER_CHECKED( false )
   , M_SAVE( 0 )
   , S0_SAVE( 0 )
   , S1_SAVE( 0 )
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: GE_SegmentPolyhedron_INT" ) ;

   plugins_map()->register_item( a_name, this ) ;
   
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_SegmentPolyhedron_INT:: GE_SegmentPolyhedron_INT(
                              PEL_Object* a_owner, size_t nb_space_dim )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , PROTO( false )
   , DIM( nb_space_dim ) 
   , INTER_CHECKED( false )
   , M_SAVE( 0 )
   , S0_SAVE( 0 )
   , S1_SAVE( 0 )
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: GE_SegmentPolyhedron_INT" ) ;
   
   PEL_CHECK( owner()==a_owner ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_SegmentPolyhedron_INT:: ~GE_SegmentPolyhedron_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: ~GE_SegmentPolyhedron_INT" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentPolyhedron_INT:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( PROTO ) ;
}

//----------------------------------------------------------------------
size_t
GE_SegmentPolyhedron_INT:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( DIM ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentPolyhedron_INT:: intersection_checked( void ) const
//----------------------------------------------------------------------
{
   return( INTER_CHECKED ) ;   
}

//----------------------------------------------------------------------
void
GE_SegmentPolyhedron_INT:: declare_intersection_checked( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: declare_intersection_checked" ) ;
   PEL_CHECK_PRE( !intersection_checked() ) ;

   INTER_CHECKED = true ;
}


//----------------------------------------------------------------------
void
GE_SegmentPolyhedron_INT:: reset( GE_Point const* S0,
                                  GE_Point const* S1,
                                  GE_Mpolyhedron const* M )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: reset" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
   PEL_CHECK_PRE( S0 != 0 ) ;
   PEL_CHECK_PRE( S0->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( S1 != 0 ) ;
   PEL_CHECK_PRE( S1->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( M != 0 ) ;
   PEL_CHECK_PRE( M->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_PRE( M->dimension() == nb_space_dimensions()-1 ) ;

   INTER_CHECKED = false ;

   M_SAVE = M ;
   S0_SAVE = S0 ;
   S1_SAVE = S1 ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !intersection_checked() ) ;
   PEL_CHECK_POST( segment_first_vertex() == S0 ) ;
   PEL_CHECK_POST( segment_second_vertex() == S1 ) ;
   PEL_CHECK_POST( target_polyhedron() == M ) ;
}

//----------------------------------------------------------------------
GE_Point const* 
GE_SegmentPolyhedron_INT:: segment_first_vertex( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: segment_first_vertex" ) ;
   PEL_CHECK_PRE( intersection_checked() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Point const* result = S0_SAVE ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_coordinates() == nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const* 
GE_SegmentPolyhedron_INT:: segment_second_vertex( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: segment_second_vertex" ) ;
   PEL_CHECK_PRE( intersection_checked() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Point const* result = S1_SAVE ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_coordinates() == nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron const* 
GE_SegmentPolyhedron_INT:: target_polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: target_polyhedron" ) ;
   PEL_CHECK_PRE( intersection_checked() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Mpolyhedron const* result = M_SAVE ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->dimension() == nb_space_dimensions()-1 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_SegmentPolyhedron_INT:: print( std::ostream& os, size_t indent_width ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegmentPolyhedron_INT:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string const space( indent_width, ' ' ) ;
   os << space << " Intersection between a segment and a polyhedron: " ;
   if( !is_a_prototype() )
   {
      os << std::endl << "Segment : " << std::endl ;
      S0_SAVE->print( os, indent_width+3 ) ;
      S1_SAVE->print( os, indent_width+3 ) ;
      os << std::endl << "Target polyhedron : " << M_SAVE->name() << std::endl ;
      M_SAVE->print( os, indent_width+3 ) ;
   }
   else
   {
      os << "prototype" << std::endl ;
   }
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron_INT:: check_intersection_PRE(
                                GE_Point const* S0,
                                GE_Point const* S1,
                                GE_Mpolyhedron const* M ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !is_a_prototype() ) ;
   PEL_ASSERT( S0 != 0 ) ;
   PEL_ASSERT( S0->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_ASSERT( S1 != 0 ) ;
   PEL_ASSERT( S1->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_ASSERT( M != 0 ) ;
   PEL_ASSERT( M->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_ASSERT( M->dimension() == nb_space_dimensions()-1 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron_INT:: check_intersection_POST(
                                GE_Point const* S0,
                                GE_Point const* S1,
                                GE_Mpolyhedron const* M ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( intersection_checked() ) ;
   PEL_ASSERT( segment_first_vertex() == S0 ) ;
   PEL_ASSERT( segment_second_vertex() == S1 ) ;
   PEL_ASSERT( target_polyhedron() == M ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron_INT:: intersection_point_PRE( GE_Point const* pt ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( intersection_checked() ) ; 
   PEL_ASSERT( one_single_intersection() ) ;
   PEL_ASSERT( pt != 0 ) ;
   PEL_ASSERT( pt->nb_coordinates() == nb_space_dimensions() ) ;
 
   return( true ) ;
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron_INT:: create_replica_PRE(
                              PEL_Object const* a_owner,
                              PEL_ModuleExplorer const* a_mod_exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool 
GE_SegmentPolyhedron_INT:: create_replica_POST(
                              GE_SegmentPolyhedron_INT const* result,
                              PEL_Object const* a_owner,
                              PEL_ModuleExplorer const* a_mod_exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
GE_SegmentPolyhedron_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_SegmentPolyhedron_INT:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "GE_SegmentPolyhedron_INT descendant" ) ;
   return( result ) ;
}
