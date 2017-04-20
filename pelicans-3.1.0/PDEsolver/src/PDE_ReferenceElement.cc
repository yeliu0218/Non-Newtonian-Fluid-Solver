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

#include <PDE_ReferenceElement.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL_Error.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <numeric>

using std::string ;

size_t PDE_ReferenceElement::NB_INSTANCES = 0 ;

//----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_ReferenceElement:: object( std::string a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement:: object" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PDE_ReferenceElement const* result =
      static_cast<PDE_ReferenceElement const*>(
                                    plugins_map()->item( a_name ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElement:: PDE_ReferenceElement(
                              std::string a_name,
                              GE_ReferencePolyhedron const* a_ref_poly )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , ONAME( a_name )
   , ID( NB_INSTANCES )
   , NB_BFS( 0 )
   , NODE_PTS( PEL_Vector::create( this, 0 ) )
   , REF_POLY( a_ref_poly )
{
   PEL_LABEL( "PDE_ReferenceElement:: PDE_ReferenceElement" ) ;

   plugins_map()->register_item( a_name, this ) ;
   ++NB_INSTANCES ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElement:: ~PDE_ReferenceElement( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_ReferenceElement:: object_with_nodes_at_vertices(
                                            GE_Mpolyhedron const* poly )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement:: object_with_nodes_at_vertices" ) ;
   PEL_CHECK_PRE( poly != 0 ) ;

   string element_name ;

   string const& name = poly->reference_polyhedron()->name() ;
   if( name=="GE_ReferenceSegment" )
   {
      element_name = "PDE_1D_P1_2nodes" ;
   }
   else if( name=="GE_ReferenceTriangle" )
   {
      element_name = "PDE_2D_P1_3nodes" ;
   }
   else if( name=="GE_ReferenceSquare" )
   {
      element_name = "PDE_2D_Q1_4nodes" ;
   }
   else if( name=="GE_ReferenceCube" )
   {
      element_name = "PDE_3D_Q1_8nodes" ;
   }
   else if( name=="GE_ReferenceTetrahedron" )
   {
      element_name = "PDE_3D_P1_4nodes" ;
   }
   else if( name=="GE_ReferencePoint" )
   {
      element_name = "PDE_0D_Q0_1node" ;
   }
   else
   {
      string mesg = name + " : unknown reference polyhedron" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
   return( PDE_ReferenceElement::object( element_name ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ReferenceElement:: nb_objects( void )
//----------------------------------------------------------------------
{
   return( NB_INSTANCES ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ReferenceElement:: id_number( void ) const
//----------------------------------------------------------------------
{
   return( ID ) ;
}

//----------------------------------------------------------------------
std::string const&
PDE_ReferenceElement:: name(void ) const
//----------------------------------------------------------------------
{
   return ONAME ;
}

//----------------------------------------------------------------------
size_t
PDE_ReferenceElement:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   
   return( REF_POLY->dimension() ) ;
}

//----------------------------------------------------------------------
GE_ReferencePolyhedron const*
PDE_ReferenceElement:: reference_polyhedron( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_ReferencePolyhedron const* result = REF_POLY ;

   PEL_CHECK_POST( reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ReferenceElement:: nb_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   
   return( NB_BFS ) ;
}

//----------------------------------------------------------------------
bool
PDE_ReferenceElement:: has_node( GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = false ;
   for( size_t iBF=0 ; iBF<NB_BFS ; iBF++ )
   {
      if( pt_ref->distance( node_location(iBF) ) <1.E-4 )
      {
         result = true ;
         break ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ReferenceElement:: local_node( GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( has_node( pt_ref ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   size_t result = NB_BFS ;

   for( size_t iBF=0 ; iBF<NB_BFS ; iBF++ )
   {
      if( pt_ref->distance( node_location(iBF) ) < 1.e-4 )
      {
         result = iBF ;
         break ;
      }
   }
   PEL_CHECK_POST( node_location( result )->distance( pt_ref ) < 1.e-4 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_ReferenceElement:: node_location( size_t node ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( node < nb_nodes() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Point const* result = static_cast<GE_Point const*>( NODE_PTS->at( node ) ) ;
   PEL_CHECK( dynamic_cast<GE_Point*>( NODE_PTS->at( node ) ) != 0 ) ;

   PEL_CHECK_POST( local_node( result ) == node ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PDE_ReferenceElement:: comparable( PEL_Object const* other ) const
//---------------------------------------------------------------------------
{
   // less restrictive than PEL_Object::comparable
   PEL_ASSERT( dynamic_cast<PDE_ReferenceElement const*>(other) != 0 ) ;

   return( true ) ;
}

//-------------------------------------------------------------------------
bool
PDE_ReferenceElement:: is_equal( PEL_Object const* other ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = ( three_way_comparison( other ) == 0 ) ;

   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
PDE_ReferenceElement:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   int result = 0 ;

   PDE_ReferenceElement const* otherElt =
                         static_cast< PDE_ReferenceElement const*> ( other ) ;

   std::string const& selfName = name() ;
   std::string const& otherName = otherElt->name() ;

   if( selfName < otherName )
   {
      result = - 1 ;
   }
   else if( selfName > otherName )
   {
      result = 1 ;
   }

   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_ReferenceElement:: hash_code( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement:: hash_code" ) ;
   std::string const& selfName = name() ;
   size_t result = std::accumulate( selfName.begin(), selfName.end(), 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_ReferenceElement:: append_node( GE_Point* node_point )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( node_point !=0 ) ;
   PEL_CHECK_PRE( node_point->nb_coordinates() == dimension() ) ;
   PEL_CHECK_PRE( reference_polyhedron()->contains( node_point ) ) ;
   PEL_CHECK_PRE( !has_node( node_point ) ) ;
   PEL_CHECK_PRE( node_point->owner() == this ) ;
   PEL_CHECK_INV( invariant() ) ;

   NODE_PTS->append( node_point ) ;
   NB_BFS++ ;

   PEL_CHECK_INV( invariant()) ;
}

//----------------------------------------------------------------------
bool
PDE_ReferenceElement:: N_local_PRE( size_t node,
                                    GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( node < nb_nodes() ) ;
   PEL_ASSERT( pt_ref!=0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( reference_polyhedron()->contains( pt_ref ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_ReferenceElement:: dN_local_PRE( size_t node,
                                     size_t a,
                                     GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( node < nb_nodes() ) ;
   PEL_ASSERT( a < dimension() ) ;
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( reference_polyhedron()->contains( pt_ref ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_ReferenceElement:: d2N_local_PRE( size_t node,
                                      size_t a,
                                      size_t b,
                                      GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( node < nb_nodes() ) ;
   PEL_ASSERT( a < dimension() ) ;
   PEL_ASSERT( b < dimension() ) ;
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( pt_ref->nb_coordinates() == dimension() ) ;
   PEL_ASSERT( reference_polyhedron()->contains( pt_ref ) ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------
bool
PDE_ReferenceElement:: reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_ReferenceElement:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( NB_BFS==NODE_PTS->count() ) ;
   PEL_ASSERT( REF_POLY!=0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PDE_ReferenceElement:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "PDE_ReferenceElement descendant" ) ;
   return( result ) ;
}
