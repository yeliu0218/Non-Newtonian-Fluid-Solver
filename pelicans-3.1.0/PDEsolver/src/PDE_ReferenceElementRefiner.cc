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

#include <PDE_ReferenceElementRefiner.hh>

#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>

#include <GE_ReferencePolyhedronRefiner.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_2D_P0_1node_RefinerA.hh>
#include <PDE_2D_P1_3nodes_RefinerA.hh>
#include <PDE_2D_P2_6nodes_RefinerA.hh>
#include <PDE_2D_Q0_1node_RefinerA.hh>
#include <PDE_2D_Q1_4nodes_RefinerA.hh>
#include <PDE_2D_Q1isoNonConfA_4nodes_RefinerA.hh>
#include <PDE_2D_Q2_9nodes_RefinerA.hh>
#include <PDE_3D_Q1_8nodes_RefinerA.hh>
#include <PDE_3D_Q2_27nodes_RefinerA.hh>

#include <iostream>
using std::cout ; using std::endl ;

PDE_ReferenceElementRefiner::OneLevelDiffRule
PDE_ReferenceElementRefiner::OLRULE = Invalid ;

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner const*
PDE_ReferenceElementRefiner:: object( std::string const& a_name )
//---------------------------------------------------------------------------
{
   PDE_ReferenceElementRefiner const* result = 0 ;
   if( a_name == "PDE_2D_Q0_1node_RefinerA" )
   {
      result = PDE_2D_Q0_1node_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_Q1_4nodes_RefinerA" )
   {
      result = PDE_2D_Q1_4nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_Q1isoNonConfA_4nodes_RefinerA" )
   {
      result = PDE_2D_Q1isoNonConfA_4nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_Q2_9nodes_RefinerA" )
   {
      result = PDE_2D_Q2_9nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_P0_1node_RefinerA" )
   {  
      result = PDE_2D_P0_1node_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_P1_3nodes_RefinerA" )
   {
      result = PDE_2D_P1_3nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_2D_P2_6nodes_RefinerA" )
   {
      result = PDE_2D_P2_6nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_3D_Q1_8nodes_RefinerA" )
   {
      result = PDE_3D_Q1_8nodes_RefinerA::object() ;
   }
   else if( a_name == "PDE_3D_Q2_27nodes_RefinerA" )
   {
      result = PDE_3D_Q2_27nodes_RefinerA::object() ;
   }
   PEL_ASSERT( result != 0 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner:: PDE_ReferenceElementRefiner(
                    std::string const& name_of_reference_element,
                    GE_ReferencePolyhedronRefiner const* rcell_refiner )
//---------------------------------------------------------------------------
   : PEL_Object( PEL_Root::object() )
   , CELL_REFINER( rcell_refiner )
   , ELM( PDE_ReferenceElement::object( name_of_reference_element ) )
   , NB_PARENTS( 0, 0 )
   , NB_CHILDS( 0, 0 )
   , PARENT( 0, 0, 0 )
   , CHILD( 0, 0, 0 )
   , LEADING_CHILD( 0, 0 )
   , REFI_COEF( 0, 0, 0 )
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: PDE_ReferenceElementRefiner" ) ;
   //???? tests de coherence CELL_REFINER et ELM

   NB_PARENTS.re_initialize( CELL_REFINER->nb_subcells(), ELM->nb_nodes() ) ;

   PARENT.re_initialize( CELL_REFINER->nb_subcells(),
                         ELM->nb_nodes(),
                         ELM->nb_nodes() ) ;
   PARENT.set( PEL::bad_index() ) ;

   NB_CHILDS.re_initialize( ELM->nb_nodes(), CELL_REFINER->nb_subcells() ) ;

   CHILD.re_initialize( ELM->nb_nodes(),
                        CELL_REFINER->nb_subcells(),
                        ELM->nb_nodes() ) ;
   CHILD.set( PEL::bad_index() ) ;

   LEADING_CHILD.re_initialize( ELM->nb_nodes(),
                                CELL_REFINER->nb_subcells() ) ;
   LEADING_CHILD.set( PEL::bad_index() ) ;

   REFI_COEF.re_initialize( CELL_REFINER->nb_subcells(),
                            ELM->nb_nodes(),
                            ELM->nb_nodes() ) ;
   REFI_COEF.set( PEL::bad_double() ) ;

   PEL_CHECK_POST( cell_refiner() == rcell_refiner ) ;
   PEL_CHECK_POST( reference_element() ==
                   PDE_ReferenceElement::object( name_of_reference_element ) ) ;
}

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner:: ~PDE_ReferenceElementRefiner( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElementRefiner:: append_child( size_t a_parent_node,
                                            size_t a_subcell,
                                            size_t a_child_node,
                                            double a_child_coef,
                                            bool is_leading_child )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: append_child" ) ;

   size_t icn = NB_CHILDS( a_parent_node, a_subcell ) ;
   CHILD( a_parent_node, a_subcell, icn ) = a_child_node ;
   ++NB_CHILDS( a_parent_node, a_subcell ) ;
   if( is_leading_child )
   {
      LEADING_CHILD( a_parent_node, a_subcell ) = a_child_node ;
   }

   size_t ipn = NB_PARENTS( a_subcell, a_child_node ) ;
   PARENT( a_subcell, a_child_node, ipn ) = a_parent_node ;
   ++NB_PARENTS( a_subcell, a_child_node ) ;

   REFI_COEF( a_subcell, a_child_node, a_parent_node ) = a_child_coef ;
}

//----------------------------------------------------------------------
void
PDE_ReferenceElementRefiner:: set_one_level_difference_rule( 
                                            OneLevelDiffRule a_rule )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: set_one_level_difference_rule" ) ;
   PEL_CHECK_PRE( a_rule != Invalid ) ;

   OLRULE = a_rule ;

   PEL_CHECK_POST( one_level_difference_rule() == a_rule ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElementRefiner::OneLevelDiffRule
PDE_ReferenceElementRefiner::one_level_difference_rule( void )
//----------------------------------------------------------------------
{
   return( OLRULE ) ;
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner const*
PDE_ReferenceElementRefiner:: cell_refiner( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: cell_refiner" ) ;

   GE_ReferencePolyhedronRefiner const* result = CELL_REFINER ;

   PEL_CHECK_POST( IMPLIES( result != 0,
                   result->subcell_reference_polyhedron() ==
                            reference_element()->reference_polyhedron() ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_ReferenceElementRefiner:: reference_element( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: reference_element" ) ;

   PDE_ReferenceElement const* result = ELM ;

   PEL_CHECK_POST( result->reference_polyhedron() ==
                   cell_refiner()->subcell_reference_polyhedron() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: nb_childs( size_t a_parent_node,
                                         size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: nb_childs" ) ;
   PEL_CHECK_PRE( a_parent_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   size_t result = NB_CHILDS( a_parent_node, ic ) ;

   PEL_CHECK_POST( result <= reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: child_node( size_t a_parent_node,
                                          size_t ic,
                                          size_t icn ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: child_node" ) ;
   PEL_CHECK_PRE( a_parent_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;
   PEL_CHECK_PRE( icn < nb_childs( a_parent_node, ic ) ) ;

   size_t result = CHILD( a_parent_node, ic, icn ) ;

   PEL_CHECK_POST( result < reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: leading_child_node( size_t a_parent_node,
                                                  size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: nb_childs" ) ;
   PEL_CHECK_PRE( a_parent_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   size_t result = LEADING_CHILD( a_parent_node, ic ) ;

   PEL_CHECK_POST( result==PEL::bad_index() ||
                   result <= reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: nb_parents( size_t a_child_node,
                                          size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: nb_parents" ) ;
   PEL_CHECK_PRE( a_child_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   size_t result = NB_PARENTS( ic, a_child_node ) ;

   PEL_CHECK_POST( result <= reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: parent_node( size_t a_child_node,
                                           size_t ic,
                                           size_t ipn ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: parent_node" ) ;
   PEL_CHECK_PRE( a_child_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;
   PEL_CHECK_PRE( ipn < nb_parents( a_child_node, ic ) ) ;

   size_t result = PARENT( ic, a_child_node, ipn ) ;

   PEL_CHECK_POST( result < reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: nb_ascendants( size_t a_fine_node,
                                             size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: nb_ascendants" ) ;
   PEL_CHECK_PRE( one_level_difference_rule() != Invalid ) ;
   PEL_CHECK_PRE( a_fine_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   size_t result = PEL::bad_index() ;
   if( OLRULE == Parents )
   {
      result = nb_parents( a_fine_node, ic ) ;
   }
   else if( OLRULE == Supports )
   {
      result = ELM->nb_nodes() ;
   }
   else if( OLRULE == No )
   {
      result = 0 ;
   }

   PEL_CHECK_POST( result != PEL::bad_index() ) ;
   PEL_CHECK_POST( 
      IMPLIES( ( one_level_difference_rule() == Parents ), 
               ( result == nb_parents( a_fine_node, ic ) ) ) ) ;
   PEL_CHECK_POST( 
      IMPLIES( ( one_level_difference_rule() == Supports ), 
               ( result == reference_element()->nb_nodes() ) ) ) ;
   PEL_CHECK_POST( 
      EQUIVALENT( ( one_level_difference_rule() == No ), 
                  ( result == 0 ) ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: ascendant_node( size_t a_fine_node,
                                              size_t ic,
                                              size_t iun ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: refinement_dep_node_down" ) ;
   PEL_CHECK_PRE( a_fine_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;
   PEL_CHECK_PRE( iun < nb_ascendants( a_fine_node, ic ) ) ;
   PEL_CHECK_PRE( OLRULE != No ) ;

   size_t result = PEL::bad_index() ;
   if( OLRULE == Parents )
   {
      result = parent_node( a_fine_node, ic, iun ) ;
   }
   else if( OLRULE == Supports )
   {
      result = iun ;
   }

   PEL_CHECK_POST( result != PEL::bad_index() ) ;
   PEL_CHECK_POST( result < reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: nb_descendants( size_t a_coarse_node,
                                              size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: nb_descendants" ) ;
   PEL_CHECK_PRE( a_coarse_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   size_t result = PEL::bad_index() ;
   if( OLRULE == Parents )
   {
      result = nb_childs( a_coarse_node, ic ) ;
   }
   else if( OLRULE == Supports )
   {
      result = ELM->nb_nodes() ;
   }
   else if( OLRULE == No )
   {
      result = 0 ;
   }

   PEL_CHECK_POST( result !=  PEL::bad_index() ) ;
   PEL_CHECK_POST( result <= reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_ReferenceElementRefiner:: descendant_node( size_t a_coarse_node,
                                               size_t ic,
                                               size_t iun ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: refinement_dep_node-up" ) ;
   PEL_CHECK_PRE( a_coarse_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;
   PEL_CHECK_PRE( iun < nb_descendants( a_coarse_node, ic ) ) ;
   PEL_CHECK_PRE( OLRULE != No ) ;

   size_t result = PEL::bad_index() ;
   if( OLRULE == Parents )
   {
    result = child_node( a_coarse_node, ic, iun ) ;
   }
   else if( OLRULE == Supports )
   {
     result = iun ;
   }

   PEL_CHECK_POST( result != PEL::bad_index() ) ;
   PEL_CHECK_POST( result < reference_element()->nb_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
PDE_ReferenceElementRefiner:: refinement_coef( size_t a_child_node,
                                               size_t a_parent_node,
                                               size_t ic ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElementRefiner:: refinement_coef" ) ;
   PEL_CHECK_PRE( a_child_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( a_parent_node < reference_element()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < cell_refiner()->nb_subcells() ) ;

   double result = REFI_COEF( ic, a_child_node, a_parent_node ) ;
   return( result ) ;
}
