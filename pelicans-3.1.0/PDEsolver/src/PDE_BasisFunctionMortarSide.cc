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

#include <PDE_BasisFunctionMortarSide.hh>

#include <PEL_Error.hh>
#include <PEL_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_MortarSideFE.hh>
#include <PDE_FaceFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_RefinementPatternProvider.hh>

#include <algorithm>
#include <iostream>
#include <list>
#include <string>

using std::cout ; using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
PDE_BasisFunctionMortarSide*
PDE_BasisFunctionMortarSide:: create( PEL_Object* a_owner )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: create" ) ;

   PDE_BasisFunctionMortarSide* result = 
                                new PDE_BasisFunctionMortarSide( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_active() ) ;
   PEL_CHECK_POST( !result->is_refined() ) ;
   PEL_CHECK_POST( !result->valid_field() ) ;
   PEL_CHECK_POST( result->refinement_level() == PEL::bad_index() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionMortarSide:: PDE_BasisFunctionMortarSide( 
                                               PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PDE_BasisFunction( a_owner )
{
}

//-------------------------------------------------------------------------
PDE_BasisFunctionMortarSide:: ~PDE_BasisFunctionMortarSide( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: refinement_level( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: refinement_level" ) ;

   size_t result = 
          SIDES.empty() ? PEL::bad_index() : SIDES[0]->refinement_level() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: nb_parents( void ) const
//-------------------------------------------------------------------------
{
   size_t result = PARENTS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionMortarSide* 
PDE_BasisFunctionMortarSide:: parent( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: parent" ) ;
   PEL_CHECK_PRE( i < nb_parents() ) ;

   PDE_BasisFunctionMortarSide* result = PARENTS[ i ] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionMortarSide:: set_child_parent_relationship( 
                                        PDE_BasisFunctionMortarSide* a_parent,
                                        double refinement_coef )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: set_child_parent_relationship" ) ;

   bool found = ( std::find( PARENTS.begin(), PARENTS.end(), a_parent ) != 
                             PARENTS.end() ) ;
   if( !found )
   {
      PARENTS.push_back( a_parent ) ;
   }

   std::vector< PDE_BasisFunctionMortarSide* >& childs = a_parent->CHILDS ;
   std::vector< double >& coefs = a_parent->REFI_COEFS ;
   found = ( std::find( childs.begin(), childs.end(), this ) != 
                        childs.end() ) ;
   if( !found )
   {
      childs.push_back( this ) ;
      coefs.push_back( refinement_coef ) ;
   }
   else
   {
      bool ok = false ;
      for( size_t i=0 ; i<childs.size() ; ++i )
      {
         if( childs[i] == this )
         {
            PEL_ASSERT( coefs[i] == refinement_coef ) ;
            ok = true ;
         }
      }
      PEL_ASSERT( ok ) ;
   }
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: nb_childs( void ) const
//-------------------------------------------------------------------------
{
   size_t result = CHILDS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionMortarSide* 
PDE_BasisFunctionMortarSide:: child( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: child" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   PDE_BasisFunctionMortarSide* result = CHILDS[ i ] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
double
PDE_BasisFunctionMortarSide:: refinement_coefficient( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: refinement_coefficient" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   double result = REFI_COEFS[ i ] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionMortarSide:: extend_pieces( PDE_MortarSideFE* a_side,
                                   size_t elm_index,
                                   size_t node_in_elm )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: extend_pieces" ) ;
   PEL_CHECK_PRE( a_side != 0 ) ;

   PEL_ASSERT( refinement_level() == PEL::bad_index() ||
               refinement_level() == a_side->refinement_level() ) ;

   bool found = false ;
   for( size_t i=0 ; !found && i<SIDES.size() ; ++i )
   {
      if( SIDES[i] == a_side )
      {
         found = true ;
         PEL_ASSERT( ELM_IDXS[i]==elm_index && ELM_NODES[i]==node_in_elm ) ;
      }
   }

   if( !found )
   {
      SIDES.push_back( a_side ) ;
      ELM_IDXS.push_back( elm_index ) ;
      ELM_NODES.push_back( node_in_elm ) ;
   }
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: nb_sides( void ) const
//-------------------------------------------------------------------------
{
   size_t result = SIDES.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_MortarSideFE*
PDE_BasisFunctionMortarSide:: side( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: bound" ) ;
   PEL_CHECK_PRE( i<nb_sides() ) ;
   
   PDE_MortarSideFE* result = SIDES[i] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: element_index_of_side( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: element_index_of_bound" ) ;
   PEL_CHECK_PRE( i<nb_sides() ) ;
   
   size_t result = ELM_IDXS[i] ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionMortarSide:: local_node_of_side( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: local_node_of_bound" ) ;
   PEL_CHECK_PRE( i<nb_sides() ) ;
   
   size_t result = ELM_NODES[i] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionMortarSide:: is_located_in_cell(
                                            PDE_CellFE const* a_cell ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: is_located_in_cell" ) ;

   PEL_Error::object()->raise_not_implemented( this, "is_located_in_cell" ) ;

   if( refinement_level() != a_cell->refinement_level() )
   {
      PEL_Error::object()->raise_internal( "refinement not handled" ) ;
   }

   bool result = false ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionMortarSide:: is_located_on_bound( 
                                            PDE_BoundFE const* a_bound ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: is_located_on_bound" ) ;

   PEL_Error::object()->raise_not_implemented( this, "is_located_on_bound" ) ;

   if( refinement_level() != a_bound->refinement_level() )
   {
      PEL_Error::object()->raise_internal( "refinement not handled" ) ;
   }

   bool result = false ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionMortarSide:: print( std::ostream& os, 
                                     size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionMortarSide:: print" ) ;

   PDE_BasisFunction::print( os, indent_width ) ;

   string bl( indent_width, ' ' ) ;
   os << bl << "nb_sides = " << nb_sides() << endl ;
   for( size_t i=0 ; i<SIDES.size() ; ++i )
   {
      os << bl << "local node " << ELM_NODES[i] 
         << " of " << SIDES[i]->reference_element( ELM_IDXS[i] ) 
         << " on mesh : " << endl ; 
      SIDES[i]->print( os, indent_width+3 ) ;
   }
}
