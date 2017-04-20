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

#include <PDE_BasisFunctionCell.hh>

#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_FaceFE.hh>
#include <PDE_MortarSideFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_RefinementPatternProvider.hh>
#include <PDE_ReferenceElementRefiner.hh>

#include <algorithm>
#include <iostream>
#include <list>
#include <string>

using std::cout ; using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_BasisFunctionCell:: create( PEL_Object* a_owner )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: create" ) ;

   PDE_BasisFunctionCell* result = new PDE_BasisFunctionCell( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_active() ) ;
   PEL_CHECK_POST( !result->is_refined() ) ;
   PEL_CHECK_POST( !result->valid_field() ) ;
   PEL_CHECK_POST( result->refinement_level() == PEL::bad_index() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell:: PDE_BasisFunctionCell( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PDE_BasisFunction( a_owner )
   , LEADING_PARENT( PEL::bad_index() )
   , REF_ELT_INDEX( PEL::bad_index() )
{
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell:: ~PDE_BasisFunctionCell( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: refinement_level( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: refinement_level" ) ;

   size_t result = ( PIECES.empty() ? 
                        PEL::bad_index() : 
                     PIECES[ 0 ].cell->refinement_level() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: nb_parents( void ) const
//-------------------------------------------------------------------------
{
   size_t result = PARENTS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_BasisFunctionCell:: parent( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: parent" ) ;
   PEL_CHECK_PRE( i < nb_parents() ) ;

   PDE_BasisFunctionCell* result = PARENTS[ i ] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_parent_of( PDE_BasisFunctionCell* a_child ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_parent_of" ) ;
   PEL_CHECK_PRE( a_child != 0 ) ;

   bool result = ( std::find( CHILDS.begin(), CHILDS.end(), a_child ) !=
                              CHILDS.end() ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: parent_is_leading( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: child_nodes_is_leading" ) ;
   PEL_CHECK_PRE( i < nb_parents() ) ;

   bool result = ( LEADING_PARENT == i ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: set_child_parent_relationship(
                                              PDE_BasisFunctionCell* a_parent,
                                              double refinement_coef,
                                              bool is_leading )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: set_child_parent_relationship" ) ;
   PEL_CHECK_PRE( a_parent != 0 )
   PEL_CHECK_PRE( refinement_coef != 0 )
   PEL_CHECK_PRE( !is_leading || refinement_coef == 1.0 )

   if( !is_child_of( a_parent ) )
   {
      PARENTS.push_back( a_parent ) ;
      if( is_leading )
      {
         PEL_ASSERT( PARENTS.size() == 1 ) ;
         LEADING_PARENT = 0 ;
      }
   }

   std::vector< PDE_BasisFunctionCell* >& childs = a_parent->CHILDS ;
   std::vector< double >& coefs = a_parent->REFI_COEFS ;
   if( !a_parent->is_parent_of( this ) )
   {
      childs.push_back( this ) ;
      coefs.push_back( refinement_coef ) ;
      if( is_leading )
      {
         PEL_ASSERT( refinement_coef == 1.0 ) ;
      }
   }
   else
   {
      //??? a supprimer en optimisé
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

   PEL_CHECK_POST( FORMAL( std::count( PARENTS.begin(),
                                       PARENTS.end(), a_parent ) == 1 ) ) ;
   PEL_CHECK_POST( FORMAL( std::count( childs.begin(),
                                       childs.end(), this ) == 1 ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: set_all_child_parent_relationship_on_cell(
                                  PDE_CellFE const* bfcell, size_t ln,
                                  size_t ee, PDE_CellFE const* pcell,
                                  size_t ic, size_t verb_level )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: set_all_child_parent_relationship_on_cell" ) ;
   PEL_CHECK_PRE( this == bfcell->basis_function( ee, ln ) )
   PEL_CHECK_PRE( pcell == bfcell->parent() )
   PEL_CHECK_PRE( bfcell == pcell->child( ic ) )

   std::string indent = "                " ;
   if( verb_level > 2 ) PEL::out() << indent
                                   << "Set parent/child relationship" << endl ;

   PDE_RefinementPatternProvider const* rpp =
                                     PDE_MeshFE::refinement_pattern_provider() ;
   PDE_ReferenceElement const* elm = bfcell->reference_element( ee ) ;
   PDE_ReferenceElementRefiner const* elrf =
                                     rpp->reference_element_refiner( elm ) ;
   //Setting with existing parents
   if( verb_level > 2 ) PEL::out() << indent << "  "
                                   << "with parents" << endl ;
   for( size_t p=0 ; p<elrf->nb_parents( ln, ic ) ; ++p )
   {
      size_t pn = elrf->parent_node( ln, ic, p ) ;
      double xx = elrf->refinement_coef( ln, pn, ic ) ;
      size_t lcn = elrf->leading_child_node( pn, ic ) ;
      bool is_leading = ( lcn == ln ) ;
      PDE_BasisFunctionCell* pbf = pcell->basis_function( ee, pn ) ;
      if( pbf != 0 )
      {
         set_child_parent_relationship( pbf, xx, is_leading ) ;
         if( verb_level > 2 ) PEL::out() << indent << "    "
                                         << "with bf" << ln
                                         << " of cell " << pcell->id_number()
                                         << " : is_leading = " << is_leading
                                         << "  coeff = " << xx << endl ;
      }
   }

   //Setting with existing children
   if( verb_level > 2 ) PEL::out() << indent << "  "
                                   << "with children" << endl ;
   for( size_t icc=0 ; icc<bfcell->nb_childs() ; ++icc )
   {
      PDE_CellFE* ccell = bfcell->child( icc ) ;
      for( size_t c=0 ; c<elrf->nb_childs( ln, icc ) ; ++c)
      {
         size_t cn = elrf->child_node( ln, icc, c ) ;
         double xx = elrf->refinement_coef( cn, ln, icc ) ;
         size_t lcn =  elrf->leading_child_node( ln, icc );
         bool is_leading = ( lcn == cn ) ;
         PDE_BasisFunctionCell* cbf = ccell->basis_function( ee, cn ) ;
         if( cbf != 0 )
         {
            cbf->set_child_parent_relationship( this, xx, is_leading ) ;
            if( verb_level > 2 ) PEL::out() << indent << "    "
                                            << "with bf" << cn
                                            << " of cell " << ccell->id_number()
                                            << " : is_leading = " << is_leading
                                            << "  coeff = " << xx << endl ;
         }
      }
   }
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: nb_childs( void ) const
//-------------------------------------------------------------------------
{
   size_t result = CHILDS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_BasisFunctionCell:: child( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: child" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   PDE_BasisFunctionCell* result = CHILDS[ i ] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_child_of( PDE_BasisFunctionCell const* a_parent ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_child_of" ) ;
   PEL_CHECK_PRE( a_parent != 0 ) ;

   bool result = ( std::find( PARENTS.begin(), PARENTS.end(), a_parent ) !=
                              PARENTS.end() ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
PDE_BasisFunctionCell:: refinement_coefficient( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: refinement_coefficient" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   double result = REFI_COEFS[ i ] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: check_child_parent_relationship( void ) const
//-------------------------------------------------------------------------
{
   PDE_RefinementPatternProvider const* rpp =
                           PDE_MeshFE::refinement_pattern_provider() ;

   bool result = true ;

   //Parcours des cellules du support de 'self'
   for( size_t ic = 0 ; ic<nb_cells() ; ic++ )
   {
      PDE_CellFE* bfcell = cell( ic ) ;
      size_t ee = element_index_of_cell( ic ) ;
      PDE_ReferenceElement const* elm = bfcell->reference_element( ee ) ;
      size_t ln = local_node_of_cell( ic ) ;
      PDE_ReferenceElementRefiner const* elrf =
                                        rpp->reference_element_refiner( elm ) ;

      //Parcours des parents existants
      PDE_CellFE* pcell = bfcell->parent() ;
      //Recherche du numero local de la cellule fille bfcell
      size_t lcc = 0 ;
      while( lcc < pcell->nb_childs( ) && !( pcell->child( lcc ) == bfcell ) )
      {
         lcc++ ;
      }
      PEL_ASSERT( bfcell == pcell->child( lcc ) ) ;

      for( size_t p=0 ; p<elrf->nb_parents( ln, lcc ) && result ; ++p )
      {
         size_t pn = elrf->parent_node( ln, lcc, p ) ;
         PDE_BasisFunctionCell* pbf = pcell->basis_function( ee, pn ) ;
         if( pbf != 0 ) result = is_child_of( pbf ) ;
      }
      if( result )
      {
         //Parcours des enfants existants
         for( size_t icc=0 ; icc<bfcell->nb_childs() ; icc++)
         {
            PDE_CellFE* ccell = bfcell->child( icc ) ;
            for( size_t c=0 ; c<elrf->nb_childs( ln, icc ) && result  ; ++c )
            {
               size_t cn = elrf->child_node( ln, icc, c ) ;
               PDE_BasisFunctionCell* cbf = ccell->basis_function( ee, cn ) ;
               if( cbf != 0 ) result = is_parent_of( cbf ) ;

            }
         }
      }
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: nb_ascendants( void ) const
//-------------------------------------------------------------------------
{
   size_t result = ASCENDS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_BasisFunctionCell:: ascendant( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: ascendant" ) ;
   PEL_CHECK_PRE( i < nb_ascendants() ) ;

   PDE_BasisFunctionCell* result = ASCENDS[ i ] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_ascendant_of(
                                  PDE_BasisFunctionCell const* a_fine ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_ascendant_of" ) ;
   PEL_CHECK_PRE( a_fine!= 0 ) ;

   bool result = ( std::find( DESCENDS.begin(), DESCENDS.end(), a_fine ) !=
                              DESCENDS.end() ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: nb_descendants( void ) const
//-------------------------------------------------------------------------
{
   size_t result = DESCENDS.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_BasisFunctionCell:: descendant( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: descendant" ) ;
   PEL_CHECK_PRE( i < nb_descendants() ) ;

   PDE_BasisFunctionCell* result = DESCENDS[ i ] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_descendant_of(
                                  PDE_BasisFunctionCell const* a_coarse ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_descendant_of" ) ;
   PEL_CHECK_PRE( a_coarse != 0 ) ;

   bool result = ( std::find( ASCENDS.begin(), ASCENDS.end(), a_coarse ) !=
                              ASCENDS.end() ) ;

   return( result ) ;
}


//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: set_ascendant_relationship(
                                            PDE_BasisFunctionCell* a_coarse )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: set_ascendant_relationship" ) ;

   if( !is_descendant_of( a_coarse ) )
   {
      ASCENDS.push_back( a_coarse ) ;
   }

   std::vector< PDE_BasisFunctionCell* >& descs = a_coarse->DESCENDS ;
   if( !a_coarse->is_ascendant_of( this ) )
   {
     descs.push_back( this ) ;
   }

   PEL_CHECK_POST( FORMAL( std::count( ASCENDS.begin(),
                                       ASCENDS.end(), a_coarse ) == 1 ) ) ;
   PEL_CHECK_POST( FORMAL( std::count( descs.begin(),
                                       descs.end(), this ) == 1 ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: set_all_ascendant_relationship_on_cell(
                                            PDE_CellFE const* bfcell, size_t ln,
                                            size_t ee, PDE_CellFE const* pcell,
                                            size_t ic, size_t verb_level )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: set_all_ascendant_relationship_on_cell" ) ;
   PEL_CHECK_PRE( this == bfcell->basis_function( ee, ln ) )
   PEL_CHECK_PRE( pcell == bfcell->parent() )
   PEL_CHECK_PRE( bfcell == pcell->child( ic ) )

   std::string indent = "                " ;
   if( verb_level > 2 ) PEL::out() << indent
                                   << "Set ref deps relationship" << endl ;

   PDE_RefinementPatternProvider const* rpp =
                                     PDE_MeshFE::refinement_pattern_provider() ;
   PDE_ReferenceElement const* elm = bfcell->reference_element( ee ) ;
   PDE_ReferenceElementRefiner const* elrf =
                                     rpp->reference_element_refiner( elm ) ;
   //Setting with coarser basis functions
   if( verb_level > 2 ) PEL::out() << indent << "  "
                                   << "with coarser bf" << endl ;
   for( size_t p=0 ; p<elrf->nb_ascendants( ln, ic ) ; ++p )
   {
      size_t pn = elrf->ascendant_node( ln, ic, p ) ;
      PDE_BasisFunctionCell* pbf = pcell->basis_function( ee, pn ) ;
      if( pbf != 0 )
      {
         set_ascendant_relationship( pbf ) ;
         if( verb_level > 2 ) PEL::out() << indent << "    "
                                         << "with bf" << pn
                                         << " of cell " << pcell->id_number()
                                         << endl ;
      }
   }
   //Setting with finer basis functions
   if( verb_level > 2 ) PEL::out() << indent << "  "
                                   << "with finer bf" << endl ;
   for( size_t icc=0 ; icc<bfcell->nb_childs() ; ++icc )
   {
      PDE_CellFE* ccell = bfcell->child( icc ) ;
      for( size_t c=0 ; c<elrf->nb_descendants( ln, icc ) ; ++c)
      {
         size_t cn = elrf->descendant_node( ln, icc, c ) ;
         PDE_BasisFunctionCell* cbf = ccell->basis_function( ee, cn ) ;
         if( cbf != 0 )
         {
            cbf->set_ascendant_relationship( this ) ;
            if( verb_level > 2 ) PEL::out() << indent << "    "
                                            << "with bf" << cn
                                            << " of cell " << ccell->id_number()
                                            << endl ;
         }
      }
   }

}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: extend_pieces( PDE_CellFE* a_cell,
                                       size_t elm_index,
                                       size_t node_in_elm )
//----------------------------- --------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: extend_pieces" ) ;
   PEL_CHECK_PRE( a_cell != 0 ) ;

   PEL_ASSERT( refinement_level() == PEL::bad_index() ||
               refinement_level() == a_cell->refinement_level() ) ;

   std::vector< Piece >::iterator where = PIECES.begin() ;
   for( ; where != PIECES.end() ; ++where )
   {
      if( where->cell == a_cell )
      {
         PEL_ASSERT( where->elm_index == elm_index ) ;
         PEL_ASSERT( where->node_in_elm == node_in_elm ) ;
         break ;
      }
   }

   if( where == PIECES.end() )
   {
      PIECES.push_back( Piece( a_cell, elm_index, node_in_elm ) ) ;
   }
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: remove_from_pieces( PDE_CellFE* a_cell,
                                            size_t elm_index,
                                            size_t node_in_elm )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: remove_from_pieces" ) ;
   PEL_CHECK_PRE( a_cell != 0 ) ;

   std::vector< Piece >::iterator where = PIECES.begin() ;
   for( ; where != PIECES.end() ; ++where )
   {
      if( where->cell == a_cell ) break ;
   }
   
   //??? these conditions should be visible to the client.
   //??? remove_from_pieces can only be called with an argument list that was
   //???    previously used in a call to extend_pieces
   PEL_ASSERT( where != PIECES.end() ) ;
   PEL_ASSERT( where->cell == a_cell ) ;
   PEL_ASSERT( where->elm_index == elm_index ) ;
   PEL_ASSERT( where->node_in_elm == node_in_elm ) ;
   
   //??? is it optimal ? the same piece may be added, removed, added...
   PIECES.erase( where ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: nb_cells( void ) const
//-------------------------------------------------------------------------
{
   size_t result = PIECES.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_CellFE*
PDE_BasisFunctionCell:: cell( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: cell" ) ;
   PEL_CHECK_PRE( i<nb_cells() ) ;

   PDE_CellFE* result = PIECES[i].cell ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: element_index_of_cell( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: element_index_of_cell" ) ;
   PEL_CHECK_PRE( i<nb_cells() ) ;

   size_t result = PIECES[i].elm_index ;

   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: local_node_of_cell( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: local_node_of_cell" ) ;
   PEL_CHECK_PRE( i<nb_cells() ) ;

   size_t result = PIECES[i].node_in_elm ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: geometrical_node( GE_Point* pt ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: geometrical_node" ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   
   PDE_CellFE* ccell = cell( 0 ) ;
   size_t ln = local_node_of_cell( 0 ) ;
   size_t ee = element_index_of_cell( 0 ) ;
  
   PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
   GE_Point const* pt_ref = elm->node_location( ln ) ;
   ccell->polyhedron()->apply_mapping( pt_ref, pt ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: ref_elts_grp_index( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: ref_elts_grp_index" ) ;
   
   size_t result = REF_ELT_INDEX ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: set_ref_elts_grp_index( size_t ind )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: set_ref_elts_grp_index" ) ;
   
   REF_ELT_INDEX = ind ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_located_in_cell( PDE_CellFE const* a_cell ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_located_in_cell" ) ;
   PEL_CHECK_PRE( is_located_in_cell_PRE( a_cell ) ) ;

   // ---> on pourrait mettre en place un système de "lazy evaluation"
   // ---> cette fonction membre risque d'etre appelée souvent,
   // ---> toujours dans les memes conditions

   bool result = check_location_in_cell( a_cell ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: is_located_on_bound( PDE_BoundFE const* a_bound ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: is_located_on_bound" ) ;
   PEL_CHECK_PRE( is_located_on_bound_PRE( a_bound ) ) ;

   // ---> on pourrait mettre en place un système de "lazy evaluation"
   // ---> cette fonction membre risque d'etre appelée souvent,
   // ---> toujours dans les memes conditions

   bool result = check_location_on_bound( a_bound ) ;
   if( !result && a_bound->has_adjacent_bound() )
   {
      result = check_location_on_bound( a_bound->adjacent_bound() ) ;
   }

   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: print" ) ;

   PDE_BasisFunction::print( os, indent_width ) ;

   string bl( indent_width, ' ' ) ;
   os << bl << "nb_cells = " << nb_cells() << endl ;
   std::vector< Piece >::const_iterator where = PIECES.begin() ;
   for( ; where != PIECES.end() ; ++where )
   {
      PDE_CellFE const* cc = where->cell ;
      os << bl << "local node " << where->node_in_elm
         << " of " << cc->reference_element( where->elm_index )->name()
         << " on mesh : " << endl ;
      cc->print( os, indent_width+3 ) ;
   }
}

//-------------------------------------------------------------------------
void
PDE_BasisFunctionCell:: print_2( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: print" ) ;

   string bl( indent_width, ' ' ) ;

   size_t i=0 ;

   os << bl << "loc bf " << PIECES[i].node_in_elm
      << " of cell " << PIECES[i].cell->id_number() << " "
      << PIECES[i].cell->reference_element( PIECES[i].elm_index )->name()
      << " (l=" << refinement_level() << ")" ;
   if( nb_cells() != 1 )
   {
      os << " (other cells:" ;
      for( i=1 ; i<nb_cells() ; ++i )
      {
         os << " " << PIECES[i].cell->id_number() ;
      }
      os << ")" ;
   }
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: check_location_in_cell(
                                          PDE_CellFE const* a_cell ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: check_location_in_cell" ) ;

   bool result = false ;

   PDE_CellFE const* cc = a_cell ;
   std::list< size_t > child_idx ;
   for( size_t i=refinement_level() ; i != a_cell->refinement_level() ; ++i )
   {
      PDE_CellFE const* cc_parent = cc->parent() ;
      size_t ichild = 0 ;
      bool ok = false ;
      for( ; ichild != cc_parent->nb_childs() ; ++ichild )
      {
         if( cc_parent->child( ichild ) == cc )
         {
            ok = true ;
            child_idx.push_back( ichild ) ;
            break ;
         }
      }
      PEL_ASSERT( ok ) ;
      cc = cc_parent ;
   }
   PEL_ASSERT( cc->refinement_level() == refinement_level() );

   std::vector< Piece >::const_iterator where = PIECES.begin() ;
   for( ; where != PIECES.end() ; ++where )
   {
      if( where->cell == cc ) break ;
   }

   if( where != PIECES.end() )
   {
      if( a_cell->refinement_level() == refinement_level() )
      {
         result = true ;
      }
      else
      {
         PDE_CellFE const* ccell = where->cell ;
         GE_ReferencePolyhedronRefiner const* cref =
             PDE_MeshFE::refinement_pattern_provider()->cell_refiner( ccell ) ;
         PDE_ReferenceElement const* elm =
                           ccell->reference_element( where->elm_index ) ;
         size_t inode = where->node_in_elm ;
         GE_Point const* pt_in_ccell = elm->node_location( inode ) ;
         GE_Point* pt = tmp_point( pt_in_ccell->nb_coordinates() ) ;
         pt->set( pt_in_ccell ) ;
         bool ok = true ;
         for( size_t i=refinement_level() ;
              ok && i!=a_cell->refinement_level() ; ++i )
         {
            size_t ichild = child_idx.back() ;
            child_idx.pop_back() ;
            cref->compute_location_in_subcell( ichild, pt, pt ) ;
            ok = cref->subcell_reference_polyhedron()->contains( pt ) ;
         }
         PEL_ASSERT( IMPLIES( ok, child_idx.empty() ) ) ;
         result = ok ;
      }
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunctionCell:: check_location_on_bound(
                                         PDE_BoundFE const* a_bound ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: check_location_on_bound" ) ;

   bool result = false ;

   PEL_ASSERT( a_bound->has_adjacent_cell() ) ;
   PDE_CellFE const* a_cell = a_bound->adjacent_cell() ;
   size_t i_face = face_index( a_bound, a_cell ) ;

   PDE_BoundFE const* bd = a_bound ;
   PDE_CellFE const*  cc = a_cell  ;

   std::list< size_t > child_idx ;
   std::list< size_t > face_idx ;
   for( size_t i=refinement_level() ; i != a_cell->refinement_level() ; ++i )
   {
      PDE_CellFE const*  cc_parent = cc->parent() ;
      PDE_BoundFE const* bd_parent = bd->parent() ;
      size_t ichild = 0 ;
      bool ok = false ;
      for( ; ichild != cc_parent->nb_childs() ; ++ichild )
      {
         if( cc_parent->child( ichild ) == cc )
         {
            ok = true ;
            child_idx.push_back( ichild ) ;
            break ;
         }
      }
      PEL_ASSERT( ok ) ;
      cc = cc_parent  ;
      bd = bd_parent ;
      face_idx.push_back( face_index( bd, cc ) ) ;
   }
   PEL_ASSERT( cc->refinement_level() == refinement_level() ) ;

   std::vector< Piece >::const_iterator where = PIECES.begin() ;
   for( ; where != PIECES.end() ; ++where )
   {
      if( where->cell == cc ) break ;
   }

   if( where != PIECES.end() )
   {
      PDE_CellFE const* ccell = where->cell ;
      PDE_ReferenceElement const* elm =
                        ccell->reference_element( where->elm_index ) ;
      size_t inode = where->node_in_elm ;
      GE_Point const* pt_in_ccell = elm->node_location( inode ) ;
      GE_ReferencePolyhedron const* poly =
                             a_cell->polyhedron()->reference_polyhedron() ;

      if( a_cell->refinement_level() == refinement_level() )
      {
         result = poly->face_contains( i_face, pt_in_ccell ) ;
      }
      else
      {
         GE_ReferencePolyhedronRefiner const* cref =
             PDE_MeshFE::refinement_pattern_provider()->cell_refiner( ccell ) ;
         GE_Point* pt = tmp_point( pt_in_ccell->nb_coordinates() ) ;
         pt->set( pt_in_ccell ) ;
         bool ok = true ;
         for( size_t i=refinement_level() ;
              ok && i!=a_cell->refinement_level() ; ++i )
         {
            size_t ichild = child_idx.back() ;
            child_idx.pop_back() ;
            i_face = face_idx.back() ;
            face_idx.pop_back() ;
            cref->compute_location_in_subcell( ichild, pt, pt ) ;
            ok = poly->face_contains( i_face, pt ) ;
         }
         PEL_ASSERT( IMPLIES( ok, child_idx.empty() && face_idx.empty() ) ) ;
         result = ok ;
      }
   }
   return( result ) ;
}

// ---> pourrait etre stocké dans PDE_BoundFE
//-------------------------------------------------------------------------
size_t
PDE_BasisFunctionCell:: face_index( PDE_BoundFE const* bound,
                                    PDE_CellFE const* a_cell )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunctionCell:: face_index" ) ;

   PEL_ASSERT( bound->refinement_level() == a_cell->refinement_level() ) ;

   size_t result = PEL::bad_index() ;

   size_t nbf = a_cell->polyhedron()->nb_faces() ;
   PEL_Vector const* faces = a_cell->faces() ;
   for( size_t i=0 ; i<nbf ; ++i )
   {
      PDE_FaceFE const* face =
                           static_cast< PDE_FaceFE* >( faces->at( i ) ) ;
      PEL_ASSERT( face->refinement_level() == bound->refinement_level() ) ;
      if( face->has_adjacent_bound() && face->adjacent_bound()==bound )
      {
         result = i ;
         break ;
      }
   }

   PEL_ASSERT( result != PEL::bad_index() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Point*
PDE_BasisFunctionCell:: tmp_point( size_t dim )
//----------------------------------------------------------------------------
{
   static GE_Point* PT1 = GE_Point::create( PEL_Root::object(), (size_t)1 ) ;
   static GE_Point* PT2 = GE_Point::create( PEL_Root::object(), (size_t)2 ) ;
   static GE_Point* PT3 = GE_Point::create( PEL_Root::object(), (size_t)3 ) ;

   if( dim == 3 )
   {
      return PT3 ;
   }
   else if( dim == 2 )
   {
      return PT2 ;
   }
   else if( dim == 1 )
   {
      return PT1 ;
   }

   PEL_Error::object()->raise_plain( "PDE_BasisFunctionCell : dim="+dim ) ;
   return( 0 ) ;
}
