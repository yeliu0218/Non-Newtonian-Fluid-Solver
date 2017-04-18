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

#include <PDE_AlgebraicCoarsener.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_GridFE.hh>
#include <PDE_SystemNumbering.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ;
using std::endl ;
using std::vector ;
using std::setw ;

//----------------------------------------------------------------------
PDE_AlgebraicCoarsener:: PDE_AlgebraicCoarsener( PEL_Object* a_owner,
                                                 PDE_GridFE* a_grid )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , GRID( a_grid )
   , NB_BLOCKS( 0 )
   , FINE_LEVEL( PEL::bad_index() )
   , LEVEL_MAX( PEL::bad_index() )
   , ROW_NMB( 0 )
   , COL_NMB( 0 )
{
}

//----------------------------------------------------------------------
PDE_AlgebraicCoarsener:: ~PDE_AlgebraicCoarsener( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: reset( void )
//----------------------------------------------------------------------
{
   LEVEL_MAX = PEL::bad_index() ;
   FINE_LEVEL = PEL::bad_index() ;

   PEL_CHECK_POST( !coarsening_is_possible() ) ;
   PEL_CHECK_POST( nb_levels() == PEL::bad_index() ) ;
   PEL_CHECK_POST( current_fine_level() == PEL::bad_index() ) ;
}

//-----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: prepare_for_coarsening(
                                          PDE_SystemNumbering const* nmb )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: prepare_for_coarsening" ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=0 ; i<nmb->nb_links() ; ++i ),
                  nmb->link( i ) != 0 ) ) ;
   PEL_CHECK_PRE( !coarsening_is_possible() ) ;

   if( ROW_NMB != 0 ) destroy_possession( ROW_NMB ) ;
   ROW_NMB = nmb->create_clone( this ) ;
   if( COL_NMB != 0 ) destroy_possession( COL_NMB ) ;
   COL_NMB = nmb->create_clone( this ) ;

   size_t old_nb_blocks = NB_BLOCKS ;
   NB_BLOCKS = ROW_NMB->nb_links() ;
   if( NB_BLOCKS > old_nb_blocks )
   {
      FF.resize( NB_BLOCKS ) ;
      BFS.resize( NB_BLOCKS ) ;
      NN_FINE.resize( NB_BLOCKS, boolVector( 0 ) ) ;
      NN_COAR.resize( NB_BLOCKS, boolVector( 0 ) ) ;
   }

   for( size_t i=0 ; i<NB_BLOCKS ; ++i )
   {
      PDE_LinkDOF2Unknown const* ll = ROW_NMB->link( i ) ;
      FF[ i ] = ll->field() ;

      size_t nb_nodes = FF[ i ]->nb_nodes() ;
      BFS[ i ].assign( nb_nodes, 0 ) ;
      NN_FINE[ i ].re_initialize( nb_nodes, false ) ;
      NN_COAR[ i ].re_initialize( nb_nodes, false ) ;
  }

   LEVEL_MAX = 0 ;

   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, GRID->cells() ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_CellFE const* cell = static_cast<PDE_CellFE*>( it->item() ) ;
      for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
      {
         PDE_DiscreteField const* ff = FF[ i_block ] ;
         std::vector< PDE_BasisFunctionCell* >& bfs = BFS[ i_block ] ;
         boolVector& nfine = NN_FINE[ i_block ] ;
         size_t ee = cell->index_of_reference_element( ff ) ;
         for( size_t i=0 ; i<cell->nb_basis_functions( ee ) ; ++i )
         {
            PDE_BasisFunctionCell* bf = cell->basis_function( ee, i ) ;
            if( bf != 0 )
            {
               size_t node = bf->node_of_DOF( ff ) ;
               if( bfs[node] == 0 )
               {
                  bfs[node] = bf ;
                  if( ff->node_is_active( node ) )
                  {
                     PEL_ASSERT( bf->is_active() ) ;
                     nfine( node ) = true ;
                     if( LEVEL_MAX<bf->refinement_level() )
                     {
                        LEVEL_MAX = bf->refinement_level() ;
                     }
                  }
               }
               else
               {
                  PEL_ASSERT( bfs[node] == bf ) ;
               }
            }
         }
      }
   }
   it->destroy() ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   LEVEL_MAX = com->max( LEVEL_MAX ) ;

   FINE_LEVEL = PEL::bad_index() ;

   PEL_CHECK_POST( nb_levels() != PEL::bad_index() ) ;
   PEL_CHECK_POST( IMPLIES( nb_levels()!=1, coarsening_is_possible() ) ) ;
   PEL_CHECK_POST( current_fine_level() == PEL::bad_index() ) ;
}

//-----------------------------------------------------------------------
size_t
PDE_AlgebraicCoarsener:: nb_levels( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: nb_levels" ) ;

   size_t result = ( LEVEL_MAX == PEL::bad_index() ?
                     PEL::bad_index() : LEVEL_MAX+1 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
bool
PDE_AlgebraicCoarsener:: coarsening_is_possible( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: coarsening_is_possible" ) ;

   bool result = ( LEVEL_MAX != PEL::bad_index() &&
                   LEVEL_MAX != 0 &&
                   ( FINE_LEVEL != 1 || FINE_LEVEL == PEL::bad_index() ) ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------
size_t
PDE_AlgebraicCoarsener:: current_fine_level( void ) const
//-----------------------------------------------------------------------
{
   return( FINE_LEVEL ) ;
}

//-----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: do_one_coarsening( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: do_one_coarsening" ) ;
   PEL_CHECK_PRE( coarsening_is_possible() ) ;
   PEL_SAVEOLD( size_t, current_fine_level, current_fine_level() ) ;

   if( FINE_LEVEL == PEL::bad_index() )
   {
      FINE_LEVEL = LEVEL_MAX ;
   }
   else
   {
      FINE_LEVEL-- ;

      PDE_SystemNumbering* dummy = ROW_NMB ;
      ROW_NMB = COL_NMB ;
      COL_NMB = dummy ;
      for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
      {
         NN_FINE[ i_block ] = NN_COAR[ i_block ] ;
         NN_COAR[ i_block ].set( false ) ;
      }
   }

   for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
   {
      PDE_DiscreteField const* ff = FF[i_block] ;
      std::vector< PDE_BasisFunctionCell* > const& bfs = BFS[ i_block ] ;
      boolVector& nfine = NN_FINE[ i_block ] ;
      boolVector& ncoar = NN_COAR[ i_block ] ;
      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         PDE_BasisFunctionCell* bf = bfs[n] ;
         PEL_ASSERT( bf != 0 ) ;
         if( nfine( n ) )
         {
            if( bf->refinement_level() == FINE_LEVEL )
            {
               for( size_t p = 0 ; p<bf->nb_parents() ; ++p )
               {
                  PDE_BasisFunctionCell* pbf = bf->parent( p ) ;
                  bool is_leading = bf->parent_is_leading( p ) ;
                  if( is_leading )
                  {
                     size_t node_parent = pbf->node_of_DOF( ff ) ;
                     ncoar( node_parent ) = true ;
                  }
               }
            }
            else if( !ncoar( n ) )
            {
               ncoar( n ) = true ;
            }
         }
      }
   }

   COL_NMB->reset( NN_COAR ) ;

   PEL_CHECK_POST( ( OLD(current_fine_level) == PEL::bad_index() &&
                     current_fine_level() == nb_levels()-1 ) ||
                   current_fine_level() == OLD(current_fine_level) - 1 ) ;
}

//-----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: build_current_prolongation_matrix(
                                       LA_Matrix* mat ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL("PDE_AlgebraicCoarsener:: build_current_prolongation_matrix") ;

   size_t nb_rows = ROW_NMB->nb_global_unknowns() ;
   size_t nb_cols = COL_NMB->nb_global_unknowns() ;
   size_t nb_local_rows = ROW_NMB->nb_unknowns_on_current_process() ;
   size_t nb_local_cols = COL_NMB->nb_unknowns_on_current_process() ;
   mat->re_initialize( nb_rows, nb_cols, nb_local_rows, nb_local_cols ) ;

   for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
   {
      PDE_DiscreteField const* ff = FF[ i_block ] ;
      PDE_LinkDOF2Unknown const* rlink = ROW_NMB->link( i_block ) ;
      PDE_LinkDOF2Unknown const* clink = COL_NMB->link( i_block ) ;
      std::vector< PDE_BasisFunctionCell* > const& bfs = BFS[ i_block ] ;
      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            if( clink->DOF_is_unknown( n, ic ) )
            {
               size_t i_col =
                      COL_NMB->global_unknown_for_DOF( n, ic, i_block ) ;
               if( rlink->DOF_is_unknown( n, ic ) )
               {
                  size_t i_row =
                         ROW_NMB->global_unknown_for_DOF( n, ic, i_block ) ;
                  mat->set_item( i_row, i_col, 1.0 ) ;
               }
               else
               {
                  PDE_BasisFunctionCell* bf = bfs[n] ;
                  for( size_t jj=0 ; jj<bf->nb_childs() ; ++jj )
                  {
                     double xx = bf->refinement_coefficient( jj ) ;
                     PDE_BasisFunctionCell* cbf = bf->child( jj ) ;
                     size_t cnode = cbf->node_of_DOF( ff ) ;
                     PEL_ASSERT( rlink->DOF_is_unknown( cnode, ic ) ) ;
                     size_t i_row =
                       ROW_NMB->global_unknown_for_DOF( cnode, ic, i_block ) ;
                     mat->set_item( i_row, i_col, xx ) ;
                  }
               }
            }
         }
      }
   }
   mat->synchronize() ;
}

//-----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: build_smoothing_lines(
                                       LA_Vector* smoothing_lines ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: build_smoothing_lines" ) ;
   PEL_CHECK_PRE( smoothing_lines != 0 ) ;

   size_t nb_rows = ROW_NMB->nb_global_unknowns() ;
   size_t nb_local_rows = ROW_NMB->nb_unknowns_on_current_process() ;
   smoothing_lines->re_initialize( nb_rows, nb_local_rows ) ;

   for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
   {
      PDE_DiscreteField const* ff = FF[ i_block ] ;
      PDE_LinkDOF2Unknown const* rlink = ROW_NMB->link( i_block ) ;
      std::vector< PDE_BasisFunctionCell* > const& bfs = BFS[ i_block ] ;
      boolVector const& nfine = NN_FINE[ i_block ] ;

      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         PDE_BasisFunctionCell* bf = bfs[n] ;
         PEL_ASSERT( bf != 0 ) ;
         if( nfine( n ) )
         {
            if( bf->refinement_level() == FINE_LEVEL )
            {
               bool to_be_smoothed = true ;
               for( size_t p = 0 ; p<bf->nb_parents() ; ++p )
               {
                  if( bf->parent_is_leading( p ) )
                  {
                     to_be_smoothed = false ;
                     break ;
                  }
               }
               if( to_be_smoothed )
               {
                  for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
                  {
                     if( rlink->DOF_is_unknown( n, ic ) )
                     {
                        size_t i_row =
                               ROW_NMB->global_unknown_for_DOF( n, ic, i_block ) ;
                        smoothing_lines->set_item( i_row, 1.0 ) ;
                     }
                  }
               }
               else
               {
                  PEL_ASSERT( bf->nb_parents() == 1 ) ;
               }
            }
         }
      }
   }
   smoothing_lines->synchronize() ;
}

//-----------------------------------------------------------------------
void
PDE_AlgebraicCoarsener:: build_current_level_unknowns(
                                    LA_Vector* current_level_unknowns ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AlgebraicCoarsener:: build_current_level_unknowns " ) ;
   PEL_CHECK_PRE( current_level_unknowns != 0 ) ;

   size_t nb_rows = ROW_NMB->nb_global_unknowns() ;
   size_t nb_local_rows = ROW_NMB->nb_unknowns_on_current_process() ;
   current_level_unknowns ->re_initialize( nb_rows, nb_local_rows ) ;

   for( size_t i_block=0 ; i_block<NB_BLOCKS ; ++i_block )
   {
      PDE_DiscreteField const* ff = FF[ i_block ] ;
      PDE_LinkDOF2Unknown const* rlink = ROW_NMB->link( i_block ) ;
      std::vector< PDE_BasisFunctionCell* > const& bfs = BFS[ i_block ] ;
      boolVector const& nfine = NN_FINE[ i_block ] ;

      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         PDE_BasisFunctionCell* bf = bfs[n] ;
         PEL_ASSERT( bf != 0 ) ;
         if( nfine( n ) )
         {
            if( bf->refinement_level() == FINE_LEVEL )
            {
               for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
               {
                  if( rlink->DOF_is_unknown( n, ic ) )
                  {
                     size_t i_row =
                        ROW_NMB->global_unknown_for_DOF( n, ic, i_block ) ;
                     current_level_unknowns->set_item( i_row, 1.0 ) ;
                  }
               }
            }
         }
      }
   }
   current_level_unknowns->synchronize() ;
}
