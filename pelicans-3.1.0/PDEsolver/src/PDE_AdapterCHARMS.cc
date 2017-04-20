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

#include <PDE_AdapterCHARMS.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_Activator.hh>
#include <PDE_AdaptationIndicator.hh>
#include <PDE_AdaptationRequestFromIndicator.hh>
#include <PDE_AlgebraicCoarsener.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_CrossProcessBFNumbering.hh>
#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DOFsInitializer.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FacesOfCellFE.hh>
#include <PDE_GridFE.hh>
#include <PDE_MeshingCoarsener.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ReferenceElementRefiner.hh>
#include <PDE_RefinementPatternProvider.hh>
#include <PDE_SetOfBasisFunctions.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <PDE_FaceFE.hh>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <list>
#include <map>
#include <set>
#include <sstream>

using std::cout ;
using std::endl ;
using std::list ;
using std::map ;
using std::set ;
using std::vector ;
using std::setw ;

struct PDE_AdapterCHARMS_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( void ) ;
   static void n3( void ) ;
} ;

//----------------------------------------------------------------------
PDE_AdapterCHARMS:: PDE_AdapterCHARMS(  PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        PDE_GridFE* a_grid,
                                        PDE_SetOfBasisFunctions* bf_set,
                                        PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , HB( false )
   , MAX_REF_LEVEL( PEL::bad_index() )
   , ACTIVATOR( 0 )
   , DOM( dom )
   , GRID( a_grid )
   , DFs( dom->set_of_discrete_fields() )
   , PT( GE_Point::create( this, a_grid->nb_space_dimensions() ) )
   , RPP( PDE_RefinementPatternProvider::object( "toto" ) )
   , INDICS( PEL_List::create( this ) )
   , INDICS_IT( 0 )
   , SOMETHING_CHANGED( false )
   , VERB_LEVEL( 0 )
   , NB_ACTIVE_CELLS( a_grid->nb_cells() )
   , NB_ACT_CELLS( 0 )
   , NB_DEACT_CELLS( 0 )
   , A_COARSENER( 0 )
   , M_COARSENER( 0 )
   , BF_NUMBERING( 0 )
   , BF_SET( bf_set )
{
   if( dom->is_distributed() )
      BF_NUMBERING = BF_SET->cross_process_numbering() ;

   if( exp->has_entry( "verbose_level" ) )
      VERB_LEVEL = exp->int_data( "verbose_level" ) ;

   PDE_MeshFE::set_refinement_pattern_provider( RPP ) ;

   std::string const& tt = exp->string_data( "type" ) ;
   if( tt == "hierarchical_basis" )
   {
      HB = true ;
   }
   else if( tt == "quasi_hierarchical_basis" )
   {
      HB = false ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value( exp, "type",
         "   \"hierarchical_basis\"\n   \"quasi_hierarchical_basis\"" ) ;
   }
   if( exp->has_entry( "highest_refinement_level" ) )
   {
      MAX_REF_LEVEL = exp->int_data( "highest_refinement_level" ) ;
   }

   PDE_ReferenceElementRefiner::OneLevelDiffRule olr =
                                PDE_ReferenceElementRefiner::Parents ;
   if( exp->has_entry( "one_level_difference_rule" ) )
   {
      std::string const& ss = exp->string_data( "one_level_difference_rule" ) ;
      if( ss == "supports" )
      {
         olr =  PDE_ReferenceElementRefiner::Supports ;
      }
      else if( ss == "parents" )
      {
         olr =  PDE_ReferenceElementRefiner::Parents ;
      }
      else if( ss == "no" )
      {
         olr =  PDE_ReferenceElementRefiner::No ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value( exp,
                   "one_level_difference_rule",
                   "- \"parents\"\n"
                   "- \"supports\"\n"
                   "- \"no\"\n" ) ;
      }
   }
   PDE_ReferenceElementRefiner::set_one_level_difference_rule( olr ) ;


   ACTIVATOR = PDE_Activator::create( this, HB, VERB_LEVEL ) ;

   if( exp->has_module( "PDE_AdaptationIndicator" ) )
   {
      if( exp->has_module( "list_of_PDE_AdaptationIndicator" ) )
         PDE_AdapterCHARMS_ERROR::n0() ;
      PEL_ModuleExplorer const* se =
                exp->create_subexplorer( 0, "PDE_AdaptationIndicator" ) ;
      PDE_AdaptationIndicator* indic =
                PDE_AdaptationIndicator::make( INDICS, dom, se, VERB_LEVEL ) ;
      se->destroy() ; se = 0 ;
      INDICS->append( indic ) ;
   }
   else if( exp->has_module( "list_of_PDE_AdaptationIndicator" ) )
   {
      if( exp->has_module( "PDE_AdaptationIndicator" ) )
         PDE_AdapterCHARMS_ERROR::n0() ;
      PEL_ModuleExplorer* se =
          exp->create_subexplorer( 0, "list_of_PDE_AdaptationIndicator" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         PEL_ModuleExplorer* sse = se->create_subexplorer( se ) ;
         PDE_AdaptationIndicator* indic =
               PDE_AdaptationIndicator::make( INDICS, dom, sse, VERB_LEVEL ) ;
         INDICS->append( indic ) ;
      }
      se->destroy() ; se = 0 ;
   }
}

//----------------------------------------------------------------------
PDE_AdapterCHARMS:: ~PDE_AdapterCHARMS( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: append_indicator( PDE_AdaptationIndicator* indic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: append_indicator" ) ;
   PEL_CHECK_PRE( indic != 0 ) ;
   PEL_CHECK_PRE( indic->owner() == this ) ;

   if( INDICS_IT != 0 ) PDE_AdapterCHARMS_ERROR::n3() ;

   change_owner( INDICS, indic ) ;
   INDICS->append( indic ) ;

   PEL_CHECK_POST( indic->is_under_ownership_of( this ) ) ;
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: add_excluded_field( PDE_DiscreteField const* ff )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: add_excluded_field" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   ACTIVATOR->add_excluded_field( ff ) ;

   PEL_CHECK_POST( is_excluded( ff ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_AdapterCHARMS:: is_excluded( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: is_excluded" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   bool result = ACTIVATOR->is_excluded( ff ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: reset( void )
//----------------------------------------------------------------------
{
   if( INDICS_IT == 0 )
   {
      if( INDICS->count() == 0 ) PDE_AdapterCHARMS_ERROR::n1() ;
      INDICS_IT = INDICS->create_iterator( INDICS ) ;
   }

   for( INDICS_IT->start() ; INDICS_IT->is_valid() ; INDICS_IT->go_next() )
   {
      PDE_AdaptationIndicator* indic =
         static_cast< PDE_AdaptationIndicator* >( INDICS_IT->item() ) ;
      indic->reset() ;
   }

   if( A_COARSENER != 0 ) A_COARSENER->reset() ;
   if( M_COARSENER != 0 ) M_COARSENER->reset() ;

   PEL_CHECK_POST( !algebraic_coarsener()->coarsening_is_possible() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->nb_levels() == PEL::bad_index() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->current_fine_level() ==
                    PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: adapt( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: adapt" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( INDICS_IT == 0 ) PDE_AdapterCHARMS_ERROR::n2() ;

   static bool const make_refined_DOFs_bad_double = true ;

   reset_counters() ;
   SOMETHING_CHANGED = false ;

   for( INDICS_IT->start() ; INDICS_IT->is_valid() ; INDICS_IT->go_next() )
   {
      PDE_AdaptationIndicator* indic =
         static_cast< PDE_AdaptationIndicator* >( INDICS_IT->item() ) ;
      indic->build() ;
   }
   PDE_AdaptationRequestFromIndicator* adap =
                 PDE_AdaptationRequestFromIndicator::create( 0, INDICS,
                                                             MAX_REF_LEVEL,
                                                             VERB_LEVEL ) ;

   adap->apply_criterion_to_infer_bfs_to_refine( GRID->cells() ) ;
   adap->apply_criterion_to_infer_bfs_to_unrefine( GRID->cells() ) ;

   std::string const indent( 10, ' ' ) ;

   // *** refinement

   adap->compute_bfs_to_refine( BF_SET ) ;

   if( adap->something_to_refine() )
   {
      SOMETHING_CHANGED = true ;

      adap->build_info_on_cells_for_refinement() ;

      split_meshes( adap ) ;

      build_new_nodes_and_bfs( adap ) ;

      ACTIVATOR->refine( adap, make_refined_DOFs_bad_double ) ;

      for( DFs->start() ; DFs->is_valid() ; DFs->go_next() )
      {
         PDE_DiscreteField* f = DFs->item() ;
         if( f->is_distributed() )
         {
            f->cross_process_numbering()->synchronize_new_nodes() ;
            f->cross_process_numbering()->synchronize_imposed_values() ;
            f->cross_process_numbering()->synchronize_valid_nodes() ;
         }
      }
      if( DOM->is_distributed() )
         BF_NUMBERING->synchronize_new_basis_functions() ;
   }

   // *** unrefinement

   adap->compute_bfs_to_unrefine( BF_SET ) ;

   if( adap->something_to_unrefine() )
   {
      SOMETHING_CHANGED = true ;

      adap->build_info_on_cells_for_unrefinement() ;

      ACTIVATOR->unrefine( adap ) ;

      ACTIVATOR->unsplit_meshes_1( adap ) ;
   }

   // ***

   if( SOMETHING_CHANGED )
   {
      bool verb = ( VERB_LEVEL > 1 ) ;
      DOM->apply_requests_of_DOFs_imposed_value_modules( verb ) ;

      GRID->change_geometric_state() ;
   }
   adap->destroy() ;

   if( A_COARSENER != 0 ) A_COARSENER->reset() ;
   if( M_COARSENER != 0 ) M_COARSENER->reset() ;

   NB_ACT_CELLS += ACTIVATOR->nb_activated_cells() ;
   NB_DEACT_CELLS += ACTIVATOR->nb_deactivated_cells() ;
   NB_ACTIVE_CELLS = NB_ACTIVE_CELLS + NB_ACT_CELLS - NB_DEACT_CELLS ;

   PEL_CHECK_POST( !algebraic_coarsener()->coarsening_is_possible() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->nb_levels() == PEL::bad_index() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->current_fine_level() ==
                    PEL::bad_index() ) ;
}

//-----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: unsplit_meshes( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: unsplit_meshes" ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_AdapterCHARMS:: unsplit_meshes ***"
                 << endl << endl ;

   reset_counters() ;
   SOMETHING_CHANGED = false ;

   list< PDE_CellFE* > cccs ;

   PEL_Vector const* cells = GRID->cells() ;
   size_t nb_cells = cells->count() ;
   for( size_t ic=0 ; ic<nb_cells ; ++ic )
   {
      PDE_CellFE* ccell = static_cast<PDE_CellFE*>( cells->at( ic ) ) ;
      if( !ccell->is_active() && ccell->nb_childs()!=0 )
      {
         possibly_unsplit_cell( ccell, cccs ) ;
      }
   }

   while( !cccs.empty() )
   {
      PDE_CellFE* ccell = cccs.front() ;
      cccs.pop_front() ;
      possibly_unsplit_cell( ccell, cccs ) ;
   }

   if( SOMETHING_CHANGED ) GRID->change_geometric_state() ;

   if( A_COARSENER != 0 ) A_COARSENER->reset() ;
   if( M_COARSENER != 0 ) M_COARSENER->reset() ;

   NB_ACT_CELLS += ACTIVATOR->nb_activated_cells() ;
   NB_DEACT_CELLS += ACTIVATOR->nb_deactivated_cells() ;
   NB_ACTIVE_CELLS = NB_ACTIVE_CELLS + NB_ACT_CELLS - NB_DEACT_CELLS ;

   PEL_CHECK_POST( !algebraic_coarsener()->coarsening_is_possible() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->nb_levels() == PEL::bad_index() ) ;
   PEL_CHECK_POST(  algebraic_coarsener()->current_fine_level() ==
                    PEL::bad_index() ) ;
}

//-----------------------------------------------------------------------
bool
PDE_AdapterCHARMS:: something_changed( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: something_changed" ) ;

   return( SOMETHING_CHANGED ) ;
}

//-----------------------------------------------------------------------
PDE_AlgebraicCoarsener*
PDE_AdapterCHARMS:: algebraic_coarsener( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: algebraic_coarsener" ) ;

   if( A_COARSENER == 0 )
   {
      A_COARSENER = new PDE_AlgebraicCoarsener(
                            const_cast< PDE_AdapterCHARMS* >( this ),
                            GRID ) ;
   }

   PDE_AlgebraicCoarsener* result = A_COARSENER ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
PDE_MeshingCoarsener*
PDE_AdapterCHARMS:: meshing_coarsener( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: meshing_coarsener" ) ;

   if( M_COARSENER == 0 )
   {
      M_COARSENER = new PDE_MeshingCoarsener(
                            const_cast< PDE_AdapterCHARMS* >( this ),
                            DOM,
                            GRID,
                            ACTIVATOR,
                            VERB_LEVEL ) ;
   }

   PDE_MeshingCoarsener* result = M_COARSENER ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: print_statistics( std::ostream& os,
                                              size_t indent_width ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: print_statistics" ) ;

   std::string const s( indent_width, ' ' ) ;

   os << s << "nb   activated cells : "
           << NB_ACT_CELLS << endl ;
   os << s << "nb deactivated cells : "
           << NB_DEACT_CELLS << endl ;
   os << s << "nb   active    cells : "
           << NB_ACTIVE_CELLS << endl ;
}

//-----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: possibly_unsplit_cell( PDE_CellFE* ccell,
                                           std::list< PDE_CellFE* >& cccs )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: possibly_unsplit_cell" ) ;

   bool do_it = true ;
   for( size_t i=0 ; i<ccell->nb_childs() ; ++i )
   {
      PDE_CellFE const* rcell = ccell->child( i ) ;
      do_it = do_it && rcell->is_active() &&
                       rcell->all_basis_functions_are_dropped() ;
   }

   if( do_it )
   {
      ACTIVATOR->activate_and_deactivate_child_cells( ccell ) ;
      SOMETHING_CHANGED = true ;

      PDE_CellFE* cpar = ccell->parent() ;
      if( cpar != 0 )
      {
         list< PDE_CellFE* >::iterator where =
                                       find( cccs.begin(), cccs.end(), cpar ) ;
         if( where == cccs.end() )
         {
            cccs.push_back( cpar ) ;
         }
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: split_meshes( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: split_meshes" ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_AdapterCHARMS:: split_meshes ***"
      << endl << endl ;

   GRID->set_verbose_level( VERB_LEVEL ) ;

   adap->start_cell_to_refine() ;
   for( ; adap->valid_cell_to_refine() ; adap->go_next_cell_to_refine() )
   {
      PDE_CellFE* ccell = adap->current_cell_to_refine() ;
      GE_ReferencePolyhedronRefiner const* crf = RPP->cell_refiner( ccell ) ;
      if( ccell->is_active() )
      {
         if( ccell->nb_childs() == 0 )
         {
            GRID->refine_cell( ccell, crf ) ;
            NB_DEACT_CELLS++ ;
            NB_ACT_CELLS += crf->nb_subcells() ;
         }
         else
         {
            //??? devrait etre dans PDE_GridFE ???
            ACTIVATOR->deactivate_and_activate_child_cells( ccell ) ;
         }
      }
   }

   if( VERB_LEVEL > 1 ) PEL::out() << endl ;
}

//-----------------------------------------------------------------------
void
PDE_AdapterCHARMS:: build_new_nodes_and_bfs( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS:: build_new_nodes_and_bfs" ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_AdapterCHARMS:: build_new_nodes_and_bfs"
      << endl << endl ;

   // STEP 0
   // ------
   //   computation and assignment of dimensions for working arrays

   size_t max_nb_ccells = 0 ;
   size_t max_nb_childs_per_ccell = 0 ;
   size_t max_nb_rcells = 0 ;
   size_t max_nb_pts_per_rcell = 0 ;
   adap->start_bf_to_refine() ;
   for( ; adap->valid_bf_to_refine() ; adap->go_next_bf_to_refine() )
   {
      PDE_BasisFunctionCell* cbf = adap->current_bf_to_refine() ;
      size_t nn = cbf->nb_cells() ;
      if( nn > max_nb_ccells ) max_nb_ccells = nn ;
      for( size_t i=0 ; i<cbf->nb_cells() ; ++i )
      {
         PDE_CellFE* ccell = cbf->cell( i ) ;

         nn = ccell->nb_childs() ;
         if( nn > max_nb_childs_per_ccell ) max_nb_childs_per_ccell = nn ;

         nn = cbf->nb_cells() * ccell->nb_childs() ;
         if( nn > max_nb_rcells ) max_nb_rcells = nn ;

         size_t ee = cbf->element_index_of_cell( i ) ;
         PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
         PDE_ReferenceElementRefiner const* elrf =
                                     RPP->reference_element_refiner( elm ) ;
         nn = elrf->reference_element()->nb_nodes() ;
         if( nn > max_nb_pts_per_rcell ) max_nb_pts_per_rcell = nn ;
      }
   }

   size_t_array2D idx_cell( max_nb_ccells, max_nb_childs_per_ccell ) ;
   size_t_array2D idx_pt_of_mesh( max_nb_pts_per_rcell, max_nb_rcells ) ;

   adap->start_bf_to_refine() ;
   for( ; adap->valid_bf_to_refine() ; adap->go_next_bf_to_refine() )
   {
      PDE_BasisFunctionCell* cbf = adap->current_bf_to_refine() ;
      PEL_ASSERT( !cbf->is_refined() ) ;

      // For all coarse basis function cbf to be refined, refined
      //    basis functions are built on the support of cbf.
      // The support of cbf is made of cells whose refinement level is the same
      //    as that of cbf, the set of all the childs of these cells forms a
      //    submeshing over which the refined basis functions are built

      // STEP 1
      // ------
      //    construction of the set of all new node points associated to the
      //    childs of the coarse basis function (these node points are
      //    recorded by PDE_GridFE)
      //
      // idx_cell( i, ic )
      // *****************
      //    index in the set of refined cells of the ic-th child of the i-th
      //    cell of the support of cbf
      //
      // idx_pt_of_mesh( in, im )
      // ************************
      //    for the im-th refined cell, index (in PDE_GridFE) of
      //    the in-th node point

      GRID->start_nodepts_recording() ;

      idx_cell.set( PEL::bad_index() ) ;
      idx_pt_of_mesh.set( PEL::bad_index() ) ;
      size_t im = 0 ;
      for( size_t i=0 ; i<cbf->nb_cells() ; ++i )
      {
         PDE_CellFE* ccell = cbf->cell( i ) ;
         size_t ee = cbf->element_index_of_cell( i ) ;
         PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
         PDE_ReferenceElementRefiner const* elrf =
                                     RPP->reference_element_refiner( elm ) ;
         size_t parent_n = cbf->local_node_of_cell( i ) ;
         if( VERB_LEVEL > 0 )
            PEL::out() << "         refined nodes from coarse node "
                       << parent_n << endl ;
         for( size_t ic=0 ; ic<ccell->nb_childs() ; ++ic )
         {
            PDE_CellFE* rcell = ccell->child( ic ) ;
            if( VERB_LEVEL > 0 ) PEL::out() << "   refined cell "
                                            << rcell->id_number() << endl ;
            GE_Mpolyhedron const* rpoly = rcell->polyhedron() ;
            for( size_t j=0 ; j<elrf->nb_childs( parent_n, ic ) ; ++j )
            {
               size_t child_n = elrf->child_node( parent_n, ic, j ) ;

               if( HB && child_n==elrf->leading_child_node( parent_n, ic ) )
                  continue ;  // <--------------------

               // the fine basis function may already exist due to the
               // refinement of an other coarse basis function
               if( rcell->basis_function( ee, child_n ) == 0 )
               {
                  size_t idx = PEL::bad_index() ;
                  bool is_new = false ;
                  GRID->record_nodept( rpoly, elm->node_location( child_n ),
                                       idx, is_new ) ;
                  PEL_ASSERT( is_new ) ;
                  idx_pt_of_mesh( child_n, im ) = idx ;
                  if( VERB_LEVEL > 0 )
                     PEL::out() << "            child node " << child_n << endl;
               }
            }
            idx_cell( i, ic ) = im ;
            ++im ;
         }
      }
      GRID->terminate_nodepts_recording() ;

      size_t nb_pts = GRID->nb_recorded_nodepts() ;
      if( nb_pts == 0 ) continue ; // <-----------------

      // STEP 2
      // ------
      //    assign an index to each new node point (two different node points
      //    may represent the same basis function, eg due to periodicity
      //    requests) and creation of fine basis functions

      size_t nb_new_nodes = GRID->nb_inferred_nodes() ;

      std::vector<PDE_BasisFunctionCell*> rbfs( nb_pts, 0 ) ;
      for( size_t ibf=0 ; ibf<nb_new_nodes ; ++ibf )
      {
         PDE_BasisFunctionCell* rbf = PDE_BasisFunctionCell::create( 0 ) ;
         BF_SET->add( rbf, cbf->ref_elts_grp_index() ) ;
         rbfs[ ibf ] = rbf ;
      }

      // STEP 3
      // ------
      //    assign a node number for each field and to each new basis function

      cbf->start_field_iterator() ;
      for( ; cbf->valid_field() ; cbf->go_next_field() )
      {
         PDE_DiscreteField* ff = cbf->field() ;

         size_t old_nb_nodes = ff->nb_nodes() ;
         ff->add_nodes( nb_new_nodes ) ;
         boolVector bf_done( nb_new_nodes ) ;
         for( size_t ip=0 ; ip<nb_pts ; ++ip )
         {
            size_t nn = GRID->inferred_node_of_recorded_nodept( ip ) ;
            PEL_ASSERT( nn != PEL::bad_index() ) ;
            size_t ibf = GRID->bf_index_of_recorded_nodept( ip ) ;
            PEL_ASSERT( ibf == nn ) ;
            if( !bf_done( ibf ) )
            {
               size_t gn = old_nb_nodes + nn ;
               rbfs[ ibf ]->append_one_field( ff, gn ) ;
               bf_done( ibf ) = true ;
            }
         }
      }

      // STEP 4
      // ------
      //    set various relationships

      for( size_t i=0 ; i<cbf->nb_cells() ; ++i )
      {
         PDE_CellFE* ccell = cbf->cell( i ) ;
         size_t ee = cbf->element_index_of_cell( i ) ;
         PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
         PDE_ReferenceElementRefiner const* elrf =
                                        RPP->reference_element_refiner( elm ) ;
         size_t parent_n = cbf->local_node_of_cell( i ) ;
         for( size_t ic=0 ; ic<ccell->nb_childs() ; ++ic )
         {
            PDE_CellFE* rcell = ccell->child( ic ) ;
            size_t i_rcell = idx_cell( i, ic ) ;
            for( size_t j=0 ; j<elrf->nb_childs( parent_n, ic ) ; ++j )
            {
               size_t child_n = elrf->child_node( parent_n, ic, j ) ;
               size_t ipt = idx_pt_of_mesh( child_n, i_rcell ) ;
               if( ipt != PEL::bad_index() )
               {
                  size_t ibf = GRID->bf_index_of_recorded_nodept( ipt ) ;
                  PDE_BasisFunctionCell* rbf = rbfs[ ibf ] ;

                  // STEP 4.1 : relationship between rbf and rcell
                  rcell->set_basis_function( ee, child_n, rbf ) ;
                  rbf->extend_pieces( rcell, ee, child_n ) ;

                  if( VERB_LEVEL > 2 ) PEL::out() << endl ;

                  // STEP 4.2 : - child/parent relationship
                  //            - ascendant relationship
                  // Note that these relationship are setting with
                  // parents, children, ascendants, descendants
                  // which are already built
                  rbf->set_all_child_parent_relationship_on_cell(
                                rcell, child_n, ee, ccell, ic, VERB_LEVEL ) ;
                  rbf->set_all_ascendant_relationship_on_cell(
                                rcell, child_n, ee, ccell, ic, VERB_LEVEL ) ;
               }
            }
         }
      }
   }
   if( VERB_LEVEL > 0 ) PEL::out() << std::endl ;
}

//----------------------------------------------------------------------------
void
PDE_AdapterCHARMS:: reset_counters( void )
//----------------------------------------------------------------------------
{
   ACTIVATOR->reset_counters() ;
   NB_ACT_CELLS = 0 ;
   NB_DEACT_CELLS = 0 ;
}

//internal--------------------------------------------------------------
void
PDE_AdapterCHARMS_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_AdapterCHARMS:" << endl ;
   msg << "    at most one of the following modules is accepted:" << endl ;
   msg << "       MODULE PDE_AdaptationIndicator" << endl ;
   msg << "       MODULE list_of_PDE_AdaptationIndicator" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_AdapterCHARMS_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_AdapterCHARMS:" << endl ;
   msg << "    at least one instance of" << endl ;
   msg << "       PDE_AdaptationIndicator" << endl ;
   msg << "    is required" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_AdapterCHARMS_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_AdapterCHARMS:" << endl ;
   msg << "    reset() should be called at least once before adapt()" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_AdapterCHARMS_ERROR:: n3( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** PDE_AdapterCHARMS:" << endl ;
   msg << "    append_indicator() can only be called before" << endl ;
   msg << "    before the first call to reset()" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

