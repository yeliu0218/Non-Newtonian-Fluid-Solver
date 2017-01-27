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

#include <PDE_AdapterHN.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>
#include <boolArray2D.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_CrossProcessBFNumbering.hh>
#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DOFsInitializer.hh>
#include <PDE_DomainAndFields.hh>
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

struct PDE_AdapterHN_ERROR
{
} ;

//----------------------------------------------------------------------
PDE_AdapterHN:: PDE_AdapterHN( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               PDE_GridFE* a_grid,
                               PDE_SetOfBasisFunctions* bf_set,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DOM( dom )
   , GRID( a_grid )
   , DFs( dom->set_of_discrete_fields() )
   , RPP( PDE_RefinementPatternProvider::object( "toto" ) )
   , VERB_LEVEL( exp->has_entry( "verbose_level" ) ? 
                    exp->int_data( "verbose_level" ) : 0 )
   , BF_NUMBERING( 0 )
   , BF_SET( bf_set )
   , CTX( 0 )
   , COORDS( 0 )
   , R_INDIC( 0 )
{
   if( dom->is_distributed() )
      BF_NUMBERING = BF_SET->cross_process_numbering() ;
   
   PEL_ContextSimple* c = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create( c, doubleVector( 0 ) ) ;
   c->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
   CTX = c ;
   R_INDIC = exp->abstract_data( this, "refinement_indicator", CTX ) ;
   if( !R_INDIC->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "refinement_indicator", R_INDIC->undefined_variables() ) ;
   }
   if( R_INDIC->data_type()!=PEL_Data::Bool )
   {
      PEL_Error::object()->raise_bad_data_type(
            exp, "refinement_indicator", PEL_Data::Bool ) ;
   }

   PDE_MeshFE::set_refinement_pattern_provider( RPP ) ;
   
   GRID->set_verbose_level( VERB_LEVEL ) ;
}

//----------------------------------------------------------------------
PDE_AdapterHN:: ~PDE_AdapterHN( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_AdapterHN:: reset( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_AdapterHN:: adapt( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN:: adapt" ) ;
   
   PEL_Vector* cells = PEL_Vector::create( 0, 0 ) ;
   find_cells_to_refine( cells ) ;
   
   PEL_VectorIterator* it = cells->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_CellFE* ccell = static_cast< PDE_CellFE* >( it->item() ) ;
      refine_cell( ccell ) ;
   }
   it->destroy() ;
   
   cells->destroy() ; cells=0 ;
   
   // ******
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
   // ******
   
   bool verb = ( VERB_LEVEL > 1 ) ;
   DOM->apply_requests_of_DOFs_imposed_value_modules( verb ) ; 
}

//----------------------------------------------------------------------
void
PDE_AdapterHN:: refine_cell( PDE_CellFE* ccell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN:: refine_cell" ) ;
   
   if( VERB_LEVEL > 0 ) PEL::out() << "ccell: " << ccell->id_number() << endl ;
   
   GE_ReferencePolyhedronRefiner const* crf = RPP->cell_refiner( ccell ) ;
   GRID->refine_cell( ccell, crf ) ;
   
   size_t nb_elms = ccell->nb_reference_elements() ;
   PEL_ASSERT( nb_elms == 1 ) ;
   for( size_t ee=0 ; ee<nb_elms ; ++ee )
   {
      build_new_nodes_and_bfs( ccell, ee ) ;
      
      PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
      PEL_ASSERT( elm->name() == "PDE_2D_Q1isoNonConfA_4nodes" ) ;
      size_t nb_bfs = ccell->nb_basis_functions( ee ) ;
      PEL_ASSERT( nb_bfs == 4 ) ;
      size_t elts_index = PEL::bad_index() ;
      for( size_t ibf=0 ; ibf<nb_bfs ; ++ibf )
      {
         PDE_BasisFunctionCell* cbf = ccell->basis_function( ee, ibf ) ;
         if( ibf == 0 )
            elts_index = cbf->ref_elts_grp_index() ;
         else 
            PEL_ASSERT( elts_index == cbf->ref_elts_grp_index() ) ;
         if( cbf->nb_cells() == 1 )
         {
            PEL_ASSERT( cbf->cell( 0 ) == ccell ) ;
            cbf->set_inactive() ;
            cbf->start_field_iterator() ;
            for( ; cbf->valid_field() ; cbf->go_next_field() )
            {
               PDE_DiscreteField* ff = cbf->field() ;
               size_t nn = cbf->node_of_field() ;
               ff->set_node_inactive( nn ) ;
               for( size_t ic=0 ; ic < ff->nb_components() ; ++ic )
               {
                  if( ff->DOF_is_constrained( nn, ic ) )
                     ff->remove_constraint_for_DOF( nn, ic ) ;
               }
            }
         }
         else
         {
            ccell->remove_basis_function( ee, ibf ) ;
            cbf->remove_from_pieces( ccell, ee, ibf ) ;
            cbf->set_refined() ;
            for( size_t j=0 ; j<cbf->nb_childs() ; ++j )
            {
               PDE_BasisFunctionCell* rbf = cbf->child( j ) ;
               double coef = cbf->refinement_coefficient( j ) ; 
               cbf->start_field_iterator() ;
               for( ; cbf->valid_field() ; cbf->go_next_field() )
               {
                  PDE_DiscreteField* ff = cbf->field() ;
                  size_t slave_n  = cbf->node_of_field() ;
                  size_t master_n = rbf->node_of_DOF( ff ) ;
                  for( size_t ic=0 ; ic < ff->nb_components() ; ++ic )
                  {
                     ff->add_constraint_for_DOF( slave_n, ic, 
                                                 master_n, ic, coef ) ;
                  }
               }
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_AdapterHN:: build_new_nodes_and_bfs( PDE_CellFE const* ccell, 
                                                   size_t ee )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN:: build_new_nodes_and_bfs" ) ;
   
   size_t dummy_idx = PEL::bad_index() ;
   bool dummy_bool = false ;
   
   PDE_ReferenceElement const* elm = ccell->reference_element( ee ) ;
   size_t nb_bfs = ccell->nb_basis_functions( ee ) ;
   PEL_ASSERT( nb_bfs == 4 ) ;
   PDE_BasisFunctionCell* cbf = ccell->basis_function( ee, 0 ) ;
   size_t elts_index = cbf->ref_elts_grp_index() ;
   
   size_t_array2D idx_pt_of_mesh( elm->nb_nodes(), ccell->nb_childs() ) ;
   idx_pt_of_mesh.set( PEL::bad_index() ) ;
   
   boolArray2D isnew_pt_of_mesh( elm->nb_nodes(), ccell->nb_childs() ) ;

   std::vector< PDE_BasisFunctionCell* > rbfs ;
   
   GRID->start_nodepts_recording() ;
   GRID->start_preexisting_nodepts_recording() ;
   for( size_t ln=0 ; ln<ccell->nb_basis_functions( ee ) ; ++ln )
   {
      PDE_BasisFunctionCell* bf = ccell->basis_function( ee, ln ) ;
      if( bf->is_refined() )
      {
         for( size_t j=0 ; j<bf->nb_childs() ; ++j )
         {
            PDE_BasisFunctionCell* rbf = bf->child( j ) ;
            GRID->record_preexisting_nodept( rbf, dummy_idx ) ;
            if( dummy_idx < rbfs.size() )
            {
               PEL_ASSERT( rbfs[ dummy_idx ] == rbf ) ;
            }
            else
            {
               PEL_ASSERT( dummy_idx == rbfs.size() ) ;
               rbfs.push_back( rbf ) ;
            }
         }
      }
   }
   GRID->terminate_preexisting_nodepts_recording() ;
   for( size_t ic=0 ; ic<ccell->nb_childs() ; ++ic )
   {
      PDE_CellFE* rcell = ccell->child( ic ) ;
      if( VERB_LEVEL > 0 ) PEL::out() << "   refined cell "
                                      << rcell->id_number() << endl ;
      GE_Mpolyhedron const* rpoly = rcell->polyhedron() ;
      for( size_t in = 0 ; in<nb_bfs ; in++ )
      {               
         GE_Point const* pt = elm->node_location( in ) ;
         GRID->record_nodept( rpoly, pt, dummy_idx, dummy_bool ) ;
         
         isnew_pt_of_mesh( in, ic ) = dummy_bool ;
         idx_pt_of_mesh( in, ic ) = dummy_idx ;
         if( VERB_LEVEL > 0 )
         {
            PEL::out() << "      local node " << in ;
            if( dummy_bool ) PEL::out() << " new" ;
            else PEL::out() << " old" ;
            PEL::out() << " (nodept: " << dummy_idx << ")" << endl ;
         }
      }
   }
   GRID->terminate_nodepts_recording() ;
   
   size_t nb_pts = GRID->nb_recorded_nodepts() ;
   
   size_t nb_new_nodes = GRID->nb_inferred_nodes() ;
   if( nb_new_nodes == 0 ) return ; // <-----------------
   
   if( VERB_LEVEL > 0 ) PEL::out() << "nb_new_nodes:" << nb_new_nodes << endl ;
   
   // build the new basis functions
   for( size_t nn=0 ; nn<nb_new_nodes ; ++nn )
   {
      PDE_BasisFunctionCell* rbf = PDE_BasisFunctionCell::create( 0 ) ;
      BF_SET->add( rbf, elts_index ) ;
      rbfs.push_back( rbf ) ;
   }
   
   // add new nodes to the fields
   // append fields to the new bfs
   size_t nb_fields = ccell->nb_discrete_fields( ee ) ;
   for( size_t i=0 ; i<nb_fields ; ++i )
   {
      PDE_DiscreteField* ff = const_cast< PDE_DiscreteField* >( 
                                          ccell->discrete_field( ee, i ) ) ;
      
      size_t old_nb_nodes = ff->nb_nodes() ;
      ff->add_nodes( nb_new_nodes ) ;

      boolVector bf_done( rbfs.size() ) ;
      for( size_t ipt=0 ; ipt<nb_pts ; ++ipt )
      {
         size_t nn = GRID->inferred_node_of_recorded_nodept( ipt ) ;
         if( VERB_LEVEL > 0 ) PEL::out() << "nodept: " << ipt ;
         if( nn != PEL::bad_index() )
         {
            size_t ibf = GRID->bf_index_of_recorded_nodept( ipt ) ;
            if( !bf_done( ibf ) )
            {
               size_t gn = old_nb_nodes + nn ;
               if( VERB_LEVEL > 0 ) PEL::out() << "   new node: " << gn ; 
               ff->set_node_active( gn ) ;
               rbfs[ ibf ]->append_one_field( ff, gn ) ;
               bf_done( ibf ) = true ;
            }
         }
         else if( VERB_LEVEL > 0 ) PEL::out() << "   old node" ;
         if( VERB_LEVEL > 0 ) PEL::out() << endl ;
      }
   }
   
   for( size_t ic=0 ; ic<ccell->nb_childs() ; ++ic )
   {
      PDE_CellFE* rcell = ccell->child( ic ) ;
      PEL_ASSERT( rcell->nb_discrete_fields( ee ) == nb_fields ) ; //???
      for( size_t in=0 ; in<elm->nb_nodes() ; ++in )
      {
         size_t ipt = idx_pt_of_mesh( in, ic ) ;
         size_t ibf = GRID->bf_index_of_recorded_nodept( ipt ) ;
         PDE_BasisFunctionCell* rbf = rbfs[ ibf ] ;
         
         // relationship between rbf and rcell
         rcell->set_basis_function( ee, in, rbf ) ;
         rbf->extend_pieces( rcell, ee, in ) ;
         
         bool is_new = isnew_pt_of_mesh( in, ic ) ;
         if( is_new )
         {
            rbf->set_all_child_parent_relationship_on_cell(
                          rcell, in, ee, ccell, ic, VERB_LEVEL ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_AdapterHN:: find_cells_to_refine( PEL_Vector* cells )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN:: find_cells_to_refine" ) ;
   PEL_CHECK( cells != 0 ) ;
   
   PEL_VectorIterator* it = GRID->cells()->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_CellFE* ccell = static_cast< PDE_CellFE* >( it->item() ) ;
      COORDS->set( ccell->polyhedron()->center()->coordinate_vector() ) ;
      if( R_INDIC->to_bool() )
      {
         cells->append( ccell ) ;
      }
   }
   it->destroy() ;
}
