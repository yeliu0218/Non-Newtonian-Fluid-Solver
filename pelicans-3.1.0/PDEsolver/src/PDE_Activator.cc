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

#include <PDE_Activator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_AdaptationRequest.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_FacesOfCellFE.hh>
#include <PDE_GridFE.hh>
#include <PDE_FaceFE.hh>
#include <PDE_ReferenceElementRefiner.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ;
using std::endl ;
using std::set ;
using std::setw ;

//----------------------------------------------------------------------
PDE_Activator*
PDE_Activator:: create( PEL_Object* a_owner,
                        bool hierarchical_basis,
                        size_t verbose_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: create" ) ;

   PDE_Activator* result = new PDE_Activator( a_owner,
                                              hierarchical_basis,
                                              verbose_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_Activator:: PDE_Activator( PEL_Object* a_owner,
                               bool hierarchical_basis,
                               size_t verbose_level )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , HB( hierarchical_basis )
   , SOMETHING_CHANGED( false )
   , VERB_LEVEL( verbose_level )
   , NB_ACT_CELLS( 0 )
   , NB_DEACT_CELLS( 0 )
{
}

//----------------------------------------------------------------------
PDE_Activator:: ~PDE_Activator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_Activator:: add_excluded_field( PDE_DiscreteField const* ff )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: add_excluded_field" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   EXCLUDED_FIELDS.insert( ff ) ;

   PEL_CHECK_POST( is_excluded( ff ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_Activator:: is_excluded( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: is_excluded" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   set< PDE_DiscreteField const* >::const_iterator it = EXCLUDED_FIELDS.find( ff ) ;
   bool result = ( it != EXCLUDED_FIELDS.end() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
void
PDE_Activator:: split_meshes( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: split_meshes" ) ;
   PEL_CHECK_PRE( adap != 0 ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_Activator:: split_meshes ***"
                 << endl << endl ;

   adap->start_cell_to_refine() ;
   for( ; adap->valid_cell_to_refine() ; adap->go_next_cell_to_refine() )
   {
      PDE_CellFE* ccell = adap->current_cell_to_refine() ;
      if( ccell->is_active() )
      {
         PEL_ASSERT( ccell->nb_childs() != 0 ) ;
         deactivate_and_activate_child_cells( ccell ) ;
      }
   }

   if( VERB_LEVEL > 0 ) PEL::out() << endl ;
}

//-----------------------------------------------------------------------
void
PDE_Activator:: refine( PDE_AdaptationRequest* adap,
                        bool make_refined_DOFs_bad_double )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: refine" ) ;
   PEL_CHECK_PRE( adap != 0 ) ;

   if( HB )
   {
      hierarchical_refine( adap, make_refined_DOFs_bad_double ) ;
   }
   else
   {
      quasi_hierarchical_refine( adap, make_refined_DOFs_bad_double ) ;
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: unrefine( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: unrefine" ) ;
   PEL_CHECK_PRE( adap != 0 ) ;

   if( HB )
   {
      hierarchical_unrefine( adap ) ;
   }
   else
   {
      quasi_hierarchical_unrefine( adap ) ;
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: unsplit_meshes_1( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: unsplit_meshes_1" ) ;
   PEL_CHECK_PRE( adap != 0 ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_Activator:: unsplit_meshes_1 ***"
                 << endl << endl ;

   adap->start_cell_to_unrefine() ;
   for( ; adap->valid_cell_to_unrefine() ; adap->go_next_cell_to_unrefine() )
   {
      PDE_CellFE* ccell = adap->current_cell_to_unrefine() ;
      PEL_ASSERT( !ccell->is_active() ) ;
      bool do_it = ccell->child( 0 )->all_basis_functions_are_dropped() ;
      for( size_t i=1 ; i<ccell->nb_childs() ; ++i )
      {
         do_it = do_it &&
         ccell->child( i )->all_basis_functions_are_dropped() ;
      }
      if( do_it )
      {
         activate_and_deactivate_child_cells( ccell ) ;
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: activate_and_deactivate_child_cells( PDE_CellFE* ccell )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: activate_and_deactivate_child_cells" ) ;
   PEL_CHECK_PRE( ccell != 0 ) ;

   for( size_t j=0 ; j<ccell->nb_childs() ; ++j )
   {
      PDE_CellFE* rcell = ccell->child( j ) ;
      PDE_FacesOfCellFE* itrs = rcell->create_faces_iterator( 0 ) ;
      for( itrs->start() ; itrs->is_valid() ; itrs->go_next() )
      {
         PDE_FaceFE* rs = itrs->item() ;
         if( rs->is_active() && !rs->has_adjacent_bound() )
         {
            //???? Modification for periodic domains
            //????    to be validated
            //???? before, instead of this test there was:
            //????    PEL_ASSERT( rs->nb_adjacent_cells() == 2 ) ;
            if( !rs->is_periodic() )
            {
               PDE_CellFE* other_cell = rs->adjacent_cell_other_than( rcell ) ;
               if( other_cell->parent() != ccell )
               {
                  size_t rlevel = rs->refinement_level() ;
                  if( other_cell->refinement_level() == rlevel )
                  {
                     rs->replace_adjacent_cell( rcell, ccell ) ;
                  }
                  else
                  {
                     PDE_FaceFE* parent_rs = rs->parent() ;
                     parent_rs->set_adjacent_cells( ccell, other_cell ) ;
                  }
               }
            }
         }
      }
      itrs->destroy() ;
   }

   if( VERB_LEVEL > 1 ) PEL::out() << "activating cell "
                                   << ccell->id_number() << endl ;
   NB_ACT_CELLS++ ;
   ccell->set_active() ;
   for( size_t j=0 ; j<ccell->nb_childs() ; ++j )
   {
      PDE_CellFE* rcell = ccell->child( j ) ;
      if( VERB_LEVEL > 1 ) PEL::out() << "   deactivating cell "
                                      << rcell->id_number() << endl ;
      PEL_ASSERT( rcell->all_basis_functions_are_dropped() ) ;
      PEL_ASSERT( rcell->is_active() ) ;
      NB_DEACT_CELLS++ ;
      rcell->set_inactive() ;
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: deactivate_and_activate_child_cells( PDE_CellFE* ccell )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: deactivate_and_activate_child_cells" ) ;
   PEL_CHECK_PRE( ccell != 0 ) ;

   for( size_t j=0 ; j<ccell->nb_childs() ; ++j )
   {
      PDE_CellFE* rcell = ccell->child( j ) ;
      size_t rlevel = rcell->refinement_level() ;
      PDE_FacesOfCellFE* itrs = rcell->create_faces_iterator( 0 ) ;
      for( itrs->start() ; itrs->is_valid() ; itrs->go_next() )
      {
         PDE_FaceFE* rs = itrs->item() ;

         if( rs->has_adjacent_bound() ) continue ; //<--------

         //???? Modification for periodic domains
         //????    to be validated
         //???? this line did not exist before
         if( rs->is_periodic() ) continue ;        //<--------

         if( rs->refinement_level() == rlevel )
         {
            PDE_CellFE* ocell = ( rs->has_adjacent_cell( rcell ) ?
                                  rs->adjacent_cell_other_than( rcell ) : 0 ) ;
            if( ocell != 0 &&
                ocell->refinement_level() == rlevel &&
                ocell->parent() == ccell )
            {
               //???? on passe plusieurs fois ici
               rs->set_adjacent_cells( rcell, ocell ) ;
            }
            else
            {
               if( rs->is_active() )
               {
                  rs->replace_adjacent_cell( ccell, rcell ) ;
               }
               else
               {
                  PDE_FaceFE* prs = rs->parent() ;
                  if( prs->is_active() )
                  {
                     ocell = prs->adjacent_cell_other_than( ccell ) ;
                     rs->set_adjacent_cells( rcell, ocell ) ;
                  }
               }
            }
         }
         else if( rs->is_active() )
         {
            rs->replace_adjacent_cell( ccell, rcell ) ;
         }
      }
      itrs->destroy() ;
   }

   if( VERB_LEVEL > 1 ) PEL::out() << "deactivating cell "
                                   << ccell->id_number() << endl ;
   NB_DEACT_CELLS++ ;
   ccell->set_inactive() ;
   for( size_t j=0 ; j<ccell->nb_childs() ; ++j )
   {
      PDE_CellFE* rcell = ccell->child( j ) ;
      if( VERB_LEVEL > 1 ) PEL::out() << "   activating cell "
                                      << rcell->id_number() << endl ;
      PEL_ASSERT( !rcell->is_active() ) ;
      NB_ACT_CELLS++ ;
      rcell->set_active() ;
   }
}

//----------------------------------------------------------------------------
void
PDE_Activator:: reset_counters( void )
//----------------------------------------------------------------------------
{
   NB_ACT_CELLS = 0 ;
   NB_DEACT_CELLS = 0 ;
}

//----------------------------------------------------------------------------
size_t
PDE_Activator:: nb_activated_cells( void ) const
//----------------------------------------------------------------------------
{
   return( NB_ACT_CELLS ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_Activator:: nb_deactivated_cells( void ) const
//----------------------------------------------------------------------------
{
   return( NB_DEACT_CELLS ) ;
}

//---------------------------------------------------------------------------
void
PDE_Activator:: quasi_hierarchical_refine(
                      PDE_AdaptationRequest* adap,
                      bool make_refined_DOFs_bad_double ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: quasi_hierarchical_refine" ) ;

   adap->start_bf_to_refine() ;
   for( ; adap->valid_bf_to_refine() ; adap->go_next_bf_to_refine() )
   {
      PDE_BasisFunctionCell* cbf = adap->current_bf_to_refine() ;

      cbf->set_inactive() ;
      cbf->set_refined() ;

      for( size_t i=0 ; i<cbf->nb_childs() ; ++i )
      {
         PDE_BasisFunctionCell* rbf = cbf->child( i ) ;

         //PEL_ASSERT( PDE_ReferenceElementRefiner::one_level_difference_rule() !=
         //                       PDE_ReferenceElementRefiner::Supports ||
         //            !rbf->is_refined()  ) ;

         if( !rbf->is_refined() )//on active uniquement
         //des fonctions de bases non raffinées
         {
            rbf->start_field_iterator() ;
            for( ; rbf->valid_field() ; rbf->go_next_field() )
            {
               PDE_DiscreteField* ff = rbf->field() ;
               if( !is_excluded( ff ) )
               {
                  size_t rn = rbf->node_of_field() ;
                  size_t cn = cbf->node_of_DOF( ff ) ;
                  if( make_refined_DOFs_bad_double )
                  {
                     for( size_t l=0 ; l<ff->storage_depth() ; ++l )
                     {
                        for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
                        {
                           //--->  if the "refinement equation" is used
                           // double xx = ff->DOF_value( l, cn, ic ) * coef ;
                           // if( rbf->is_active() )
                           // {
                           //    xx += ff->DOF_value( l, rn, ic ) ;
                           // }
                           // ff->set_DOF_value( l, rn, xx, ic ) ;
                           //---
                           ff->set_DOF_value( l, rn, PEL::bad_double(), ic ) ;
                        }
                     }
                  }
                  ff->set_node_active( rn ) ;
                  ff->set_node_inactive( cn ) ;
               }
            }
            PEL_ASSERT( rbf->check_child_parent_relationship() ) ;
            rbf->set_active() ;
         }
      }
   }
}

//---------------------------------------------------------------------------
void
PDE_Activator:: hierarchical_refine(
                             PDE_AdaptationRequest* adap,
                             bool make_refined_DOFs_bad_double ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: hierarchical_refine" ) ;

   adap->start_bf_to_refine() ;
   for( ; adap->valid_bf_to_refine() ; adap->go_next_bf_to_refine() )
   {
      PDE_BasisFunctionCell* cbf = adap->current_bf_to_refine() ;

      cbf->set_refined() ;

      for( size_t i=0 ; i<cbf->nb_childs() ; ++i )
      {
         PDE_BasisFunctionCell* rbf = cbf->child( i ) ;

         rbf->start_field_iterator() ;
         for( ; rbf->valid_field() ; rbf->go_next_field() )
         {
            PDE_DiscreteField* ff = rbf->field() ;
            if( !is_excluded( ff ) )
            {
               size_t rn = rbf->node_of_field() ;
               if( make_refined_DOFs_bad_double )
               {
                  for( size_t l=0 ; l<ff->storage_depth() ; ++l )
                  {
                     for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
                     {
                        //---> if the "refinement equation" is used
                        // ff->set_DOF_value( l, rn, 0.0, ic ) ;
                        //---
                        ff->set_DOF_value( l, rn, PEL::bad_double(), ic ) ;
                     }
                  }
               }
               ff->set_node_active( rn ) ;
            }
         }

         rbf->set_active() ;
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: quasi_hierarchical_unrefine( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: quasi_hierarchical_unrefine" ) ;

   if( VERB_LEVEL > 1 )
      PEL::out() << "*** PDE_Activator::"
                 << " quasi_hierarchical_unrefine ***" << endl << endl ;

   adap->start_bf_to_unrefine() ;
   for( ; adap->valid_bf_to_unrefine() ; adap->go_next_bf_to_unrefine() )
   {
      PDE_BasisFunctionCell* bf = adap->current_bf_to_unrefine() ;
      if( VERB_LEVEL > 1 )
      {
         PEL::out() << "activating: " ;
         print_bf( PEL::out(), bf ) ;
         PEL::out() << endl ;
      }

      PEL_ASSERT( bf->is_refined() ) ;
      PEL_ASSERT( !bf->is_active() ) ;
      bf->set_unrefined() ;
      bf->set_active() ;

      if( VERB_LEVEL > 1 ) PEL::out() << "   validating..." ;
      size_t nn = 0 ;
      bf->start_field_iterator() ;
      for( ; bf->valid_field() ; bf->go_next_field() )
      {
         PDE_DiscreteField* ff = bf->field() ;
         if( !is_excluded( ff ) )
         {
            if( VERB_LEVEL > 1 )
            {
               PEL::out() << "  \"" << ff->name() << "\":"
                          << bf->node_of_field() ;
            }
            ff->set_node_active( bf->node_of_field() ) ;
            nn++ ;
         }
      }
      if( VERB_LEVEL > 1 )
      {
         if( nn == 0 ) PEL::out() << "  nothing" ;
         PEL::out() << endl ;
      }

      for( size_t i=0 ; i<bf->nb_childs() ; ++i )
      {
         PDE_BasisFunctionCell* rbf = bf->child( i ) ;
         PEL_ASSERT( !rbf->is_refined() ) ;
         bool to_deactivate = true ;
         if( VERB_LEVEL > 1 )
         {
            PEL::out() << "   attempt to deactivate: " ;
            print_bf( PEL::out(), rbf ) ;
            PEL::out() << endl ;
         }
         for( size_t j=0 ; j<rbf->nb_parents() ; ++j )
         {
            PDE_BasisFunctionCell* p_rbf = rbf->parent( j ) ;
            if( VERB_LEVEL > 1 )
            {
               PEL::out() << "      checking parent: " ;
               print_bf( PEL::out(), p_rbf ) ;
               PEL::out() << endl ;
            }
//            if( p_rbf!=0 && p_rbf->is_refined() )
            if( p_rbf->is_refined() )
            {
               to_deactivate = false ;
               break ;
            }
         }
         if( to_deactivate )
         {
            if( VERB_LEVEL > 1 )
            {
               PEL::out() << "         desactivating: " ;
               print_bf( PEL::out(), rbf ) ;
               PEL::out() << endl ;
            }
            rbf->set_inactive() ;

            if( VERB_LEVEL > 1 ) PEL::out() << "         invalidating..." ;
            nn = 0 ;
            rbf->start_field_iterator() ;
            for( ; rbf->valid_field() ; rbf->go_next_field() )
            {
               PDE_DiscreteField* ff = rbf->field() ;
               if( !is_excluded( ff ) )
               {
                  if( VERB_LEVEL > 1 )
                  {
                     PEL::out() << "  \"" << ff->name() << "\":"
                                << rbf->node_of_field() ;
                  }
                  ff->set_node_inactive( rbf->node_of_field() ) ;
                  nn++ ;
               }
            }
            if( VERB_LEVEL > 1 )
            {
               if( nn == 0 ) PEL::out() << "  nothing" ;
               PEL::out() << endl ;
            }
         }
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_Activator:: hierarchical_unrefine( PDE_AdaptationRequest* adap )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: hierarchical_unrefine" ) ;

   if( VERB_LEVEL > 1 )
      PEL::out() << "*** PDE_Activator::"
                 << " hierarchical_unrefine ***" << endl << endl ;

   adap->start_bf_to_unrefine() ;
   for( ; adap->valid_bf_to_unrefine() ; adap->go_next_bf_to_unrefine() )
   {
      PDE_BasisFunctionCell* bf = adap->current_bf_to_unrefine() ;
      PEL_ASSERT( bf->is_refined() ) ;
      PEL_ASSERT( bf->is_active() ) ;
      bf->set_unrefined() ;

      for( size_t i=0 ; i<bf->nb_childs() ; ++i )
      {
         PDE_BasisFunctionCell* rbf = bf->child( i ) ;
         bool to_deactivate = true ;
         for( size_t j=0 ; j<rbf->nb_parents() ; ++j )
         {
            PDE_BasisFunctionCell* p_rbf = rbf->parent( j ) ;
            if( p_rbf!=0 && p_rbf->is_refined() )
            {
               to_deactivate = false ;
               break ;
            }
         }
         if( to_deactivate )
         {
            rbf->set_inactive() ;
            rbf->start_field_iterator() ;
            for( ; rbf->valid_field() ; rbf->go_next_field() )
            {
               PDE_DiscreteField* ff = rbf->field() ;
               if( !is_excluded( ff ) )
               {
                  ff->set_node_inactive( rbf->node_of_field() ) ;
               }
            }
         }
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_Activator:: print_bf( std::ostream& os,
                                      PDE_BasisFunctionCell* bf )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Activator:: print_bf" ) ;

   if( bf->is_active() )
      os << "  A" ;
   else
      os << "  I" ;
   if( bf->is_refined() )
      os << "R" ;
   else
      os << " " ;

   bf->start_field_iterator() ;
   PEL_ASSERT( bf->valid_field() ) ;
   os << "  \"" << bf->field()->name() << "\":" << bf->node_of_field() ;
   if( bf->field()->node_is_active( bf->node_of_field() ) )
      os << "(V)" ;
   else
      os << "(I)" ;
   for( bf->go_next_field() ; bf->valid_field() ; bf->go_next_field() )
   {
      os << "  \"" << bf->field()->name() << "\":" << bf->node_of_field() ;
      if( bf->field()->node_is_active( bf->node_of_field() ) )
         os << "(V)" ;
      else
         os << "(I)" ;
   }

   os << "  cells:" << setw(3) << bf->cell(0)->id_number() ;
   for( size_t i=1 ; i<bf->nb_cells() ; ++i )
   {
      os << "," << setw(3) << bf->cell(i)->id_number() ;
   }
}

