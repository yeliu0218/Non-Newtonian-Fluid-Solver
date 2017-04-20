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

#include <PDE_AdaptationRequest.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <PDE_AdaptationIndicator.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_CrossProcessBFNumbering.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ReferenceElementRefiner.hh>
#include <PDE_SetOfBasisFunctions.hh>

#include <algorithm>
#include <functional>
#include <iostream>

using std::vector ;
using std::endl ;

//---------------------------------------------------------------------------
PDE_AdaptationRequest:: PDE_AdaptationRequest( PEL_Object* a_owner,
                                               size_t highest_refinement_level,
                                               size_t verbose_level )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , iM_R( PEL::bad_index() )
   , iL_R( PEL::bad_index() )
   , iM_U( PEL::bad_index() )
   , iL_U( PEL::bad_index() )
   , MAX_REF_LEVEL( highest_refinement_level )
   , VERB_LEVEL( verbose_level )
   , BFS_DUM_R( PEL_ListIdentity::create( this ) )
   , BFS_DUM_U( PEL_ListIdentity::create( this ) )
   , IT_R_CRITERIA( PEL_ListIterator::create( this ) )
   , BFS_R_CRITERIA( 0 )
   , IT_U_CRITERIA( PEL_ListIterator::create( this ) )
   , BFS_U_CRITERIA( 0 )
   , IT_R( PEL_ListIterator::create( this ) )
   , BFS_R( 0 )
   , IT_U( PEL_ListIterator::create( this ) )
   , BFS_U( 0 )
{
}

//---------------------------------------------------------------------------
PDE_AdaptationRequest:: ~PDE_AdaptationRequest( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: apply_criterion_to_infer_bfs_to_refine(
                                                    PEL_Vector const* cells )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: apply_criterion_to_infer_bfs_to_refine" ) ;
   PEL_CHECK_PRE( cells != 0 ) ;

   PEL_ASSERT( BFS_R_CRITERIA == 0 ) ;//tjrs apres un create
   BFS_R_CRITERIA = PEL_ListIdentity::create( this ) ;

   if( VERB_LEVEL > 0 )
   {
      PEL::out() << "*** PDE_AdaptationRequest:: "
                    "apply_criterion_to_infer_bfs_to_refine"
                 << std::endl << std::endl ;
   }
      
   size_t nb_cells = cells->count() ;

   for( size_t ic=0 ; ic<nb_cells ; ++ic )
   {
      PDE_CellFE const* cell = static_cast<PDE_CellFE*>( cells->at( ic ) ) ;

      size_t level = cell->refinement_level() ;
      PEL_ASSERT( level <= MAX_REF_LEVEL ) ;  //?????????
      if( level == MAX_REF_LEVEL ) continue ; // <-------

      bool print_cell = true ;
      for( size_t ee=0 ; ee<cell->nb_reference_elements() ; ++ee )
      {
         PDE_ReferenceElement const* elm = cell->reference_element( ee ) ;
         for( size_t ln=0 ; ln<elm->nb_nodes() ; ++ln )
         {
            PDE_BasisFunctionCell* bf = cell->basis_function( ee, ln ) ;
            if( bf != 0 )
            {
               if( bf->is_active() && !bf->is_refined() &&
                   to_be_refined( cell, bf, elm, ln ) )
               {
                  if( VERB_LEVEL > 0 )
                  {
                     if( print_cell )
                     {
                        PEL::out() << "cell " << cell->id_number() ;
                        PDE_CellFE const* pcell = cell->parent() ;
                        if( pcell!= 0 )
                           PEL::out() << " (parent : " << pcell->id_number()
                                      << ")" ;
                        PEL::out() << std::endl ;
                        print_cell = false ;
                     }
                     PEL::out() << "   local bf " << ln << " must be refined"
                                << " (" << elm->name() << ")" << std::endl ;
                  }

                  PEL_ASSERT( bf->refinement_level() == level ) ;
                  BFS_R_CRITERIA->extend( bf ) ;
               }
            }
         }
      }
   }
   IT_R_CRITERIA->re_initialize( BFS_R_CRITERIA ) ;
   
   if( VERB_LEVEL > 0 && BFS_R_CRITERIA->index_limit() != 0 ) 
   {
      PEL::out() << std::endl ;
   }
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: something_to_refine( void ) const
//---------------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( true ) ;

   bool something_to_ref = ( BFS_R != 0 ) &&
                           ( BFS_R->index_limit() ) ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   bool result = com->boolean_or( something_to_ref ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: extend_bfs_to_refine( PEL_ListIdentity* bfs_to_refine,
                                              PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: extend_bfs_to_refine" ) ;
   PEL_CHECK_COLLECTIVE( bf_set != 0 && bf_set->is_distributed() ) ;
   PEL_CHECK_PRE( bfs_to_refine != 0 ) ;
   PEL_CHECK_PRE( bfs_to_refine->index_limit() == 0 ) ;

   if( VERB_LEVEL > 0 )
   {
      PEL::out() << "***  PDE_AdaptationRequest:: extend_bfs_to_refine"
                 << std::endl << std::endl ;
   }

   IT_R_CRITERIA->start() ;
   for( ; IT_R_CRITERIA->is_valid() ; IT_R_CRITERIA->go_next() )
   {
      PDE_BasisFunctionCell* bf =
               static_cast< PDE_BasisFunctionCell* >( IT_R_CRITERIA->item() ) ;
      extend_with_ascendants( bf, bfs_to_refine, 6 ) ;
   }

   if( VERB_LEVEL > 0 && BFS_R_CRITERIA->index_limit() != 0 ) 
   {
      PEL::out() << std::endl ;
   }
   
   if( bf_set != 0 && bf_set->is_distributed() )
      synchronize_selected_bfs_to_refine( bfs_to_refine, bf_set ) ;

   if( VERB_LEVEL > 0 ) PEL::out() << std::endl ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: compute_bfs_to_refine( PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: compute_next_bfs_to_refine" ) ;
   PEL_CHECK_COLLECTIVE( bf_set != 0 && bf_set->is_distributed() ) ;

   BFS_DUM_R->clear() ;

   extend_bfs_to_refine( BFS_DUM_R, bf_set ) ;

   read_bfs_to_refine( BFS_DUM_R ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: synchronize_selected_bfs_to_refine(
                                              PEL_ListIdentity* bfs,
                                              PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: synchronize_selected_bfs_to_refine" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( bf_set != 0 ) ;
   PEL_CHECK_PRE( bfs != 0 ) ;

   // La strategie de communication est surement a ameliorer
   // Objectif : on ajoute une fonction de base dans la liste de chaque process
   // qui la voit des que cette fonction de base est dans la liste
   // d'un des processes.

   if( VERB_LEVEL > 0 )
   {
      PEL::out() << "***  PDE_AdaptationRequest::"
                 << " synchronize_selected_bfs_to_refine"
                 << std::endl << std::endl ;
   }
   
   PDE_CrossProcessBFNumbering* bf_numbering =
                                bf_set->cross_process_numbering() ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;

   size_t_vector bfs_global_number( 0 ) ;
   IT_R->re_initialize( bfs ) ;
   for( ; IT_R->is_valid() ; IT_R->go_next() )
   {
      PDE_BasisFunctionCell const* bf =
                       static_cast< PDE_BasisFunctionCell* >( IT_R->item() ) ;
      size_t id = bf->id_number() ;
      size_t e = bf->ref_elts_grp_index() ;
      bfs_global_number.append(
                        bf_numbering->global_basis_function_index( id, e ) ) ;
   }

   size_t_vector bfs_number_dum( 0 ) ;
   for( size_t r=0 ; r < com->nb_ranks() ; ++r )
   {
      if( com->rank() == r )
      {
         bfs_number_dum = bfs_global_number ;
      }
      com->broadcast( bfs_number_dum, r ) ;
      for( size_t i = 0 ; i < bfs_number_dum.size() ; ++i )
      {
         size_t local_e =  PEL::bad_index() ;
         size_t local_id =
                bf_numbering->local_basis_function_index( bfs_number_dum(i),
                                                          local_e ) ;
         if( local_id != PEL::bad_index() )
         {
            PDE_BasisFunctionCell* bf = 
               static_cast< PDE_BasisFunctionCell* >( 
                                        bf_set->item( local_id, local_e ) ) ;
            if( !bfs->has( bf ) )
            {
               bfs->append( bf ) ;
               if (VERB_LEVEL > 0 ) 
               { 
                  bf->print_2( PEL::out(), 3 ) ;
                  PEL::out() << " added" << std::endl ;
               }
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: read_bfs_to_refine(
                                    PEL_ListIdentity const* bfs_to_refine )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: read_bfs_to_refine" ) ;
   PEL_CHECK_PRE( bfs_to_refine != 0 ) ;

   BFS_R = bfs_to_refine ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: apply_criterion_to_infer_bfs_to_unrefine(
                                                  PEL_Vector const* cells )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: apply_criterion_to_infer_bfs_to_unrefine" ) ;
   PEL_CHECK_PRE( cells != 0 ) ;

   PEL_ASSERT( BFS_U_CRITERIA == 0 ) ;//tjrs apres un create
   BFS_U_CRITERIA = PEL_ListIdentity::create( this ) ;

   if( VERB_LEVEL > 0 )
      PEL::out() << "*** PDE_AdaptationRequest:: "
                    "apply_criterion_to_infer_bfs_to_unrefine"
                 << std::endl << std::endl ;

   size_t nb_cells = cells->count() ;

   for( size_t ic=0 ; ic<nb_cells ; ++ic )
   {
      PDE_CellFE const* cell = static_cast<PDE_CellFE*>( cells->at( ic ) ) ;
      bool print_cell = true ;

      for( size_t ee=0 ; ee<cell->nb_reference_elements() ; ++ee )
      {
         PDE_ReferenceElement const* elm = cell->reference_element( ee ) ;
         for( size_t ln=0 ; ln<elm->nb_nodes() ; ++ln )
         {
            PDE_BasisFunctionCell* bf = cell->basis_function( ee, ln ) ;
            if( bf != 0 )
            {
               if( bf->is_refined() &&
                   to_be_unrefined( cell, bf, elm, ln ) )
               {
                  PEL_ASSERT( !cell->is_active() ) ;
                  if( VERB_LEVEL > 0 )
                  {
                     if( print_cell )
                     {
                        PEL::out() << "cell " << cell->id_number() ;
                        PDE_CellFE const* pcell = cell->parent() ;
                        if( pcell!= 0 )
                           PEL::out() << " (parent : " << pcell->id_number()
                                      << ")" ;
                        PEL::out() << endl ;
                        print_cell = false ;
                     }
                     PEL::out() << "   local bf " << ln
                                << " may be unrefined ("
                                << elm->name() << ")" << std::endl ; ;
                  }
                  BFS_U_CRITERIA->extend( bf ) ;
               }
            }
         }
      }
   }
   IT_U_CRITERIA->re_initialize( BFS_U_CRITERIA ) ;
   
   if( VERB_LEVEL > 0 && BFS_U_CRITERIA->index_limit() != 0 ) 
   {
      PEL::out() << std::endl ;
   }
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: extend_bfs_to_unrefine(
                                            PEL_ListIdentity* bfs_to_unrefine,
                                            PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: extend_bfs_to_unrefine" ) ;
   PEL_CHECK_COLLECTIVE( bf_set != 0 && bf_set->is_distributed() ) ;
   PEL_CHECK_PRE( bfs_to_unrefine != 0 ) ;
   PEL_CHECK_PRE( bfs_to_unrefine->index_limit() == 0 ) ;

   if( VERB_LEVEL > 0 )
   {
      PEL::out() << "***  PDE_AdaptationRequest:: extend_bfs_to_unrefine"
                 << std::endl << std::endl ;
   }

   IT_U_CRITERIA->start() ;
   for( ; IT_U_CRITERIA->is_valid() ; IT_U_CRITERIA->go_next() )
   {
      PDE_BasisFunctionCell* bf =
               static_cast< PDE_BasisFunctionCell* >( IT_U_CRITERIA->item() ) ;
      consider_unrefining_bf( bf, bfs_to_unrefine ) ;
   }

   if( VERB_LEVEL > 0 && BFS_U_CRITERIA->index_limit() != 0 ) 
   {
      PEL::out() << endl ;
   }
   
   if( bf_set != 0 && bf_set->is_distributed() )
         synchronize_selected_bfs_to_unrefine( bfs_to_unrefine, bf_set ) ;

   if( VERB_LEVEL > 0 ) PEL::out() << endl ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: compute_bfs_to_unrefine(
                                              PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: compute_bfs_to_unrefine" ) ;
   PEL_CHECK_COLLECTIVE( bf_set != 0 && bf_set->is_distributed() ) ;

   BFS_DUM_U->clear() ;

   extend_bfs_to_unrefine( BFS_DUM_U, bf_set ) ;
 
   read_bfs_to_unrefine( BFS_DUM_U ) ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: something_to_unrefine( void ) const
//---------------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( true ) ;

   bool something_to_unref = ( BFS_U != 0 ) &&
                             ( BFS_U->index_limit() != 0 ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   bool result = com->boolean_or( something_to_unref ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: synchronize_selected_bfs_to_unrefine(
                                              PEL_ListIdentity* bfs,
                                              PDE_SetOfBasisFunctions* bf_set )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: synchronize_selected_bfs_to_unrefine" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;
   PEL_CHECK_PRE( bf_set != 0 ) ;
   PEL_CHECK_PRE( bfs != 0 ) ;

   // La strategie de communication est surement a ameliorer
   // Objectif : on enleve une fonction de base des listes de chaque process
   // si il existe un process voyant cette fonction de base mais qui ne la
   // contient pas dans sa liste

   if( VERB_LEVEL > 0 )
   {
      PEL::out() << "***  PDE_AdaptationRequest:: "
                 << " synchronize_selected_bfs_to_unrefine"
                 << std::endl << std::endl ;
   }
   
   PDE_CrossProcessBFNumbering* bf_numbering =
                                bf_set->cross_process_numbering() ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;

   size_t_vector bfs_global_number( 0 ) ;
   IT_R->re_initialize( bfs ) ;
   for( ; IT_R->is_valid() ; IT_R->go_next() )
   {
      PDE_BasisFunctionCell const* bf =
                       static_cast< PDE_BasisFunctionCell* >( IT_R->item() ) ;
      size_t id = bf->id_number() ;
      size_t e = bf->ref_elts_grp_index() ;
      bfs_global_number.append(
                        bf_numbering->global_basis_function_index( id, e ) ) ;
   }

   size_t_vector bfs_number_dum( 0 ) ;
   for( size_t r=0 ; r < com->nb_ranks() ; ++r )
   {
      if( com->rank() == r )
      {
         bfs_number_dum = bfs_global_number ;
      }
      com->broadcast( bfs_number_dum, r ) ;

      for( size_t i = 0 ; i < bfs_number_dum.size() ; ++i )
      {
         size_t local_e = PEL::bad_index() ;
         size_t local_id =
               bf_numbering->local_basis_function_index( bfs_number_dum(i),
                                                         local_e ) ;
         PDE_BasisFunctionCell const* bf = 0 ;
         bool on_process = ( local_id != PEL::bad_index() ) ;
         bool is_not_in_list = false ;
         if( on_process )
         {
            bf = static_cast< PDE_BasisFunctionCell* >(  
                                        bf_set->item( local_id, local_e ) ) ;
            is_not_in_list = !( bfs->has( bf ) ) ;
         }
         bool to_remove = com->boolean_or( is_not_in_list ) ;
         if( to_remove && !is_not_in_list && on_process )
         {
            PEL_ASSERT( bf != 0 ) ;
            bfs->remove( bf ) ;
            if (VERB_LEVEL > 0 ) 
            {
               bf->print_2( PEL::out(), 3 ) ;
               PEL::out() << " removed" << std::endl ;
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: read_bfs_to_unrefine(
                                    PEL_ListIdentity const* bfs_to_unrefine )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: read_bfs_to_refine" ) ;
   PEL_CHECK_PRE( bfs_to_unrefine != 0 ) ;

   BFS_U = bfs_to_unrefine ;
}

//----------------------------------------------------------------------------
void
PDE_AdaptationRequest:: build_info_on_cells_for_refinement( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: build_info_on_cells_for_refinement" ) ;

   PEL_ASSERT( BFS_R != 0 ) ;
   for( size_t l=0 ; l<CELLS_TO_REFINE.size() ; ++l )
   {
      CELLS_TO_REFINE[ l ].clear() ;
   }
   IT_R->re_initialize( BFS_R ) ;
   for( ; IT_R->is_valid() ; IT_R->go_next() )
   {
      PEL_CHECK( dynamic_cast<PDE_BasisFunctionCell*>(IT_R->item())!= 0) ;
      PDE_BasisFunctionCell const* bf =
                  static_cast< PDE_BasisFunctionCell* >( IT_R->item() ) ;
      for( size_t i=0 ; i<bf->nb_cells() ; ++i )
      {
         PDE_CellFE* cell = bf->cell( i ) ;
         size_t ll = cell->refinement_level() ;
         if( ll+1 > CELLS_TO_REFINE.size() ) CELLS_TO_REFINE.resize( ll+1 ) ;
         bool found = ( std::find( CELLS_TO_REFINE[ ll ].begin(),
                                   CELLS_TO_REFINE[ ll ].end(), cell )
                        != CELLS_TO_REFINE[ ll ].end() ) ;
         if( !found ) CELLS_TO_REFINE[ ll ].push_back( cell ) ;
      }
   }

   if( VERB_LEVEL > 0 ) PEL::out() << endl ;
}

//----------------------------------------------------------------------------
void
PDE_AdaptationRequest:: build_info_on_cells_for_unrefinement( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: build_unrefinement_info_on_cells" ) ;

   PEL_ASSERT( BFS_U != 0 ) ;

   for( size_t l=0 ; l<CELLS_TO_UNREFINE.size() ; ++l )
   {
      CELLS_TO_UNREFINE[ l ].clear() ;
   }
   IT_U->re_initialize( BFS_U ) ;
   for( ; IT_U->is_valid() ; IT_U->go_next() )
   {
      PEL_CHECK( dynamic_cast<PDE_BasisFunctionCell*>(IT_U->item())!= 0) ;
      PDE_BasisFunctionCell const* bf =
         static_cast< PDE_BasisFunctionCell* >( IT_U->item() ) ;
      for( size_t i=0 ; i<bf->nb_cells() ; ++i )
      {
         PDE_CellFE* cell = bf->cell( i ) ;
         size_t ll = cell->refinement_level() ;
         if( ll+1 > CELLS_TO_UNREFINE.size() )
         {
            CELLS_TO_UNREFINE.resize( ll+1 ) ;
         }
         bool found = ( std::find( CELLS_TO_UNREFINE[ ll ].begin(),
                                   CELLS_TO_UNREFINE[ ll ].end(), cell )
                        != CELLS_TO_UNREFINE[ ll ].end() ) ;
         if( !found ) CELLS_TO_UNREFINE[ ll ].push_back( cell ) ;
      }
   }

   if( VERB_LEVEL > 0 ) PEL::out() << endl ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: start_bf_to_refine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: start_bf_to_refine" ) ;
   if( BFS_R != 0 )
   {
      IT_R->re_initialize( BFS_R ) ;
   }
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: valid_bf_to_refine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: valid_bf_to_refine" ) ;

   bool result = ( (BFS_R != 0) && IT_R->is_valid() ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: go_next_bf_to_refine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: go_next_bf_to_unrefine" ) ;
   PEL_CHECK_PRE( valid_bf_to_refine() ) ;

   IT_R->go_next() ;
}

//---------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_AdaptationRequest:: current_bf_to_refine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: current_bf_to_refine" ) ;
   PEL_CHECK_PRE( valid_bf_to_refine() ) ;

   PEL_CHECK( dynamic_cast<PDE_BasisFunctionCell*>( IT_R->item() ) != 0 ) ;
   PDE_BasisFunctionCell* result =
            static_cast< PDE_BasisFunctionCell* >( IT_R->item() ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: start_bf_to_unrefine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: start_bf_to_unrefine" ) ;
   if( BFS_U != 0 )
   {
      IT_U->re_initialize( BFS_U ) ;
   }
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: valid_bf_to_unrefine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: valid_bf_to_unrefine" ) ;

   bool result = ( (BFS_U != 0) && IT_U->is_valid() ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: go_next_bf_to_unrefine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: go_next_bf_to_unrefine" ) ;
   PEL_CHECK_PRE( valid_bf_to_unrefine() ) ;

   IT_U->go_next() ;
}

//---------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_AdaptationRequest:: current_bf_to_unrefine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: current_bf_to_unrefine" ) ;
   PEL_CHECK_PRE( valid_bf_to_unrefine() ) ;

   PEL_CHECK( dynamic_cast<PDE_BasisFunctionCell*>( IT_U->item() ) != 0 ) ;
   PDE_BasisFunctionCell* result =
            static_cast< PDE_BasisFunctionCell* >( IT_U->item() ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: start_cell_to_refine( void )
//---------------------------------------------------------------------------
{
   iL_R = 0 ;
   while( iL_R <  CELLS_TO_REFINE.size() &&
          CELLS_TO_REFINE[ iL_R ].size() == 0 )
   {
      ++iL_R ;
   }
   iM_R = 0 ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: valid_cell_to_refine( void ) const
//---------------------------------------------------------------------------
{
   bool result = ( iL_R < CELLS_TO_REFINE.size() ) ;
   PEL_ASSERT( IMPLIES( iL_R < CELLS_TO_REFINE.size(),
                        iM_R < CELLS_TO_REFINE[ iL_R ].size() ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: go_next_cell_to_refine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: go_next_cell_to_refine" ) ;
   PEL_CHECK_PRE( valid_cell_to_refine() ) ;

   ++iM_R ;
   if( iM_R ==  CELLS_TO_REFINE[ iL_R ].size() )
   {
      do
      {
         ++iL_R ;
      } while( iL_R < CELLS_TO_REFINE.size() &&
               CELLS_TO_REFINE[ iL_R ].size() == 0 ) ;
      iM_R = 0 ;
   }
   PEL_ASSERT( IMPLIES( iL_R < CELLS_TO_REFINE.size(),
                        CELLS_TO_REFINE[ iL_R ].size() != 0 ) ) ;
   PEL_ASSERT(  IMPLIES( iL_R < CELLS_TO_REFINE.size(),
                         iM_R < CELLS_TO_REFINE[ iL_R ].size() ) ) ;
}

//---------------------------------------------------------------------------
PDE_CellFE*
PDE_AdaptationRequest:: current_cell_to_refine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: current_cell_to_refine" ) ;
   PEL_CHECK_PRE( valid_cell_to_refine() ) ;

   PDE_CellFE* result = CELLS_TO_REFINE[iL_R][iM_R] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: start_cell_to_unrefine( void )
//---------------------------------------------------------------------------
{
   iL_U = 0 ;
   while( iL_U < CELLS_TO_UNREFINE.size() &&
          CELLS_TO_UNREFINE[ iL_U ].size() == 0 )
   {
      ++iL_U ;
   }
   iM_U = 0 ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: valid_cell_to_unrefine( void ) const
//---------------------------------------------------------------------------
{
   bool result = ( iL_U < CELLS_TO_UNREFINE.size() ) ;
   PEL_ASSERT( IMPLIES( iL_U < CELLS_TO_UNREFINE.size(),
                        iM_U < CELLS_TO_UNREFINE[ iL_U ].size() ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdaptationRequest:: go_next_cell_to_unrefine( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: go_next_cell_to_unrefine" ) ;
   PEL_CHECK_PRE( valid_cell_to_unrefine() ) ;

   ++iM_U ;
   if( iM_U ==  CELLS_TO_UNREFINE[ iL_U ].size() )
   {
      do
      {
         ++iL_U ;
      } while( iL_U < CELLS_TO_UNREFINE.size() &&
               CELLS_TO_UNREFINE[ iL_U ].size() == 0 ) ;
      iM_U = 0 ;
   }
   PEL_ASSERT( IMPLIES( iL_U < CELLS_TO_UNREFINE.size(),
                        CELLS_TO_UNREFINE[ iL_U ].size() != 0 ) ) ;
   PEL_ASSERT(  IMPLIES( iL_U < CELLS_TO_UNREFINE.size(),
                         iM_U < CELLS_TO_UNREFINE[ iL_U ].size() ) ) ;
}

//---------------------------------------------------------------------------
PDE_CellFE*
PDE_AdaptationRequest:: current_cell_to_unrefine( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: current_cell_to_unrefine" ) ;
   PEL_CHECK_PRE( valid_cell_to_unrefine() ) ;

   PDE_CellFE* result = CELLS_TO_UNREFINE[iL_U][iM_U] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_AdaptationRequest:: extend_with_ascendants(
                                    PDE_BasisFunctionCell* bf,
                                    PEL_ListIdentity* bfs_to_refine,
                                    size_t indent ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: extend_with_parents" ) ;

   bfs_to_refine->extend( bf ) ;

   std::string bl( indent, ' ' ) ;
   if( VERB_LEVEL > 0 )
   {
      PEL::out() << bl << "checking ascendants of: " ;
      bf->print_2( PEL::out(), 0 ) ; PEL::out() << endl ;
   }

   for( size_t ip=0 ; ip<bf->nb_ascendants() ; ++ip )
   {
      PDE_BasisFunctionCell* pbf = bf->ascendant( ip ) ;
      if( pbf->is_active() && !pbf->is_refined() )
      {
         if( ( VERB_LEVEL > 0 ) && ( !bfs_to_refine->has( pbf ) ) )
         {
            pbf->print_2( PEL::out(), indent ) ;
            PEL::out() << " ok" << std::endl ;
         }
         PEL_ASSERT( !pbf->is_refined() ) ;
         extend_with_ascendants( pbf, bfs_to_refine, indent+3 ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_AdaptationRequest:: consider_unrefining_bf(
                                 PDE_BasisFunctionCell* bf,
                                 PEL_ListIdentity* bfs_to_unrefine )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequest:: consider_unrefining_bf" ) ;

   if (VERB_LEVEL > 0 ) bf->print_2( PEL::out(), 3 ) ;

   bool may_be_unrefined = true ;
   for( size_t ir=0 ; ir<bf->nb_childs() ; ++ir )
   {
      PDE_BasisFunctionCell* rbf = bf->child( ir ) ;
      if( rbf->is_refined() )
      {
         if( VERB_LEVEL > 0 ) PEL::out() << " no" ;
         may_be_unrefined = false ;
         break ;
      }
   }

   if( may_be_unrefined )
   {
      //??? si les descendants contiennent les enfants, la partie du dessus
      //??? est inutiles (faite deux fois)
      for( size_t ir=0 ; ir<bf->nb_descendants() ; ++ir )
      {
         PDE_BasisFunctionCell* rbf = bf->descendant( ir ) ;
         if( rbf->is_refined() )
         {
            if (VERB_LEVEL > 0) PEL::out() << " no" ;
            may_be_unrefined = false ;
            break ;
         }

         if( BFS_R!=0 && BFS_R->has( rbf ) )
         {
            if (VERB_LEVEL > 0) PEL::out() << " no2" ;
            may_be_unrefined = false ;
            break ;
         }
      }
   }

   if( may_be_unrefined )
   {
      bfs_to_unrefine->extend( bf ) ;
      if( VERB_LEVEL > 0 ) PEL::out() << " ok" ;
   }
   if( VERB_LEVEL > 0 ) PEL::out() << std::endl ;
}

//----------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: to_be_refined_PRE( PDE_CellFE const* cell,
                                           PDE_BasisFunctionCell const* bf,
                                           PDE_ReferenceElement const* elm,
                                           size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( cell != 0 ) ;
   PEL_ASSERT( bf != 0 ) ;
   PEL_ASSERT( bf->is_active() ) ;
   PEL_ASSERT( !bf->is_refined() ) ;
   PEL_ASSERT( elm != 0 ) ;
   PEL_ASSERT( local_node < elm->nb_nodes() ) ;
   PEL_ASSERT( bf == cell->basis_function(
                           cell->index_of_element( elm ), local_node ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
PDE_AdaptationRequest:: to_be_unrefined_PRE( PDE_CellFE const* cell,
                                             PDE_BasisFunctionCell const* bf,
                                             PDE_ReferenceElement const* elm,
                                             size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( cell != 0 ) ;
   PEL_ASSERT( bf != 0 ) ;
   PEL_ASSERT( bf->is_refined() ) ;
   PEL_ASSERT( elm != 0 ) ;
   PEL_ASSERT( local_node < elm->nb_nodes() ) ;
   PEL_ASSERT( bf == cell->basis_function(
                           cell->index_of_element( elm ), local_node ) ) ;
   return( true ) ;
}
