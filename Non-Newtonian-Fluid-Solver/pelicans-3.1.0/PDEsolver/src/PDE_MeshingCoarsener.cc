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

#include <PDE_MeshingCoarsener.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_Vector.hh>

#include <PDE_Activator.hh>
#include <PDE_AdaptationRequestFromLevels.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridFE.hh>
#include <PDE_LocalFEcell.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ;
using std::endl ;
using std::setw ;

//----------------------------------------------------------------------
PDE_MeshingCoarsener:: PDE_MeshingCoarsener( PEL_Object* a_owner,
                                             PDE_DomainAndFields const* dom,
                                             PDE_GridFE* a_grid,
                                             PDE_Activator* activator,
                                             size_t verbose_level )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DOM( dom )
   , GRID( a_grid )
   , ACTIVATOR( activator )
   , BF_LISTS( 0 )
   , LEVEL_MAX( PEL::bad_index() )
   , FINE_LEVEL( PEL::bad_index() )
   , VERB_LEVEL( verbose_level )
{
}

//----------------------------------------------------------------------
PDE_MeshingCoarsener:: ~PDE_MeshingCoarsener( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_MeshingCoarsener:: reset( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshingCoarsener:: reset" ) ;

   LEVEL_MAX = PEL::bad_index() ;
   FINE_LEVEL = PEL::bad_index() ;
   destroy_possession( BF_LISTS ) ;
   BF_LISTS = 0 ;

   PEL_CHECK_POST( nb_levels() == PEL::bad_index() ) ;
   PEL_CHECK_POST( current_fine_level() == PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void
PDE_MeshingCoarsener:: prepare_for_coarsening( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshingCoarsener:: prepare_for_coarsening" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   if( BF_LISTS != 0 )
   {
      std::ostringstream mesg ;
      mesg << "*** PDE_MeshingCoarsener error" << endl ;
      mesg << "    the member function \"reset\"" << endl
           << "       should be called priot to" << endl
           << "    the member function \"prepare_for_coarsening\"" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   //????? la determination de LEVEL_MAX devrait se faire à partir
   //????? des fonctions de base actives
   LEVEL_MAX = 0 ;
   PDE_LocalFEcell* cFE = DOM->create_LocalFEcell( 0 ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      if( cFE->refinement_level() > LEVEL_MAX )
         LEVEL_MAX = cFE->refinement_level() ;
   }
   cFE->destroy() ;
   
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   LEVEL_MAX = com->max( LEVEL_MAX ) ;

   FINE_LEVEL = LEVEL_MAX ;

   BF_LISTS = PEL_Vector::create( this, LEVEL_MAX ) ;

   PEL_CHECK_POST( nb_levels() != PEL::bad_index() ) ;
   PEL_CHECK_POST( current_fine_level() == nb_levels()-1 ) ;
}

//-----------------------------------------------------------------------
size_t
PDE_MeshingCoarsener:: nb_levels( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshingCoarsener:: nb_levels" ) ;

   size_t result = ( LEVEL_MAX == PEL::bad_index() ?
                     PEL::bad_index() : LEVEL_MAX+1 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
size_t
PDE_MeshingCoarsener:: current_fine_level( void ) const
//-----------------------------------------------------------------------
{
   return( FINE_LEVEL ) ;
}

//-----------------------------------------------------------------------
void
PDE_MeshingCoarsener:: do_one_coarsening( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshingCoarsener:: do_one_coarsening" ) ;
   PEL_CHECK_PRE( current_fine_level() != PEL::bad_index() ) ;
   PEL_CHECK_PRE( current_fine_level() != 0 ) ;
   PEL_SAVEOLD( size_t, current_fine_level, current_fine_level() ) ;

   if( VERB_LEVEL > 1 )
   {
      PEL::out() << "PDE_MeshingCoarsener - coarsening : " << FINE_LEVEL
                 << "->" << FINE_LEVEL-1 << std::endl << std::endl ;
   }

   PDE_AdaptationRequestFromLevels* adap =
   PDE_AdaptationRequestFromLevels::create_for_unrefinement( 0,
                                                             FINE_LEVEL,
                                                             VERB_LEVEL ) ;

   PEL_ListIdentity* bfs =
           static_cast< PEL_ListIdentity* >( BF_LISTS->at( FINE_LEVEL-1 ) ) ;
   if( bfs == 0 )
   {
      bfs = PEL_ListIdentity::create( this ) ;

      adap->apply_criterion_to_infer_bfs_to_unrefine( GRID->cells() ) ;
      adap->extend_bfs_to_unrefine( bfs ) ;
      adap->read_bfs_to_unrefine( bfs ) ;
      adap->build_info_on_cells_for_unrefinement() ;

      BF_LISTS->set_at( FINE_LEVEL-1, bfs ) ;
   }
   else
   {
      adap->read_bfs_to_unrefine( bfs ) ;
      adap->build_info_on_cells_for_unrefinement() ;
   }

   PEL_ASSERT( adap->something_to_unrefine() ) ;

   ACTIVATOR->unrefine( adap ) ;
   ACTIVATOR->unsplit_meshes_1( adap ) ;

   adap->destroy() ;

   FINE_LEVEL-- ;

   PEL_CHECK_POST( current_fine_level() == OLD(current_fine_level) - 1 ) ;
}

//-----------------------------------------------------------------------
void
PDE_MeshingCoarsener:: do_one_uncoarsening( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshingCoarsener:: do_one_uncoarsening" ) ;
   PEL_CHECK_PRE( current_fine_level() != nb_levels()-1 ) ;
   PEL_SAVEOLD( size_t, current_fine_level, current_fine_level() ) ;

   if( VERB_LEVEL > 1 )
   {
      PEL::out() << "PDE_MeshingCoarsener - uncoarsening : " << FINE_LEVEL
                 << "->" << FINE_LEVEL+1 << std::endl << std::endl ;
   }

   PEL_ListIdentity const* bfs =
           static_cast< PEL_ListIdentity* >( BF_LISTS->at( FINE_LEVEL ) ) ;
   PEL_ASSERT( bfs != 0 ) ;

   PDE_AdaptationRequestFromLevels* adap =
   PDE_AdaptationRequestFromLevels::create_for_refinement( 0,
                                                           FINE_LEVEL,
                                                           VERB_LEVEL ) ;
   adap->read_bfs_to_refine( bfs ) ;
   adap->build_info_on_cells_for_refinement() ;

   PEL_ASSERT( adap->something_to_refine() ) ;

   ACTIVATOR->split_meshes( adap ) ;
   bool make_refined_DOFs_bad_double = false ;
   ACTIVATOR->refine( adap, make_refined_DOFs_bad_double ) ;
   ACTIVATOR->unsplit_meshes_1( adap ) ;

   adap->destroy() ;

   FINE_LEVEL++ ;

   PEL_CHECK_POST( current_fine_level() == OLD(current_fine_level) + 1 ) ;
}
