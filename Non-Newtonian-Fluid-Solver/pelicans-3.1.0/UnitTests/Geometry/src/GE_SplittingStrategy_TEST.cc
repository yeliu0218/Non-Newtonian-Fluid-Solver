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

#include <GE_SplittingStrategy_TEST.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <boolVector.hh>
#include <size_t_array2D.hh>

#include <GE_Meshing.hh>
#include <GE_SplittingStrategy.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::ostringstream ;

//---------------------------------------------------------------------------
GE_SplittingStrategy_TEST*
GE_SplittingStrategy_TEST:: REGISTRATOR = new GE_SplittingStrategy_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
GE_SplittingStrategy_TEST:: GE_SplittingStrategy_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_SplittingStrategy", 
                     "GE_SplittingStrategy_TEST" )
{
}

//---------------------------------------------------------------------------
GE_SplittingStrategy_TEST:: ~GE_SplittingStrategy_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_SplittingStrategy_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy_TEST:: process_one_test" ) ;
   
   out() << "|" << endl << "| " << exp->name() << endl ;

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "GE_Meshing" ) ;
   GE_Meshing* mm = GE_Meshing::create( 0, ee, 
                                exp->int_data( "nb_space_dimensions" ) ) ;   
   ee->destroy() ; ee = 0 ;
   ee = exp->create_subexplorer( 0, "splitting_strategy" ) ;
   
   size_t_array2D face2cell( 0, 0 ) ;
   size_t_array2D cell2face( 0, 0 ) ;
   build_connectivities( mm, cell2face, face2cell ) ;
   
   size_t nb_ranks = PEL::bad_index() ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( exp->has_entry( "nb_ranks" ) )
   {
      nb_ranks = exp->int_data( "nb_ranks" ) ;
      if( com->nb_ranks() != 1 )
      {
         ostringstream mesg ;
         mesg << "*** GE_SplittingStrategy_TEST error: " << endl ;
         mesg << "    the entry of keyword \"nb_ranks\" is allowed" << endl ;
         mesg << "    only for non-distributed calculations" ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }
   }
   else
   {
      nb_ranks = com->nb_ranks() ;
   }
   
   if( com->nb_ranks() == 1 )
   {
      for( size_t rank=0 ; rank<nb_ranks ; ++rank )
      {
         out() << "| ... meshing of rank " << rank << endl ;
         GE_SplittingStrategy* strategy = 
            GE_SplittingStrategy::create( 0, ee, mm, nb_ranks, rank ) ;

         test_meshing( cell2face, face2cell, strategy ) ;

         strategy->destroy() ; strategy = 0 ;
      }
   }
   else
   {
      GE_SplittingStrategy* strategy = 
         GE_SplittingStrategy::create( 0, ee, mm, com ) ;

      test_meshing( cell2face, face2cell, strategy ) ;

      strategy->destroy() ; strategy = 0 ;
   }
   
   ee->destroy() ; ee = 0 ;
   mm->destroy() ; mm = 0 ;
}

//---------------------------------------------------------------------------
void
GE_SplittingStrategy_TEST:: test_meshing( size_t_array2D const& cell2face, 
                                          size_t_array2D const& face2cell,
                                          GE_SplittingStrategy const* strategy )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy_TEST:: test_meshing" ) ;
   
   size_t nb_cells = cell2face.index_bound( 0 ) ;
   
   boolVector ok_cell( nb_cells ) ;
   ok_cell.set( false ) ;

   size_t icell_init = 0 ;
   while( icell_init < nb_cells &&
          strategy->cell_rank( icell_init ) != strategy->rank() )
   {
      ++icell_init ;
   }
   bool ok = ( icell_init != nb_cells ) ;
   notify_one_test_result( "not empty", ok ) ;
   
   if( ok )
   {
      append_neigh( icell_init, cell2face, face2cell, strategy, ok_cell ) ;

      size_t nb_connected_cells = 0 ;
      for( size_t i=0 ; i<nb_cells ; ++i )
      {
         if( strategy->cell_rank( i ) == strategy->rank() )
         {
            if( !ok_cell( i ) )
            {
               ok = false ;
            }
            else
            {
               ++nb_connected_cells ;
            }
         }
      }
      {
         ostringstream mesg ;
         mesg << "connected (" << nb_connected_cells 
              << " cells traversed)" ;
         notify_one_test_result( mesg.str(), ok ) ;            
      }
   }
}

//---------------------------------------------------------------------------
void
GE_SplittingStrategy_TEST:: append_neigh( size_t icell,
                                          size_t_array2D const& cell2face, 
                                          size_t_array2D const& face2cell,
                                          GE_SplittingStrategy const* strategy,
                                          boolVector& ok_cell )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy_TEST:: append_neigh" ) ;
   
   PEL_ASSERT( strategy->cell_rank( icell ) == strategy->rank() ) ;
   
   ok_cell( icell ) = true ;
   
   for( size_t i=0 ; i<cell2face.index_bound( 1 ) ; ++i )
   {
      size_t iface = cell2face( icell, i ) ;
      if( iface == PEL::bad_index() ) break ;
      
      size_t icell_other = PEL::bad_index() ;
      if( face2cell( iface, 0 ) == icell )
      {
         icell_other = face2cell( iface, 1 ) ;
      }
      else
      {
         PEL_ASSERT( face2cell( iface, 1 ) == icell ) ;
         icell_other = face2cell( iface, 0 ) ;
      }
      
      if( ( icell_other != PEL::bad_index() ) && 
            !ok_cell( icell_other ) && 
          ( strategy->cell_rank( icell_other ) == strategy->rank() ) )
         append_neigh( icell_other, cell2face, face2cell, strategy, ok_cell ) ;
   }
}

//---------------------------------------------------------------------------
void
GE_SplittingStrategy_TEST:: build_connectivities( GE_Meshing* mm, 
                                                  size_t_array2D& cell2face, 
                                                  size_t_array2D& face2cell )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy_TEST:: build_connectivities" ) ;
   
   // maximum 8 faces per cell (for cubes)
   size_t const max_nb_faces_per_cell = 8 ;
   
   face2cell.re_initialize( mm->nb_faces(), 2 ) ;
   face2cell.set( PEL::bad_index() ) ;
   
   cell2face.re_initialize( mm->nb_cells(), max_nb_faces_per_cell ) ;
   cell2face.set( PEL::bad_index() ) ;
   
   size_t icell = 0 ;
   mm->start_cell_iterator() ;
   for( ; mm->valid_cell() ; mm->go_next_cell(), ++icell )
   {
      size_t_vector const& cf = mm->cell_faces() ;
      PEL_ASSERT( cf.size() < max_nb_faces_per_cell ) ;
      for( size_t i=0 ; i<cf.size() ; ++i )
      {
         size_t iface = cf( i ) ;
         
         PEL_ASSERT( cell2face( icell, i ) == PEL::bad_index() ) ;
         cell2face( icell, i ) = iface ;
         
         PEL_ASSERT( face2cell( iface, 1 ) == PEL::bad_index() ) ;
         if( face2cell( iface, 0 ) == PEL::bad_index() )
         {
            face2cell( iface, 0 ) = icell ; 
         }
         else
         {
            face2cell( iface, 1 ) = icell ;
         }
      }
   }
}

