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

#include <PDE_CrossProcessBFNumbering.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_VectorIterator.hh>
#include <PEL_assertions.hh>
#include <boolArray2D.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_GridFE.hh>
#include <PDE_MeshFE.hh>
#include <PDE_CellFE.hh>

#include <PDE_ReferenceElement.hh>
#include <PDE_SetOfBasisFunctions.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using std::endl ;

//----------------------------------------------------------------------
PDE_CrossProcessBFNumbering:: PDE_CrossProcessBFNumbering(
                                             PEL_Object* a_owner,
                                             PDE_DomainBuilder const* a_dom,
                                             PEL_Communicator const* a_com )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DIM( a_dom->nb_space_dimensions() )
   , DOMB( a_dom )
   , PT( 0 )
   , COMM( a_com )
   , NB_GLOBAL_BF( 0 )
   , BFS( a_dom->set_of_basis_functions() )
   , IS_SYNC( false )
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: PDE_CrossProcessBFNumbering" ) ;

   PT = GE_Point::create( this, DIM ) ;

   NB_GLOBAL_BF = size_t_vector( BFS->nb_ref_elts_grps() ) ;
   size_t_vector const dummy( 0 ) ;
   LOCAL_BF = std::vector<size_t_vector>( BFS->nb_ref_elts_grps(),
                                          dummy ) ;
   GLOBAL_BF = std::vector<size_t_vector>( BFS->nb_ref_elts_grps(),
                                          dummy ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( BFS != 0 ) ;
}

//----------------------------------------------------------------------
PDE_CrossProcessBFNumbering:: ~PDE_CrossProcessBFNumbering( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: ~PDE_CrossProcessBFNumbering" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: resize_internals( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: resize_internals" ) ;

   NB_GLOBAL_BF.resize( BFS->nb_ref_elts_grps() ) ;
   size_t_vector const dummy( 0 ) ;
   LOCAL_BF.resize( BFS->nb_ref_elts_grps(), dummy ) ;
   GLOBAL_BF.resize( BFS->nb_ref_elts_grps(), dummy ) ;

   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PDE_SetOfBasisFunctions const*
PDE_CrossProcessBFNumbering:: set_of_basis_functions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: set_of_basis_functions" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_SetOfBasisFunctions* result = BFS ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Communicator const*
PDE_CrossProcessBFNumbering:: communicator( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: communicator" ) ;

   PEL_Communicator const* result = COMM ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessBFNumbering:: nb_global_basis_functions( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: nb_global_basis_functions" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( e < set_of_basis_functions()->nb_ref_elts_grps() ) ;

   return( NB_GLOBAL_BF( e ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessBFNumbering:: nb_global_basis_functions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: nb_global_basis_functions" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;

   return( NB_GLOBAL_BF.sum() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessBFNumbering:: global_basis_function_index( size_t local_ind,
                                                           size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: global_basis_function_index" ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( e < set_of_basis_functions()->nb_ref_elts_grps() ) ;
   PEL_CHECK_PRE( local_ind < 
                        set_of_basis_functions()->nb_basis_functions( e ) ) ;
   size_t shift = 0 ;
   for( size_t i = 0 ; i < e ; ++i )
   {
      shift += nb_global_basis_functions( i ) ;
   }

   size_t result = shift + GLOBAL_BF[e]( local_ind ) ;
   PEL_ASSERT( result < shift + nb_global_basis_functions( e ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessBFNumbering:: local_basis_function_index( size_t global_ind,
                                                          size_t& local_e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: local_basis_function_index" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_synchronized() ) ;
   PEL_CHECK_PRE( global_ind <  nb_global_basis_functions() ) ;

   size_t e = 0 ;
   size_t shift = 0 ;
   size_t shiftp = nb_global_basis_functions( e ) ;
   while( global_ind + 1 > shiftp )
   {
      e++ ;
      shift = shiftp ;
      shiftp += nb_global_basis_functions( e ) ;
   }
   local_e = e ;
   size_t result = LOCAL_BF[e]( global_ind - shift ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: synchronize_new_basis_functions( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: synchronize_new_basis_functions" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t e = 0 ; e < BFS->nb_ref_elts_grps() ; ++e )
   {
      size_t old_nb_bfs = GLOBAL_BF[e].size() ;
      size_t old_nb_global_bfs = NB_GLOBAL_BF( e ) ;
      size_t nb_new_bfs = BFS->nb_basis_functions( e ) - old_nb_bfs ;

      bool go_on = COMM->boolean_or( nb_new_bfs > 0 ) ;
      if( go_on )
      {
         doubleArray2D coord_glob( DIM, BFS->nb_basis_functions( e ) ) ;
         recover_coordinates( BFS, e, coord_glob ) ;

         // Only take new basis functions:
         doubleArray2D coord( DIM, nb_new_bfs ) ;
         for( size_t i=0 ; i<nb_new_bfs ; ++i )
         {
            for( size_t ic=0 ; ic<DIM ; ++ic )
            {
               coord( ic, i ) = coord_glob( ic, i+old_nb_bfs ) ;
            }
         }

         // Recover new basis functions global numbering:
         size_t_vector global_bfs( 0 ) ;
         size_t nb_glob_bf ;
         recover_global_numbering( coord, global_bfs, nb_glob_bf, false ) ;

         // Update numberings:
         GLOBAL_BF[e].resize( old_nb_bfs + global_bfs.size() ) ;
         for( size_t i=0 ; i<global_bfs.size() ; ++i )
         {
            GLOBAL_BF[e]( old_nb_bfs+i ) = global_bfs(i) + old_nb_global_bfs ;
         }
         NB_GLOBAL_BF(e) += nb_glob_bf ;

         set_global_to_local( GLOBAL_BF[e], e, LOCAL_BF[e] ) ;
      }
   }
   IS_SYNC = true ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: globalize(  void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: globalize" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t e= 0 ; e < BFS->nb_ref_elts_grps() ; ++e )
   {
      doubleArray2D coord( DIM, BFS->nb_basis_functions( e ) ) ;

      recover_coordinates( BFS, e,  coord ) ;
      recover_global_numbering( coord, GLOBAL_BF[e], NB_GLOBAL_BF(e), true ) ;
      set_global_to_local( GLOBAL_BF[e], e,  LOCAL_BF[e] ) ;
   }
   IS_SYNC = true ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: recover_coordinates(
                                            PDE_SetOfBasisFunctions const* bf_set,
                                            size_t e,
                                            doubleArray2D& coord ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessBFNumbering:: recover_coordinates" ) ;
   PEL_CHECK( coord.index_bound(0) == DIM ) ;
   PEL_CHECK( coord.index_bound(1) == bf_set->nb_basis_functions( e ) ) ;

   //--For subsequent check
   boolVector is_ok( bf_set->nb_basis_functions( e ) ) ;
   is_ok.set( false ) ;
   //--

   bf_set->start( e ) ;
   for(  ; bf_set->is_valid( e ) ; bf_set->go_next( e ) )
   {
      PDE_BasisFunctionCell* bf =
                    static_cast<PDE_BasisFunctionCell*>( bf_set->item( e ) ) ;
      bf->geometrical_node( PT ) ;
      size_t n = bf->id_number() ;
      for( size_t d=0 ; d<DIM ; ++d )
      {
         coord( d, n ) = PT->coordinate( d ) ;
      }
      //--For subsequent check
      is_ok( n ) = true ;
      //--
   }

   //--For subsequent check
   for( size_t n=0 ; n<bf_set->nb_basis_functions( e ); ++n )
   {
      PEL_ASSERT( is_ok( n ) ) ;
   }
   //--
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: recover_global_numbering(
                                               doubleArray2D& coord,
                                               size_t_vector& global_bf_id,
                                               size_t& nb_glob_bf,
                                               bool verbose ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: recover_global_numbering" ) ;
   PEL_CHECK( coord.index_bound(0) == DIM ) ;

   size_t nb_local_bf = coord.index_bound(1) ;

   global_bf_id.re_initialize( nb_local_bf ) ;
   global_bf_id.set( PEL::bad_index() ) ;

   COMM->merge( GE_Point::coordinates_comparator(), coord, global_bf_id ) ;

   nb_glob_bf = coord.index_bound( 1 ) ;
   COMM->broadcast( nb_glob_bf ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: set_global_to_local(
                                    size_t_vector const& global_bf_id,
                                    size_t e,
                                    size_t_vector& local_bf_id ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: set_global_to_local" ) ;
   PEL_CHECK(
      FORALL( ( size_t i=0 ; i<global_bf_id.size() ; i++ ),
                global_bf_id(i)<NB_GLOBAL_BF(e) ) ) ;

   local_bf_id.re_initialize( NB_GLOBAL_BF(e) ) ;
   local_bf_id.set( PEL::bad_index() ) ;
   for( size_t i=0 ; i<global_bf_id.size() ; i++ )
   {
      size_t ig = global_bf_id( i ) ;
      local_bf_id( ig ) = i ;
   }
}

//----------------------------------------------------------------------
bool
PDE_CrossProcessBFNumbering:: is_synchronized( void ) const
//----------------------------------------------------------------------
{
   return( IS_SYNC ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessBFNumbering:: set_unsynchronized_state( void )
//----------------------------------------------------------------------
{
   IS_SYNC = false ;
   PEL_CHECK_POST( !is_synchronized() ) ;
}

//----------------------------------------------------------------------
bool
PDE_CrossProcessBFNumbering:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( DOMB->facade()->is_distributed() ) ;
   PEL_ASSERT( BFS != 0 ) ;
   bool size_ok = true ;
   size_ok = size_ok&& ( NB_GLOBAL_BF.size() == BFS->nb_ref_elts_grps() ) ;
   size_ok = size_ok&& ( LOCAL_BF.size() == BFS->nb_ref_elts_grps() ) ;
   size_ok = size_ok&& ( GLOBAL_BF.size() == BFS->nb_ref_elts_grps() ) ;
   if( !size_ok )
      PEL_Error::object()->raise_plain(
                 "Number of reference elements groups changed "
                 "\"resize_internal()\" has to be called before"
                 " usage of PDE_CrossProcessBFNumbering" ) ;
   return( true ) ;
}


