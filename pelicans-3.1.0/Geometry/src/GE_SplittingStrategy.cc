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

#include <GE_SplittingStrategy.hh>

#include <GE_Color.hh>
#include <GE_Meshing.hh>

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <iostream>
#include <string>

//-------------------------------------------------------------------------
GE_SplittingStrategy*
GE_SplittingStrategy:: create( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp,
                               GE_Meshing* meshing,
                               PEL_Communicator const* com )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy:: create(com)" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( meshing != 0 ) ;
   PEL_CHECK_PRE( com != 0 ) ;
   
   std::string const& name = exp->string_data( "concrete_name" ) ;
   GE_SplittingStrategy const* proto =
      static_cast<GE_SplittingStrategy const*>( plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
   
   GE_SplittingStrategy* result =
                 proto->create_replica( a_owner, exp, meshing, com ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_ranks() == com->nb_ranks() ) ;
   PEL_CHECK_POST( result->rank() == com->rank() ) ;
   PEL_CHECK_POST( result->nb_cells() == meshing->nb_cells() ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_SplittingStrategy*
GE_SplittingStrategy:: create( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp,
                               GE_Meshing* meshing,
                               size_t nb_rks, size_t rk )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplittingStrategy:: create(nb_rks,rk)" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( meshing != 0 ) ;
   PEL_CHECK_PRE( rk<nb_rks ) ;
   
   std::string const& name = exp->string_data( "concrete_name" ) ;
   GE_SplittingStrategy const* proto =
      static_cast<GE_SplittingStrategy const*>( plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
   
   GE_SplittingStrategy* result =
      proto->create_replica( a_owner, exp, meshing, nb_rks, rk ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_ranks() == nb_rks ) ;
   PEL_CHECK_POST( result->rank() == rk ) ;
   PEL_CHECK_POST( result->nb_cells() == meshing->nb_cells() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_SplittingStrategy:: GE_SplittingStrategy( PEL_Object* a_owner,
					     PEL_ModuleExplorer const* exp,
					     GE_Meshing* meshing,
					     PEL_Communicator const* com )
//------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NB_CELLS( meshing->nb_cells() )
   , NB_RANKS( com->nb_ranks() )
   , RANK( com->rank() )
{
   PEL_LABEL( "GE_SplittingStrategy:: GE_SplittingStrategy" ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
GE_SplittingStrategy:: GE_SplittingStrategy( PEL_Object* a_owner,
					     PEL_ModuleExplorer const* exp,
					     GE_Meshing* meshing,
					     size_t nb_rks, size_t rk  )
//------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NB_CELLS( meshing->nb_cells() )
   , NB_RANKS( nb_rks )
   , RANK( rk )
{
   PEL_LABEL( "GE_SplittingStrategy:: GE_SplittingStrategy" ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
GE_SplittingStrategy:: GE_SplittingStrategy( std::string const& name )
//------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , NB_CELLS( PEL::bad_index() )
   , NB_RANKS( PEL::bad_index() )
   , RANK( PEL::bad_index() )
{
   PEL_LABEL( "GE_SplittingStrategy:: GE_SplittingStrategy" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//------------------------------------------------------------------------
GE_SplittingStrategy:: ~GE_SplittingStrategy( void )
//------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
bool
GE_SplittingStrategy:: is_a_prototype( void ) const
//------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//-------------------------------------------------------------------------
size_t
GE_SplittingStrategy:: nb_cells( void ) const
//-------------------------------------------------------------------------
{
   return( NB_CELLS ) ;
}

//-------------------------------------------------------------------------
size_t
GE_SplittingStrategy:: nb_ranks( void ) const
//-------------------------------------------------------------------------
{
   return( NB_RANKS ) ;
}

//-------------------------------------------------------------------------
size_t
GE_SplittingStrategy:: rank( void ) const
//-------------------------------------------------------------------------
{
   return( RANK ) ;
}

//-------------------------------------------------------------------------
bool
GE_SplittingStrategy:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
GE_SplittingStrategy:: cell_rank_PRE( size_t mesh_id ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( mesh_id<nb_cells() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
GE_SplittingStrategy:: cell_rank_POST( size_t mesh_id, size_t result ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( result<nb_ranks() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool 
GE_SplittingStrategy:: create_replica_PRE( PEL_Object* a_owner, 
                                           PEL_ModuleExplorer const* exp,
                                           GE_Meshing* meshing, 
                                           PEL_Communicator const* com ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( exp!=0 ) ;
   PEL_ASSERT( meshing!=0 ) ;
   PEL_ASSERT( com!=0 ) ;
   return( true ) ;
}


//-------------------------------------------------------------------------
bool 
GE_SplittingStrategy:: create_replica_POST( PEL_Object* a_owner, 
                                            PEL_ModuleExplorer const* exp,
                                            GE_Meshing* meshing, 
                                            PEL_Communicator const* com, 
                                            GE_SplittingStrategy const* result ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->owner()==a_owner ) ;
   PEL_ASSERT( result->nb_ranks() == com->nb_ranks() ) ;
   PEL_ASSERT( result->rank() == com->rank() ) ;
   PEL_ASSERT( result->nb_cells() == meshing->nb_cells() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool 
GE_SplittingStrategy:: create_replica_PRE( PEL_Object* a_owner, 
                                           PEL_ModuleExplorer const* exp,
                                           GE_Meshing* meshing, 
                                           size_t nb_rks, size_t rk ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( exp!=0 ) ;
   PEL_ASSERT( meshing!=0 ) ;
   PEL_ASSERT( nb_rks>0 ) ;
   PEL_ASSERT( rk<nb_rks ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool 
GE_SplittingStrategy:: create_replica_POST( PEL_Object* a_owner, 
                                            PEL_ModuleExplorer const* exp,
                                            GE_Meshing* meshing, 
                                            size_t nb_rks, size_t rk, 
                                            GE_SplittingStrategy const* result ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->owner()==a_owner ) ;
   PEL_ASSERT( result->nb_ranks() == nb_rks ) ;
   PEL_ASSERT( result->rank() == rk ) ;
   PEL_ASSERT( result->nb_cells() == meshing->nb_cells() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_SplittingStrategy:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
				  "GE_SplittingStrategy descendant" ) ;
   return( result ) ;
}
