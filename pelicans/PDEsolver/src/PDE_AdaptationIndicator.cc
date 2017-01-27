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

#include <PDE_AdaptationIndicator.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <GE_Mpolyhedron.hh>

#include <PDE_ReferenceElement.hh>

#include <iostream>
#include <string>

using std::string ;
using std::cout ; using std::endl ;

//---------------------------------------------------------------------------
PDE_AdaptationIndicator*
PDE_AdaptationIndicator:: make( PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                PEL_ModuleExplorer const* exp,
                                size_t a_verbose_level  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationIndicator:: make" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   string name = exp->string_data( "concrete_name" ) ;
   PDE_AdaptationIndicator const* proto =
   static_cast<PDE_AdaptationIndicator const*>( plugins_map()->item( name ) ) ;
      
   PDE_AdaptationIndicator* result = proto->create_replica( a_owner, 
                                                            dom,
                                                            exp,
                                                            a_verbose_level ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationIndicator:: PDE_AdaptationIndicator( PEL_Object* a_owner,
                                                   size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , VERB( a_verbose_level )
{
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
PDE_AdaptationIndicator:: PDE_AdaptationIndicator( std::string const& a_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , VERB( PEL::bad_index() )
{
   PEL_LABEL( "PDE_AdaptationIndicator:: PDE_AdaptationIndicator" ) ;

   if( !a_name.empty() )
   {
      plugins_map()->register_item( a_name, this ) ;
   }

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationIndicator:: ~PDE_AdaptationIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationIndicator:: is_a_prototype( void ) const
//---------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//---------------------------------------------------------------------------
size_t
PDE_AdaptationIndicator:: verbose_level( void ) const
//---------------------------------------------------------------------------
{
   return( VERB ) ;
}

//---------------------------------------------------------------------------
bool 
PDE_AdaptationIndicator:: to_be_refined_PRE( double bf_indicator,
                                             GE_Mpolyhedron const* poly,
                                             PDE_ReferenceElement const* elm,
                                             size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( poly != 0 ) ;
   PEL_ASSERT( elm != 0 ) ;
   PEL_ASSERT( elm->reference_polyhedron() == poly->reference_polyhedron() ) ;
   PEL_ASSERT( local_node < elm->nb_nodes() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool 
PDE_AdaptationIndicator:: to_be_unrefined_PRE( double bf_indicator,
                                               GE_Mpolyhedron const* poly,
                                               PDE_ReferenceElement const* elm,
                                               size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( poly != 0 ) ;
   PEL_ASSERT( elm != 0 ) ;
   PEL_ASSERT( elm->reference_polyhedron() == poly->reference_polyhedron() ) ;
   PEL_ASSERT( local_node < elm->nb_nodes() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationIndicator:: create_replica_PRE( PEL_Object* a_owner,
                                              PDE_DomainAndFields const* dom,
                                              PEL_ModuleExplorer const* exp,
                                              size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationIndicator:: create_replica_POST(  
                                       PDE_AdaptationIndicator const* result,
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       PEL_ModuleExplorer const* exp,
                                       size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->verbose_level() == a_verbose_level ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
PEL_ObjectRegister*
PDE_AdaptationIndicator:: plugins_map( void )
//---------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "PDE_AdaptationIndicator descendant" ) ;
   return( result ) ;
}
