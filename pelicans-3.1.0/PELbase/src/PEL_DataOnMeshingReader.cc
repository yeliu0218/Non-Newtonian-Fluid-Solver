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

#include <PEL_DataOnMeshingReader.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_DataOnMeshingReader*
PEL_DataOnMeshingReader:: make( PEL_Object* a_owner,
		                std::string const& name,
                                PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataOnMeshingReader:: make" ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PEL_DataOnMeshingReader const* proto =
      static_cast<PEL_DataOnMeshingReader const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PEL_DataOnMeshingReader* result = proto->create_replica( a_owner, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DataOnMeshingReader:: PEL_DataOnMeshingReader( std::string const& name )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "PEL_DataOnMeshingReader:: PEL_DataOnMeshingReader" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_DataOnMeshingReader:: PEL_DataOnMeshingReader( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_DataOnMeshingReader:: ~PEL_DataOnMeshingReader( void  )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: meshing_POST(
                                PEL_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, result->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->name()=="meshing" ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->has_entry( "nb_sp_dims" ) ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->has_entry( "vertices" ) ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->has_entry( "cell2vertex" ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: fields_POST(
                                PEL_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, result->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->name()=="fields" ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: integration_domain_POST(
                                PEL_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, result->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->name()=="integration_domain" ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->has_entry( "inner_boundary" ) ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->has_entry( "polygon" ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: variables_POST(
                                PEL_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, result->owner()==this ) ) ;
   PEL_ASSERT( IMPLIES( result!=0, result->name()=="variables" ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: create_replica_PRE( 
                                        PEL_Object const* a_owner,
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingReader:: create_replica_POST(
                                    PEL_DataOnMeshingReader const* result,
                                    PEL_Object const* a_owner,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_DataOnMeshingReader:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "PEL_DataOnMeshingReader descendant" ) ;
   return( result ) ;
}
