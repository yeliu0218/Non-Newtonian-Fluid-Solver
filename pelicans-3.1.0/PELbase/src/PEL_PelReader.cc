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

#include <PEL_PelReader.hh>

#include <PEL.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <sstream>

PEL_PelReader const* PEL_PelReader::PROTOTYPE = new PEL_PelReader() ;

//----------------------------------------------------------------------
PEL_PelReader:: PEL_PelReader( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingReader( "PEL_PelReader" )
   , MESHING_EXP( 0 )
   , FIELDS_EXP( 0 )
   , I_DOM_EXP( 0 )
   , VAR_EXP( 0 )
{
}

//----------------------------------------------------------------------
PEL_PelReader*
PEL_PelReader:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_PelReader* result = new PEL_PelReader( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_PelReader:: PEL_PelReader( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingReader( a_owner )
   , MESHING_EXP( 0 )
   , FIELDS_EXP( 0 )
   , I_DOM_EXP( 0 )
   , VAR_EXP( 0 )
{
   PEL_LABEL( "PEL_PelReader:: PEL_PelReader" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t i_cycle = exp->int_data( "cycle" ) ;

   std::string PEL_file = exp->string_data( "files_basename" )+".pel" ;

   PEL_Module* mod = PEL_Module::create( this, "MAIN", PEL_file ) ;

   MESHING_EXP = restore_cycle( mod, i_cycle, "meshing" ) ;
   if( MESHING_EXP==0 )
   {
      MESHING_EXP = restore_cycle( mod, 1, "meshing" ) ;
   }
   FIELDS_EXP = restore_cycle( mod, i_cycle, "fields" ) ;
   I_DOM_EXP = restore_cycle( mod, i_cycle, "integration_domain" ) ;
   if( I_DOM_EXP==0 )
   {
      I_DOM_EXP = restore_cycle( mod, 1, "integration_domain" ) ;
   }
   VAR_EXP = restore_cycle( mod, 1, "variables" ) ;
      
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_PelReader:: ~PEL_PelReader( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: ~PEL_PelReader" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_PelReader:: meshing( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: meshing" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = MESHING_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( meshing_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_PelReader:: fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: fields" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = FIELDS_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( fields_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_PelReader:: integration_domain( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: integration_domain" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = I_DOM_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( integration_domain_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_PelReader:: variables( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: variables" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = VAR_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( variables_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_PelReader:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_DataOnMeshingReader::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_PelReader:: restore_cycle( PEL_Module const* m,
                               size_t i_cycle,
                               std::string const& name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelReader:: restore_cycle" ) ;
   PEL_CHECK( m!=0 ) ;
   PEL_CHECK( !name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_ModuleExplorer* result = 0 ;

   std::ostringstream os ;
   os <<  "cycle_" << (int) i_cycle ;

   if( m->has_module( os.str() ) )
   {
      PEL_Module const* n = m->module( os.str() ) ;
      if( n->has_module( name ) )
      {
         result = PEL_ModuleExplorer::create( this, n->module( name ) ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->owner()==this ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->name()==name ) ) ;
   return( result ) ;
}
