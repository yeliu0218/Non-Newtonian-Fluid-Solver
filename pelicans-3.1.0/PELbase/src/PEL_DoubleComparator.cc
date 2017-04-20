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

#include <PEL_DoubleComparator.hh>

#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

//----------------------------------------------------------------------
PEL_DoubleComparator const*
PEL_DoubleComparator:: make( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleComparator:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   PEL_DoubleComparator const* proto =
      static_cast<PEL_DoubleComparator const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PEL_DoubleComparator const* result = proto->create_replica( a_owner, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_DoubleComparator:: PEL_DoubleComparator( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparator:: PEL_DoubleComparator( std::string const& name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "PEL_DoubleComparator:: PEL_DoubleComparator" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_under_ownership_of( plugins_map() ) ) ;
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparator:: ~PEL_DoubleComparator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
PEL_DoubleComparator:: create_replica_PRE(
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DoubleComparator:: create_replica_POST(
                                   PEL_DoubleComparator const* result,
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}
//----------------------------------------------------------------------
bool
PEL_DoubleComparator:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DoubleComparator:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_DoubleComparator:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
            PEL_ObjectRegister::create( PEL_Root::object(),
                                        "PEL_DoubleComparator descendant" ) ;
   return( result ) ;
}
