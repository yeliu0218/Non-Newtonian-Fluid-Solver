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

#include <GE_PointPoint_INT.hh>

#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <GE_Point.hh>

//----------------------------------------------------------------------
GE_PointPoint_INT*
GE_PointPoint_INT:: create( PEL_Object* a_owner,
                            std::string const& a_name,
                            PEL_ModuleExplorer const* a_mod_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint_INT:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_mod_exp!=0 ) ;

   GE_PointPoint_INT const* proto =
      static_cast<GE_PointPoint_INT const*>(
                                    plugins_map()->item( a_name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   GE_PointPoint_INT* result =  proto->create_replica( a_owner, a_mod_exp ) ;
   PEL_CHECK( !result->is_a_prototype() ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointPoint_INT:: GE_PointPoint_INT( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , PROTO( true )
{
   PEL_LABEL( "GE_PointPoint_INT:: GE_PointPoint_INT" ) ;

   plugins_map()->register_item( a_name, this ) ;
   
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointPoint_INT:: GE_PointPoint_INT( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , PROTO( false )
{
   PEL_LABEL( "GE_PointPoint_INT:: GE_PointPoint_INT" ) ;

   PEL_CHECK( owner()==a_owner ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointPoint_INT:: ~GE_PointPoint_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint_INT:: ~GE_PointPoint_INT" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_PointPoint_INT:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( PROTO ) ;
}

//----------------------------------------------------------------------
bool
GE_PointPoint_INT:: points_are_close_PRE( GE_Point const* P1,
                                          GE_Point const* P2 ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( P1!=0 ) ;
   PEL_ASSERT( P2!=0 ) ;
   PEL_ASSERT( P1->nb_coordinates()==P2->nb_coordinates() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
GE_PointPoint_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_PointPoint_INT:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "GE_PointPoint_INT descendant" ) ;
   return( result ) ;
}

