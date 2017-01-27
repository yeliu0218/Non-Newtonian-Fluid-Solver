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

#include <GE_Transform.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_Translation.hh>

#include <iostream>

using std::string ;

//----------------------------------------------------------------------
GE_Transform*
GE_Transform:: create( PEL_Object* a_owner,
                      size_t nb_sp_dims,
                      PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Transform:: create" ) ;
   
   GE_Transform* result = 0 ;
   string nn = exp->string_data( "concrete_name" ) ;
   if( nn == "GE_Translation"  )
   {
      result = GE_Translation::create( a_owner, nb_sp_dims, exp ) ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value( exp, "concrete_name",
                                                 "   \"GE_Translation\"" ) ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == nb_sp_dims ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Transform:: GE_Transform( PEL_Object* a_owner,
                             size_t nb_sp_dims,
                             PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_DIMS( nb_sp_dims )
   , S_COLOR( GE_Color::object( exp->string_data( "source_color" ) ) )
   , T_COLOR( GE_Color::object( exp->string_data( "target_color" ) ) )
{
}

//----------------------------------------------------------------------
GE_Transform:: GE_Transform( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//----------------------------------------------------------------------
GE_Transform:: ~GE_Transform( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
GE_Transform:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( NB_DIMS ) ;
}

//----------------------------------------------------------------------
GE_Color const*
GE_Transform:: source_color( void ) const
//----------------------------------------------------------------------
{
   return( S_COLOR ) ;
}

//----------------------------------------------------------------------
GE_Color const*
GE_Transform:: target_color( void ) const
//----------------------------------------------------------------------
{
   return( T_COLOR ) ;
}

//----------------------------------------------------------------------
void
GE_Transform:: initialize_as_inverse( GE_Transform const* other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Transform:: initialize_as_inverse" ) ;
   
   NB_DIMS = other->NB_DIMS ;
   S_COLOR = other->T_COLOR ;
   T_COLOR = other->S_COLOR ;
   
   PEL_CHECK_POST( nb_space_dimensions() == other->nb_space_dimensions() ) ;
   PEL_CHECK_POST( source_color() == other->target_color() ) ;
   PEL_CHECK_POST( target_color() == other->source_color() ) ;
}

//----------------------------------------------------------------------
bool
GE_Transform:: apply_PRE( GE_Point* pt ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( pt!=0 ) ;
   PEL_ASSERT( pt->nb_coordinates()==nb_space_dimensions() ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
GE_Transform:: inverse_POST( GE_Transform const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == owner() ) ;
   PEL_ASSERT( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   PEL_ASSERT( result->source_color() == target_color() ) ;
   PEL_ASSERT( result->target_color() == source_color() ) ;
   return true ;
}
