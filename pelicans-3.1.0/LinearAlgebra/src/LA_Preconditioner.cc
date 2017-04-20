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

#include <LA_Preconditioner.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <iostream>

using std::string ;
using std::endl ;

//----------------------------------------------------------------------
LA_Preconditioner*
LA_Preconditioner:: make( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Preconditioner:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   LA_Preconditioner const* proto =
      static_cast<LA_Preconditioner const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   LA_Preconditioner* result = proto->create_replica( a_owner, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
LA_Preconditioner:: LA_Preconditioner( PEL_Object* a_owner )
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
   PEL_LABEL( "LA_Preconditioner:: LA_Preconditioner" ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------------
LA_Preconditioner:: LA_Preconditioner( std::string const& class_name )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "LA_Preconditioner:: LA_Preconditioner" ) ;
   
   plugins_map()->register_item( class_name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//-----------------------------------------------------------------------------
LA_Preconditioner:: ~LA_Preconditioner( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Preconditioner:: ~LA_Preconditioner" ) ;
}

//----------------------------------------------------------------------
void
LA_Preconditioner:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Preconditioner:: print" ) ;
   std::string s( indent_width, ' ') ;
   
   os << s << "LA_Preconditioner: \"" << type_name() << "\"" << endl ;
   print_more( os, indent_width+3 ) ;
}

//----------------------------------------------------------------------
void
LA_Preconditioner:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: is_a_prototype( void ) const
//-----------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: create_clone_POST( LA_Preconditioner* result,
                                       PEL_Object* a_owner ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::create_clone_POST( result, a_owner ) ) ;
   PEL_ASSERT( !result->is_valid() ) ;
   PEL_ASSERT( !result->successful_solve() ) ;
   return( true ) ;
}


//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: dimension_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: build_PRE( LA_Matrix const* mat ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( mat!=0 ) ;
   PEL_ASSERT( mat->nb_rows() == mat->nb_cols() ) ;
   PEL_ASSERT( mat->nb_rows() > 0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: build_POST( LA_Matrix const* mat ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_valid(), dimension() == mat->nb_rows() ) ) ;
   PEL_ASSERT( !successful_solve() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: unbuild_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: unbuild_POST( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( !is_valid() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: solve_PRE( LA_Vector const* rhs,
                               LA_Vector const* sol ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( rhs->is_synchronized() ) ;
   PEL_ASSERT( rhs->nb_rows() == dimension() ) ;
   PEL_ASSERT( sol->nb_rows() == dimension() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_Preconditioner:: solve_POST( LA_Vector const* rhs,
                                LA_Vector const* sol ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( sol->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Preconditioner:: create_replica_PRE( PEL_Object const* a_owner,
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Preconditioner:: create_replica_POST( 
                                    LA_Preconditioner const* result,
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
LA_Preconditioner:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "LA_Preconditioner descendant" ) ;
   return( result ) ;
}
