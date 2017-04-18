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

#include <LA_UzawaPreconditioner.hh>

#include <LA_CahouetChabard_UP.hh>
#include <LA_Identity_UP.hh>
#include <LA_Jacobi_UP.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>

//-----------------------------------------------------------------------------
LA_UzawaPreconditioner*
LA_UzawaPreconditioner:: create( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaPreconditioner:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   LA_UzawaPreconditioner* result = 0 ;

   std::string const& name = exp->string_data( "concrete_name" ) ;
   if( name == "LA_CahouetChabard_UP" )
   {
      result = LA_CahouetChabard_UP::create( a_owner, exp ) ;
   }
   else if( name == "LA_Jacobi_UP" )
   {
      result = LA_Jacobi_UP::create( a_owner, exp ) ;
   }
   else if( name == "LA_Identity_UP" )
   {
      result = LA_Identity_UP::create( a_owner ) ;
   }
   else
   {
      PEL_Error::object()->raise_plain( 
                      "invalid LA_UzawaPreconditioner descendant : " + name ) ;
   }
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
LA_UzawaPreconditioner:: LA_UzawaPreconditioner( PEL_Object* a_owner  )
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//-----------------------------------------------------------------------------
LA_UzawaPreconditioner:: ~LA_UzawaPreconditioner( void )
//-----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_UzawaPreconditioner:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_UzawaPreconditioner:: print" ) ;
   
   std::string s( indent_width, ' ' ) ;
   os << s << "Uzawa preconditioner: \"" << type_name() << "\"" << std::endl ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: create_clone_POST( LA_UzawaPreconditioner* result,
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
LA_UzawaPreconditioner:: dimension_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: implementation_PRE( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: build_PRE( LA_Matrix const* A, 
                                    LA_Matrix const* B,
                                    LA_Matrix const* C ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( B != 0 ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( A->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == B->nb_cols() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( B->implementation() == A->implementation() ) ;
   PEL_ASSERT( IMPLIES( C!=0, C->nb_rows()==C->nb_cols() ) ) ;
   PEL_ASSERT( IMPLIES( C!=0, C->nb_rows()==B->nb_rows() ) ) ;
   PEL_ASSERT( IMPLIES( C!=0, C->is_synchronized() ) ) ;

   PEL_ASSERT( IMPLIES( C!=0, C->is_synchronized() ) ) ;
   PEL_ASSERT( IMPLIES( C!=0, C->implementation() == A->implementation() ) ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: build_POST( LA_Matrix const* A, 
                                     LA_Matrix const* B,
                                     LA_Matrix const* C ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_valid(), dimension() == B->nb_rows() ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), implementation() == A->implementation() ) ) ;
   PEL_ASSERT( !successful_solve() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: solve_PRE( LA_Vector const* rhs, 
                                    LA_Vector const* sol ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( rhs != 0 ) ;
   PEL_ASSERT( rhs->is_synchronized() ) ;
   PEL_ASSERT( rhs->implementation() == implementation() ) ;
   PEL_ASSERT( sol != 0 ) ;
   PEL_ASSERT( sol->is_synchronized() ) ;
   PEL_ASSERT( sol->implementation() == implementation() ) ;
   PEL_ASSERT( rhs->nb_rows() == dimension() ) ;
   PEL_ASSERT( sol->nb_rows() == dimension() ) ;
   PEL_ASSERT( rhs->is_synchronized() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
LA_UzawaPreconditioner:: solve_POST( LA_Vector const* rhs,
                                     LA_Vector const* sol ) const
//-----------------------------------------------------------------------------
{
   return( true ) ;
}
