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

#include <LA_Identity_PC.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

#include <iostream>

LA_Identity_PC const* LA_Identity_PC::PROTOTYPE = new LA_Identity_PC() ;

//----------------------------------------------------------------------
LA_Identity_PC:: LA_Identity_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_Identity_PC" )
   , SIZE( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_Identity_PC*
LA_Identity_PC:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_Identity_PC* result = new LA_Identity_PC( a_owner ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Identity_PC:: LA_Identity_PC( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , SIZE( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_Identity_PC:: ~LA_Identity_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Identity_PC* 
LA_Identity_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: create_clone" ) ;

   LA_Identity_PC* result = new LA_Identity_PC( a_owner ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_Identity_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: is_valid" ) ;
   return( SIZE != 0 ) ;
}

//----------------------------------------------------------------------
size_t
LA_Identity_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( SIZE ) ;
}

//----------------------------------------------------------------------
void
LA_Identity_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_Identity_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   SIZE = mat->nb_rows() ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_Identity_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_Identity_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   SIZE=0 ;
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Identity_PC:: solve( LA_Vector const* rhs,
                        LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   sol->set( rhs ) ;
   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_Identity_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}
