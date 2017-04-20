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

#include <LA_Identity_UP.hh>

#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Vector.hh>


//-------------------------------------------------------------------------
LA_Identity_UP* 
LA_Identity_UP:: create( PEL_Object* a_owner )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: create" ) ;

   LA_Identity_UP* result = new LA_Identity_UP( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
LA_Identity_UP:: LA_Identity_UP( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : LA_UzawaPreconditioner( a_owner )
   , DIM( 0 )
   , IMPL( 0 )
   , SOLVE_OK( false )
{
}

//-------------------------------------------------------------------------
LA_Identity_UP:: ~LA_Identity_UP( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
LA_Identity_UP*
LA_Identity_UP:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: create_clone" ) ;

   LA_Identity_UP* result = new LA_Identity_UP( a_owner ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
LA_Identity_UP:: is_valid( void ) const
//-------------------------------------------------------------------------
{
   return( DIM != 0 ) ;
}

//-------------------------------------------------------------------------
size_t
LA_Identity_UP:: dimension( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( DIM ) ;
}

//-------------------------------------------------------------------------
LA_Implementation const*
LA_Identity_UP:: implementation( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: implementation" ) ;
   PEL_CHECK_PRE( implementation_PRE() ) ;

   return( IMPL ) ;
}

//-------------------------------------------------------------------------
void
LA_Identity_UP:: build( LA_Matrix const* A,
                        LA_Matrix const* B,
                        LA_Matrix const* C )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: build" ) ;
   PEL_CHECK_PRE( build_PRE( A, B, C ) ) ;

   DIM = B->nb_rows() ;
   IMPL = A->implementation() ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( A, B, C ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//-------------------------------------------------------------------------
void
LA_Identity_UP:: solve( LA_Vector const* rhs,
                        LA_Vector* sol )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   sol->set( rhs ) ;
   SOLVE_OK = true ;
   
   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_Identity_UP:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_UP:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}



