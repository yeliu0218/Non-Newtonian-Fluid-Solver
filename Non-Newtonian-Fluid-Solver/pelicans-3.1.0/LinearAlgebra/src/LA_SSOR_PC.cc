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

#include <LA_SSOR_PC.hh>

#include <LA_SOR_PC.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

LA_SSOR_PC const* LA_SSOR_PC::PROTOTYPE = new LA_SSOR_PC() ;

//----------------------------------------------------------------------
LA_SSOR_PC:: LA_SSOR_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_SSOR_PC" )
   , SOR_PC( 0 )
{
}

//----------------------------------------------------------------------
LA_SSOR_PC*
LA_SSOR_PC:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double const smallest_inverted_item =
                         exp->double_data( "smallest_inverted_item" ) ;
   exp->test_data( "smallest_inverted_item", "smallest_inverted_item>0." ) ;
   
   double const omega = exp->double_data( "omega" ) ;
   exp->test_data( "omega", "omega>0. && omega<2." ) ;

   size_t const nb_iters = (size_t) exp->int_data( "nb_inner_iterations" ) ;
   exp->test_data( "nb_inner_iterations", "nb_inner_iterations>0" ) ;
   
   LA_SSOR_PC* result = new LA_SSOR_PC( a_owner,
                                        omega,
                                        nb_iters,
                                        smallest_inverted_item ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SSOR_PC*
LA_SSOR_PC:: create( PEL_Object* a_owner, 
                     double omega,
                     size_t nb_iters,
                     double smallest_inverted_item )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: create" ) ;
   PEL_CHECK_PRE( omega>0 && omega<2 ) ;
   PEL_CHECK_PRE( nb_iters>0 ) ;
   PEL_CHECK_PRE( smallest_inverted_item > 0 ) ;

   LA_SSOR_PC* result = new LA_SSOR_PC( a_owner, 
                                        omega,
                                        nb_iters, 
                                        smallest_inverted_item ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SSOR_PC:: LA_SSOR_PC( PEL_Object* a_owner, 
                         double omega,
                         size_t nb_iters, 
                         double smallest_inverted_item )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , SOR_PC( LA_SOR_PC::create( this, omega, "symmetric",
                                nb_iters, smallest_inverted_item ) )
{
}

//----------------------------------------------------------------------
LA_SSOR_PC:: ~LA_SSOR_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_SSOR_PC*
LA_SSOR_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: create_clone" ) ;

   LA_SSOR_PC* result =new LA_SSOR_PC( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SSOR_PC:: LA_SSOR_PC( PEL_Object* a_owner, LA_SSOR_PC const* other )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , SOR_PC( other->SOR_PC->create_clone( a_owner ) )
{
}

//----------------------------------------------------------------------
bool
LA_SSOR_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( SOR_PC->is_valid() ) ;
}

//----------------------------------------------------------------------
size_t
LA_SSOR_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( SOR_PC->dimension() ) ;
}

//----------------------------------------------------------------------
void
LA_SSOR_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   SOR_PC->build( mat ) ;
   
   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_SSOR_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_SSOR_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;

   SOR_PC->unbuild() ;
   
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SSOR_PC:: solve( LA_Vector const* rhs, LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SSOR_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   SOR_PC->solve( rhs, sol ) ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_SSOR_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   return( SOR_PC->successful_solve() ) ;
}

//----------------------------------------------------------------------
void
LA_SSOR_PC:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   SOR_PC->print_more( os, indent_width ) ;
}
