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

#include <LA_ILU0_PC.hh>

#include <LA_CRSmatrix.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <sstream>


LA_ILU0_PC const* LA_ILU0_PC::PROTOTYPE = new LA_ILU0_PC() ;

//----------------------------------------------------------------------
LA_ILU0_PC:: LA_ILU0_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_ILU0_PC" )
   , MODIFIED( false )
   , PIV_MIN( -PEL::max_double() )
   , LU( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILU0_PC*
LA_ILU0_PC:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   bool const diagonal_compensation =
                                 exp->bool_data( "diagonal_compensation" ) ;
   double const smallest_nonzero_pivot =
                              exp->double_data( "smallest_nonzero_pivot" ) ;
   exp->test_data( "smallest_nonzero_pivot", "smallest_nonzero_pivot>0." ) ;

   LA_ILU0_PC* result = new LA_ILU0_PC( a_owner,
                                        diagonal_compensation,
                                        smallest_nonzero_pivot ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ILU0_PC*
LA_ILU0_PC:: create( PEL_Object* a_owner,
                     bool diagonal_compensation,
                     double smallest_nonzero_pivot )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: create" ) ;
   PEL_CHECK_PRE( smallest_nonzero_pivot > 0. ) ;

   LA_ILU0_PC* result = new LA_ILU0_PC( a_owner, 
                                        diagonal_compensation,
                                        smallest_nonzero_pivot ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ILU0_PC:: LA_ILU0_PC( PEL_Object* a_owner,
                         bool diagonal_compensation,
                         double smallest_nonzero_pivot )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , MODIFIED( diagonal_compensation )
   , PIV_MIN( smallest_nonzero_pivot )
   , LU( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_ILU0_PC:: ~LA_ILU0_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_ILU0_PC*
LA_ILU0_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: create_clone" ) ;

   LA_ILU0_PC* result = new LA_ILU0_PC( a_owner, MODIFIED, PIV_MIN ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_ILU0_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: is_valid" ) ;
   return( LU != 0 ) ;
}

//----------------------------------------------------------------------
size_t
LA_ILU0_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( LU->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_ILU0_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   LA_SeqMatrix const* smat = dynamic_cast<LA_SeqMatrix const*>(mat) ;
   if( smat == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_ILU0_PC error:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected" ) ;
   }
   
   LU = LA_CRSmatrix::create( this, smat ) ;
   
   LU->factorize_MILU0( MODIFIED, PIV_MIN ) ;

   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_ILU0_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_ILU0_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   
   destroy_possession( LU ) ;
   LU = 0 ;
   
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_ILU0_PC:: solve( LA_Vector const* rhs, LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( rhs ) != 0 ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( sol ) != 0 ) ;

   LA_SeqVector const* brhs = static_cast<LA_SeqVector const*>( rhs ) ;
   LA_SeqVector* bsol = static_cast<LA_SeqVector*>( sol ) ;
   
   LU->solve_LU( brhs, bsol ) ;
   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_ILU0_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
void
LA_ILU0_PC:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ILU0_PC:: print_more" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string s( indent_width, ' ') ;

   os << s << "smallest_nonzero_pivot : " << PIV_MIN << std::endl ;
   os << s << "diagonal_compensation  : "
      << ( MODIFIED ? "true" : "false" ) << std::endl ;
}
