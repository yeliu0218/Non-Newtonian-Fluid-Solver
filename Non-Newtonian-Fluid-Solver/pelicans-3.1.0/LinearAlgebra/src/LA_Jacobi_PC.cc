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

#include <LA_Jacobi_PC.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL.hh>

#include <sstream>

struct LA_Jacobi_PC_ERROR
{
   static void n0( int i ) ;
} ;

LA_Jacobi_PC const* LA_Jacobi_PC::PROTOTYPE = new LA_Jacobi_PC() ;

//----------------------------------------------------------------------
LA_Jacobi_PC:: LA_Jacobi_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_Jacobi_PC" )
   , MIN_DIAG( -PEL::max_double() )
   , INV_DIAG( 0 )
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_Jacobi_PC*
LA_Jacobi_PC:: create_replica( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double a_smallest_item = exp->double_data( "smallest_inverted_item" ) ;
   if( a_smallest_item<=0. )
   {
      PEL_Error::object()->raise_bad_data_value( exp,
                                                 "smallest_inverted_item",
                                                 "greater than zero" ) ;
   }
   LA_Jacobi_PC* result = LA_Jacobi_PC::create( a_owner, a_smallest_item ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Jacobi_PC*
LA_Jacobi_PC:: create( PEL_Object* a_owner,
                       double smallest_inverted_item )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: create" ) ;
   PEL_CHECK_PRE( smallest_inverted_item > 0.0 ) ;

   LA_Jacobi_PC* result = new LA_Jacobi_PC( a_owner,
                                            smallest_inverted_item ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->successful_solve() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Jacobi_PC:: LA_Jacobi_PC( PEL_Object* a_owner,
                             double smallest_inverted_item )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , MIN_DIAG( smallest_inverted_item )
   , INV_DIAG(0)
   , BUILD_OK( false )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_Jacobi_PC:: ~LA_Jacobi_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Jacobi_PC*
LA_Jacobi_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: create_clone" ) ;

   LA_Jacobi_PC* result = new LA_Jacobi_PC( a_owner, MIN_DIAG ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_Jacobi_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: is_valid" ) ;
   return( BUILD_OK ) ;
}

//----------------------------------------------------------------------
size_t
LA_Jacobi_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( INV_DIAG->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_Jacobi_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   if( INV_DIAG==0 ) INV_DIAG = mat->create_vector( this ) ;
   INV_DIAG->re_initialize( mat->nb_cols(),
                            mat->nb_local_cols() ) ;

   mat->extract_diag( INV_DIAG ) ;
   INV_DIAG->set_as_reciprocal( INV_DIAG, MIN_DIAG, 1.0 ) ;

   BUILD_OK = true ;
   SOLVE_OK = false ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_Jacobi_PC:: unbuild( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   BUILD_OK = false ;
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Jacobi_PC:: solve( LA_Vector const* rhs,
                      LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   sol->set_as_v_product( rhs, INV_DIAG ) ;
   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_Jacobi_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Jacobi_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//internal--------------------------------------------------------------
void
LA_Jacobi_PC_ERROR:: n0( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** Jacobi preconditioner :" << std::endl ;
   mesg << "***   vanishing diagonal term at line " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
