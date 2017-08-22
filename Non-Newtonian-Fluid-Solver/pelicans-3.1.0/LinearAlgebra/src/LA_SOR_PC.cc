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

#include <LA_SOR_PC.hh>

#include <LA_MatrixIterator.hh>
#include <LA_SeqImplementation.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <stringVector.hh>

#include <sstream>

struct LA_SOR_PC_ERROR
{
   static void n0( int i ) ;
} ;

LA_SOR_PC const* LA_SOR_PC::PROTOTYPE = new LA_SOR_PC() ;

//----------------------------------------------------------------------
LA_SOR_PC:: LA_SOR_PC( void )
//----------------------------------------------------------------------
   : LA_Preconditioner( "LA_SOR_PC" )
   , MIN_DIAG( -PEL::max_double() )
   , OMEGA( -PEL::max_double() )
   , FORWARD( false )
   , BACKWARD( false )
   , NB_ITERS( PEL::bad_index() )
   , OMEGA_INV_DIAG( 0 )
   , MAT( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_SOR_PC*
LA_SOR_PC:: create_replica( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   double const smallest_inverted_item =
                         exp->double_data( "smallest_inverted_item" ) ;
   exp->test_data( "smallest_inverted_item", "smallest_inverted_item>0." ) ;
   
   double const omega = exp->double_data( "omega" ) ;
   exp->test_data( "omega", "omega>0. && omega<2." ) ;

   size_t const nb_iters = (size_t) exp->int_data( "nb_inner_iterations" ) ;
   exp->test_data( "nb_inner_iterations", "nb_inner_iterations>0" ) ;

   std::string const& sweep = exp->string_data( "sweep" ) ;
   exp->test_data_in( "sweep", "symmetric,forward,backward" ) ;
   
   LA_SOR_PC* result =
      new LA_SOR_PC( a_owner,
                     smallest_inverted_item,
                     omega,
                     ( sweep == "forward" || sweep == "symmetric" ),
                     ( sweep == "backward" || sweep == "symmetric" ),
                     nb_iters ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SOR_PC*
LA_SOR_PC:: create( PEL_Object* a_owner,
                    double omega,
                    std::string const& sweep,
                    size_t nb_iters,
                    double smallest_inverted_item )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: create" ) ;
   PEL_CHECK_PRE( omega>0 && omega<2 ) ;
   PEL_CHECK_PRE( sweep == "forward" ||
                  sweep == "backward" ||
                  sweep == "symmetric" ) ;
   PEL_CHECK_PRE( nb_iters>0 ) ;
   PEL_CHECK_PRE( smallest_inverted_item > 0 ) ;

   LA_SOR_PC* result =
      new LA_SOR_PC( a_owner,
                     smallest_inverted_item,
                     omega,
                     ( sweep == "forward" || sweep == "symmetric" ),
                     ( sweep == "backward" || sweep == "symmetric" ),
                     nb_iters ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SOR_PC:: LA_SOR_PC( PEL_Object* a_owner,
                       double smallest_inverted_item,
                       double omega,
                       bool forward,
                 bool backward,
                       size_t nb_iters )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , MIN_DIAG( smallest_inverted_item )
   , OMEGA( omega )
   , FORWARD( forward )
   , BACKWARD( backward )
   , NB_ITERS( nb_iters )
   , OMEGA_INV_DIAG( LA_SeqVector::create( this, 0 ) )
   , MAT( 0 )
   , SOLVE_OK( false )
{
}

//----------------------------------------------------------------------
LA_SOR_PC:: ~LA_SOR_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_SOR_PC*
LA_SOR_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: create_clone" ) ;

   LA_SOR_PC* result = new LA_SOR_PC(
               a_owner, MIN_DIAG, OMEGA, FORWARD, BACKWARD, NB_ITERS ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_SOR_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( MAT != 0 ) ;
}

//----------------------------------------------------------------------
size_t
LA_SOR_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
void
LA_SOR_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;
   
   MAT = dynamic_cast<LA_SeqMatrix const*>(mat) ;
   if( MAT == 0 )
   {
      PEL_Error::object()->raise_internal(
         "*** LA_SOR_PC error:\n"
         "    a matrix of type \"LA_SeqMatrix\" is expected" ) ;
   }

   size_t const n = MAT->nb_rows() ;

   if( OMEGA_INV_DIAG->nb_rows() != n )
   {
      OMEGA_INV_DIAG->re_initialize( n ) ;
   }
   MAT->extract_diag( OMEGA_INV_DIAG ) ;
   double* ptr_omega_inv_diag = OMEGA_INV_DIAG->data() ;
   for( size_t i=0 ; i<n ; ++i )
   {
      double d = ptr_omega_inv_diag[i] ;
      if( PEL::abs(d) <= MIN_DIAG )
      {
         d = 1. ;
         LA_SOR_PC_ERROR::n0( i ) ;
      }
      ptr_omega_inv_diag[i] = OMEGA/d ;
   }
   
   SOLVE_OK = false ;
   
   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_SOR_PC:: unbuild( void )
//----------------------------------------------------------------------
{   
   PEL_LABEL( "LA_SOR_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;

   MAT = 0 ;
   
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_SOR_PC:: solve( LA_Vector const* rhs, LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SOR_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( rhs ) != 0 ) ;
   LA_SeqVector const* seq_rhs = static_cast<LA_SeqVector const*>( rhs ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( sol ) != 0 ) ;
   LA_SeqVector* seq_sol = static_cast<LA_SeqVector*>( sol ) ;

   LA_SeqMatrix::relaxation_mode const m =
      ( FORWARD && BACKWARD ? LA_SeqMatrix::symmetric
                            : ( FORWARD ? LA_SeqMatrix::forward
                                        : LA_SeqMatrix::backward ) ) ;
   
   seq_sol->nullify() ;
   for( size_t iter=0 ; iter<NB_ITERS ; ++iter )
   {
      MAT->relax( OMEGA, m, OMEGA_INV_DIAG, seq_rhs, seq_sol ) ;
   }
   
   SOLVE_OK = true ;

   PEL_CHECK_POST( solve_POST( rhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;
}

//----------------------------------------------------------------------
bool
LA_SOR_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
void
LA_SOR_PC:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ') ;

   os << s << "omega: " << OMEGA << std::endl ;
   os << s << "smallest_inverted_item: " << MIN_DIAG << std::endl ;
   os << s << "nb_inner_iterations: " << NB_ITERS << std::endl ;
   os << s << "sweed: "
      << ( FORWARD && BACKWARD ? "\"symmetric\""
                               : ( FORWARD ? "\"forward\"" : "\"backward\"" ) )
      << std::endl ;
}

//internal--------------------------------------------------------------
void 
LA_SOR_PC_ERROR:: n0( int i )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** SOR preconditioner :" << std::endl ;
   mesg << "***   vanishing diagonal term at line " << i << std::endl ;
   PEL_Error::object()->display_info( mesg.str() ) ;
}
