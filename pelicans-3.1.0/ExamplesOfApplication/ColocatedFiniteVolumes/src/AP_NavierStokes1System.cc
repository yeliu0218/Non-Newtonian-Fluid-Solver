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

#include <AP_NavierStokes1System.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>
#include <LA_TwoBlocksMethod.hh>
#include <LA_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_SystemNumbering.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

//----------------------------------------------------------------------
AP_NavierStokes1System*
AP_NavierStokes1System:: create( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp,
                                 PDE_LinkDOF2Unknown* uu_link,
                                 PDE_LinkDOF2Unknown* pp_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( uu_link != 0 ) ;
   PEL_CHECK_PRE( uu_link->components_table().size() > 1 ) ;
   PEL_CHECK_PRE( uu_link->owner() == 0 ) ;
   PEL_CHECK_PRE( pp_link != 0 ) ;
   PEL_CHECK_PRE( pp_link->components_table().size() == 1 ) ;
   PEL_CHECK_PRE( pp_link->owner() == 0 ) ;

   AP_NavierStokes1System* result = 
                new AP_NavierStokes1System( a_owner, exp, uu_link, pp_link ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_initialized() ) ; //??? oui ou non
   PEL_CHECK_POST( result->linkDOF2Unknown_U() == uu_link ) ;
   PEL_CHECK_POST( result->linkDOF2Unknown_P() == pp_link ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_NavierStokes1System:: AP_NavierStokes1System(
                                         PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp,
                                         PDE_LinkDOF2Unknown* uu_link,
                                         PDE_LinkDOF2Unknown* pp_link )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , METH( invalid )
   , HAS_INIT_U( false )
   , HAS_INIT_P( false )
   , BoverDT( PEL::bad_double() )
   , RR( 0.0 )
   , CONVERGED( false )
   , TOL_VELO( -PEL::max_double() )
   , TOL_DIV( -PEL::max_double() )
   , VERBOSE( exp->has_entry( "verbose_level" ) ? 
                 exp->int_data( "verbose_level" ) : 0 )
   , INDENT( "" )
   , NMB( 0 )
   , NMB_U( 0 )
   , NMB_P( 0 )
   , idx_U( 0 )
   , idx_P( 1 )
   , P( 0 )
   , A( 0 )
   , F( 0 )
   , B( 0 )
   , G( 0 )
   , S( 0 )
   , C( 0 )
   , U( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
   , P_LOC( LA_SeqVector::create( this, 0 ) )
   , SOLVER( 0 )
{
   //??? for the sake of clarity, do not use the same A for all cases
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "method" ) ;
   string const& mtype = se->string_data( "type" ) ;
   if( mtype == "monolithic" )
   {
      METH = MONO ;
      
      PEL_Vector* vec = PEL_Vector::create( 0, 2 ) ;
      vec->set_at( idx_U, uu_link ) ;
      vec->set_at( idx_P, pp_link ) ;
      std::string ordering = "sequence_of_the_discrete_fields" ;
      NMB = PDE_SystemNumbering::create( this, vec, ordering ) ;
      vec->destroy() ; vec = 0 ;

      F = A->create_vector( this ) ;
      U = P = A->create_vector( this ) ;   //???
      
      ee = se->create_subexplorer( 0, "LA_Solver") ;
      SOLVER= LA_Solver::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   else if( mtype == "augmented_Lagrangian" )
   {
      METH = AL ;
      std::string ordering = "sequence_of_the_discrete_fields" ; //??????
      NMB_U = PDE_SystemNumbering::create( this, uu_link ) ;
      NMB_P = PDE_SystemNumbering::create( this, pp_link ) ;

      U = A->create_vector( this ) ;
      P = A->create_vector( this ) ;
      
      B = A->create_matrix( this ) ; 
      C = A->create_matrix( this ) ; 
      S = A->create_vector( this ) ; 
      F = A->create_vector( this ) ; 
      G = A->create_vector( this ) ; 
      
      ee = se->create_subexplorer( 0, "LA_TwoBlocksMethod" ) ;
      SOLVER_TB = LA_TwoBlocksMethod::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
      
      SOLVER_TB->set_matrix_prototype( A ) ;
   } 
   else
   {
      PEL_Error::object()->raise_bad_data_value( se, "type", 
           "\"augmented_Lagrangian\" \"monolithic\""  ) ;
   }
   se->destroy() ; se = 0 ;
}

//----------------------------------------------------------------------
AP_NavierStokes1System:: ~AP_NavierStokes1System( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: re_initialize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: re_initialize" ) ;

   CONVERGED = false ;
   HAS_INIT_U = false ;
   HAS_INIT_P = false ;
   
   if( METH == MONO )
   {
      NMB->reset() ;
      
      size_t n_glob = NMB_U->nb_global_unknowns() ;
      size_t n_loc  = NMB_U->nb_unknowns_on_current_process() ;
      
      A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
      F->re_initialize( n_glob, n_loc ) ;
      U->re_initialize( n_glob, n_loc ) ;
      
      U_LOC->re_initialize( NMB->link( idx_U )->unknown_vector_size() ) ;
      P_LOC->re_initialize( NMB->link( idx_P )->unknown_vector_size() ) ;
      
      NMB->define_scatters( U ) ;
   }
   else
   {
      NMB_U->reset() ;
      NMB_P->reset() ;
      
      size_t nv_glob = NMB_U->nb_global_unknowns() ;
      size_t np_glob = NMB_P->nb_global_unknowns() ;

      size_t nv_loc = NMB_U->nb_unknowns_on_current_process() ;
      size_t np_loc = NMB_P->nb_unknowns_on_current_process() ;

      A->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
      U->re_initialize( nv_glob, nv_loc ) ;
      P->re_initialize( np_glob, np_loc ) ; 
      
      B->re_initialize( np_glob, nv_glob, np_loc, nv_loc ) ;
      C->re_initialize( np_glob, np_glob, np_loc, np_loc ) ;
      S->re_initialize( np_glob, np_loc ) ;
      F->re_initialize( nv_glob, nv_loc ) ;
      G->re_initialize( np_glob, np_loc ) ;
      
      SOLVER_TB->re_initialize_internals( nv_glob, np_glob, nv_loc, np_loc ) ;
      
      U_LOC->re_initialize( NMB_U->link()->unknown_vector_size() ) ;
      P_LOC->re_initialize( NMB_P->link()->unknown_vector_size() ) ;
      
      NMB_U->define_scatters( U ) ;
      NMB_P->define_scatters( P ) ;
   }
   
   INIT = true ;
   
   PEL_CHECK_POST( unknown_vector_U()->nb_rows() == 
                   linkDOF2Unknown_U()->unknown_vector_size() ) ;
   PEL_CHECK_POST( unknown_vector_P()->nb_rows() == 
                   linkDOF2Unknown_P()->unknown_vector_size() ) ;
}

//----------------------------------------------------------------------
bool
AP_NavierStokes1System:: is_initialized( void ) const
//----------------------------------------------------------------------
{
   return( INIT ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
AP_NavierStokes1System:: linkDOF2Unknown_U( void ) const
//----------------------------------------------------------------------
{
   PDE_LinkDOF2Unknown const* result = 
                       ( METH == MONO ?  NMB->link( idx_U ) : NMB_U->link() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
AP_NavierStokes1System:: linkDOF2Unknown_P( void ) const
//----------------------------------------------------------------------
{
   PDE_LinkDOF2Unknown const* result = 
                       ( METH == MONO ? NMB->link( idx_P ) :  NMB_P->link() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector const*
AP_NavierStokes1System:: unknown_vector_U( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: unknown_vector_U" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   if( METH == MONO )
   {
      NMB->scatter( idx_U )->get( U, U_LOC ) ;
   }
   else
   {
      NMB_U->scatter()->get( U, U_LOC ) ;
   }
   LA_SeqVector const* result = U_LOC ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_rows() == 
                   linkDOF2Unknown_U()->unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector const*
AP_NavierStokes1System:: unknown_vector_P( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: unknown_vector_P" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   if( METH == MONO )
   {
      NMB->scatter( idx_P )->get( P, P_LOC ) ;
   }
   else
   {
      NMB_P->scatter()->get( P, P_LOC ) ;
   }
   LA_SeqVector const* result = P_LOC ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_rows() == 
                   linkDOF2Unknown_P()->unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
AP_NavierStokes1System:: global_unknown_for_DOF_of_U( size_t n,
                                                      size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: global_unknown_for_DOF_of_U" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   size_t result = PEL::bad_index() ;
   if( METH == MONO )
   {
      result = NMB->global_unknown_for_DOF( n, ic, idx_U ) ;
   }
   else
   {
      result = NMB_U->global_unknown_for_DOF( n, ic ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
AP_NavierStokes1System:: global_unknown_for_DOF_of_P( size_t n ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: global_unknown_for_DOF_of_P" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   size_t result = PEL::bad_index() ;
   if( METH == MONO )
   {
      result = NMB->global_unknown_for_DOF( n, 0, idx_P ) ;
   }
   else
   {
      result = NMB_P->global_unknown_for_DOF( n, 0 ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: set_initial_guess_U( 
                                           PDE_DiscreteField const* uu,
                                           size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: set_initial_guess_U" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( uu != 0 ) ;
   PEL_CHECK_PRE( uu == linkDOF2Unknown_U()->field() ) ;
   PEL_CHECK_PRE( level < uu->storage_depth() ) ;
   
   if( METH == MONO )
   {
      uu->extract_unknown_DOFs_value( level, U_LOC, NMB->link( idx_U ) ) ;
      NMB->scatter( idx_U )->set( U_LOC, U ) ;
   }
   else
   {
      uu->extract_unknown_DOFs_value( level, U_LOC, NMB_U->link() ) ;
      NMB_U->scatter()->set( U_LOC, U ) ;
   }
   HAS_INIT_U = true ;
}

//----------------------------------------------------------------------
bool
AP_NavierStokes1System:: initial_guess_U_is_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_INIT_U ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: set_initial_guess_P( PDE_DiscreteField const* pp,
                                              size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: set_initial_guess_P" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( pp != 0 ) ;
   PEL_CHECK_PRE( pp == linkDOF2Unknown_P()->field() ) ;
   PEL_CHECK_PRE( level < pp->storage_depth() ) ;

   if( METH == MONO )
   {
      pp->extract_unknown_DOFs_value( level, P_LOC, NMB->link( idx_P ) ) ;
      // NMB->scatter( idx_P )->set( P_LOC, P0 ) ;
      NMB->scatter( idx_P )->set( P_LOC, P ) ;
   }
   else
   {
      pp->extract_unknown_DOFs_value( level, P_LOC, NMB_P->link() ) ;
//      NMB_P->scatter()->set( P_LOC, P0 ) ;
      NMB_P->scatter()->set( P_LOC, P ) ;
   }
   HAS_INIT_P = true ;
}

//----------------------------------------------------------------------
bool
AP_NavierStokes1System:: initial_guess_P_is_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_INIT_P ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: set_leading_BDF_over_dt( double value )
//----------------------------------------------------------------------
{
   BoverDT = value ; //??? never used 
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_A_item( size_t i_row, size_t j_col, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_A_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   A->add_to_item( i_row, j_col, xx ) ;

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_F_item( size_t i_row, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_F_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   F->add_to_item( i_row, xx ) ;

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_B_item( size_t i_row, size_t j_col, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_B_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   if( METH == MONO )
   {
      A->add_to_item( i_row, j_col, xx ) ;
      A->add_to_item( j_col, i_row, xx ) ;
   }
   else
   {
      B->add_to_item( i_row, j_col, xx ) ;
   }

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_G_item( size_t i_row, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_G_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   if( METH == MONO )
   {
      F->add_to_item( i_row, xx ) ;
   }
   else
   {
      G->add_to_item( i_row, xx ) ;
   }

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
bool
AP_NavierStokes1System:: MPl_is_required( void ) const
//----------------------------------------------------------------------
{
   bool result = ( METH == AL ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_MPl_item( size_t i_row, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_MPl_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   S->add_to_item( i_row, xx ) ;

//   MP_INVERTED = false ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: add_to_C_item( size_t i_row, size_t j_col, double xx )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: add_to_C_item" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   //??? more preconditions

   CONVERGED = false ;

   if( METH == MONO )
   {
      A->add_to_item( i_row, j_col, xx ) ;
   }
   else
   {
      C->add_to_item( i_row, j_col, xx ) ;
   }

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: estimate_unknowns( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: estimate_unknowns" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   switch( METH )
   {
      case AL :
         estimate_unknowns_AL() ;
         break ;
      case MONO :
         estimate_unknowns_MONO() ;
         break ;
      case invalid :
         break ;
   }
   HAS_INIT_U = false ;
   HAS_INIT_P = false ;
   
   PEL_CHECK_POST( !initial_guess_U_is_set() ) ;
   PEL_CHECK_POST( !initial_guess_P_is_set() ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: estimate_unknowns_AL( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: estimate_unknowns_AL" ) ;
   PEL_ASSERT( METH == AL ) ;
   
   if( S != 0 )
   {
      S->synchronize() ;
      S->set_as_reciprocal( S ) ;
      SOLVER_TB->set_S( S ) ;
   }

   A->synchronize() ;
   B->synchronize() ;
   C->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   
   SOLVER_TB->set_system( A, B, F, G, C ) ;
   SOLVER_TB->estimate_unknowns( HAS_INIT_U, U, HAS_INIT_P, P ) ;
   SOLVER_TB->unset_system() ;
   
   CONVERGED = SOLVER_TB->successful_estimation() ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: estimate_unknowns_MONO( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: estimate_unknowns_MONO" ) ;
   PEL_ASSERT( METH == MONO ) ;
   
   SOLVER->set_initial_guess_nonzero( HAS_INIT_U && HAS_INIT_P ) ;
   A->synchronize() ;
   F->synchronize() ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, U ) ;
   SOLVER->unset_matrix() ;

   CONVERGED = SOLVER->solution_is_achieved() ;
}

//----------------------------------------------------------------------
bool
AP_NavierStokes1System:: unknowns_are_solution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: unknowns_are_solution" ) ;

   return( CONVERGED ) ;
}

//----------------------------------------------------------------------
void
AP_NavierStokes1System:: set_indent( std::string const& indent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_NavierStokes1System:: set_indent" ) ;

   INDENT = indent ;
   if( SOLVER_TB != 0 ) SOLVER_TB->set_indent( indent ) ;
}
