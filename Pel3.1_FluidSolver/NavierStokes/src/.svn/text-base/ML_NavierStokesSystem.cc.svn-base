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

#include <ML_NavierStokesSystem.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <PDE.hh>
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

#include <fstream>

//----------------------------------------------------------------------
ML_NavierStokesSystem*
ML_NavierStokesSystem:: create( PEL_Object* a_owner, 
                                   PEL_ModuleExplorer const* exp,
                                   PDE_LinkDOF2Unknown* uu_link,
                                   PDE_LinkDOF2Unknown* pp_link,
                                   PDE_LinkDOF2Unknown* tt_link,
                                   PDE_LinkDOF2Unknown* gamma_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( uu_link != 0 ) ;
   PEL_CHECK_PRE( uu_link->components_table().size() > 1 ) ;
   PEL_CHECK_PRE( uu_link->owner() == 0 ) ;
   PEL_CHECK_PRE( pp_link != 0 ) ;
   PEL_CHECK_PRE( pp_link->components_table().size() == 1 ) ;
   PEL_CHECK_PRE( pp_link->owner() == 0 ) ;

   PEL_CHECK_PRE( tt_link != 0 ) ;
   PEL_CHECK_PRE( tt_link->components_table().size() == 3 || tt_link->components_table().size() == 6 ) ;
   PEL_CHECK_PRE( tt_link->owner() == 0 ) ;
   PEL_CHECK_PRE( gamma_link != 0 ) ;
   PEL_CHECK_PRE( gamma_link->components_table().size() == 3 || gamma_link->components_table().size() == 6 ) ;
   PEL_CHECK_PRE( gamma_link->owner() == 0 ) ;
   
   ML_NavierStokesSystem* result = 
                new ML_NavierStokesSystem( a_owner, exp, uu_link, pp_link, tt_link, gamma_link ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_initialized() ) ; //??? oui ou non
   PEL_CHECK_POST( result->system_numbering_U()->link() == uu_link ) ;
   PEL_CHECK_POST( result->system_numbering_P()->link() == pp_link ) ;
   PEL_CHECK_POST( result->system_numbering_T()->link() == tt_link ) ;
   PEL_CHECK_POST( result->system_numbering_GAMMA()->link() == gamma_link ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
ML_NavierStokesSystem:: ML_NavierStokesSystem(
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp,
                                    PDE_LinkDOF2Unknown* uu_link,
                                    PDE_LinkDOF2Unknown* pp_link,
                                    PDE_LinkDOF2Unknown* tt_link,
                                    PDE_LinkDOF2Unknown* gamma_link )
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
   , INIT( false )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
   , P_LOC( LA_SeqVector::create( this, 0 ) )
   , T_LOC( LA_SeqVector::create( this, 0 ) )
   , SOLVER_A( 0 )
   , SOLVER_L( 0 )
   , SOLVER_Mv( 0 )
   , MP_INVERTED( true )
   , LHS_flag( true )
{
   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "method" ) ;
   string const& mtype = se->string_data( "type" ) ;
   if( mtype == "augmented_Lagrangian" )
   {
      METH = AL ;
      RR = se->double_data( "augmentation_parameter" ) ;
      TOL_VELO = se->double_data( "tolerance_on_velocity_increment" ) ;
      TOL_DIV = se->double_data( "tolerance_on_divergence" ) ;

   } 
   else if( mtype == "Yosida" )
   {
      METH = YOS ;

   }
   else if( mtype == "penalty_projection" )
   {
      METH = PP ;
      RR = se->double_data( "augmentation_parameter" ) ;
      se->test_data( "augmentation_parameter",
                     "augmentation_parameter>=0." ) ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value( se, "type", 
           "\"augmented_Lagrangian\" \"Yosida\" \"augmented_projection\""  ) ;
   }
   if( METH == YOS || METH == PP ) 
   {
      PEL_ModuleExplorer const* ee = se->create_subexplorer( 0, "solver_L") ;
      SOLVER_L = LA_Solver::make( this, ee ) ;
      ee->destroy() ;
      ee = se->create_subexplorer( 0, "solver_Mv") ;
      SOLVER_Mv = LA_Solver::make( this, ee ) ;
      ee->destroy() ;
   }
   se->destroy() ; se = 0 ;

   NMB_U = PDE_SystemNumbering::create( this, uu_link ) ;
   NMB_P = PDE_SystemNumbering::create( this, pp_link ) ;
   NMB_T = PDE_SystemNumbering::create( this, tt_link ) ;
   NMB_GAMMA = PDE_SystemNumbering::create( this, gamma_link ) ;
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ; //??? LHS( V_UNK, V_UNK, 0 )
   ee->destroy() ; ee = 0 ;
   
   A_explicit = A->create_matrix( this ) ; //??? LHS( V_UNK, V_UNK, 2 ) ; 
   M = A->create_matrix( this ) ; //??? LHS( V_UNK, V_UNK, 1 ) ;
   BtMpInv = A->create_matrix( this ) ; //??? SYS->LHS( V_UNK, P_UNK, 0 ) ;
   B = A->create_matrix( this ) ; //??? SYS->LHS( P_UNK, V_UNK, 0 ) ;
   L = A->create_matrix( this ) ; //??? SYS->LHS( P_UNK, P_UNK, 0 ) ;
   
   U  = A->create_vector( this ) ; //??? SYS->unknown( V_UNK, 0 ) ;
   DU = A->create_vector( this ) ; //??? SYS->unknown( V_UNK, 1 ) ;
   P  = A->create_vector( this ) ; //??? SYS->unknown( P_UNK, 0 ) ;
   P0 = A->create_vector( this ) ; //??? SYS->unknown( P_UNK, 1 ) ;
   
   F    = A->create_vector( this ) ; //??? SYS->RHS( V_UNK, 0 ) ;
   F0   = A->create_vector( this ) ; //??? SYS->RHS( V_UNK, 1 ) ;
   F_explicit = A->create_vector( this ) ; //??? SYS->RHS( V_UNK, 2 ) ;
   G    = A->create_vector( this ) ; //??? SYS->RHS( P_UNK, 0 ) ;
   Mpl  = A->create_vector( this ) ; //??? SYS->RHS( P_UNK, 1 ) ;
   PRES = A->create_vector( this ) ; //??? SYS->RHS( P_UNK, 2 ) ;
   
   A_Outer = A->create_matrix( this ) ;
   A_Inner = A->create_matrix( this ) ; 
   F_Outer = A->create_vector( this ) ;
   F_Inner = A->create_vector( this ) ;
   F_Flowrate = A->create_vector( this ) ;
   F_Stress = A->create_vector( this ) ;
      
   // visco-plastic problem
   F_Visco = A->create_vector( this ) ;
   T  = A->create_vector( this ) ;
   A_TtimesGAMMADOTv = A->create_matrix( this ) ;
   TminusRGAMMA = A->create_vector( this ) ;
   
   ee = exp->create_subexplorer( 0, "solver_A" ) ;
   SOLVER_A = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
ML_NavierStokesSystem:: ~ML_NavierStokesSystem( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: re_initialize" ) ;

   CONVERGED = false ;
   HAS_INIT_U = false ;
   HAS_INIT_P = false ;

   NMB_U->reset() ;
   NMB_P->reset() ;
   NMB_T->reset() ;
   NMB_GAMMA->reset() ;

   size_t nv_glob = NMB_U->nb_global_unknowns() ;
   size_t np_glob = NMB_P->nb_global_unknowns() ;
   size_t nt_glob = NMB_T->nb_global_unknowns() ;

   size_t nv_loc = NMB_U->nb_unknowns_on_current_process() ;
   size_t np_loc = NMB_P->nb_unknowns_on_current_process() ;
   size_t nt_loc = NMB_T->nb_unknowns_on_current_process() ;

   re_initialize_matrices_and_vectors( nv_glob, np_glob, nt_glob, nv_loc, np_loc, nt_loc ) ;

   U_LOC->re_initialize( NMB_U->link()->unknown_vector_size() ) ;
   P_LOC->re_initialize( NMB_P->link()->unknown_vector_size() ) ;
   T_LOC->re_initialize( NMB_T->link()->unknown_vector_size() ) ;

   NMB_U->define_scatters( U ) ;
   NMB_P->define_scatters( P ) ;
   NMB_T->define_scatters( T ) ;
   
   INIT = true ;

   PEL_CHECK_POST( is_initialized() ) ;
   PEL_CHECK_POST( unknown_vector_U()->nb_rows() == 
                   system_numbering_U()->link()->unknown_vector_size() ) ;
   PEL_CHECK_POST( unknown_vector_P()->nb_rows() == 
                   system_numbering_P()->link()->unknown_vector_size() ) ;
   PEL_CHECK_POST( unknown_vector_T()->nb_rows() == 
                   system_numbering_T()->link()->unknown_vector_size() ) ;                   
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: is_initialized( void ) const
//----------------------------------------------------------------------
{
   return( INIT ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering const*
ML_NavierStokesSystem:: system_numbering_U( void ) const
//----------------------------------------------------------------------
{
   PDE_SystemNumbering const* result = NMB_U ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering const*
ML_NavierStokesSystem:: system_numbering_P( void ) const
//----------------------------------------------------------------------
{
   PDE_SystemNumbering const* result = NMB_P ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering const*
ML_NavierStokesSystem:: system_numbering_T( void ) const
//----------------------------------------------------------------------
{
   PDE_SystemNumbering const* result = NMB_T ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering const*
ML_NavierStokesSystem:: system_numbering_GAMMA( void ) const
//----------------------------------------------------------------------
{
   PDE_SystemNumbering const* result = NMB_GAMMA ;
   return( result ) ;
}


//----------------------------------------------------------------------
LA_SeqVector const*
ML_NavierStokesSystem:: unknown_vector_U( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: unknown_vector_U" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   NMB_U->scatter()->get( U, U_LOC ) ;
   LA_SeqVector const* result = U_LOC ;

   PEL_CHECK_POST( result->nb_rows() == 
                   system_numbering_U()->link()->unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector const*
ML_NavierStokesSystem:: unknown_vector_P( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: unknown_vector_P" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   NMB_P->scatter()->get( P, P_LOC ) ;
   LA_SeqVector const* result = P_LOC ;

   PEL_CHECK_POST( result->nb_rows() == 
                   system_numbering_P()->link()->unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector const*
ML_NavierStokesSystem:: unknown_vector_T( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: unknown_vector_T" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   NMB_T->scatter()->get( T, T_LOC ) ;
   LA_SeqVector const* result = T_LOC ;

   PEL_CHECK_POST( result->nb_rows() == 
                   system_numbering_T()->link()->unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: nullify_LHS_and_RHS( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: nullify_LHS_and_RHS" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

// NOT used so far

   CONVERGED = false ;
   HAS_INIT_U = false ;
   HAS_INIT_P = false ;
   
   //??? some of these nullify may be useless
   A->nullify() ;
   M->nullify() ;
   BtMpInv->nullify() ;
   B->nullify() ;
   L->nullify() ;

   U->nullify() ;
   DU->nullify() ;
   P->nullify() ;
   P0->nullify() ;
   
   F->nullify() ;
   F0->nullify() ;
   F_explicit->nullify() ;
   G->nullify() ;
   Mpl->nullify() ;
   PRES->nullify() ;
   
   A_Outer->nullify() ;
   F_Outer->nullify() ;
   A_Inner->nullify() ;
   F_Inner->nullify() ;
   
   F_Flowrate->nullify() ;
   F_Stress->nullify() ;
   F_Visco->nullify() ;
   
   T->nullify() ;
   A_TtimesGAMMADOTv->nullify() ;
   TminusRGAMMA->nullify() ;
   
   PEL_CHECK_POST( !initial_guess_U_is_set() ) ;
   PEL_CHECK_POST( !initial_guess_P_is_set() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: nullify_for_new_internal_iteration( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: nullify_for_new_internal_iteration" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   CONVERGED = false ;
   HAS_INIT_U = false ;
   HAS_INIT_P = false ;
      
   A->nullify() ;
   //BtMpInv->nullify() ; see augment_system() and set_LHS  
   M->nullify() ;
   
   F->nullify() ;
   F0->nullify() ;

   PRES->nullify() ;
   
   U->nullify() ;
   DU->nullify() ;
   P->nullify() ;
   P0->nullify() ;
   
   T->nullify() ;
   F_Visco->nullify() ;
   TminusRGAMMA->nullify() ;
    
   PEL_CHECK_POST( !initial_guess_U_is_set() ) ;
   PEL_CHECK_POST( !initial_guess_P_is_set() ) ;
}


//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: nullify_A_F_Inner( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: nullify_A_F_Inner" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
      
   A_Inner->nullify() ; 
   F_Inner->nullify() ;
}


//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_explicit_terms( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: copy_explicit_terms" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   A_explicit->synchronize() ;
   A->add_Mat( A_explicit, 1.0 ) ;

   F_explicit->synchronize() ;
   F->sum( F_explicit, 1.0 ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_Outer_terms( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: add_Outer_terms" ) ;
   
   A_Outer->synchronize() ;
   A->add_Mat( A_Outer, 1.0 ) ;

   F_Outer->synchronize() ;   
   F->sum( F_Outer, 1.0) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_Inner_terms( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: add_Outer_terms" ) ;
   
   A_Inner->synchronize() ;
   A->add_Mat( A_Inner, 1.0 ) ;
   
   F_Inner->synchronize() ;   
   F->sum( F_Inner, 1.0 ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_Flowrate_terms( double inner_rhs_scale )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: add_Flowrate_terms" ) ;
   
   F_Flowrate->synchronize() ;   
   F->sum( F_Flowrate, inner_rhs_scale ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_BCstress_terms( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: add_BCstress_terms" ) ;
   
   F_Stress->synchronize() ; 
   F->sum( F_Stress, 1.0 ) ;
}

/** @brief add visco-plastic terms
 * \f$ [T-r \gamma] \dot\gamma \f$
 */
//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: add_Visco_terms( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: add_Visco_terms" ) ;
    
   // multiply A_TtimesGAMMADOTv with TminusRGAMMA and write into RHS F_Visco
    A_TtimesGAMMADOTv->multiply_vec_then_add( TminusRGAMMA, F_Visco, 1.0, 1.0) ;    
    F_Visco->synchronize() ;
    F->sum( F_Visco, 1.0 ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: set_initial_guess_U( 
                                           PDE_DiscreteField const* uu,
                                           size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: set_initial_guess_U" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( uu != 0 ) ;
   PEL_CHECK_PRE( uu == system_numbering_U()->link()->field() ) ;
   PEL_CHECK_PRE( level < uu->storage_depth() ) ;
   
   uu->extract_unknown_DOFs_value( level, U_LOC, NMB_U->link() ) ;
   NMB_U->scatter()->set( U_LOC, U ) ;
   HAS_INIT_U = true ;

   PEL_CHECK_POST( initial_guess_U_is_set() ) ;
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: initial_guess_U_is_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_INIT_U ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: set_initial_guess_P( 
                                           PDE_DiscreteField const* pp,
                                           size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: set_initial_guess_P" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( pp != 0 ) ;
   PEL_CHECK_PRE( pp == system_numbering_P()->link()->field() ) ;
   PEL_CHECK_PRE( level < pp->storage_depth() ) ;

   pp->extract_unknown_DOFs_value( level, P_LOC, NMB_P->link() ) ;
   NMB_P->scatter()->set( P_LOC, P0 ) ;
   HAS_INIT_P = true ;

   PEL_CHECK_POST( initial_guess_P_is_set() ) ;
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: initial_guess_P_is_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_INIT_P ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: set_leading_BDF_over_dt( double value )
//----------------------------------------------------------------------
{
   BoverDT = value ;
}

/** @brief Set TminusRGAMMA to 
 * \f$ T-r \gamma \f$. Needed in step one of the Augmented Lagrangian algorithm.
 */
//----------------------------------------------------------------------
void
ML_NavierStokesSystem::set_TminusRGAMMA( double Lagrange_parameter)
//----------------------------------------------------------------------
{
    PEL_LABEL( "ML_NavierStokesSystem:: set_TminusRGAMMA" ) ;
    
    size_t Lambda_size = system_numbering_T()->link()->unknown_vector_size();
    
    LA_SeqVector* VL_tmp = LA_SeqVector::create( 0, Lambda_size ) ;       
    PDE_DiscreteField const* vl_field = system_numbering_T()->link()->field() ;    
    vl_field->extract_unknown_DOFs_value( 0, VL_tmp, NMB_T->link() ) ;
        
    LA_SeqVector* d_tmp = LA_SeqVector::create( 0, Lambda_size ) ;
    PDE_DiscreteField const* d_field = NMB_GAMMA->link()->field() ;
    d_field->extract_unknown_DOFs_value( 0, d_tmp, NMB_GAMMA->link() ) ;
        
    LA_SeqVector* TminusRGAMMA_tmp = LA_SeqVector::create( 0, Lambda_size ) ;
    for (size_t i=0; i < Lambda_size; ++i)
        TminusRGAMMA_tmp->set_item( i, VL_tmp->item(i)-Lagrange_parameter
                                        * d_tmp->item(i));
                                        
    NMB_T->scatter()->set( TminusRGAMMA_tmp, TminusRGAMMA ) ;
                                        
    VL_tmp->destroy() ; 
    d_tmp->destroy() ;  
    TminusRGAMMA_tmp->destroy() ;                             
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_A_F( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_A_F" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_vector_0( A, F, leq, NMB_U ) ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}


//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_A_F_explicit( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_A_F_explicit" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_vector_0( A_explicit, F_explicit, leq, NMB_U ) ;
  
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_A_F_Outer( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_A_F_Outer" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_vector_0( A_Outer, F_Outer, leq, NMB_U ) ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_A_F_Inner( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_A_F_Inner" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_vector_0( A_Inner, F_Inner, leq, NMB_U ) ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_F_Flowrate( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_F_Flowrate" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==  
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;
   
   PDE::assemble_in_vector_1( F_Flowrate, leq, NMB_U ) ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_F_Stress( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_F_Stress" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==  
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;
   
   PDE::assemble_in_vector_1( F_Stress, leq, NMB_U ) ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_B_G( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_B_G" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_vector_0( B, G, leq, NMB_P, NMB_U ) ;

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_A_LagrangeMultDivMatrix( 
                                            PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
    PEL_LABEL( "ML_NavierStokesSystem:: assemble_A_LagrangeMultDivMatrix" ) ;
    PEL_CHECK_PRE( leq != 0 ) ;
    PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                   system_numbering_U()->link()->components_table().size() ) ;
    PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                   system_numbering_T()->link()->components_table().size() ) ;

    PDE::assemble_in_matrix_0( A_TtimesGAMMADOTv, leq, NMB_U, NMB_T ) ;
    PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: MPl_is_required( void ) const
//----------------------------------------------------------------------
{
   bool result = ( METH==AL || METH==PP ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_MPl( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_MPl" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 1 ) ;

   CONVERGED = false ;

   // the local matrix is already lumped
   for( size_t il=0 ; il<leq->nb_rows() ; il++ )
   {
      size_t nn = leq->row_node( il ) ;
      if( NMB_P->link()->DOF_is_unknown( nn ) )
      {
         int ig = NMB_P->global_unknown_for_DOF( nn, 0, 0 ) ;
         double aii = leq->matrix_item( il, il ) ;
         Mpl->add_to_item( ig, aii ) ;
      }
   }     
   MP_INVERTED = false ;
   
   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: L_is_required( void ) const
//----------------------------------------------------------------------
{
   bool result = ( ( METH == YOS ) || ( METH == PP ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_L( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_L" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 1 ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_0( L, leq, NMB_P ) ;

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: MV_is_required( void ) const
//----------------------------------------------------------------------
{
   bool result = ( METH == PP ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: assemble_MV( PDE_LocalEquation const* leq )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: assemble_MV" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  system_numbering_U()->link()->components_table().size() ) ;

   CONVERGED = false ;

   PDE::assemble_in_matrix_0( M, leq, NMB_U ) ;

   PEL_CHECK_POST( !unknowns_are_solution() ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: estimate_unknowns( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: estimate_unknowns" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   switch( METH )
   {
      case AL :
         estimate_unknowns_AL() ;
         break ;
      case YOS :
         estimate_unknowns_YOS() ;
         break ;
      case PP :
         estimate_unknowns_PP() ;
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
ML_NavierStokesSystem:: augment_system( double augmentation_parameter, bool set_LHS )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: augment_system" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   
   if( augmentation_parameter > 0.0 )
   {
      if( set_LHS ) 
      {       
        inverseP() ;

        // Bt = transpose( B )
        BtMpInv->nullify() ; 
        BtMpInv->add_tMat( B ) ;
        BtMpInv->synchronize() ;
      
        // BtMPInv = Bt.INVP
        BtMpInv->scale_as_mat_diag_mat( Mpl ) ;
   
        // BtMPInv = r.Bt.INVP
        BtMpInv->scale( augmentation_parameter ) ;
        

        if( VERBOSE > 1) PEL::out() << "----- LHS done -----" << std::endl;       
      }
       
      // A += BBT.B
      A->add_Mat_Mat( BtMpInv, B, 1.0 ) ;     
      A->synchronize() ;
      
      // F += BtMPInv.G
      BtMpInv->multiply_vec_then_add( G, F, 1.0, 1.0 ) ;
   }
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: inverseP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: inverseP" ) ;
   if( !MP_INVERTED )
   {
      Mpl->synchronize() ;
      Mpl->set_as_reciprocal( Mpl ) ;
      MP_INVERTED = true ;
   }
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: estimate_unknowns_AL( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: estimate_unknowns_AL" ) ;
   
   A->synchronize() ;
   B->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   BtMpInv->synchronize() ;
   
   double uvar = PEL::max_double() ;
   double delta_P = PEL::max_double() ;

   augment_system( RR, LHS_flag ) ;
   
   // F0 = F
   F0->set( F ) ;
   
   CONVERGED = false ;
   size_t nb_it = PEL::bad_index() ;
   size_t n = 0 ;
   bool failed = false ;
   do
   {
      ++n ;      
      // rhsA = F
      F->set( F0 ) ;
      
      // rhsA -= BT.P
      B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      
      double u_old_norm = U->two_norm() ;

      // DU = U
      DU->set( U ) ;
      
      // U = A-1.F
      SOLVER_A->set_initial_guess_nonzero( HAS_INIT_U ) ;
      SOLVER_A->set_matrix( A ) ;
      SOLVER_A->solve( F, U ) ;
      SOLVER_A->unset_matrix() ;

      failed = ( ! SOLVER_A->solution_is_achieved() ) ;

      if( !failed )
      {
         if( SOLVER_A->is_iterative() ) nb_it = SOLVER_A->nb_iterations_achieved() ;
      
         // DU = U - DU
         DU->scale( -1.0 ) ;
         DU->sum( U ) ;
      
         uvar = DU->two_norm() ;

         if( u_old_norm > 1.e-12 ) uvar /= u_old_norm ;
         bool cv1 = ( uvar < TOL_VELO ) ;

         // RES = -G
         PRES->set(G) ;

         // RES = B.U - G
         B->multiply_vec_then_add( U, PRES, 1.0, -1.0 ) ;
                                     
         // RES = INVP*RES
         PRES->set_as_v_product( PRES, Mpl ) ;
      
         // P0 += r.RES
         if( RR != 0.0 )
            P0->sum( PRES, RR ) ;

         delta_P = PRES->two_norm() ;
      
         bool cv2 = ( delta_P < TOL_DIV ) ;
         CONVERGED = ( cv1 && cv2 ) ;

         if( VERBOSE > 1 )
            print_errors_AL( INDENT, n, delta_P, uvar, nb_it ) ;
      }

   } while( !failed && !CONVERGED ) ;

   if( CONVERGED )
   {
      // P = P0
      P->set( P0 ) ;
      
      if( VERBOSE==1 ) print_end_errors_AL( INDENT, n, delta_P, uvar ) ;
   }
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: estimate_unknowns_YOS( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: estimate_unknowns_YOS" ) ;

   bool failed = false ;
   A->synchronize() ;
   B->synchronize() ;
   L->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   U->synchronize() ;

   if( ! failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 1 " ;
      // F0 = F
      F0->set( F ) ;
      
      // F -= BT.P0
      B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      
      // U = A-1.F
      SOLVER_A->set_initial_guess_nonzero( HAS_INIT_U ) ;
      SOLVER_A->set_matrix( A ) ;
      SOLVER_A->solve( F, U ) ;
      
      if( ! SOLVER_A->solution_is_achieved( ) )
      {
         failed = true ;
      }
      else
      { 
         if( VERBOSE > 0 && SOLVER_A->is_iterative() )
            PEL::out() << ": " << SOLVER_A->nb_iterations_achieved() 
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
      }
   }

   if( ! failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 2 " ;

      // G = BU - G
      B->multiply_vec_then_add( U, G, 1.0, -1.0 ) ;
      
      // P = L-1.G
      SOLVER_L->set_initial_guess_nonzero( HAS_INIT_P ) ;
      SOLVER_L->set_matrix( L ) ;
      SOLVER_L->solve( G, P ) ;
      SOLVER_L->unset_matrix() ;
      if( ! SOLVER_L->solution_is_achieved() )
      {
         failed = true ;
      }
      else
      {
         // P = P0 + beta*P/DT
         P->scale( BoverDT ) ;
         P->sum( P0 ) ;
         
         if( VERBOSE > 0 && SOLVER_L->is_iterative()  )
            PEL::out() << ": " << SOLVER_L->nb_iterations_achieved() 
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
      }
   }

   if( ! failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 3 " ;
      // F = F0
      F->set( F0 ) ;
      
      // F -= BT.P
      B->tr_multiply_vec_then_add( P, F, -1.0, 1.0 ) ;
      
      // U = A-1.F
      SOLVER_A->solve( F, U ) ;
      if( ! SOLVER_A->solution_is_achieved() )
      {
         failed = true ;
      }
      else
      {
         if( VERBOSE > 0 && SOLVER_A->is_iterative() )
            PEL::out() << ": " << SOLVER_A->nb_iterations_achieved()
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
         CONVERGED = true ;        
      }
   }
   SOLVER_A->unset_matrix() ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: estimate_unknowns_PP( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: estimate_unknowns_PP" ) ;

   bool failed = false ;
   A->synchronize() ;
   B->synchronize() ;
   M->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   L->synchronize() ;
   U->synchronize() ;

   if( ! failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 1 : " ;
      augment_system( RR, true ) ;
      // F -= BT.P0
      B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      
      // U = A-1.F
      SOLVER_A->set_initial_guess_nonzero( HAS_INIT_U ) ;
      SOLVER_A->set_matrix( A ) ;
      SOLVER_A->solve( F, U ) ;
      SOLVER_A->unset_matrix() ;
      if( ! SOLVER_A->solution_is_achieved() )
      {
         failed = true ;
      }
      else
      {
         if( VERBOSE > 0 && SOLVER_A->is_iterative() )
            PEL::out() << ": " << SOLVER_A->nb_iterations_achieved()
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
      }
   }

   if( ! failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 2 " ;
      // G = B.U -G
      B->multiply_vec_then_add( U, G, 1.0, -1.0 ) ;
      
      // P = L-1.G
      SOLVER_L->set_initial_guess_nonzero( false ) ;
      SOLVER_L->set_matrix( L ) ;
      SOLVER_L->solve( G, P ) ;
      SOLVER_L->unset_matrix() ;
      if( ! SOLVER_L->solution_is_achieved() )
      {
         failed = true ;
      }
      else
      {
         if( VERBOSE > 0 && SOLVER_L->is_iterative() )
            PEL::out() << ": " << SOLVER_L->nb_iterations_achieved()
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
      }
   }

   if( !failed )
   {
      if( VERBOSE > 0 ) PEL::out() << INDENT << "   step 3 " ;
      // F0 = M.U
      M->multiply_vec_then_add( U, F0, 1.0, 0.0 ) ;
      
      // F0 -= BT.PHI
      B->tr_multiply_vec_then_add( P, F0, -1.0, 1.0 ) ;
      
      // U = M-1*F0
      SOLVER_Mv->set_initial_guess_nonzero( true ) ;
      SOLVER_Mv->set_matrix( M ) ;
      SOLVER_Mv->solve( F0, U ) ;
      SOLVER_Mv->unset_matrix() ;
      
      if( ! SOLVER_Mv->solution_is_achieved() )
      {
         failed = true ;
      }
      else
      {
         // P = P0 + beta*PHI/DT
         P->scale( BoverDT ) ;
         P->sum(P0) ;
         
         if( RR != 0.0 )
         {
            // G = INVP*G
            G->set_as_v_product( G, Mpl ) ;
            
            // P = P + RR*G
            P->sum( G, RR ) ;
         }

         if( VERBOSE > 0 && SOLVER_Mv->is_iterative() )
            PEL::out() << ": " << SOLVER_Mv->nb_iterations_achieved()
                       << " iterations" ;
         if( VERBOSE > 0 ) PEL::out() << std::endl ;
         CONVERGED = true ;
      }
   }
}

//----------------------------------------------------------------------
bool
ML_NavierStokesSystem:: unknowns_are_solution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: unknowns_are_solution" ) ;

   return( CONVERGED ) ;
}

//----------------------------------------------------------------------
void
ML_NavierStokesSystem:: set_indent( std::string const& indent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: set_indent" ) ;

   INDENT = indent ;
}

//-----------------------------------------------------------------------
void
ML_NavierStokesSystem:: re_initialize_matrices_and_vectors( size_t nv_glob,
                                                               size_t np_glob,
                                                               size_t nt_glob,                                                               
                                                               size_t nv_loc,
                                                               size_t np_loc, 
                                                               size_t nt_loc)
//-----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: re_initialize_matrices_and_vectors" ) ;

   A->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   A_explicit->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   A_Outer->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   A_Inner->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   
   B->re_initialize( np_glob, nv_glob, np_loc, nv_loc ) ;
   BtMpInv->re_initialize( nv_glob, np_glob, nv_loc, np_loc ) ;

   U->re_initialize( nv_glob, nv_loc ) ;
   DU->re_initialize( nv_glob, nv_loc ) ;
   P->re_initialize( np_glob, np_loc ) ;
   P0->re_initialize( np_glob, np_loc ) ;

   F->re_initialize( nv_glob, nv_loc ) ;
   F0->re_initialize( nv_glob, nv_loc ) ;
   F_explicit->re_initialize( nv_glob, nv_loc ) ;
   F_Outer->re_initialize( nv_glob, nv_loc ) ;
   F_Inner->re_initialize( nv_glob, nv_loc ) ;
   F_Flowrate->re_initialize( nv_glob, nv_loc ) ;
   F_Stress->re_initialize( nv_glob, nv_loc ) ;
   
   G->re_initialize( np_glob, np_loc ) ;
   Mpl->re_initialize( np_glob, np_loc ) ;
   PRES->re_initialize( np_glob, np_loc ) ;

   L->re_initialize( np_glob, np_glob, np_loc, np_loc ) ;
   M->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
 
   // viscoplastic 
   T->re_initialize( nt_glob, nt_loc ) ;
   A_TtimesGAMMADOTv->re_initialize( nv_glob, nt_glob, nv_loc, nt_loc ) ;
   TminusRGAMMA->re_initialize( nt_glob, nt_loc ) ;
   F_Visco->re_initialize( nv_glob, nv_loc ) ;
}

//-----------------------------------------------------------------------
void 
ML_NavierStokesSystem:: print_errors_AL( std::string const& indent, size_t n,
                                         double err1, double err2,
                                         size_t nb_it )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( n==1 )
   {
      PEL::out() << indent 
           << setw( 3+15 ) << "[BU-G]/Mp" 
           << setw( 15 )   << "DU/U" ;
      if( nb_it != PEL::bad_index() ) PEL::out() << setw( 9 ) << "AM_its" ;
      PEL::out() << endl ;
   }
   PEL::out() << indent 
        << setw( 3 ) << n 
        << setprecision( 6 ) << setw( 15 ) << err1
        << setprecision( 6 ) << setw( 15 ) << err2 ;
   if( nb_it != PEL::bad_index() ) PEL::out() << setw( 9 ) << nb_it ;
   PEL::out() << endl ;

   PEL::out().flags( original_flags ) ;
}

//-----------------------------------------------------------------------
void 
ML_NavierStokesSystem:: print_end_errors_AL( 
                                      std::string const& indent, size_t n,
                                      double err1, double err2 )
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   PEL::out() << indent 
        << n << " iterations" 
        << "   [BU-G]/Mp = " << setprecision( 6 ) << err1
        << "   DU/U = "      << setprecision( 6 ) << err2 
        << endl ;

   PEL::out().flags( original_flags ) ;
}

//-----------------------------------------------------------------------
void
ML_NavierStokesSystem::  set_calculate_LHS( bool flag )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "ML_NavierStokesSystem:: set_calculate_LHS" ) ;
   
   LHS_flag = flag;
}   
//-----------------------------------------------------------------------