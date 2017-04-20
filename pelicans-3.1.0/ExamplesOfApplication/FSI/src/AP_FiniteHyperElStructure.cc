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

#include <AP_FiniteHyperElStructure.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <AP.hh>
#include <AP_ConstitutiveLaw.hh>
#include <AP_KinematicState.hh>
#include <AP_LoadCalculator.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_FiniteHyperElStructure const* 
AP_FiniteHyperElStructure::PROTOTYPE = new AP_FiniteHyperElStructure() ;

//---------------------------------------------------------------------------
AP_FiniteHyperElStructure:: AP_FiniteHyperElStructure( void )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( "AP_FiniteHyperElStructure" )
{
}

//---------------------------------------------------------------------------
AP_FiniteHyperElStructure*
AP_FiniteHyperElStructure:: create_replica( PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_FiniteHyperElStructure* result = 
                    new AP_FiniteHyperElStructure( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_FiniteHyperElStructure:: AP_FiniteHyperElStructure( 
                               PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( a_owner, dom, exp )
   , DISP( dom->set_of_discrete_fields()->item( 
                                   exp->string_data( "displacement" ) ) )
   , VELO( 0 )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , DYNAMIC( false )
   , SMALL_DEF( false )
   , ST( AP_KinematicState::create( this, dom->nb_space_dimensions() ) ) 
   , LAW( 0 )
   , RHO( PEL::bad_double() )
   , BCs( dom->set_of_boundary_conditions() )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , LC( 0 )
   , LOAD( 0 )
   , EXT_LOAD( 0 ) 
   , RHS( prms->item( exp->string_data( "param_source" ) ) )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , ITER( 1 )
   , MAX_ITER( exp->int_data( "max_nb_iterations" ) )
   , EPS_DISP( exp->double_data( "disp_tolerance" ) )
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: AP_FiniteHyperElStructure" ) ;

   cFE->require_field_calculation( DISP, PDE_LocalFE::N )  ;
   cFE->require_field_calculation( DISP, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( DISP, PDE_LocalFE::N )  ;
 
   PEL_ModuleExplorer* se = 0 ;
   if( exp->has_module( "dynamic" ) )
   {
      se = exp->create_subexplorer( 0, "dynamic" ) ;
      RHO = se->double_data( "density" ) ;
      VELO = dom->set_of_discrete_fields()->item( 
                                            se->string_data( "velocity" ) ) ;
      cFE->require_field_calculation( VELO, PDE_LocalFE::N )  ;
      DYNAMIC = true ;
      se->destroy() ; se = 0 ;
   }

   if( exp->has_entry( "small_deformations" ) )
   {
      SMALL_DEF = exp->bool_data( "small_deformations" ) ;
   }

   if( exp->has_entry( "param_external_load" ) ) //???? coherence avec CL
   {
      EXT_LOAD = prms->item( exp->string_data( "param_external_load" ) ) ; 
      check_param_nb_components( EXT_LOAD, "param_external_load", 
                                 dom->nb_space_dimensions() ) ;
      EXT_LOAD->transfer_bound_calculation_requirements( bFE, 
                                                       FE_Parameter::Val ) ;
   }

   check_param_nb_components( RHS, "param_source", 
                              dom->nb_space_dimensions() ) ;
   RHS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val );

   se = exp->create_subexplorer( 0, "AP_ConstitutiveLaw" ) ;
   LAW = AP_ConstitutiveLaw::create( this, se ) ;
   se->destroy() ;
   
   PDE_LinkDOF2Unknown* lnk = PDE_LinkDOF2Unknown::create( 0, 
                                          DISP, 
                                          "sequence_of_the_components", 
                                          true ) ;
   NMB = PDE_SystemNumbering::create( this, lnk ) ;

   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   if( exp->has_entry( "load_calculator" ) )
   {
      LC = AP_LoadCalculator::object( exp->string_data( "load_calculator" ) ) ;
      LOAD = LA_SeqVector::create( this, 0 ) ;
   }
}

//---------------------------------------------------------------------------
AP_FiniteHyperElStructure:: ~AP_FiniteHyperElStructure( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
size_t
AP_FiniteHyperElStructure:: nb_unknowns( void ) const
//---------------------------------------------------------------------------
{
   return( 1 ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField*
AP_FiniteHyperElStructure:: field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: field" ) ;
   PEL_CHECK_PRE( field_PRE( i_unk ) ) ;

   PDE_DiscreteField* result = DISP ;

   PEL_CHECK_POST( field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
AP_FiniteHyperElStructure:: level_of_field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: level_of_field" ) ;
   PEL_CHECK_PRE( level_of_field_PRE( i_unk ) ) ;

   size_t result = L_UPDATE ;

   PEL_CHECK_POST( level_of_field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
AP_FiniteHyperElStructure:: link_DOF_2_unknown( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( link_DOF_2_unknown_PRE( i_unk ) ) ;

   PDE_LinkDOF2Unknown const* result = NMB->link() ;

   PEL_CHECK_POST( link_DOF_2_unknown_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_FiniteHyperElStructure:: build_function_and_jacobian( 
                                           FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: build_function_and_jacobian" ) ;

   start_assembling_timer() ;
   // -------------------
   
   if( ITER == 1 )
   {
      NMB->reset() ;

      size_t nb_global_unk = NMB->nb_global_unknowns() ;
      size_t nb_local_unk  = NMB->nb_unknowns_on_current_process() ;

      A->re_initialize( nb_global_unk, nb_global_unk,
                        nb_local_unk, nb_local_unk ) ;
      F->re_initialize( nb_global_unk, nb_local_unk ) ;
      X->re_initialize( nb_global_unk, nb_local_unk ) ;

      NMB->define_scatters( X ) ;

      size_t nn = NMB->link()->unknown_vector_size() ;
      X_LOC->re_initialize( nn ) ;
      if( LOAD != 0 ) LOAD->re_initialize( nn ) ;
   }
   else
   {
      A->nullify() ;
      F->nullify() ; 
   }

   loop_on_cells( t_it ) ;
   loop_on_bounds( t_it ) ;
   
   if( LC != 0 )
   {
      LC->update_load( NMB->link(), LOAD ) ;
      for( size_t i=0 ; i<LOAD->nb_rows() ; ++i )
      {
//         std::cout << "LOAD->item( i )=" << LOAD->item( i ) << endl ;
         F->add_to_item( i, LOAD->item( i ) ) ;
//         GLOBAL_EQ->add_to_RHS_subvector_item( i, LOAD->item( i ) ) ;
      }
   }

   stop_assembling_timer() ;
   // ------------------
}

//---------------------------------------------------------------------------
LA_SeqVector const*
AP_FiniteHyperElStructure:: create_function( PEL_Object* a_owner,
                                             size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: create_function" ) ;
   PEL_CHECK_PRE( create_function_PRE( a_owner, i_unk ) ) ;
   
   F->synchronize() ; //??????? pas fait dans CH_CahnHilliard ???
   size_t nn = NMB->link( i_unk )->unknown_vector_size() ;
   LA_SeqVector* result = LA_SeqVector::create( a_owner, nn ) ;
   NMB->scatter( i_unk )->get( F, result ) ;

   PEL_CHECK_POST( create_function_POST( result, a_owner, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SeqMatrix const*
AP_FiniteHyperElStructure:: create_jacobian( PEL_Object* a_owner,
                                             size_t i_eq, 
                                             size_t j_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: create_jacobian" ) ;
   PEL_CHECK_PRE( create_jacobian_PRE( a_owner, i_eq, j_unk ) ) ;

//   LA_SeqMatrix const* result = GLOBAL_EQ->create_block_matrix( a_owner, i_eq, j_unk, 0 ) ;
   LA_SeqMatrix const* result =
                PDE::create_extracted_block( a_owner, A, NMB, i_eq, j_unk ) ;

   PEL_CHECK_POST( create_jacobian_POST( result, a_owner, i_eq, j_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
AP_FiniteHyperElStructure:: do_one_inner_iteration( 
                                                 FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_FiniteHyperElStructure:: do_one_inner_iteration" ) ;
   // --------------

   ITER = 0 ;
   size_t nint = PEL::bad_index() ;
   bool conv = false ;
   do
   {
      ITER++ ;
      if( ITER >= MAX_ITER )
         PEL_Error::object()->raise_plain( "Newton convergence failure" ) ;

      build_function_and_jacobian( t_it ) ;

      start_solving_timer() ;
      // ----------------
      
      A->synchronize() ;
      F->synchronize() ;

      //??? cf CH_CahnHillard
      SOLVER->set_matrix( A ) ;
      SOLVER->solve( F, X ) ;
      SOLVER->unset_matrix() ;
      
      if( SOLVER->is_iterative() ) nint = SOLVER->nb_iterations_achieved() ;

      stop_solving_timer() ;
      // ---------------

      LA_Scatter const* sca = NMB->scatter() ;
      PDE_LinkDOF2Unknown const* lnk = NMB->link() ;
      sca->get( X, X_LOC ) ;
      DISP->add_to_free_DOFs_value( L_UPDATE, X_LOC, lnk, 1.0 ) ;

      double delta_DISP = X->max_norm() ;
      conv = ( delta_DISP < EPS_DISP ) ; //??????? à ajuster ????

      if( verbose_level() >= 1 )
         print_errors( ITER, delta_DISP, F->two_norm(), nint ) ;

   } while( !conv ) ;

   if( DYNAMIC )
   {
      size_t nbc = VELO->nb_components() ;
      for( size_t n=0 ; n<VELO->nb_nodes() ; ++n )
      {
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            double xx = 
               2.0/t_it->time_step() * ( DISP->DOF_value( L_UPDATE, n, ic ) - 
                                         DISP->DOF_value( L_EXPLICIT, n, ic ) )
               - VELO->DOF_value( L_EXPLICIT, n, ic ) ;
            VELO->set_DOF_value( L_UPDATE, n, xx, ic ) ;
         }
      }
   }

   stop_total_timer() ;
   // -------------
}

//------------------------------------------------------------------------
void
AP_FiniteHyperElStructure:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: loop_on_cells" ) ;

   size_t nbc = DISP->nb_components() ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector coef( nbc ) ;
   double dt = t_it->time_step() ;

   doubleArray2D S( 3, 3 ) ;
   doubleArray4D C( 3, 3, 3, 3 ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( DISP, DISP ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         if( DYNAMIC )
         {
            // pourrait ne pas etre assemblé à chaque pas de Newton
            double rddt = 2.0 * RHO /dt/dt ;
            FE::add_row_col_S( ELEMENT_EQ, cFE, rddt ) ;

            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               coef( d ) = RHS->cell_value_at_IP( t_it, cFE, d )
                         - rddt * ( cFE->value_at_IP( DISP, L_UPDATE, d ) -
                                    cFE->value_at_IP( DISP, L_EXPLICIT, d ) )
                     + 2.0*RHO/dt * cFE->value_at_IP( VELO, L_EXPLICIT, d ) ;
            }
            FE::add_row( ELEMENT_EQ, cFE, coef ) ;
         }

         ST->set_state( DISP, L_UPDATE, cFE ) ;
         LAW->update_S_dSdE( ST, S, C ) ;
   
         double ss = 1.0 ;
         if( DYNAMIC ) ss = 0.5 ;
         if( SMALL_DEF )
         {
            AP::add_C_gradsym_row_gradsym_col( ELEMENT_EQ, cFE, C, ss ) ;
            AP::add_S_gradsym_row( ELEMENT_EQ, cFE, S, -ss ) ;
         }
         else
         {
            AP::add_C_gradGL_row_gradGL_col( ELEMENT_EQ, cFE, 
                                             ST->grad_disp(), C, ss ) ;
            AP::add_S_graddGL_row_col( ELEMENT_EQ, cFE, S, ss ) ;
            AP::add_S_gradGL_row( ELEMENT_EQ, cFE, 
                                  ST->grad_disp(), S, -ss ) ;
         }

         if( DYNAMIC )
         {
            ST->set_state( DISP, L_EXPLICIT, cFE ) ;
            LAW->update_S_dSdE( ST, S, C ) ;
            if( SMALL_DEF )
            {
               AP::add_S_gradsym_row( ELEMENT_EQ, cFE, S, -0.5 ) ;
            }
            else
            {
               AP::add_S_gradGL_row( ELEMENT_EQ, cFE,
                                      ST->grad_disp(), S, -0.5 ) ;
            }
         }
      }
      PDE::assemble_in_vector_1( F, ELEMENT_EQ, NMB ) ;
      PDE::assemble_in_matrix_0( A, ELEMENT_EQ, NMB ) ;
   }
}

//---------------------------------------------------------------------------
void
AP_FiniteHyperElStructure:: loop_on_bounds( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FiniteHyperElStructure:: loop_on_bounds" ) ;

   size_t nbc = DISP->nb_components() ;
   size_t nb_dims = bFE->nb_space_dimensions() ;

   boolVector done( DISP->nb_nodes() ) ;
   doubleVector coef( nb_dims ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      bFE->set_row_and_col_fields( DISP, DISP ) ;

      ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), nbc,
                              bFE->col_field_node_connectivity(), nbc ) ;

      GE_Color const* color = bFE->color() ;
      if( BCs->has_BC( color, DISP ) )
      {
         PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, DISP ) ;
         string const& type = ee->string_data( "type" ) ;
         if( type == "Neumann" )
         {
            bFE->start_IP_iterator( QRP ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {  
               for( size_t i=0 ; i<bFE->nb_basis_functions( row ) ; ++i )
               {
                  for( size_t ic=0 ; ic<nb_dims ; ++ic )
                  {     
                     double tr = EXT_LOAD->bound_value_at_IP( t_it, bFE, ic ) ;
                     double xx = bFE->weight_of_IP() * tr
                               * bFE->N_at_IP( row, i ) ;
                     ELEMENT_EQ->add_to_vector( xx, i, ic ) ;
                  }
               }
            }   
            PDE::assemble_in_vector_1( F, ELEMENT_EQ, NMB ) ; //??? nothing in A
//            GLOBAL_EQ->assemble_in_RHS( ELEMENT_EQ ) ; 
         }
      }
   }
}

//---------------------------------------------------------------------------
void
AP_FiniteHyperElStructure:: print_errors( size_t iter, 
                                          double delta_DISP, 
                                          double L2_norm_RHS, 
                                          size_t nint ) const
//---------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( iter==1 )
   {
      PEL::out() << indent()
                 << setw( 3+15 ) << "max(d_disp)" 
                 << setw( 15 )   << "L2(rhs)" ;
      if( nint != PEL::bad_index() ) PEL::out() << setw( 9 ) << "A_its" ;
      PEL::out() << endl ;
   }   
   PEL::out() << indent() 
              << setw( 3 ) << iter
              << setprecision( 6 ) << setw( 15 ) << delta_DISP
              << setprecision( 6 ) << setw( 15 ) << L2_norm_RHS ;
   if( nint != PEL::bad_index() )  PEL::out() << setw( 9 ) << nint ;
   PEL::out() << endl ;

   PEL::out().flags( original_flags ) ;
}
