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

#include <AP_QIncompHyperelasticity.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
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

AP_QIncompHyperelasticity const* 
AP_QIncompHyperelasticity::PROTOTYPE = new AP_QIncompHyperelasticity() ;

//---------------------------------------------------------------------------
AP_QIncompHyperelasticity:: AP_QIncompHyperelasticity( void )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( "AP_QIncompHyperelasticity" )
   , LAMBDA( 0 )
{
}

//---------------------------------------------------------------------------
AP_QIncompHyperelasticity*
AP_QIncompHyperelasticity:: create_replica( PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_QIncompHyperelasticity* result =
                  new AP_QIncompHyperelasticity( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_QIncompHyperelasticity:: AP_QIncompHyperelasticity( 
                                      PEL_Object* a_owner,
                                      PDE_DomainAndFields const* dom,
                                      FE_SetOfParameters const* prms,
                                      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIterationOpen( a_owner, dom, exp )
   , DISP( dom->set_of_discrete_fields()->item( 
                                   exp->string_data( "displacement" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , ST( AP_KinematicState::create( this, dom->nb_space_dimensions() ) ) 
   , LAW( 0 )
   , BCs( dom->set_of_boundary_conditions() )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , LEQ_2( PDE_LocalEquation::create( this ) )
   , LEQ_3( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , LC( 0 )
   , LOAD( 0 )
   , EXT_LOAD( 0 )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , PENALTY( exp->double_data( "penalty_parameter" ) )
   , ITER( 1 )
   , MAX_ITER( exp->int_data( "max_nb_iterations" ) )
   , EPS_DISP( exp->double_data( "disp_tolerance" ) )
   , TOL_AL( exp->double_data( "multiplier_tolerance" ) )
   , LAMBDA( 0 ) 
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: AP_QIncompHyperelasticity" ) ;

   cFE->require_field_calculation( DISP, PDE_LocalFE::N )  ;
   cFE->require_field_calculation( DISP, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( DISP, PDE_LocalFE::N )  ;
 
   if( exp->has_entry( "param_external_load" ) )
   {
      EXT_LOAD = prms->item( exp->string_data( "param_external_load" ) ) ; 
      check_param_nb_components( EXT_LOAD, "param_external_load", 
                                 dom->nb_space_dimensions() ) ;
      EXT_LOAD->transfer_bound_calculation_requirements( bFE, 
                                                         FE_Parameter::Val ) ;
   }

   PEL_ModuleExplorer const* ee =
                      exp->create_subexplorer( 0, "AP_ConstitutiveLaw" ) ;
   LAW = AP_ConstitutiveLaw::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   PDE_LinkDOF2Unknown* lnk = PDE_LinkDOF2Unknown::create( 0, 
                                          DISP, 
                                          "sequence_of_the_components", 
                                          true ) ;
   NMB = PDE_SystemNumbering::create( this, lnk ) ;

   ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
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
   
   LAMBDA.re_initialize( cFE->nb_meshes() ) ;
   if( exp->has_entry( "multiplier" ) )
   {
      LAMBDA.set( exp->double_data( "multiplier" ) ) ;
   }
}

//---------------------------------------------------------------------------
AP_QIncompHyperelasticity:: ~AP_QIncompHyperelasticity( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
size_t
AP_QIncompHyperelasticity:: nb_unknowns( void ) const
//---------------------------------------------------------------------------
{
   return( 1 ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField*
AP_QIncompHyperelasticity:: field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: field" ) ;
   PEL_CHECK_PRE( field_PRE( i_unk ) ) ;

   PDE_DiscreteField* result = DISP ;

   PEL_CHECK_POST( field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
AP_QIncompHyperelasticity:: level_of_field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: level_of_field" ) ;
   PEL_CHECK_PRE( level_of_field_PRE( i_unk ) ) ;

   size_t result = L_UPDATE ;

   PEL_CHECK_POST( level_of_field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
AP_QIncompHyperelasticity:: link_DOF_2_unknown( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( link_DOF_2_unknown_PRE( i_unk ) ) ;

   PDE_LinkDOF2Unknown const* result = NMB->link() ;

   PEL_CHECK_POST( link_DOF_2_unknown_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_QIncompHyperelasticity:: build_function_and_jacobian( 
                                           FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: build_function_and_jacobian" ) ;

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
         F->add_to_item( i, LOAD->item( i ) ) ;
//         GLOBAL_EQ->add_to_RHS_subvector_item( i, LOAD->item( i ) ) ;
      }
   }

   stop_assembling_timer() ;
   // ------------------
}

//---------------------------------------------------------------------------
LA_SeqVector const*
AP_QIncompHyperelasticity:: create_function( PEL_Object* a_owner,
                                             size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: create_function" ) ;
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
AP_QIncompHyperelasticity:: create_jacobian( PEL_Object* a_owner,
                                             size_t i_eq, 
                                             size_t j_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: create_jacobian" ) ;
   PEL_CHECK_PRE( create_jacobian_PRE( a_owner, i_eq, j_unk ) ) ;

//   LA_SeqMatrix const* result = GLOBAL_EQ->create_block_matrix( a_owner, i_eq, j_unk, 0 ) ;
   LA_SeqMatrix const* result =
                PDE::create_extracted_block( a_owner, A, NMB, i_eq, j_unk ) ;

   PEL_CHECK_POST( create_jacobian_POST( result, a_owner, i_eq, j_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
AP_QIncompHyperelasticity:: do_one_inner_iteration( 
                                                FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_QIncompHyperelasticity:: do_one_inner_iteration" ) ;
   // --------------

   size_t iterAL = 0 ;
   bool convAL = false ;
   LAMBDA.set( 0.0 ) ;

   do
   {
      iterAL++ ;
      
      bool convNewton = false ;
      ITER = 0 ;
      size_t nint = PEL::bad_index() ;
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
         convNewton = ( delta_DISP < EPS_DISP ) ; //??????? à ajuster ????

         PEL::out() << "Newton it : " << ITER << "  " << delta_DISP << "  "
                    << F->two_norm()  << endl ;

      } while( !convNewton ) ;

      size_t e = 0 ;
      double h_max = 0.0 ;
      double vol = 0.0 ;
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         double theta_e = cell_dilatation() ;
         double h_e = hh( theta_e ) ;
         LAMBDA( e ) += PENALTY * h_e ;

         if( PEL::abs( h_e ) > h_max ) h_max = PEL::abs( h_e ) ;

         vol += theta_e * cFE->polyhedron()->measure() ;
         ++e ;
      }
      convAL = ( h_max < TOL_AL ) ;

      PEL::out() << " AL it : " << iterAL << endl ;
      PEL::out() << "       h_max = " << h_max << "  vol = " << vol << endl ;
      PEL::out() << "-------------- " << endl ;

   } while( !convAL ) ;

   PEL::out() << " convergence of AL and Newton achieved. " << endl ;

   stop_total_timer() ;
   // ---------------------
}

//------------------------------------------------------------------------
void
AP_QIncompHyperelasticity:: loop_on_cells( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: loop_on_cells" ) ;

   size_t nbc = DISP->nb_components() ;
   size_t nb_dims = cFE->nb_space_dimensions() ;

   doubleArray2D S( 3, 3 ) ;
   doubleArray2D SS( 3, 3 ) ;
   doubleArray4D C( 3, 3, 3, 3 ) ;

   size_t e = 0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( DISP, DISP ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;

      LEQ_2->initialize( cFE->row_field_node_connectivity(), nbc,
                         cFE->col_field_node_connectivity(), nbc ) ;

      LEQ_3->initialize( cFE->row_field_node_connectivity(), nbc,
                         cFE->col_field_node_connectivity(), nbc ) ;

      //??? boucle sur les IPs que l'on doit pouvoir éviter car 
      //??? theta_e, p_e, dp_e sont constants sur la maille

      double theta_e = cell_dilatation() ;

      double p_e = LAMBDA( e ) * d_hh( theta_e ) + 
                   PENALTY * d_gamma( theta_e ) ;
      
      double dp_e = LAMBDA( e ) * d2_hh( theta_e ) + 
                    PENALTY * d2_gamma( theta_e ) ;
                         
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         ST->set_state( DISP, L_UPDATE, cFE ) ;
         for( size_t a=0 ; a<3 ; ++a )
            for( size_t b=0 ; b<3 ; ++b )
               SS( a, b ) = ST->detF() * ST->inv_rCG()(a,b) ;

         LAW->update_S_dSdE( ST, S, C ) ;
         add_volumetric_part( ST, nb_dims, p_e, S, C ) ;

         AP::add_C_gradGL_row_gradGL_col( ELEMENT_EQ, cFE, 
                                          ST->grad_disp(), C, 1.0 ) ;
         AP::add_S_graddGL_row_col( ELEMENT_EQ, cFE, S, 1.0 ) ;
         AP::add_S_gradGL_row( ELEMENT_EQ, cFE, ST->grad_disp(), S, -1.0 ) ;

         AP::add_S_gradGL_row( LEQ_2, cFE, ST->grad_disp(), SS,
                               1.0/cFE->polyhedron()->measure() ) ;
         AP::add_S_gradGL_row( LEQ_3, cFE, ST->grad_disp(), SS, 1.0 ) ;
      }

      size_t nb_nodes = cFE->nb_basis_functions( row ) ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            for( size_t j=0 ; j<nb_nodes ; ++j )
            {
               for( size_t jc=0 ; jc<nbc ; ++jc )
               {
                  double xx = LEQ_2->vector_item( j, jc ) *
                              LEQ_3->vector_item( i, ic ) ;
                  ELEMENT_EQ->add_to_matrix( xx*dp_e, i, j, ic, jc ) ;
               }
            }
         }
      }

      PDE::assemble_in_vector_1( F, ELEMENT_EQ, NMB ) ;
      PDE::assemble_in_matrix_0( A, ELEMENT_EQ, NMB ) ;
//      GLOBAL_EQ->assemble_in_RHS( ELEMENT_EQ ) ;
//      GLOBAL_EQ->assemble_in_LHS( ELEMENT_EQ ) ;

      ++e ;
   }
}

//---------------------------------------------------------------------------
void
AP_QIncompHyperelasticity:: loop_on_bounds( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: loop_on_bounds" ) ;

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
                     double tr = EXT_LOAD->bound_value_at_IP( t_it, bFE, ic )  ;
                     double xx = bFE->weight_of_IP() * tr
                                                     * bFE->N_at_IP( row, i ) ;
                    
                     ELEMENT_EQ->add_to_vector( xx, i, ic ) ;
                  }
               }
            }  
            PDE::assemble_in_vector_1( F, ELEMENT_EQ, NMB ) ; //??? nothing in A
//            GLOBAL_EQ->assemble_in_RHS( ELEMENT_EQ ) ; //???? rien dans A
         }
      }
   }
}

//---------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: cell_dilatation( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: cell_dilatation" ) ;

   double result = 0.0 ;
   cFE->start_IP_iterator( QRP ) ;
   for( ; cFE->valid_IP() ; cFE->go_next_IP() )
   {
      ST->set_state( DISP, L_UPDATE, cFE ) ;
      result += ST->detF() * cFE->weight_of_IP() ;
   }

   result /= cFE->polyhedron()->measure() ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_QIncompHyperelasticity:: add_volumetric_part( AP_KinematicState const* st,
                                                 size_t nb_dims, 
                                                 double pp,
                                                 doubleArray2D& Pio2,
                                                 doubleArray4D& dPio2dGL )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_QIncompHyperelasticity:: add_volumetric_part" ) ;

   double jac = st->detF() ;
   doubleArray2D const& inv_rCG = st->inv_rCG() ;

   for( size_t a=0 ; a<nb_dims ; ++a )
   {
      for( size_t b=0 ; b<nb_dims ; ++b )
      {
         double xx = jac * pp * inv_rCG( a, b ) ;

         for( size_t k=0 ; k<nb_dims ; ++k )
         {
            for( size_t l=0 ; l<nb_dims ; ++l )
            {
               double yy =  pp ;
               yy *= inv_rCG( a, b ) * inv_rCG( k, l ) ;
               yy -= pp * ( inv_rCG( a, k ) * inv_rCG( b, l ) 
                          + inv_rCG( a, l ) * inv_rCG( b, k ) ) ; 
               yy *= jac ; 

               dPio2dGL( a, b, k, l ) += yy ;
            }
         }
         Pio2( a, b ) += xx ;
      }
   }
}

//----------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: d_gamma( double x ) const
//----------------------------------------------------------------------------
{
   double result = x - 1.0/x ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: d2_gamma( double x ) const
//----------------------------------------------------------------------------
{
   double result = 1.0 + 1.0/x/x ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: hh( double x ) const
//----------------------------------------------------------------------------
{
   double result = x - 1.0/x ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: d_hh( double x ) const
//----------------------------------------------------------------------------
{
   double result = 1.0 + 1.0/x/x ;
   return( result ) ;
}

//----------------------------------------------------------------------------
double
AP_QIncompHyperelasticity:: d2_hh( double x ) const
//----------------------------------------------------------------------------
{
   double result = -2.0/x/x/x ;
   return( result ) ;
}
