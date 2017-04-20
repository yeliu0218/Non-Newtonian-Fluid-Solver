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

#include <CH_CahnHilliard.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_MemoryTracer.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_QRprovider.hh>
#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <CH_BulkChemicalPotential.hh>

#include <fstream>
#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ; using std::cout ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

struct CH_CahnHilliard_ERROR
{
   static void n0( std::string const& fname ) ;
} ;

CH_CahnHilliard const*
CH_CahnHilliard::PROTOTYPE = new CH_CahnHilliard() ;

//----------------------------------------------------------------------
CH_CahnHilliard:: CH_CahnHilliard( void )
//----------------------------------------------------------------------
   : FE_OneStepIterationOpen( "CH_CahnHilliard" )
   , PRINTED_RES( 0 )
{
}

//----------------------------------------------------------------------
CH_CahnHilliard*
CH_CahnHilliard:: create_replica( PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   CH_CahnHilliard* result =
                      new CH_CahnHilliard( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
CH_CahnHilliard:: CH_CahnHilliard( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIterationOpen( a_owner, dom, exp )
   , C1( dom->set_of_discrete_fields()->item(
                     exp->string_data( "phase_field_1" ) ) )
   , C1_EXP( dom->set_of_discrete_fields()->item(
                     exp->string_data( "phase_field_1_explicit" ) ) )
   , M1( dom->set_of_discrete_fields()->item(
                     exp->string_data( "generalized_potential_1" ) ) )
   , C2( dom->set_of_discrete_fields()->item(
                     exp->string_data( "phase_field_2" ) ) )
   , C2_EXP( dom->set_of_discrete_fields()->item(
                     exp->string_data( "phase_field_2_explicit" ) ) )
   , M2( dom->set_of_discrete_fields()->item(
                     exp->string_data( "generalized_potential_2" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , L_EXPLICIT( exp->int_data( "level_of_explicit" ) )
   , AA( 0 )
   , L_AA( PEL::bad_index() )
   , EPS( exp->double_data( "thickness" ) )
   , CAP1( 0.0 )
   , CAP2( 0.0 )
   , THETA( 1.0 )
   , NB_PRE_EULER_STEP( 1 )//valeur par defaut à 0
   , MOB_CST( exp->double_data( "mobility_cst" ) )
   , MOB_DEG( exp->double_data( "mobility_deg" ) )
   , MOB_EXP( exp->bool_data( "explicit_mobility" ) )
   , MOB_MOD( false )
   , MOB_MAX( 0.1 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object(
                      exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , VOL1( PEL::bad_double() )
   , VOL2( PEL::bad_double() )
   , ITER( 1 )
   , NB_ITER_MAX( exp->int_data( "nb_iterations_max" ) )
   , TOL( exp->double_data( "newton_tolerance" ) )
   , idx_C1( 0 )
   , idx_M1( 1 )
   , idx_C2( 2 )
   , idx_M2( 3 )
   , LHS_cst( 0 )
   , LHS( 0 )
   , RHS( 0 )
   , UNK( 0 )
   , UNK_LOC( LA_SeqVector::create( this, 0 ) )
   , SOLVER( 0 )
   , PRINTED_RES( 4 )
{
   cFE->require_field_calculation( C1, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( C1, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( C1_EXP, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( C1_EXP, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( C2, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( C2, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( C2_EXP, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( C2_EXP, PDE_LocalFE::dN ) ;

   cFE->require_field_calculation( M1, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( M1, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( M2, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( M2, PDE_LocalFE::dN ) ;

   if( exp->has_module( "advection" ) )
   {
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "advection" ) ;
      AA = dom->set_of_discrete_fields()->item( e->string_data( "field" ) ) ;
      check_field_nb_components( AA, dom->nb_space_dimensions() ) ;
      L_AA = e->int_data( "level" ) ;
      e->destroy() ;
      cFE->require_field_calculation( AA, PDE_LocalFE::N ) ;
   }

   PEL_ModuleExplorer const* ee =
                      exp->create_subexplorer( 0, "CH_BulkChemicalPotential" ) ;
   BULK_MU = CH_BulkChemicalPotential::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   CAP1 = 3.*EPS*BULK_MU->Sigma1()/4. ;
   CAP2 = 3.*EPS*BULK_MU->Sigma2()/4. ;

   if( exp->has_entry( "theta_coef" ) )
   {
      THETA = exp->double_data( "theta_coef") ;
   }
   if( exp->has_entry( "nb_pre_euler_step" ) )
   {
      NB_PRE_EULER_STEP = exp->int_data( "nb_pre_euler_step") ;
   }
   if( exp->has_entry( "mobility_modification" ) )
   {
      MOB_MOD = exp->bool_data( "mobility_modification") ;
      if( exp->has_entry( "mobility_max" ) )
         MOB_MAX =  exp->double_data( "mobility_max") ;
      M1_EXP = dom->set_of_discrete_fields()->item(
                exp->string_data( "generalized_potential_1_explicit" ) ) ;
      M2_EXP = dom->set_of_discrete_fields()->item(
                     exp->string_data( "generalized_potential_2_explicit" ) ) ;
      cFE->require_field_calculation( M1_EXP, PDE_LocalFE::dN ) ;
      cFE->require_field_calculation( M2_EXP, PDE_LocalFE::dN ) ;
   }

   // -- Global Equation
   PEL_Vector* vec = PEL_Vector::create( 0, 4 ) ;
   vec->set_at( idx_C1, PDE_LinkDOF2Unknown::create( 0, C1, true ) ) ;
   vec->set_at( idx_M1, PDE_LinkDOF2Unknown::create( 0, M1, true ) ) ;
   vec->set_at( idx_C2, PDE_LinkDOF2Unknown::create( 0, C2, true ) ) ;
   vec->set_at( idx_M2, PDE_LinkDOF2Unknown::create( 0, M2, true ) ) ;

   //std::string ordering = "sequence_of_the_unknowns" ;
   std::string ordering = "sequence_of_the_discrete_fields" ;
   NMB = PDE_SystemNumbering::create( this, vec, ordering ) ;

   ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   LHS = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   LHS_cst = LHS->create_matrix( this ) ;
   RHS = LHS->create_vector( this ) ;
   UNK = LHS->create_vector( this ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   std::string nn = "multilevel_preconditioner_name" ;
   configure_multilevel_preconditioner( exp, nn, dom, NMB ) ;
   vec->destroy() ; vec = 0 ;
}

//---------------------------------------------------------------------------
CH_CahnHilliard:: ~CH_CahnHilliard( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
size_t
CH_CahnHilliard:: nb_unknowns( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: nb_unknowns" ) ;

   size_t result = 4 ;

   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField*
CH_CahnHilliard:: field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: field" ) ;
   PEL_CHECK_PRE( field_PRE( i_unk ) ) ;

   PDE_DiscreteField* result = 0 ;
   if( i_unk == idx_C1 ) result = C1 ;
   if( i_unk == idx_M1 ) result = M1 ;
   if( i_unk == idx_C2 ) result = C2 ;
   if( i_unk == idx_M2 ) result = M2 ;

   PEL_CHECK_POST( field_POST( result , i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
CH_CahnHilliard:: level_of_field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: level_of_field" ) ;
   PEL_CHECK_PRE( level_of_field_PRE( i_unk ) ) ;

   size_t result = L_UPDATE ;

   PEL_CHECK_POST( level_of_field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
CH_CahnHilliard:: link_DOF_2_unknown( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( link_DOF_2_unknown_PRE( i_unk ) ) ;

   PDE_LinkDOF2Unknown const* result = NMB->link( i_unk ) ;

   PEL_CHECK_POST( link_DOF_2_unknown_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
CH_CahnHilliard:: build_function_and_jacobian( FE_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: build_function_and_jacobian" ) ;

   start_assembling_timer() ;

   if( ITER == 1 )
   {
      NMB->reset() ;

      size_t nb_global_unk = NMB->nb_global_unknowns() ;
      size_t nb_local_unk  = NMB->nb_unknowns_on_current_process() ;

      LHS->re_initialize( nb_global_unk, nb_global_unk,
                          nb_local_unk, nb_local_unk ) ;
      LHS_cst->re_initialize( nb_global_unk, nb_global_unk,
                              nb_local_unk, nb_local_unk ) ;
      RHS->re_initialize( nb_global_unk, nb_local_unk ) ;
      UNK->re_initialize( nb_global_unk, nb_local_unk ) ;

      NMB->define_scatters( UNK ) ;

      size_t nn = NMB->link( 0 )->unknown_vector_size() ;
      for( size_t i=1 ; i<NMB->nb_links() ; ++i )
      {
         size_t mm = NMB->link( i )->unknown_vector_size() ;
         if( mm > nn ) nn = mm ;
      }
      UNK_LOC->re_initialize( nn ) ;
   }
   else
   {
      LHS->nullify() ;
      RHS->nullify() ;
   }

   loop_on_cells( t_it ) ;

   LHS->synchronize() ;
   if( !LHS_cst->is_synchronized() ) LHS_cst->synchronize() ;
   LHS->add_Mat( LHS_cst ) ;

   RHS->synchronize() ;

   stop_assembling_timer() ;
}

//---------------------------------------------------------------------------
LA_SeqVector const*
CH_CahnHilliard:: create_function( PEL_Object* a_owner, size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: create_function" ) ;
   PEL_CHECK_PRE( create_function_PRE( a_owner, i_unk ) ) ;

   size_t nn = NMB->link( i_unk )->unknown_vector_size() ;
   LA_SeqVector* result = LA_SeqVector::create( a_owner, nn ) ;
   NMB->scatter( i_unk )->get( RHS, result ) ;

   PEL_CHECK_POST( create_function_POST( result, a_owner, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SeqMatrix const*
CH_CahnHilliard:: create_jacobian( PEL_Object* a_owner,
                                   size_t i_eq, size_t j_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: create_jacobian" ) ;
   PEL_CHECK_PRE( create_jacobian_PRE( a_owner, i_eq, j_unk ) ) ;

   LA_SeqMatrix const* result =
                PDE::create_extracted_block( a_owner, LHS, NMB, i_eq, j_unk ) ;

   PEL_CHECK_POST( create_jacobian_POST( result, a_owner, i_eq,  j_unk ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "CH_CahnHilliard:: do_one_inner_iteration" ) ;

   ITER = 0 ;
   VOL1 = 0.0 ;
   VOL2 = 0.0 ;
   do
   {
      ITER++ ;
      if( ITER >= NB_ITER_MAX )
         PEL_Error::object()->raise_plain( "Newton convergence failure" ) ;

      build_function_and_jacobian( t_it ) ;

      start_solving_timer() ;

      SOLVER->set_stop_on_error( false ) ;

      if( ( verbose_level() >= 1 ) && ( ITER == 1 ) )
      {
         PEL::out() << indent() << "   nb unknowns: "
                                << UNK->nb_rows() << endl ;
         PEL::out() << indent() << "   nb stored items: "
                                << LHS->nb_stored_items() << endl ;
         PEL::out() << indent() << "   set_matrix: " ;
         PEL_MemoryTracer::display_memory( PEL::out(),
                                           PEL_MemoryTracer::used_memory() ) ;
         PEL::out() << " -> " ;
      }
      //??? voir AP_NavierStokes2CFV
      SOLVER->set_matrix( LHS ) ;
      if( ( verbose_level() >= 1 ) && ( ITER == 1 ) )
      {
         PEL_MemoryTracer::display_memory( PEL::out(),
                                           PEL_MemoryTracer::used_memory() ) ;
         PEL::out() << std::endl ;
      }
      SOLVER->solve( RHS, UNK ) ;
      if( !SOLVER->solution_is_achieved() )
      {
         PEL::out() << endl << "writing rhs in vec.mtx" << endl ;
         RHS->write( "vec.mtx" ) ;
         PEL::out() << endl << "writing lhs in mat.mtx" << endl ;
         LHS->writeMM( "mat.mtx" ) ;
         PEL_Error::object()->raise_plain(
            "*** CH_CahnHilliard:\n"
            "    Linear system resolution failure" ) ;
      }

      //??? on devrait pouvoir conserver le profil, qui ne change
      //??? pas d'une iteration de Newton à l'autre
      SOLVER->unset_matrix() ;

      stop_solving_timer() ;

      for( size_t i=0 ; i<NMB->nb_links() ; ++i )
      {
         LA_Scatter const* sca = NMB->scatter( i ) ;
         PDE_LinkDOF2Unknown const* link = NMB->link( i ) ;
         PEL_ASSERT( field( i ) == link->field() ) ;
         sca->get( UNK, UNK_LOC ) ;
         field( i )->add_to_free_DOFs_value( L_UPDATE, UNK_LOC, link, 1.0 ) ;

         if( verbose_level() >= 1 )
         {
            PRINTED_RES( i ) = UNK_LOC->max_norm() ;
         }
      }

   } while( !convergence_achieved() ) ;

   // memory release
   LHS->re_initialize( 0, 0, 0, 0 ) ;
   LHS_cst->re_initialize( 0, 0, 0, 0 ) ;
   RHS->re_initialize( 0, 0 ) ;
   UNK->re_initialize( 0, 0 ) ;

   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: save_other_than_time_and_fields(
                                                 FE_TimeIterator const* t_it,
                                                 PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Communicator const* com = communicator() ;

   if( VOL1 == PEL::bad_double() )
   {
      PEL_ASSERT( VOL2 == PEL::bad_double() ) ;

      if( com->rank() == 0 )
      {
         std::ofstream ofs( "volumes.txt", std::ios::out | std::ios::trunc ) ;
         if( !ofs ) CH_CahnHilliard_ERROR::n0( "volumes.txt" ) ;
         ofs << "#" << endl ;
         ofs << "# CH_CahnHilliard generated file" << endl ;
         ofs << "# meaning of the columns:" << endl ;
         ofs << "#    TIME : current time" << endl ;
         ofs << "#    VOL1 : explicit volume of phase 1" << endl ;
         ofs << "#    VOL2 : explicit volume of phase 2" << endl ;
         ofs << "#" << endl ;
         ofs << "#       TIME         VOL1         VOL2" << endl ;
         ofs.close() ;
      }
   }
   else
   {
      VOL1 = com->sum( VOL1 ) ;
      VOL2 = com->sum( VOL2 ) ;

      if( com->rank() == 0 )
      {
         std::ofstream ofs( "volumes.txt", std::ios::out | std::ios::app ) ;
         ofs.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
         ofs << std::setprecision( 5 ) ;

         ofs << setw( 12 ) << t_it->time() << " " ;
         ofs << setw( 12 ) << VOL1 << " " ;
         ofs << setw( 12 ) << VOL2 << " " ;
         ofs << endl ;

         ofs.close() ;
      }
   }
}

//----------------------------------------------------------------------
bool
CH_CahnHilliard:: convergence_achieved( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: convergence_achieved" ) ;

   double rhs = RHS->two_norm() ;
   bool result =  ( rhs < TOL )  ;

   if( verbose_level() >=1 )
   {
      ios_base::fmtflags original_flags = PEL::out().flags() ;
      PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;
      bool iterative = SOLVER->is_iterative() ;

      if( ITER==1 )
      {
         PEL::out() << indent()
                    << setw( 3+12 ) << "max(d_C1)"
                    << setw( 12 )   << "max(d_M1)"
                    << setw( 12 )   << "max(d_C2)"
                    << setw( 12 )   << "max(d_M2)"
                    << setw( 12 )   << "L2(rhs)" ;
         if( iterative )
            PEL::out() << setw( 6 ) << "its" ;
         PEL::out() << endl ;
      }

      PEL::out() << indent()
                 << setw( 3 ) << ITER
                 << setprecision( 3 ) << setw( 12 )
                 << PRINTED_RES( idx_C1 )
                 << setprecision( 3 ) << setw( 12 )
                 << PRINTED_RES( idx_M1 )
                 << setprecision( 3 ) << setw( 12 )
                 << PRINTED_RES( idx_C2 )
                 << setprecision( 3 ) << setw( 12 )
                 << PRINTED_RES( idx_M2 )
                 << setprecision( 3 ) << setw( 12 )
                 << rhs ;
      if( iterative )
      {
         PEL::out() << setw( 6 ) << SOLVER->nb_iterations_achieved() ;
      }
      PEL::out() << endl ;
      PEL::out().flags( original_flags ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_on_cells( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_on_cells" ) ;
   double dt = t_it->time_step()  ;

   double theta = ( t_it->iteration_number() < NB_PRE_EULER_STEP+1 ? 1.0
                                                                   : THETA ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( M1, C1 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Mi_Ci( C1, C1_EXP, M1, VOL1, dt ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_M1 ) ;
      if( ITER == 1 )
         PDE::assemble_in_matrix_0( LHS_cst, ELEMENT_EQ, NMB, idx_M1, idx_C1 ) ;

      cFE->set_row_and_col_fields( M1, M1 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         if( !MOB_MOD )
         {
            loop_Mi_Mi( M1, BULK_MU->Sigma1() ) ;
         }
         else
         {
            loop_Mi_Mi_MobMod( M1, M1_EXP, BULK_MU->Sigma1() ) ;
         }
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_M1 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_M1, idx_M1 ) ;

      cFE->set_row_and_col_fields( C1, M1 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Mi( M1 ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C1 ) ;
      if( ITER == 1 )
         PDE::assemble_in_matrix_0( LHS_cst, ELEMENT_EQ, NMB, idx_C1, idx_M1 ) ;

      cFE->set_row_and_col_fields( C1, C1 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Ci( 1, CAP1, theta ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C1 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_C1, idx_C1 ) ;

      cFE->set_row_and_col_fields( M2, C2 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Mi_Ci( C2, C2_EXP, M2, VOL2, dt ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_M2 ) ;
      if( ITER == 1 )
         PDE::assemble_in_matrix_0( LHS_cst, ELEMENT_EQ, NMB, idx_M2, idx_C2 ) ;

      cFE->set_row_and_col_fields( M2, M2 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         if( !MOB_MOD )
         {
            loop_Mi_Mi( M2, BULK_MU->Sigma2() ) ;
         }
         else
         {
            loop_Mi_Mi_MobMod( M2, M2_EXP, BULK_MU->Sigma2() ) ;
         }
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_M2 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_M2, idx_M2 ) ;

      cFE->set_row_and_col_fields( C2, M2 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Mi( M2 ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C2 ) ;
      if( ITER == 1 )
         PDE::assemble_in_matrix_0( LHS_cst, ELEMENT_EQ, NMB, idx_C2, idx_M2 ) ;

      cFE->set_row_and_col_fields( C2, C2 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Ci( 2, CAP2, theta ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C2 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_C2, idx_C2 ) ;

      cFE->set_row_and_col_fields( C1, C2 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Cj( 1, 2 ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C1 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_C1, idx_C2 ) ;

      cFE->set_row_and_col_fields( C2, C1 ) ;
//    ------------------------------------
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         loop_Ci_Cj( 2, 1 ) ;
      }
      PDE::assemble_in_vector_1( RHS, ELEMENT_EQ, NMB, idx_C2 ) ;
      PDE::assemble_in_matrix_0( LHS, ELEMENT_EQ, NMB, idx_C2, idx_C1 ) ;
   }
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Mi_Ci( PDE_DiscreteField const* ci,
                              PDE_DiscreteField const* ci_exp,
                              PDE_DiscreteField const* mi,
                              double& vol,
                              double dt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Mi_Ci" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector vv( nb_dims ) ;
   doubleVector cv( nb_dims ) ;
   doubleVector aa( nb_dims ) ;

   double val_exp = cFE->value_at_IP( ci_exp, L_EXPLICIT ) ;
   double xx = ( cFE->value_at_IP( ci, L_UPDATE ) - val_exp ) / dt ;
   if( AA != 0 )
   {
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         vv( d ) = cFE->value_at_IP( AA, L_AA, d ) ;
         xx += vv( d ) * cFE->gradient_at_IP( ci, L_UPDATE, d ) ;
      }
   }
   FE::add_row( ELEMENT_EQ, cFE, -xx ) ;

   if( ITER==1 )
   {
      FE::add_row_col_NS( ELEMENT_EQ, cFE, 1.0/dt ) ;
      if( AA != 0 )
      {
         FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, vv, 1.0 ) ;
      }

      double ww = cFE->weight_of_IP() ;
      if( FE::geometry() == FE::axisymmetrical )
      {
         ww *= 2.0 * PEL::pi() * cFE->coordinates_of_IP()->coordinate( 0 ) ;
      }
      vol += val_exp * ww ;
   }
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Mi_Mi( PDE_DiscreteField const* mi,
                              double sigma )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Mi_Mi" ) ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector vv( nb_dims ) ;

   double x = 0.0 ;
   double y = 0.0 ;
   if( MOB_EXP )
   {
      x =  cFE->value_at_IP( C1_EXP, L_EXPLICIT ) ;
      y =  cFE->value_at_IP( C2_EXP, L_EXPLICIT ) ;
   }
   else
   {
      x =  cFE->value_at_IP( C1, L_UPDATE ) ;
      y =  cFE->value_at_IP( C2, L_UPDATE ) ;
   }
   double g = (1.-x)*(1.-x)*(1.-y)*(1.-y)*(x+y)*(x+y) ;

   double mo = (MOB_CST + MOB_DEG * g )/sigma ;

   for( size_t d=0 ; d<nb_dims ; ++d )
   {
      vv( d ) = - mo * cFE->gradient_at_IP( mi, L_UPDATE, d ) ;
   }
   FE::add_grad_row( ELEMENT_EQ, cFE, vv ) ;

   //????? on devrait pouvoir utiliser la version _S : ce sont
   //????? des champs différents, mais avec la meme discretisation
   //
   //????? si on n'est pas à l'itération 1 et MOX_EXP ???
   FE::add_grad_row_grad_col_NS( ELEMENT_EQ, cFE, mo ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Mi_Mi_MobMod( PDE_DiscreteField const* mi,
                                   PDE_DiscreteField const* mi_exp,
                                   double sigma )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Mi_Mi_MobMod" ) ;
   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector vv( nb_dims ) ;

   double x = 0.0 ;
   double y = 0.0 ;
   if( MOB_EXP )
   {
      x =  cFE->value_at_IP( C1_EXP, L_EXPLICIT ) ;
      y =  cFE->value_at_IP( C2_EXP, L_EXPLICIT ) ;
   }
   else
   {
      x =  cFE->value_at_IP( C1, L_UPDATE ) ;
      y =  cFE->value_at_IP( C2, L_UPDATE ) ;
   }
   double g = (1.-x)*(1.-x)*(1.-y)*(1.-y)*(x+y)*(x+y) ;

   double mo = ( MOB_CST + MOB_DEG * g )/sigma ;
   double moc = ( MOB_CST + MOB_DEG * MOB_MAX )/sigma ;

   PEL_ASSERT( moc - mo > 0 ) ;

   for( size_t d=0 ; d<nb_dims ; ++d )
   {
      vv( d ) = - moc * cFE->gradient_at_IP( mi, L_UPDATE, d )
                + ( moc - mo ) * cFE->gradient_at_IP( mi_exp, L_UPDATE, d );
   }
   FE::add_grad_row( ELEMENT_EQ, cFE, vv ) ;

   //????? on devrait pouvoir utiliser la version _S : ce sont
   //????? des champs différents, mais avec la meme discretisation
   //
   //????? si on n'est pas à l'itération 1 et MOX_EXP ???
   FE::add_grad_row_grad_col_NS( ELEMENT_EQ, cFE, moc ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Ci_Ci( size_t i, double cap, double theta )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Ci_Ci" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector vv( nb_dims ) ;

   double x = cFE->value_at_IP( C1, L_UPDATE ) ;
   double y = cFE->value_at_IP( C2, L_UPDATE ) ;
   double s = cFE->value_at_IP( C1_EXP, L_EXPLICIT ) ;
   double t = cFE->value_at_IP( C2_EXP, L_EXPLICIT ) ;

   double dF = BULK_MU->DDiF( x, y, s, t, i, EPS ) ;

   FE::add_row( ELEMENT_EQ, cFE, dF ) ;

   for( size_t d=0 ; d<nb_dims ; ++d )
   {
      if( i==1 )
      {
         vv( d ) = theta * cap * cFE->gradient_at_IP( C1, L_UPDATE, d ) ;
         if( theta != 1.0 )
         {
            vv( d ) += ( 1.0 - theta ) * cap *
                       cFE->gradient_at_IP( C1_EXP, L_EXPLICIT, d ) ;
         }
      }
      else if( i==2 )
      {
         vv( d ) =  theta * cap * cFE->gradient_at_IP( C2, L_UPDATE, d ) ;
         if( theta != 1.0 )
         {
            vv( d ) += ( 1.0 - theta ) * cap *
                       cFE->gradient_at_IP( C2_EXP, L_EXPLICIT, d ) ;
         }
      }
   }
   FE::add_grad_row( ELEMENT_EQ, cFE, vv ) ;

   FE::add_grad_row_grad_col_NS( ELEMENT_EQ, cFE, -cap*theta ) ;

   double d2F = BULK_MU->dj_DDiF( x, y, s, t, i, i, EPS ) ;

   FE::add_row_col_NS( ELEMENT_EQ, cFE, -d2F ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Ci_Cj( size_t i, size_t j )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Ci_Cj" ) ;

   double x = cFE->value_at_IP( C1, L_UPDATE ) ;
   double y = cFE->value_at_IP( C2, L_UPDATE ) ;

   double s = cFE->value_at_IP( C1_EXP, L_EXPLICIT ) ;
   double t = cFE->value_at_IP( C2_EXP, L_EXPLICIT ) ;

   double d2F = BULK_MU->dj_DDiF( x, y, s, t, i, j, EPS );

   FE::add_row_col_NS( ELEMENT_EQ, cFE, -d2F ) ;
}

//----------------------------------------------------------------------
void
CH_CahnHilliard:: loop_Ci_Mi( PDE_DiscreteField const* mi )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: loop_Ci_Mi" ) ;

   FE::add_row( ELEMENT_EQ, cFE, -cFE->value_at_IP( mi, L_UPDATE ) ) ;
   if( ITER==1 ) FE::add_row_col_NS( ELEMENT_EQ, cFE, 1.0 ) ;
}

//---------------------------------------------------------------------------
void
CH_CahnHilliard:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CahnHilliard:: print" ) ;

   FE_OneStepIteration:: print( os, indent_width ) ;

   std::string const s( indent_width+3, ' ' ) ;

   os << s << "Mobility: " << std::endl ;
   os << s << "  explicit: " ;
   if( MOB_EXP )
   {
      os << "yes" << std::endl ;
   }
   else
   {
      os << "no" << std::endl ;
   }
   os << s << "  constant   part: " << MOB_CST << std::endl ;
   os << s << "  degenerate part: " << MOB_DEG << std::endl ;
   os << s << "  discretization trick: " ;
   if( MOB_MOD )
   {
      os << "yes" << std::endl ;
   }
   else
   {
      os << "no" << std::endl ;
   }
}

//internal--------------------------------------------------------------
void
CH_CahnHilliard_ERROR:: n0( std::string const& fname  )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** CH_CahnHilliard:" << endl << endl ;
   msg << "    unable to open file \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

