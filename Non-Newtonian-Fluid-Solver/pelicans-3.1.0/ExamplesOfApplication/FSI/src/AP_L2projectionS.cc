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

#include <AP_L2projectionS.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <LA_Matrix.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_TimeIterator.hh>

#include <AP_ConstitutiveLaw.hh>
#include <AP_KinematicState.hh>

#include <iostream>

using std::cout ; using std::endl ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

AP_L2projectionS const* 
AP_L2projectionS::PROTOTYPE = new AP_L2projectionS() ;

//---------------------------------------------------------------------------
AP_L2projectionS:: AP_L2projectionS( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_L2projectionS" )
{
}

//---------------------------------------------------------------------------
AP_L2projectionS*
AP_L2projectionS:: create_replica( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms, 
                                   PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_L2projectionS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_L2projectionS* result = new AP_L2projectionS( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_L2projectionS:: AP_L2projectionS( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , SS( dom->set_of_discrete_fields()->item( 
                               exp->string_data( "structure_stress" ) ) )
   , L_SS( exp->int_data( "level_of_stress" ) )
   , SD( dom->set_of_discrete_fields()->item( 
                               exp->string_data( "structure_displacement" ) ) )
   , L_SD( exp->int_data( "level_of_displacement" ) )
   , NB_COMP_S( PEL::bad_index() )
   , NB_COMP_D( PEL::bad_index() )
   , ST( AP_KinematicState::create( this, dom->nb_space_dimensions() ) ) 
   , LAW( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "AP_L2projectionS:: AP_L2projectionS" ) ;

   NB_COMP_S = SS->nb_components() ;
   NB_COMP_D = SD->nb_components() ;

   cFE->require_field_calculation( SS, PDE_LocalFE::N )  ;

   cFE->require_field_calculation( SD, PDE_LocalFE::N )  ;
   cFE->require_field_calculation( SD, PDE_LocalFE::dN ) ;
 
   PEL_ModuleExplorer const* ee =
                      exp->create_subexplorer( 0, "AP_ConstitutiveLaw" ) ;
   LAW = AP_ConstitutiveLaw::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   PDE_LinkDOF2Unknown* lnk = PDE_LinkDOF2Unknown::create( 0, SS, 
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
}

//---------------------------------------------------------------------------
AP_L2projectionS:: ~AP_L2projectionS( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_L2projectionS:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_L2projectionS:: save_other_than_time_and_fields( 
                                                  FE_TimeIterator const* t_it,
                                                  PDE_ResultSaver* rs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_L2projectionS:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   start_total_timer( "AP_L2projectionS:: save_other_than_time_and_fields" ) ;
   // --------------

   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   start_assembling_timer() ;
   // -------------------

   doubleVector s( NB_COMP_S ) ;
   doubleArray2D S( 3, 3 ) ;
   doubleArray4D C( 3, 3, 3, 3 ) ;
 
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( SS, SS ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), NB_COMP_S,
                              cFE->col_field_node_connectivity(), NB_COMP_S ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         FE::add_row_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
  
         ST->set_state( SD, L_SD, cFE ) ;
         LAW->update_S_dSdE( ST, S, C ) ;
         
         stress_to_smooth( S, s ) ;

         FE::add_row( ELEMENT_EQ, cFE, s ) ;
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }

   stop_assembling_timer() ;
   start_solving_timer() ;
   // ----------------

   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
  
   stop_solving_timer() ;
   // ---------------
   
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, X_LOC ) ;
   SS->update_free_DOFs_value( L_SS, X_LOC, NMB->link() ) ;

   stop_total_timer() ;
   // -------------
}
 
//---------------------------------------------------------------------------
void
AP_L2projectionS:: stress_to_smooth( doubleArray2D const& S,
                                     doubleVector& s ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_L2projectionS:: stress_to_smooth" ) ;

   s( 0 ) = S( 0, 0 ) ;
   s( 1 ) = S( 1, 1 ) ;
   s( 2 ) = S( 2, 2 ) ;
   s( 3 ) = S( 0, 1 ) ;
   s( 4 ) = S( 0, 2 ) ;
   s( 5 ) = S( 1, 2 ) ;  
}
