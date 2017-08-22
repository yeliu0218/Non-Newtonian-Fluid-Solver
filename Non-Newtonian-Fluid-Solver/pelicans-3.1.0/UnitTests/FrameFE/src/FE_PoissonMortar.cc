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

#include <FE_PoissonMortar.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
using std::cout ;
using std::endl ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE_PoissonMortar const* 
FE_PoissonMortar:: PROTOTYPE = new FE_PoissonMortar() ;

//-------------------------------------------------------------------------
FE_PoissonMortar:: FE_PoissonMortar( void )
//-------------------------------------------------------------------------
   : FE_OneStepIterationOpen( "FE_PoissonMortar" ) 
{
}

//---------------------------------------------------------------------------
FE_PoissonMortar*
FE_PoissonMortar:: create_replica( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_PoissonMortar* result = 
                     new FE_PoissonMortar( a_owner, dom, prms, exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_PoissonMortar:: FE_PoissonMortar( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_OneStepIterationOpen( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( 
                                      exp->string_data( "current_field" ) ) )
   , UU_link( 0 )
   , L_UU( exp->int_data( "level_of_current" ) )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , ORIGIN( GE_Point::create( this, 0.0, 0.0 ) )
   , NMB( 0 )
   , LHS( 0 )
   , RHS( 0 )
   , UNK( 0 )
   , UNK_LOC( LA_SeqVector::create( this, 0 ) )
   , SOLVER( 0 )
{
   PEL_LABEL( "FE_PoissonMortar:: FE_PoissonMortar" ) ;
   
   UU_link = PDE_LinkDOF2Unknown::create( this, UU, 
                                          "sequence_of_the_components", 
                                          true ) ;
   
   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   
   if( exp->has_module( "discrete_system" ) )
   {
      PDE_LinkDOF2Unknown* lnk = PDE_LinkDOF2Unknown::create( 0, UU, 
                                     "sequence_of_the_components", true ) ;
      NMB = PDE_SystemNumbering::create( this, lnk ) ;
      
      PEL_ModuleExplorer const* ee = 
                         exp->create_subexplorer( 0, "discrete_system" ) ;
      PEL_ModuleExplorer const* eee = 
                         ee->create_subexplorer( 0, "LA_Matrix" ) ;
      LHS = LA_Matrix::make( this, eee ) ;
      eee->destroy() ; eee = 0 ;
      
      RHS = LHS->create_vector( this ) ;
      UNK = LHS->create_vector( this ) ;
      
      eee = ee->create_subexplorer( 0, "LA_Solver" ) ;
      SOLVER = LA_Solver::make( this, eee ) ;
      eee->destroy() ;
      ee->destroy() ;
   }
}

//-------------------------------------------------------------------------
FE_PoissonMortar:: ~FE_PoissonMortar( void )
//-------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
size_t
FE_PoissonMortar:: nb_unknowns( void ) const
//------------------------------------------------------------------------
{
   return( 1 ) ;
}

//------------------------------------------------------------------------
PDE_DiscreteField*
FE_PoissonMortar:: field( size_t i_unk ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: field" ) ;
   PEL_CHECK_PRE( field_PRE( i_unk ) ) ;
   
   PDE_DiscreteField* result = UU ;
   
   PEL_CHECK_POST( field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
FE_PoissonMortar:: level_of_field( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: level_of_field" ) ;
   PEL_CHECK_PRE( level_of_field_PRE( i_unk ) ) ;

   size_t result = L_UU ;

   PEL_CHECK_POST( level_of_field_POST( result, i_unk ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
PDE_LinkDOF2Unknown const* 
FE_PoissonMortar:: link_DOF_2_unknown( size_t i_unk ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( link_DOF_2_unknown_PRE( i_unk ) ) ;

   PDE_LinkDOF2Unknown const* result = UU_link ;

   PEL_CHECK_POST( link_DOF_2_unknown_POST( result, i_unk ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_PoissonMortar:: assemble_contribution( 
                                   FE_TimeIterator const* t_it,
                                   LA_Matrix* matrix,
                                   LA_Vector* vector,
                                   PDE_SystemNumbering const* nmb,
                                   size_t i_link ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: assemble_contribution" ) ;
   PEL_CHECK_PRE( assemble_contribution_PRE( t_it, matrix, vector, 
                                             nmb, i_link ) ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
      }
      if( cFE->polyhedron()->contains( ORIGIN )  ) 
      { 
         add_row_at_pt( ELEMENT_EQ, cFE, ORIGIN ) ;
      }
      PDE::assemble_in_matrix_vector_0( matrix, vector, ELEMENT_EQ, 
                                        nmb, i_link, i_link ) ;
   }
}

//--------------------------------------------------------------------------
void
FE_PoissonMortar:: update_DOFs( LA_Vector* vector,
                                PDE_SystemNumbering const* nmb,
                                size_t i_link ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: update_DOFs" ) ;
   PEL_CHECK_PRE( update_DOFs_PRE( vector, nmb, i_link ) ) ;

   LA_Scatter const* sca = nmb->scatter( i_link ) ;
   PDE_LinkDOF2Unknown const* link = nmb->link( i_link ) ;
   UNK_LOC->re_initialize( link->unknown_vector_size() ) ;
   sca->get( vector, UNK_LOC ) ;
   UU->update_free_DOFs_value( L_UU, UNK_LOC, link ) ;
}

//-------------------------------------------------------------------------
void
FE_PoissonMortar:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   
   start_total_timer( "FE_PoissonMortar:: do_one_inner_iteration" ) ;
   // --------------
 
   check_assembled_system( NMB, "do_one_inner_iteration" ) ;
   
   NMB->reset() ;
   size_t dim = NMB->nb_global_unknowns() ;
   LHS->re_initialize( dim, dim ) ;
   RHS->re_initialize( dim ) ;
   UNK->re_initialize( dim ) ;
   
   NMB->define_scatters( RHS ) ;
   
   start_assembling_timer() ;
   // -------------------
   
   assemble_contribution( t_it, LHS, RHS, NMB, 0 ) ;
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ----------------
   
   SOLVER->set_matrix( LHS ) ;
   SOLVER->solve( RHS, UNK ) ;
   SOLVER->unset_matrix() ;
   
   stop_solving_timer() ;
   // ---------------
   
   update_DOFs( UNK, NMB, 0 ) ;
   
   stop_total_timer() ;
   // -------------
}

//---------------------------------------------------------------------------
void
FE_PoissonMortar:: add_row_at_pt( PDE_LocalEquation* leq,
                                  PDE_LocalFEcell* fe,
                                  GE_Point const* pt ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_PoissonMortar:: add_row_at_pt" ) ;

   fe->set_calculation_point( pt ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = fe->N_at_pt( fe->field(row), i ) ;
      leq->add_to_vector( xx, i ) ;
   }
}

