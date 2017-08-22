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

#include <FE_Vorticity.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_QRprovider.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <PDE.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>

#include <iostream>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE_Vorticity const* FE_Vorticity::PROTOTYPE = new FE_Vorticity() ;

//----------------------------------------------------------------------
FE_Vorticity:: FE_Vorticity( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_Vorticity" )
{
}

//----------------------------------------------------------------------
FE_Vorticity*
FE_Vorticity:: create_replica( PEL_Object* a_owner, 
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Vorticity:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   if( dom->nb_space_dimensions() != 2 )
   {
      PEL_Error::object()->raise_plain( "FE_Vorticity requires 2D" ) ;
   }

   FE_Vorticity* result = new FE_Vorticity( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_Vorticity:: FE_Vorticity( PEL_Object* a_owner,
			     PDE_DomainAndFields const* dom,
			     FE_SetOfParameters const* prms,
 			     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , VORT( dom->set_of_discrete_fields()->item( 
                                    exp->string_data( "vorticity" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , VV( dom->set_of_discrete_fields()->item( 
                                    exp->string_data( "velocity" ) ) )
   , L_VV( exp->int_data( "level_of_velocity" ) )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( "GE_QRprovider_5" ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "FE_Vorticity:: FE_Vorticity" ) ;

   check_field_storage_depth( VORT, L_UPDATE ) ;
   check_field_storage_depth( VV, L_VV ) ;

   cFE->require_field_calculation( VORT, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( VV, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( VV, PDE_LocalFE::dN ) ;

   PDE_LinkDOF2Unknown* ll = PDE_LinkDOF2Unknown::create( 0, VORT, true ) ;
   NMB = PDE_SystemNumbering::create( this, ll ) ;
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
FE_Vorticity:: ~FE_Vorticity( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_Vorticity:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_Vorticity:: save_other_than_time_and_fields( FE_TimeIterator const* t_it,
                                                PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Vorticity:: save_other_than_time_and_fields" ) ;

   start_total_timer( "FE_Vorticity:: save_other_than_time_and_fields" ) ;

   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   start_assembling_timer() ;
   // ---------------------
   
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( VORT, VORT ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() ) 
      {
         double vort = cFE->gradient_at_IP( VV, L_VV, 0, 1 )
 	                  - cFE->gradient_at_IP( VV, L_VV, 1, 0 ) ;
         FE::add_row( ELEMENT_EQ, cFE, vort ) ;

         FE::add_row_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }
   
   stop_assembling_timer() ;
   start_solving_timer() ;
   // ------------------
   
   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
   stop_solving_timer() ;
   
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, X_LOC ) ;
   VORT->update_free_DOFs_value( L_UPDATE, X_LOC, NMB->link() ) ;
   
   stop_total_timer() ;
}

