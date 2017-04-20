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

#include <AP_TutorialG.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <string>
using std::string ;

AP_TutorialG const* AP_TutorialG::PROTOTYPE = new AP_TutorialG() ;

//----------------------------------------------------------------------------
AP_TutorialG:: AP_TutorialG( void )
//----------------------------------------------------------------------------
   : PEL_Application( "AP_TutorialG" ) 
{
}

//----------------------------------------------------------------------------
AP_TutorialG:: ~AP_TutorialG( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
AP_TutorialG* 
AP_TutorialG:: create_replica( PEL_Object* a_owner, 
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialG:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_TutorialG* result = new AP_TutorialG( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
AP_TutorialG:: AP_TutorialG( PEL_Object* a_owner, 
                             PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
{
   PEL_LABEL( "AP_TutorialG:: AP_TutorialG" ) ;

   PEL_ModuleExplorer const* ee = 
                      exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   TT = dom->set_of_discrete_fields()->item( "temperature" ) ;

   ELEMENT_EQ = PDE_LocalEquation::create( this ) ;

   cFE = dom->create_LocalFEcell( this ) ;
   cFE->require_field_calculation( TT, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( TT, PDE_LocalFE::dN ) ;

   bFE = dom->create_LocalFEbound( this ) ;
   bFE->require_field_calculation( TT, PDE_LocalFE::N ) ;

   BCs = dom->set_of_boundary_conditions() ;

   ee = exp->create_subexplorer( 0, "AP_TutorialG" ) ;
   CONDUCTIVITY = ee->double_data( "conductivity" ) ;

   PDE_LinkDOF2Unknown* tt_link = PDE_LinkDOF2Unknown::create( 0, TT, true ) ;
   NMB = PDE_SystemNumbering::create( this, tt_link ) ;
   
   PEL_ModuleExplorer const* eee = ee->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, eee ) ;
   eee->destroy() ; eee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC = LA_SeqVector::create( this, NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;

   eee = ee->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, eee ) ;
   eee->destroy() ; eee = 0 ;
   
   ee->destroy() ; ee = 0 ;

   SAVER = dom->result_saver() ;
}

//----------------------------------------------------------------------------
void AP_TutorialG:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialG:: run" ) ;

   loop_on_cells() ;

   loop_on_bounds() ;

   A->synchronize() ;
   F->synchronize() ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
   
   NMB->scatter()->get( X, X_LOC ) ;
   TT->update_free_DOFs_value( 0, X_LOC, NMB->link() ) ;
   
   SAVER->start_cycle() ;
   SAVER->save_grid() ;
   SAVER->save_fields( 0 ) ;
   SAVER->terminate_cycle() ;
}

//----------------------------------------------------------------------------
void AP_TutorialG:: loop_on_cells( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialG:: loop_on_cells" ) ;

   PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
   PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( TT, TT ) ;

      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                              cFE->col_field_node_connectivity(), 1 ) ;

      cFE->start_IP_iterator( qrp ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         for( size_t i=0 ; i<cFE->nb_basis_functions( row ) ; ++i )
         {
            for( size_t j=0 ; j<cFE->nb_basis_functions( row ) ; ++j )
            {
               double xx = 0.0 ;
               for( size_t d=0 ; d<nb_dims ; ++d )
               {
                  xx += cFE->dN_at_IP( col, j, d ) * 
                        cFE->dN_at_IP( row, i, d ) ;
               }
               xx = xx * cFE->weight_of_IP() * CONDUCTIVITY ;
               ELEMENT_EQ->add_to_matrix( xx, i, j ) ;
            }
         }
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }
}

//----------------------------------------------------------------------------
void AP_TutorialG:: loop_on_bounds( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialG:: loop_on_bounds" ) ;

   PDE_LocalFE::field_id const row = PDE_LocalFE::row ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      bFE->set_row_and_col_fields( TT, TT ) ;

      ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 1,
                              bFE->col_field_node_connectivity(), 1 ) ;

      GE_Color const* color = bFE->color() ;
      if( BCs->has_BC( color, TT ) )
      {
         PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, TT ) ;
         string bc_type = ee->string_data( "type" ) ;
         if( bc_type=="imposed_flux" )
         {
            double flux = ee->double_data( "flux_value" ) ;

            bFE->start_IP_iterator( qrp ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {
               for( size_t i=0 ; i<bFE->nb_basis_functions( row ) ; ++i )
               {
                  double xx = bFE->weight_of_IP() * 
                              flux * bFE->N_at_IP( row, i ) ;
                  ELEMENT_EQ->add_to_vector( xx, i ) ;
               }
            }
         }
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }
}
