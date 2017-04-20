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

#include <AP_TutorialFV.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_Mpolyhedron.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
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

AP_TutorialFV const* AP_TutorialFV::PROTOTYPE = new AP_TutorialFV() ;

//---------------------------------------------------------------------------
AP_TutorialFV:: AP_TutorialFV( void )
//---------------------------------------------------------------------------
   : PEL_Application( "AP_TutorialFV" )
{
}

//---------------------------------------------------------------------------
AP_TutorialFV:: ~AP_TutorialFV( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
AP_TutorialFV*
AP_TutorialFV:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialFV:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_TutorialFV* result = new AP_TutorialFV( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_TutorialFV:: AP_TutorialFV( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
{
   PEL_LABEL( "AP_TutorialFV:: AP_TutorialFV" ) ;

   PEL_ModuleExplorer const* ee = 
                      exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   TT = dom->set_of_discrete_fields()->item( "temperature" ) ;

   sFE = dom->create_CursorFEside( this ) ;
   sFE->require_field_calculation( TT, PDE_LocalFE::node ) ;

   bFE = dom->create_LocalFEbound( this ) ;
   bFE->require_field_calculation( TT, PDE_LocalFE::node ) ;

   BCs = dom->set_of_boundary_conditions() ;

   ee = exp->create_subexplorer( 0, "AP_TutorialFV" ) ;
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

//---------------------------------------------------------------------------
void 
AP_TutorialFV:: run( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialFV:: run" ) ;

   loop_on_sides() ;

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

//------------------------------------------------------------------------
void
AP_TutorialFV:: loop_on_sides( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialFV:: loop_on_sides" ) ;

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      size_t n_K = fe_K->global_node( TT, 0 ) ;
      size_t n_L = fe_L->global_node( TT, 0 ) ;

      double d_KL = sFE->distance_to_adjacent_finite_volume_center( 0 ) +
                    sFE->distance_to_adjacent_finite_volume_center( 1 ) ;

      double xx = CONDUCTIVITY * sFE->polyhedron()->measure() / d_KL ;

      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;
      size_t u_L = NMB->global_unknown_for_DOF( n_L, 0 ) ;

      A->add_to_item( u_K, u_K,  xx ) ;
      A->add_to_item( u_L, u_K, -xx ) ;

      A->add_to_item( u_K, u_L, -xx ) ;
      A->add_to_item( u_L, u_L,  xx ) ;
   }
}

//------------------------------------------------------------------------
void
AP_TutorialFV:: loop_on_bounds( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "AP_TutorialFV:: loop_on_bounds" ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Color const* color = bFE->color() ;
      if( BCs->has_BC( color, TT ) )
      {
         PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, TT ) ;

         size_t n_K = bFE->global_node( TT, 0 ) ;
         size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;

         string const& type = ee->string_data( "type" ) ;
         if( type == "Dirichlet" )
         {
            double val = ee->double_data( "value" ) ;
            double xx = CONDUCTIVITY * bFE->polyhedron()->measure() / 
                        bFE->distance_to_adjacent_finite_volume_center() ;
            A->add_to_item( u_K, u_K, xx ) ;
            F->add_to_item( u_K, xx * val ) ; 
         }
         else if( type == "imposed_flux" )
         {
            double xx = bFE->polyhedron()->measure() * 
                        ee->double_data( "flux_value" ) ;
            F->add_to_item( u_K, xx ) ;
         }
      }
   }
}
