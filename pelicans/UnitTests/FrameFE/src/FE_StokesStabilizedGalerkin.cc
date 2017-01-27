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

#include <FE_StokesStabilizedGalerkin.hh>

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_TwoBlocksMethod.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>

using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE_StokesStabilizedGalerkin const* 
FE_StokesStabilizedGalerkin:: PROTOTYPE = new FE_StokesStabilizedGalerkin() ;

//----------------------------------------------------------------------
FE_StokesStabilizedGalerkin:: FE_StokesStabilizedGalerkin( void )
//----------------------------------------------------------------------
   : PEL_Application( "Stokes_StabilizedGalerkin" )
   , UU( 0 )
   , PP( 0 )
   , SAVER( 0 )
   , ELEMENT_EQ( 0 )
   , cFE( 0 )
   , STABILIZATION( none )
   , ALPHA( 0.0 )
{
}

//----------------------------------------------------------------------
FE_StokesStabilizedGalerkin*
FE_StokesStabilizedGalerkin:: create_replica( PEL_Object* a_owner, 
                            PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StokesStabilizedGalerkin:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   FE_StokesStabilizedGalerkin* result = new FE_StokesStabilizedGalerkin( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_StokesStabilizedGalerkin:: FE_StokesStabilizedGalerkin( 
                                       PEL_Object* a_owner, 
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , UU( 0 )
   , PP( 0 )
   , SAVER( 0 )
   , ELEMENT_EQ( PDE_LocalEquation::create( this )  )
   , cFE( 0 )
   , STABILIZATION( none )
   , ALPHA( 0.0 )
   , S( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
   , P_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_ModuleExplorer const* sexp = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = 
       PDE_DomainAndFields::create( a_owner, sexp, PEL_Exec::communicator() ) ;
   sexp->destroy() ; sexp = 0 ;
   
   UU = dom->set_of_discrete_fields()->item( "velocity" ) ;
   PP = dom->set_of_discrete_fields()->item( "pressure" ) ;

   SAVER = dom->result_saver() ;

   cFE = dom->create_LocalFEcell( this ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::N  ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::N ) ;

   sexp = exp->create_subexplorer( 0, "FE_StokesStabilizedGalerkin" ) ;

   string const& type = sexp->string_data( "type" ) ;
   if( type == "Galerkin" ) 
   {
      STABILIZATION = none ;
   }
   else if( type == "Galerkin_Least_Square" )
   {
      STABILIZATION = GLS ;
      ALPHA = sexp->double_data( "stabilization_constant" ) ;
      cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value( sexp, "type",
                               "\"Galerkin\" or \"Galerkin_Least_Square\"" ) ;
   }

   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU, 
                                          "sequence_of_the_components",
                                          true ) ;
   PDE_LinkDOF2Unknown* pp_link = PDE_LinkDOF2Unknown::create( 0, PP, true ) ;

   NMB_U = PDE_SystemNumbering::create( this, uu_link ) ;
   NMB_P = PDE_SystemNumbering::create( this, pp_link ) ;
   
   PEL_ModuleExplorer const* ee = 
                      sexp->create_subexplorer( 0, "LA_TwoBlocksMethod" ) ;
   SOLVER = LA_TwoBlocksMethod::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = sexp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   B = A->create_matrix( this ) ;
   C = A->create_matrix( this ) ;
   
   U = A->create_vector( this ) ;
   P = A->create_vector( this ) ;
   if( SOLVER->S_is_required() ) S = A->create_vector( this ) ;
   
   F = A->create_vector( this ) ; //??? RHS(0,0) ;
   G = A->create_vector( this ) ; //??? RHS(1,0) ;
   
   SOLVER->set_matrix_prototype( A ) ;
   
   FE::set_geometry( FE::cartesian ) ;
   
   sexp->destroy() ; sexp = 0 ;
}

//----------------------------------------------------------------------
FE_StokesStabilizedGalerkin:: ~FE_StokesStabilizedGalerkin( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_StokesStabilizedGalerkin:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StokesStabilizedGalerkin:: run" ) ;

   NMB_U->reset() ;
   NMB_P->reset() ;
   
   size_t nv_glob = NMB_U->nb_global_unknowns() ;
   size_t np_glob = NMB_P->nb_global_unknowns() ;

   size_t nv_loc = NMB_U->nb_unknowns_on_current_process() ;
   size_t np_loc = NMB_P->nb_unknowns_on_current_process() ;

   re_initialize_matrices_and_vectors( nv_glob, np_glob, nv_loc, np_loc ) ;

   SOLVER->re_initialize_internals( nv_glob, np_glob, nv_loc, np_loc ) ;

   U_LOC->re_initialize( NMB_U->link()->unknown_vector_size() ) ;
   P_LOC->re_initialize( NMB_P->link()->unknown_vector_size() ) ;
   
   NMB_U->define_scatters( U ) ;
   NMB_P->define_scatters( P ) ;
   
   //??? to be given by the data deck
   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_5" ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( UU, UU ) ;

      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 
                              UU->nb_components(),
                              cFE->col_field_node_connectivity(), 
                              UU->nb_components() ) ;

      cFE->start_IP_iterator( qrp ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() ) 
      {
         FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB_U ) ;

      cFE->set_row_and_col_fields( PP, UU ) ;

      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 
                              1,
                              cFE->col_field_node_connectivity(), 
                              UU->nb_components() ) ;

      cFE->start_IP_iterator( qrp ) ;      
      for(  ; cFE->valid_IP() ; cFE->go_next_IP() ) 
      {
         FE::add_row_div_col( ELEMENT_EQ, cFE, -1.0 ) ;
      }
      PDE::assemble_in_matrix_vector_0( B, G, ELEMENT_EQ, NMB_P, NMB_U ) ;

      if( STABILIZATION == GLS )
      {
         cFE->set_row_and_col_fields( PP, PP ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), 1,
                                 cFE->col_field_node_connectivity(), 1 ) ;

         double coef = cFE->polyhedron()->inter_vertices_maximum_distance() ;
         coef = - ALPHA * coef*coef ;

         cFE->start_IP_iterator( qrp ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() ) 
         {
            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, coef ) ;
         }
         PDE::assemble_in_matrix_vector_0( C, G, ELEMENT_EQ, NMB_P ) ;
      }
   }

   estimate_unknowns() ;
   
   LA_Scatter const* sca = NMB_U->scatter() ;
   sca->get( U, U_LOC ) ;
   UU->update_free_DOFs_value( 0, U_LOC, NMB_U->link() ) ;

   sca = NMB_P->scatter() ;
   sca->get( P, P_LOC ) ;
   PP->update_free_DOFs_value( 0, P_LOC, NMB_P->link() ) ; 

   SAVER->start_cycle() ;
   SAVER->save_grid() ;
   SAVER->save_fields( 0 ) ;
   SAVER->terminate_cycle() ;
}

//----------------------------------------------------------------------
void
FE_StokesStabilizedGalerkin:: estimate_unknowns( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_StokesStabilizedGalerkin:: estimate_unknowns" ) ;
   
   A->synchronize() ;
   B->synchronize() ;
   C->synchronize() ;
   F->synchronize() ;
   G->synchronize() ;
   
   if( S != 0 )
   {
      S->set( 1.0 ) ;
      SOLVER->set_S( S ) ;
   }
   SOLVER->set_system( A, B, F, G, C ) ;
   
   SOLVER->estimate_unknowns( false, U, false, P ) ;
   
   if( ! SOLVER->successful_estimation() )
      PEL_Error::object()->raise_plain( 
                           "FE_StokesStabilizedGalerkin:: solution failure!") ;
   
   SOLVER->unset_system() ;
}

//-----------------------------------------------------------------------
void
FE_StokesStabilizedGalerkin:: re_initialize_matrices_and_vectors( 
                                                 size_t nv_glob,
                                                 size_t np_glob,
                                                 size_t nv_loc,
                                                 size_t np_loc )
//-----------------------------------------------------------------------
{
   PEL_LABEL( 
       "FE_StokesStabilizedGalerkin:: re_initialize_matrices_and_vectors" ) ;

   A->re_initialize( nv_glob, nv_glob, nv_loc, nv_loc ) ;
   B->re_initialize( np_glob, nv_glob, np_loc, nv_loc ) ;
   C->re_initialize( np_glob, np_glob, np_loc, nv_loc ) ;

   U->re_initialize( nv_glob, nv_loc ) ;
   P->re_initialize( np_glob, np_loc ) ;
   S->re_initialize( np_glob, np_loc ) ; //????? si S non nul ???

   F->re_initialize( nv_glob, nv_loc ) ;
   G->re_initialize( np_glob, np_loc ) ;
}


