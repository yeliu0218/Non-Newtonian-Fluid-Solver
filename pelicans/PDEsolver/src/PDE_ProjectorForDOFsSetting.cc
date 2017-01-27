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

#include <PDE_ProjectorForDOFsSetting.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleVector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SystemNumbering.hh>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

//----------------------------------------------------------------------
PDE_ProjectorForDOFsSetting:: PDE_ProjectorForDOFsSetting( 
                                         PEL_Object* a_owner,
                                         PDE_DiscreteField* a_field,
                                         size_t a_field_level,
                                         PDE_DomainAndFields const* a_dom,
                                         PEL_ModuleExplorer* a_exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FIELD( a_field )
   , FIELD_LEVEL( a_field_level )
   , QRP( GE_QRprovider::object( a_exp->string_data( "QRprovider_name" ) ) )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , cFE( 0 )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "PDE_ProjectorForDOFsSetting:: PDE_ProjectorForDOFsSetting" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   cFE = a_dom->create_LocalFEcell( this ) ;
   cFE->require_field_calculation( FIELD, PDE_LocalFE::N ) ;

   PDE_LinkDOF2Unknown* ll = PDE_LinkDOF2Unknown::create( 0, FIELD,
                                                  "sequence_of_the_nodes",
                                                  true ) ;
   NMB = PDE_SystemNumbering::create( this, ll ) ;
   
   PEL_ModuleExplorer const* ee = a_exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   ee = a_exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_ProjectorForDOFsSetting:: ~PDE_ProjectorForDOFsSetting( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ProjectorForDOFsSetting:: ~PDE_ProjectorForDOFsSetting" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_ProjectorForDOFsSetting:: field( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ProjectorForDOFsSetting:: field" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PDE_DiscreteField const* result = FIELD ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;  
}

//----------------------------------------------------------------------
void
PDE_ProjectorForDOFsSetting:: project_and_update_field( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ProjectorForDOFsSetting:: project_and_update_field" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   NMB->reset() ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   size_t n_loc  = NMB->nb_unknowns_on_current_process() ;

   A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
   F->re_initialize( n_glob, n_loc ) ;
   X->re_initialize( n_glob, n_loc ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   size_t nbc = FIELD->nb_components() ;
   doubleVector val( nbc ) ;
   
   // Assemble :
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( FIELD, FIELD ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(),
                              FIELD->nb_components(),
                              cFE->col_field_node_connectivity(),
                              FIELD->nb_components() ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         doubleVector const& N = cFE->Ns_at_IP( row ) ;
         compute_value_at_IP( cFE, val ) ;
         size_t const nb_nodes = cFE->nb_basis_functions( row ) ;
         double const w = cFE->weight_of_IP() ;
         for( size_t i=0 ; i<nb_nodes ; ++i )
         {
            double const w_Ni = w*N(i) ;
            for( size_t j=i ; j<nb_nodes ; ++j )
            {
               double const xx = w_Ni*N(j) ;
               for( size_t i_comp=0 ; i_comp<nbc ; ++i_comp )
               {
                  ELEMENT_EQ->add_to_matrix( xx, i, j, i_comp, i_comp ) ;
                  if( i!=j )
                  {
                     ELEMENT_EQ->add_to_matrix( xx, j, i, i_comp, i_comp ) ;
                  }
               }
            }
            for( size_t i_comp=0 ; i_comp<nbc ; ++i_comp )
            {
               ELEMENT_EQ->add_to_vector( w_Ni*val(i_comp), i, i_comp ) ;
            }
         }
      }
      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMB ) ;
   }

   // Solve :
   A->synchronize() ;
   F->synchronize() ;
   
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;

   // Update :
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, X_LOC ) ;
   FIELD->update_free_DOFs_value( FIELD_LEVEL, X_LOC, NMB->link() ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_ProjectorForDOFsSetting:: add_field_requirement_on_cells(
                                           PDE_DiscreteField const* ff,
                                           int derivation_order )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ProjectorForDOFsSetting:: add_field_requirement_on_cells" ) ;
   PEL_CHECK( ff != 0 ) ;
   PEL_CHECK( derivation_order == PDE_LocalFE::N  ||
              derivation_order == PDE_LocalFE::dN ||
              derivation_order == PDE_LocalFE::d2N ) ;
   PEL_CHECK_INV( invariant() ) ;

   cFE->require_field_calculation( ff, derivation_order ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
PDE_ProjectorForDOFsSetting:: compute_value_at_IP_PRE( 
                                           PDE_LocalFEcell const* fe,
                                           doubleVector& result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( fe!=0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( result.size()==field()->nb_components() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_ProjectorForDOFsSetting:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

