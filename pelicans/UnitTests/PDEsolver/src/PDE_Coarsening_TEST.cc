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

#include <PDE_Coarsening_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_BlockSeqMatrix.hh>
#include <LA_MatrixIterator.hh>

#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_AlgebraicCoarsener.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>
#include <PDE_SystemNumbering.hh>

// --> This dependency upon the FrameFE library could be avoided by
// implementing here the methods of FE that are used
// --> Nevertheless, this test is really a test of the PDEsolver library
#include <FE.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::ostringstream ;
using std::setw ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

struct PDE_Coarsening_TEST_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
} ;

//---------------------------------------------------------------------------
PDE_Coarsening_TEST*
PDE_Coarsening_TEST:: REGISTRATOR = new PDE_Coarsening_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_Coarsening_TEST:: PDE_Coarsening_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_AdapterCHARMS",
                     "PDE_Coarsening_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_Coarsening_TEST:: ~PDE_Coarsening_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_Coarsening_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Coarsening_TEST:: process_one_test" ) ;

   TEST_NAME = exp->name() ;
   MY_EPS = exp->double_data( "dbl_epsilon" ) ;
   MY_MIN = exp->double_data( "dbl_minimum" ) ;

   PEL_ModuleExplorer* se =
      exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   DOM = PDE_DomainAndFields::create( 0, se, PEL_Exec::communicator() ) ;
   se->destroy() ; se = 0 ;
   cFE = DOM->create_LocalFEcell( 0 ) ;
   DA  = DOM->adapter_CHARMS() ;
   if( DA == 0 ) PDE_Coarsening_TEST_ERROR::n0() ;

   ELEMENT_EQ = PDE_LocalEquation::create( 0 ) ;

   QRP = GE_QRprovider::object(
                        exp->string_data( "quadrature_rule_provider" ) ) ;


   se = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   MAT_PROTO = LA_Matrix::make( 0, se ) ;
   se->destroy() ;

   FE::set_geometry( FE::cartesian ) ;

   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;

   size_t nb_fields = sdf->nb_fields() ;

   PEL_Vector* cmats = PEL_Vector::create( 0, nb_fields+1 ) ;
   PEL_Vector* links = PEL_Vector::create( 0, nb_fields ) ;
   PEL_Vector* nmbs = PEL_Vector::create( 0, nb_fields+1 ) ;
   for( size_t ii=0 ; ii<nb_fields ; ++ii )
   {
      PDE_DiscreteField const* ff = sdf->item( ii ) ;

      cFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
      cFE->require_field_calculation( ff, PDE_LocalFE::dN ) ;

      PDE_LinkDOF2Unknown* ff_link = PDE_LinkDOF2Unknown::create( 0, ff,
                                           "sequence_of_the_components",
                                           true ) ;
      links->set_at( ii, ff_link ) ;


      PDE_LinkDOF2Unknown* link = ff_link->create_clone( 0 ) ; 
      PEL_Vector* vec = PEL_Vector::create( 0, 1 ) ;
      vec->set_at( 0, link ) ;
      PDE_SystemNumbering* ff_nmb = PDE_SystemNumbering::create( this, vec,
                                          "sequence_of_the_discrete_fields" ) ;
      nmbs->set_at( ii, ff_nmb ) ;

      LA_Matrix* mm = create_assembled_matrix( cmats, ff_nmb ) ;
      vec->destroy() ; vec=0 ;

      cmats->set_at( ii, mm ) ;
   }
   std::string ordering = "sequence_of_the_discrete_fields" ;
   PDE_SystemNumbering* nmb =
                    PDE_SystemNumbering::create( this, links, ordering ) ;
   LA_Matrix* mm = create_assembled_matrix( cmats, nmb ) ;
   cmats->set_at( nb_fields, mm ) ;
   nmbs->set_at( nb_fields, nmb ) ;


   DA->reset() ;
   bool keep_adapting = false ;
   do
   {
      DA->adapt() ;
      keep_adapting = DA->something_changed() ;
      if( keep_adapting )
      {
         DOM->apply_requests_of_DOFs_values_modules( true ) ;
      }
   } while( keep_adapting ) ;

   for( size_t ii=0 ; ii<nb_fields ; ++ii )
   {
      LA_Matrix const* coarsest_mat =
         static_cast< LA_Matrix* >( cmats->at( ii ) ) ;
      PDE_SystemNumbering* ff_nmb =
         static_cast< PDE_SystemNumbering* >( nmbs->at( ii ) ) ;

      check_coarsening( ff_nmb, coarsest_mat ) ;
   }

   LA_Matrix const* coarsest_mat =
                   static_cast< LA_Matrix* >( cmats->at( nb_fields ) ) ;
   PDE_SystemNumbering* system_nmb =
                   static_cast< PDE_SystemNumbering* >( nmbs->at( nb_fields ) ) ;
   check_coarsening( system_nmb, coarsest_mat ) ;

   cmats->destroy() ;
   links->destroy() ;
   nmbs->destroy() ;

   ELEMENT_EQ->destroy() ;
   cFE->destroy() ;
   MAT_PROTO->destroy() ;
   DOM->destroy() ;
}

//---------------------------------------------------------------------------
void
PDE_Coarsening_TEST:: check_coarsening( PDE_SystemNumbering* nmb,
                                        LA_Matrix const* coarsest_mat )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Coarsening_TEST:: check_coarsening" ) ;

   PDE_AlgebraicCoarsener* coar = DA->algebraic_coarsener() ;

   nmb->reset() ;

   string fnames( "|" ) ;
   for( size_t i=0 ; i<nmb->nb_links() ; ++i )
   {
      PDE_LinkDOF2Unknown const* link = nmb->link( i ) ;
      fnames += link->field()->name() + "|" ;
   }

   LA_Matrix const* finest_mat = create_assembled_matrix( 0, nmb ) ;

   coar->prepare_for_coarsening( nmb ) ;

   size_t nb_levels = coar->nb_levels() ;
   if( !( nb_levels > 1 ) ) PDE_Coarsening_TEST_ERROR::n1() ;

   std::vector< LA_Matrix const* > coarse_to_fine( nb_levels-1, 0 ) ;
   std::vector< LA_Matrix* > aa( nb_levels-1, 0 ) ;

   while( coar->coarsening_is_possible() )
   {
      LA_Matrix* pr = MAT_PROTO->create_matrix( 0 ) ;

      coar->do_one_coarsening() ;
      coar->build_current_prolongation_matrix( pr ) ;

      size_t fine_l   = coar->current_fine_level() ;
      size_t coarse_l = fine_l - 1 ;
      
      coarse_to_fine[coarse_l] = pr ;

      LA_Matrix* mat = MAT_PROTO->create_matrix( 0 ) ;
      size_t nb_unks = pr->nb_cols() ;
      mat->re_initialize( nb_unks, nb_unks ) ;
      aa[coarse_l] = mat ;
   }
   
   LA_Matrix* dummy_0 = MAT_PROTO->create_matrix( 0 ) ;
   LA_Matrix* dummy_1 = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t level=nb_levels-1 ; level!=0 ; --level )
   {
      LA_Matrix const* fine_A = ( level==(nb_levels-1) ?
                                        finest_mat : aa[level] ) ;
      LA_Matrix* coar_A = aa[level-1] ;

      dummy_0->re_initialize( coar_A->nb_rows(),
                                                fine_A->nb_rows() ) ;

      LA_Matrix const* pr = coarse_to_fine[ level-1 ] ;

      dummy_1->re_initialize( pr->nb_cols(), pr->nb_rows() ) ;
      dummy_1->add_tMat( pr ) ;
      dummy_1->synchronize() ;

      dummy_0->add_Mat_Mat( dummy_1, fine_A ) ;
      dummy_0->synchronize() ;
      coar_A->add_Mat_Mat( dummy_0, pr ) ;    // coar_A = dummy * Pr
      coar_A->synchronize() ;
    }
   dummy_0->destroy() ;
   dummy_1->destroy() ;

   compare_matrices( fnames, coarsest_mat, aa[0] ) ;

   for( size_t i=0 ; i<nb_levels-1 ; ++i )
   {
      coarse_to_fine[i]->destroy() ;
      aa[i]->destroy() ;
   }
   finest_mat->destroy() ;
}

//----------------------------------------------------------------------------
LA_Matrix*
PDE_Coarsening_TEST:: create_assembled_matrix(
                                       PEL_Object* a_owner,
                                       PDE_SystemNumbering* nmb )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Coarsening_TEST:: create_assembled_matrix" ) ;

   size_t nb_dims = cFE->nb_space_dimensions() ;
   doubleVector aa( nb_dims ) ;

   LA_Matrix* result = MAT_PROTO->create_matrix( a_owner ) ;

   result->re_initialize( nmb->nb_global_unknowns(), 
                          nmb->nb_global_unknowns() ) ;

   // terms are only assembled on diagonal block
   // --> could be extended
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      for( size_t i=0 ; i<nmb->nb_links() ; ++i )
      {
        PDE_LinkDOF2Unknown const* link = nmb->link( i ) ;
         PDE_DiscreteField const* ff = link->field() ;
         size_t nb_c = ff->nb_components() ;

         cFE->set_row_and_col_fields( ff, ff ) ;
         ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nb_c,
                                 cFE->col_field_node_connectivity(), nb_c ) ;

         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               aa( d ) = (double)d + 1.0 ;
            }
            FE::add_row_vvgrad_col( ELEMENT_EQ, cFE, aa, 1.0 ) ;

            FE::add_grad_row_grad_col_S( ELEMENT_EQ, cFE, 1.0 ) ;
         }
         PDE::assemble_in_matrix_0( result, ELEMENT_EQ, nmb, i, i ) ;
      }
   }

   result->synchronize() ;

   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_Coarsening_TEST:: compare_matrices( std::string const& fnames,
                                        LA_Matrix const* m1,
                                        LA_Matrix const* m2 )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Coarsening_TEST:: compare_matrices" ) ;

   bool ok = true ;
   LA_MatrixIterator* it1 = m1->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* it2 = m2->create_stored_item_iterator( 0 ) ;

   it1->start_all_items() ;
   it2->start_all_items() ;
   while( it1->is_valid() && it2->is_valid() && ok )
   {
      bool same = ( it1->row() == it2->row() )
               && ( it1->col() == it2->col() ) ;
      double xx1 = it1->item() ;
      double xx2 = it2->item() ;
      bool eq = PEL::double_equality( xx1, xx2, MY_EPS, MY_MIN  ) ;
      ok = ok && same && eq ;
      if( !ok ) display_error( it1->row(), it1->col(), xx1,
                               it2->row(), it2->col(), xx2 ) ;
      it1->go_next() ;
      it2->go_next() ;
   }
   ok = ok && ( !it1->is_valid() && !it2->is_valid() ) ;
   notify_one_test_result( TEST_NAME+"("+fnames+")", ok ) ;
   it1->destroy() ;
   it2->destroy() ;
}

//-----------------------------------------------------------------------
void
PDE_Coarsening_TEST:: display_error(
                              size_t row1, size_t col1, double xx_1,
                              size_t row2, size_t col2, double xx_2 ) const
//-----------------------------------------------------------------------
{
   std::ios_base::fmtflags original_flags = out().flags() ;
   out().setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;
   out() << setw( 5 ) << row1 << setw( 5 ) << col1 << "   :  "
         << std::setw( 20 ) << xx_1
         << endl
         << setw( 5 ) << row2 << setw( 5 ) << col2 << "   :  "
         << std::setw( 20 ) << xx_2
         << endl ;
   out().flags( original_flags ) ;
}

//internal--------------------------------------------------------------
void
PDE_Coarsening_TEST_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PDE_Coarsening_TEST: " << endl ;
   msg << "    a discretization adapter is requested" ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_Coarsening_TEST_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PDE_Coarsening_TEST: " << endl ;
   msg << "    the number of levels should be > 1" ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

