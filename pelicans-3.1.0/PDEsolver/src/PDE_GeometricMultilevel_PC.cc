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

#include <PDE_GeometricMultilevel_PC.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <PDE_AlgebraicCoarsener.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <LA_MatrixIterator.hh>
#include <LA_Matrix.hh>
#include <LA_DistMatrix.hh>
#include <LA_DistScatter.hh>

#include <LA_Vector.hh>

#include <size_t_vector.hh>
#include <intVector.hh>

#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::set ;
using std::string ; using std::ostringstream ;

struct PDE_GeometricMultilevel_PC_ERROR
{
   static void n0( std::string const& a_name ) ;
   static void n1( std::string const& a_name ) ;
   static void n2( void ) ;
   static void n3( std::string const& tname ) ;
} ;

std::map< std::string, PDE_GeometricMultilevel_PC* >
PDE_GeometricMultilevel_PC::OBJS ;

//----------------------------------------------------------------------
PDE_GeometricMultilevel_PC:: PDE_GeometricMultilevel_PC(
                                               std::string const& class_name )
//----------------------------------------------------------------------
   : LA_Preconditioner( class_name )
{
}

//----------------------------------------------------------------------
PDE_GeometricMultilevel_PC:: PDE_GeometricMultilevel_PC(
                                           PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
   , NAME( exp->has_entry( "name") ? exp->string_data( "name" ) : "unknown" )
   , NMB( 0 )
   , COAR( 0 )
   , MAT_PROTO( 0 )
   , AA( PEL_Vector::create( this, 0 ) )
   , FINEST_A( 0 )
   , B( PEL_Vector::create( this, 0 ) )
   , X( PEL_Vector::create( this, 0 ) )
   , RES( PEL_Vector::create( this, 0 ) )
   , COARSE_TO_FINE( PEL_Vector::create( this, 0 ) )
   , SMOO_LINE( PEL_Vector::create( this, 0 ) )
   , UNKNOWNS_LEVEL( PEL_Vector::create( this, 0 ) )
   , NB_LEVELS( 1 )
   , BUILD_OK( false )
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: PDE_GeometricMultilevel_PC" ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   MAT_PROTO = LA_Matrix::make( this, se ) ;
   se->destroy() ;

   if( OBJS.count( NAME ) != 0 ) PDE_GeometricMultilevel_PC_ERROR::n0( NAME ) ;
   OBJS[ NAME ] = this ;
}

//----------------------------------------------------------------------
PDE_GeometricMultilevel_PC:: ~PDE_GeometricMultilevel_PC( void )
//----------------------------------------------------------------------
{
   OBJS.erase( NAME ) ;
}

//----------------------------------------------------------------------
PDE_GeometricMultilevel_PC:: PDE_GeometricMultilevel_PC(
                                 PEL_Object* a_owner,
                                 PDE_GeometricMultilevel_PC const* other )
//----------------------------------------------------------------------
   : LA_Preconditioner( a_owner )
{
}

//----------------------------------------------------------------------
PDE_GeometricMultilevel_PC*
PDE_GeometricMultilevel_PC:: object( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: object" ) ;
   std::map<std::string,PDE_GeometricMultilevel_PC*>::const_iterator it =
                                                         OBJS.find( a_name ) ;
   if( it == OBJS.end() ) PDE_GeometricMultilevel_PC_ERROR::n1( a_name ) ;

   PDE_GeometricMultilevel_PC* result = (*it).second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PDE_GeometricMultilevel_PC:: name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: name" ) ;

   return( NAME ) ;
}

//----------------------------------------------------------------------
bool
PDE_GeometricMultilevel_PC:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: is_valid" ) ;
   return( BUILD_OK ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GeometricMultilevel_PC:: dimension( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: dimension" ) ;
   PEL_CHECK_PRE( dimension_PRE() ) ;

   size_t result = FINEST_A->nb_rows() ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: set_discretization_scene(
                                 PDE_DomainAndFields const* dom,
                                 PDE_SystemNumbering const* nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: set_discretization_scene" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE(
    FORALL( ( size_t i=0 ; i<nmb->nb_links() ; ++i ),
            nmb->link(i) !=0  &&
            dom->set_of_discrete_fields()->has( nmb->link(i)->field()->name() ) ) ) ;

   COAR = dom->adapter_CHARMS()->algebraic_coarsener() ;
   NMB = nmb ;

   PEL_CHECK_POST( algebraic_coarsener() ==
                   dom->adapter_CHARMS()->algebraic_coarsener() ) ;
}


//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: build( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( mat ) ) ;

   if( algebraic_coarsener() == 0 || mat->nb_rows() != nb_fine_unknowns() )
      PDE_GeometricMultilevel_PC_ERROR::n3( type_name() ) ;

   //??? ne le faire que si le discretization_adapter n'a pas changé
   get_prolongation_matrices( algebraic_coarsener() ) ;

   FINEST_A = mat ;

   LA_Matrix* dummy_0 = MAT_PROTO->create_matrix( 0 ) ;
   LA_Matrix* dummy_1 = MAT_PROTO->create_matrix( 0 ) ;
   for( size_t level=NB_LEVELS-1 ; level!=0 ; --level )
   {
      LA_Matrix const* fine_A =
         ( level==(NB_LEVELS-1) ? FINEST_A :
           static_cast< LA_Matrix* >( AA->at( level ) ) ) ;
      PEL_ASSERT( fine_A->nb_rows() == fine_A->nb_cols() ) ;
      size_t nb_fine = fine_A->nb_rows() ;
      size_t nb_local_fine = fine_A->nb_local_rows() ;
      PEL_ASSERT( fine_A->row_distribution()->global_number() == nb_fine ) ;
      PEL_ASSERT( fine_A->col_distribution()->global_number() == nb_fine ) ;
      PEL_ASSERT( fine_A->nb_local_rows() == fine_A->nb_local_cols() ) ;

      LA_Matrix* coar_A =
           static_cast< LA_Matrix* >( AA->at( level-1 ) ) ;
      PEL_ASSERT( coar_A->nb_rows() == coar_A->nb_cols() ) ;
      size_t nb_coar = coar_A->nb_rows() ;
      size_t nb_local_coar = coar_A->nb_local_rows() ;
      PEL_ASSERT( coar_A->row_distribution()->global_number() == nb_coar ) ;
      PEL_ASSERT( coar_A->col_distribution()->global_number() == nb_coar ) ;
      PEL_ASSERT( coar_A->nb_local_rows() ==  coar_A->nb_local_cols() ) ;

      LA_Matrix* pr =
         static_cast< LA_Matrix* >( COARSE_TO_FINE->at( level-1 ) ) ;
      PEL_ASSERT( pr->nb_rows() == nb_fine ) ;
      PEL_ASSERT( pr->nb_cols() == nb_coar ) ;
      PEL_ASSERT( pr->row_distribution()->global_number() == nb_fine ) ;
      PEL_ASSERT( pr->col_distribution()->global_number() == nb_coar ) ;
      PEL_ASSERT( pr->nb_local_rows() == nb_local_fine ) ;
      PEL_ASSERT( pr->nb_local_cols()== nb_local_coar ) ;

      dummy_1->re_initialize( nb_coar, nb_fine,
                              nb_local_coar, nb_local_fine ) ;
      dummy_1->add_tMat( pr ) ;
      dummy_1->synchronize() ;

      dummy_0->re_initialize( nb_coar, nb_fine,
                              nb_local_coar, nb_local_fine ) ;
      dummy_0->add_Mat_Mat( dummy_1, fine_A ) ;
      dummy_0->synchronize() ;

      //dummy->add_tMat_Mat( pr, fine_A ) ;   // dummy = tr(pr) * fine_A
      coar_A->add_Mat_Mat( dummy_0, pr ) ;    // coar_A = dummy * Pr
      coar_A->synchronize() ;
   }
   dummy_0->destroy() ;
   dummy_1->destroy() ;

   BUILD_OK = true  ;

   PEL_CHECK_POST( build_POST( mat ) ) ;
   PEL_CHECK_POST( is_valid() ) ;
   PEL_CHECK_POST( finest_mat() == mat ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: unbuild( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: unbuild" ) ;
   PEL_CHECK_PRE( unbuild_PRE() ) ;
   BUILD_OK = false ;
   PEL_CHECK_POST( unbuild_POST() ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering const*
PDE_GeometricMultilevel_PC:: system_numbering( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: system_numbering" ) ;

   PDE_SystemNumbering const* result = NMB ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GeometricMultilevel_PC:: nb_fine_unknowns( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: nb_fine_unknowns" ) ;

   size_t result = 0 ;
   if( NMB != 0 )
   {
      result = NMB->nb_global_unknowns() ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GeometricMultilevel_PC:: nb_local_fine_unknowns( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: nb_local_fine_unknowns" ) ;

   size_t result = 0 ;
   if( NMB != 0 )
   {
      result = NMB->nb_unknowns_on_current_process() ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_AlgebraicCoarsener*
PDE_GeometricMultilevel_PC:: algebraic_coarsener( void ) const
//----------------------------------------------------------------------
{
   return( COAR ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: smooth_GaussSeidel( size_t nb_steps,
                                                 LA_Matrix const* mat,
                                                 LA_Vector const* rhs,
                                                 LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: smooth_GaussSeidel" ) ;
   PEL_CHECK_PRE( rhs->is_synchronized() ) ;
   PEL_CHECK_PRE( mat->is_synchronized() ) ;

   //A faire avec le parallelisme : voir l'autre GaussSeidel ci-dessous
   
   LA_MatrixIterator* it = mat->create_stored_item_iterator( 0 ) ;
   PEL_ASSERT( it->nb_rows() == mat->nb_rows() ) ;

   for( size_t is=0 ; is<nb_steps ; ++is )
   {
      for( size_t i=0 ; i<mat->nb_rows() ; ++i )
      {
         if( i < mat->row_distribution()->local_index_limit() &&
             i >= mat->row_distribution()->first_local_index() )
         {
            double xx = 0.0 ;
            double diag_term = PEL::bad_double() ;
            it->start_row_items( i ) ;
            for( ; it->is_valid() ; it->go_next() )
            {
               size_t j = it->col() ;
               if( j < mat->row_distribution()->local_index_limit() &&
                   j >= mat->row_distribution()->first_local_index() )
               {
                  xx += it->item() * sol->item( j ) ;//?? C'est pas bon
                  if( i == j ) diag_term = it->item() ;
               }
            }
            double rr = rhs->item( i ) - xx ;
            sol->add_to_item( i, rr/diag_term ) ;
         }
         sol->synchronize() ;//Ca fait beaucoup de synchronisation
      }
   }
   it->destroy() ; it=0 ;

   PEL_CHECK_POST( sol->is_synchronized() ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: extra_columns(
                                    LA_Matrix const* mat,
                                    intVector& ext_cols )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: extra_columns" ) ;
   PEL_CHECK_PRE( mat->row_distribution()->is_compatible(
                                              mat->col_distribution() ) ) ;
   PEL_CHECK_PRE( mat->nb_rows() == mat->nb_cols() ) ;

   size_t_vector receive_index( 0 ) ;

   LA_MatrixIterator* it = mat->create_stored_item_iterator( 0 ) ;

   size_t first = mat->row_distribution()->first_local_index() ;
   size_t last  = mat->row_distribution()->local_index_limit() ;

   boolVector is_ext( mat->nb_rows(), false ) ;
   size_t nb_ext = 0 ;

   for( size_t i=first ; i<last ; ++i )
   {
      it->start_row_items( i ) ;
      for( ; it->is_valid() ; it->go_next() )
      {
         size_t j = it->col() ;
         if( j<first || j>=last )
         {
            if( !is_ext( j ) )
            {
               nb_ext++ ;
               is_ext( j ) = true ;
            }
         }
      }
   }
   it->destroy() ; it=0 ;

   ext_cols.resize( nb_ext ) ;
   size_t idx = 0 ;
   for( size_t j=0 ; j<is_ext.size() ; j++ )
   {
      if( is_ext( j ) )
      {
         ext_cols( idx ) = (int) j ;
         ++idx ;
      }
   }
   PEL_ASSERT( idx == nb_ext ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: priority_rows(
                                    PEL_DistributedPartition const* row_dist,
                                    intVector const& rows_to_send,
                                    size_t_vector& priority_rows )
//----------------------------------------------------------------------
{
   PEL_Communicator const* com = PEL_Exec::communicator() ;

   intVector nb_to_send_to_proc( com->nb_ranks() ) ;
   intVector nb_priority_from_proc( com->nb_ranks() ) ;

   size_t rank = com->rank() ;

   for( size_t i=0 ; i<rows_to_send.size() ; i++ )
   {
      size_t idx = (size_t) rows_to_send( i ) ;
      size_t row_owner = row_dist->rank_of( idx ) ;
      nb_to_send_to_proc( row_owner )++ ;
   }

   com->all_to_all( nb_to_send_to_proc, nb_priority_from_proc ) ;
   nb_priority_from_proc( rank ) = 0 ;

   size_t nb_total_priority = nb_priority_from_proc.sum() ;

   priority_rows.resize( nb_total_priority ) ;

   //Communications

   size_t nb_sent = 0 ;
   size_t nb_received = 0 ;

   for( size_t r=0 ; r<com->rank() ; r++ )
   {
      size_t nb_from_r = (size_t) nb_priority_from_proc( r ) ;

      if( nb_from_r > 0 )
      {
         intVector received( nb_from_r ) ;
         com->receive( r, const_cast<int*>(received.data()), nb_from_r ) ;
         for( size_t j=0 ; j<nb_from_r ; j++ )
         {
            priority_rows( nb_received + j ) = (size_t) received( j ) ;
         }
         nb_received += nb_from_r ;
      }
   }

   for( size_t r=0 ; r<com->nb_ranks() ; r++ )
   {
      size_t nb_to_r = nb_to_send_to_proc( r ) ;

      if( nb_to_r > 0 )
      {
         if( r != com->rank() )
         {
            intVector sent( nb_to_r ) ;
            for( size_t j=0 ; j<nb_to_r ; j++ )
            {
               sent( j ) = rows_to_send( j + nb_sent ) ;
            }
            com->send( r, sent.data(), nb_to_r ) ;
         }
         nb_sent += nb_to_r ;
      }
   }

   for( size_t r=com->rank()+1 ; r<com->nb_ranks() ; r++ )
   {
      size_t nb_from_r = nb_priority_from_proc( r ) ;
      if( nb_from_r > 0 )
      {
         intVector received( nb_from_r ) ;
         com->receive( r, const_cast<int*>(received.data()), nb_from_r ) ;
         for( size_t j=0 ; j<nb_from_r ; j++ )
         {
            priority_rows( nb_received + j ) =  (size_t) received( j ) ;
         }
         nb_received += nb_from_r ;
      }
   }
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<priority_rows.size() ; ++i ),
                        priority_rows( i ) >= row_dist->first_local_index() &&
                        priority_rows( i ) < row_dist->local_index_limit() ) ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: smooth_GaussSeidel(
                                    size_t nb_steps,
                                    LA_Vector const* to_be_smoothed,
                                    LA_Matrix const* mat,
                                    LA_Vector const* rhs,
                                    LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: smooth_GaussSeidel" ) ;
   PEL_CHECK_PRE( to_be_smoothed != 0 ) ;
   PEL_CHECK_PRE( mat != 0 ) ;
   PEL_CHECK_PRE( rhs != 0 ) ;
   PEL_CHECK_PRE( sol != 0 ) ;
   PEL_CHECK_PRE( rhs->is_synchronized() ) ;
   PEL_CHECK_PRE( to_be_smoothed->is_synchronized() ) ;
   PEL_CHECK_PRE( mat->is_synchronized() ) ;
   PEL_CHECK_PRE( mat->distribution_strategy() ==
                  rhs->distribution_strategy() ) ;
   PEL_CHECK_PRE( mat->distribution_strategy() ==
                  sol->distribution_strategy() ) ;
//   PEL_CHECK_PRE( EQUIVALENT( mat->is_distributed(),
//                              rhs->is_distributed() ) ) ;
//   PEL_CHECK_PRE( EQUIVALENT( mat->is_distributed(),
//                              sol->is_distributed() ) ) ;
   PEL_CHECK_PRE( IMPLIES( mat->distribution_strategy() != LA::NoDistribution,
                  sol->row_distribution()->is_compatible(
                                              mat->col_distribution() ) ) ) ;
   PEL_CHECK_PRE( IMPLIES( mat->distribution_strategy() != LA::NoDistribution,
                  rhs->row_distribution()->is_compatible(
                                              mat->row_distribution() ) ) ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   boolVector is_smoothed( mat->nb_rows(), false ) ;

   //A attacher à la matrice
   LA_MatrixIterator* it = mat->create_stored_item_iterator( 0 ) ;
   PEL_ASSERT( it->nb_rows() == mat->nb_rows() ) ;

   LA_SeqVector* glo_vec = LA_SeqVector::create( 0, sol->nb_rows() ) ;

   //Ca pourrait etre attache à la matrice
   //et fait une fois pour toute
   size_t_vector  priority( 0 ) ;
   intVector ext_cols( 0 ) ;
   if( mat->distribution_strategy() != LA::NoDistribution )
   {
      extra_columns( mat, ext_cols ) ;
      priority_rows( mat->row_distribution(), ext_cols, priority ) ;
   }

   LA_Scatter* scatter = sol->create_scatter( 0, ext_cols, ext_cols ) ;

   PEL_ASSERT( mat->row_distribution()->is_compatible(
                                                  scatter->distribution() ) ) ;
   size_t first = mat->row_distribution()->first_local_index() ;
   size_t last  = mat->row_distribution()->local_index_limit() ;

   sol->start_local_modifs() ;
   for( size_t is=0 ; is<nb_steps ; ++is )
   {
      is_smoothed.set( false ) ;
      for( size_t r=0 ; r<com->nb_ranks() ; ++r )
      {
         scatter->get( sol, glo_vec ) ;
         if( r == com->rank() )
         {
            for( size_t i=0 ; i<priority.size() ; ++i )
            {
               size_t i_row = priority( i ) ;
               PEL_ASSERT(  i_row >= first && i_row < last ) ;
               if( to_be_smoothed->item( i_row ) > 0.9 )
               {
                  is_smoothed( i_row ) = true ;
                  double xx = 0.0 ;
                  double diag_term = PEL::bad_double() ;
                  it->start_row_items( i_row ) ;
                  for( ; it->is_valid() ; it->go_next() )
                  {
                     size_t j = it->col() ;
                     double sol_j = ( j >= first && j < last ) ?
                           sol->item( j ) : glo_vec->item( j ) ;
                           xx += it->item() * sol_j ;
                           if( i_row == j ) diag_term = it->item() ;
                  }
                  double rr = rhs->item( i_row ) - xx ;
                  sol->add_to_item( i_row, rr/diag_term ) ;
               }
            }
         }
      }
      scatter->get( sol, glo_vec ) ;
      for( size_t i_row=first ; i_row<last ; ++i_row )
      {
         if( !is_smoothed( i_row ) && to_be_smoothed->item( i_row ) > 0.9 )
         {
            is_smoothed( i_row ) = true ;
            double xx = 0.0 ;
            double diag_term = PEL::bad_double() ;
            it->start_row_items( i_row ) ;
            for( ; it->is_valid() ; it->go_next() )
            {
               size_t j = it->col() ;
               double sol_j = ( j >= first && j < last ) ?
                     sol->item( j ) : glo_vec->item( j ) ;
                     xx += it->item() * sol_j ;
                     if( i_row == j ) diag_term = it->item() ;
            }
            double rr = rhs->item( i_row ) - xx ;
            sol->add_to_item( i_row, rr/diag_term ) ;
         }
      }
   }
   sol->stop_local_modifs() ;
   it->destroy() ; it=0 ;
   glo_vec->destroy() ; glo_vec = 0 ;
   scatter->destroy() ; scatter = 0 ;

   PEL_CHECK_POST( sol->is_synchronized() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GeometricMultilevel_PC:: nb_levels( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: nb_levels" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   size_t result = NB_LEVELS ;

   return( result ) ;
}

//----------------------------------------------------------------------
LA_Matrix*
PDE_GeometricMultilevel_PC:: mat_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: mat_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels()-1 ) ;

   LA_Matrix* result = static_cast< LA_Matrix* >( AA->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Matrix const*
PDE_GeometricMultilevel_PC:: finest_mat( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: finest_mat" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   LA_Matrix const* result = FINEST_A ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Vector*
PDE_GeometricMultilevel_PC:: rhs_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: rhs_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels()-1 ) ;

   LA_Vector* result = static_cast< LA_Vector* >( B->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Vector*
PDE_GeometricMultilevel_PC:: res_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: res_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels() ) ;

   LA_Vector* result = static_cast< LA_Vector* >( RES->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Vector*
PDE_GeometricMultilevel_PC:: sol_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: sol_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels()-1 ) ;

   LA_Vector* result = static_cast< LA_Vector* >( X->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Matrix const*
PDE_GeometricMultilevel_PC:: coarse_to_fine( size_t coarse_level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: coarse_to_fine" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( coarse_level < nb_levels()-1 ) ;

   LA_Matrix const* result =
      static_cast< LA_Matrix* >( COARSE_TO_FINE->at( coarse_level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Vector const*
PDE_GeometricMultilevel_PC:: smoothing_lines_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: smoothing_lines_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels()-1 ) ;

   LA_Vector const* result =
                    static_cast< LA_Vector* >( SMOO_LINE->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Vector const*
PDE_GeometricMultilevel_PC:: unknowns_of_level( size_t level ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: unknowns_of_level" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( level < nb_levels()-1 ) ;

   LA_Vector const* result =
                    static_cast< LA_Vector* >( UNKNOWNS_LEVEL->at( level ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: compute_residual( LA_Matrix const* mat,
                                               LA_Vector const* rhs,
                                               LA_Vector const* sol,
                                               LA_Vector* residual )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: compute_residual" ) ;
   PEL_CHECK_PRE( mat->is_synchronized() ) ;
   PEL_CHECK_PRE( rhs->is_synchronized() ) ;
   PEL_CHECK_PRE( sol->is_synchronized() ) ;
   PEL_CHECK_PRE( mat->distribution_strategy() == 
                  rhs->distribution_strategy() ) ;
   PEL_CHECK_PRE( mat->distribution_strategy() == 
                  sol->distribution_strategy() ) ;
   PEL_CHECK_PRE( mat->distribution_strategy() == 
                  residual->distribution_strategy() ) ;
//   PEL_CHECK_PRE( EQUIVALENT( mat->is_distributed(),
//                              rhs->is_distributed() ) ) ;
//   PEL_CHECK_PRE( EQUIVALENT( mat->is_distributed(),
//                              sol->is_distributed() ) ) ;
//   PEL_CHECK_PRE( EQUIVALENT( mat->is_distributed(),
//                              residual->is_distributed() ) ) ;
   PEL_CHECK_PRE( IMPLIES( mat->distribution_strategy() != LA::NoDistribution,
                  sol->row_distribution()->is_compatible(
                                              mat->col_distribution() ) ) ) ;
   PEL_CHECK_PRE( IMPLIES( mat->distribution_strategy() != LA::NoDistribution,
                  rhs->row_distribution()->is_compatible(
                                              mat->row_distribution() ) ) ) ;
   PEL_CHECK_PRE( IMPLIES( mat->distribution_strategy() != LA::NoDistribution,
                  residual->row_distribution()->is_compatible(
                                              mat->row_distribution() ) ) ) ;

   // residual = rhs
   // residual = -mat*rhs + residual
   residual->set( rhs ) ;
   mat->multiply_vec_then_add( sol, residual, -1.0, 1.0 ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: print_residuals(
                                      std::string const& indent, size_t n,
                                      double norm_res,
                                      double norm_res_0 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: print_residuals" ) ;
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = PEL::out().precision() ;
   PEL::out() << setprecision( 6 ) ;

   if( n==1 )
   {
      PEL::out() << indent
           << setw( 3+15 ) << "|r|_2"
           << setw( 15 )   << "|r|_2/|r0|_2" ;
      PEL::out() << endl ;
   }
   PEL::out() << indent
        << setw( 3 ) << n
        << setw( 15 ) << norm_res
        << setw( 15 ) << norm_res/norm_res_0 ;
   PEL::out() << endl ;

   PEL::out() << std::setprecision(p) ;
   PEL::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC:: get_prolongation_matrices(
                                              PDE_AlgebraicCoarsener* coar )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricMultilevel_PC:: get_prolongation_matrices" ) ;

   coar->prepare_for_coarsening( system_numbering() ) ;

   NB_LEVELS = coar->nb_levels() ;
   if( NB_LEVELS == 1 ) PDE_GeometricMultilevel_PC_ERROR::n2() ;

   size_t old_size = COARSE_TO_FINE->index_limit() ;
   size_t new_size = NB_LEVELS - 1 ;
   if( new_size > old_size )
   {
      COARSE_TO_FINE->resize( new_size ) ;
      SMOO_LINE->resize( new_size ) ;
      UNKNOWNS_LEVEL->resize( new_size ) ;
      AA->resize( new_size ) ;
      B->resize( new_size ) ;
      X->resize( new_size ) ;
      RES->resize( NB_LEVELS ) ;
      RES->set_at( NB_LEVELS-1, MAT_PROTO->create_vector( this ) ) ;
   }
   for( size_t i=old_size ; i<new_size ; ++i )
   {
      COARSE_TO_FINE->set_at( i, MAT_PROTO->create_matrix( this ) ) ;
      SMOO_LINE->set_at( i, MAT_PROTO->create_vector( this ) ) ;
      UNKNOWNS_LEVEL->set_at( i, MAT_PROTO->create_vector( this ) ) ;
      AA->set_at( i, MAT_PROTO->create_matrix( this ) ) ;
      B->set_at( i, MAT_PROTO->create_vector( this ) ) ;
      X->set_at( i, MAT_PROTO->create_vector( this ) ) ;
      RES->set_at( i, MAT_PROTO->create_vector( this ) ) ;
   }

   //??? si on ne change que la matrice, mais pas les niveaux,
   //??? on va exécuter cette boucle inutilement
   size_t coarse_l = PEL::bad_index() ;
   while( coar->coarsening_is_possible() )
   {
      coar->do_one_coarsening() ;
      size_t fine_l = coar->current_fine_level() ;
      coarse_l = fine_l - 1 ;

      LA_Matrix* mm =
                 static_cast< LA_Matrix* >( COARSE_TO_FINE->at( coarse_l ) )  ;
      coar->build_current_prolongation_matrix( mm ) ;
      size_t nb_coar = mm->nb_cols() ;
      size_t nb_local_coar = mm->nb_local_cols() ;

      LA_Vector* vv = static_cast< LA_Vector* >( SMOO_LINE->at( coarse_l ) ) ;
      coar->build_smoothing_lines( vv ) ;

      vv= static_cast< LA_Vector* >( UNKNOWNS_LEVEL->at( coarse_l ) ) ;
      coar->build_current_level_unknowns( vv ) ;

      mm = static_cast< LA_Matrix* >( AA->at( coarse_l ) ) ;
      mm->re_initialize( nb_coar, nb_coar, nb_local_coar, nb_local_coar ) ;

      vv = static_cast< LA_Vector* >( B->at( coarse_l ) ) ;
      vv->re_initialize( nb_coar, nb_local_coar ) ;

      vv = static_cast< LA_Vector* >( X->at( coarse_l ) ) ;
      vv->re_initialize( nb_coar, nb_local_coar ) ;

      vv = static_cast< LA_Vector* >( RES->at( coarse_l ) ) ;
      vv->re_initialize( nb_coar, nb_local_coar ) ;
   }
   PEL_ASSERT( coarse_l == 0 ) ;
   LA_Vector* vv = static_cast< LA_Vector* >( RES->at( NB_LEVELS-1 ) ) ;
   vv->re_initialize( nb_fine_unknowns(),
                      NMB->nb_unknowns_on_current_process() ) ;
}

//internal-------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC_ERROR:: n0( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PDE_GeometricMultilevel_PC error:" << endl ;
   mesg << "    attempt to create two instances with the" << endl ;
   mesg << "    same name: \"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC_ERROR:: n1( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PDE_GeometricMultilevel_PC error:" << endl ;
   mesg << "    there is no recorded instance of name: \"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC_ERROR:: n2( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PDE_GeometricMultilevel_PC error:" << endl ;
   mesg << "    a multilevel hierarchy is requested" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
PDE_GeometricMultilevel_PC_ERROR:: n3( std::string const& tname )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** PDE_GeometricMultilevel_PC error:" << endl << endl ;
   mesg << "    invalid use of an object of type \"" << tname << "\"" << endl ;
   mesg << "    the member function" << endl ;
   mesg << "       \"set_discretization_scene\"" << endl ;
   mesg << "    may not have been called properly for that object" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

