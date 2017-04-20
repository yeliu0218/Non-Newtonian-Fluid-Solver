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

#include <PDE.hh>

#include <LA_MatrixIterator.hh>
#include <LA_PelMatrix.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DOFconstraintsIterator.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_SystemNumbering.hh>

#include <iostream>

using std::endl ;
using std::string ;

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_vector_0( LA_Matrix* matrix,
                                   LA_Vector* vector,
                                   PDE_LocalEquation const* leq,
                                   PDE_SystemNumbering const* r_nmb,
                                   size_t i_link,
                                   PDE_SystemNumbering const* c_nmb,
                                   size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_vector_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( r_nmb != 0 ) ;
   PEL_CHECK_PRE( c_nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == c_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( i_link < r_nmb->nb_links() ) ;
   PEL_CHECK_PRE( j_link < c_nmb->nb_links() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       r_nmb->link( i_link )->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       c_nmb->link( j_link )->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
   
   PDE_LinkDOF2Unknown const* r_lnk = r_nmb->link( i_link ) ;
   size_t_vector const& r_cmps = r_lnk->components_table() ;
   PDE_DiscreteField const* r_ff = r_lnk->field() ;
   
   for( size_t il=0 ; il<leq->nb_rows() ; ++il )
   {
      size_t nn = leq->row_node( il ) ;
      for( size_t i=0 ; i<r_cmps.size() ; ++i )
      {
         size_t ic = r_cmps( i ) ;
         if( r_lnk->DOF_is_unknown( nn, ic ) )
         {
            size_t ig = r_nmb->global_unknown_for_DOF( nn, ic, i_link ) ;

            add_line( il, i, ig, matrix, vector, leq, 1.0, c_nmb, j_link ) ;
         }
         else if( r_ff->DOF_is_constrained( nn, ic ) )
         {
            PDE_DOFconstraintsIterator* ct = 
                              r_ff->create_constraints_iterator( 0 ) ;
            for( ct->start( nn, ic ) ; ct->is_valid() ; ct->go_next() )
            {
               size_t ig = r_nmb->global_unknown_for_DOF( 
                              ct->node_of_constraining_DOF(),
                              ct->component_of_constraining_DOF(),
                              i_link ) ;
               double coef = ct->constraint_coefficient() ;
               add_line( il, i, ig, matrix, vector, leq, coef, c_nmb, j_link ) ;
            }
            ct->destroy() ; //??? each time ???????????????
         }
      }
   }

   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
void
PDE:: add_line( size_t il, size_t i, size_t ig,
                LA_Matrix* matrix,
                LA_Vector* vector,
                PDE_LocalEquation const* leq,
                double coef,
                PDE_SystemNumbering const* c_nmb,
                size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: add_line" ) ;
   
   PDE_LinkDOF2Unknown const* c_lnk = c_nmb->link( j_link ) ;
   size_t_vector const& c_cmps = c_lnk->components_table() ;
   PDE_DiscreteField const* c_ff = c_lnk->field() ;
   
   vector->add_to_item( ig, coef * leq->vector_item( il, i ) ) ;
   
   for( size_t jl=0 ; jl<leq->nb_columns() ; ++jl )
   {
      size_t mm = leq->column_node( jl ) ;
      for( size_t j=0 ; j<c_cmps.size() ; ++j )
      {
         size_t jc = c_cmps( j ) ;
         if( leq->matrix_item_is_set( il, jl, i, j ) )
         {
            double aij = leq->matrix_item( il, jl, i, j ) * coef ;
            if( c_lnk->DOF_is_unknown( mm, jc ) )
            {
               size_t jg = 
                      c_nmb->global_unknown_for_DOF( mm, jc, j_link ) ;
               matrix->add_to_item( ig, jg, aij ) ;
            }
            else if( c_ff->DOF_is_constrained( mm, jc ) )
            {
               PDE_DOFconstraintsIterator* ct = 
                                 c_ff->create_constraints_iterator( 0 ) ;
               for( ct->start( mm, jc ) ; ct->is_valid() ; ct->go_next() )
               {
                  size_t jg = c_nmb->global_unknown_for_DOF( 
                                 ct->node_of_constraining_DOF(),
                                 ct->component_of_constraining_DOF(), 
                                 j_link ) ;
                  double x = aij * ct->constraint_coefficient() ;
                  matrix->add_to_item( ig, jg, x ) ;                           
               }
               ct->destroy() ; //??? each time ???
            }
            else
            {
               double x = - aij * c_ff->DOF_imposed_value( mm, jc ) ;
               vector->add_to_item( ig, x ) ;
            }
         }
      }
   }   
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_vector_0( LA_Matrix* matrix,
                                   LA_Vector* vector,
                                   PDE_LocalEquation const* leq,
                                   PDE_SystemNumbering const* r_nmb,
                                   PDE_SystemNumbering const* c_nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_vector_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( r_nmb != 0 ) ;
   PEL_CHECK_PRE( c_nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == c_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( r_nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( c_nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       r_nmb->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       c_nmb->link()->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
   
   assemble_in_matrix_vector_0( matrix, vector, 
                                leq, 
                                r_nmb, 0, c_nmb, 0 ) ;

   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_vector_0( LA_Matrix* matrix,
                                   LA_Vector* vector,
                                   PDE_LocalEquation const* leq,
                                   PDE_SystemNumbering const* nmb,
                                   size_t i_link,
                                   size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_vector_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( i_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( j_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link( i_link )->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       nmb->link( j_link )->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
   
   assemble_in_matrix_vector_0( matrix, vector, 
                                leq, 
                                nmb, i_link, nmb, j_link ) ;

   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_vector_0( LA_Matrix* matrix,
                                   LA_Vector* vector,
                                   PDE_LocalEquation const* leq,
                                   PDE_SystemNumbering const* nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_vector_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       nmb->link()->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
   
   assemble_in_matrix_vector_0( matrix, vector, 
                                leq, 
                                nmb, 0, nmb, 0 ) ;

   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_0( LA_Matrix* matrix,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* r_nmb,
                            size_t i_link, 
                            PDE_SystemNumbering const* c_nmb,
                            size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( r_nmb != 0 ) ;
   PEL_CHECK_PRE( c_nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == c_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( i_link < r_nmb->nb_links() ) ;
   PEL_CHECK_PRE( j_link < c_nmb->nb_links() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       r_nmb->link( i_link )->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       c_nmb->link( j_link )->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   
   PDE_LinkDOF2Unknown const* r_lnk = r_nmb->link( i_link ) ;
   PDE_LinkDOF2Unknown const* c_lnk = c_nmb->link( j_link ) ;

   size_t_vector const& r_cmps = r_lnk->components_table() ;
   size_t_vector const& c_cmps = c_lnk->components_table() ;

   for( size_t il=0 ; il<leq->nb_rows() ; ++il )
   {
      size_t nn = leq->row_node( il ) ;
      for( size_t i=0 ; i<r_cmps.size() ; ++i )
      {
         size_t ic = r_cmps( i ) ;
         if( r_lnk->DOF_is_unknown( nn, ic ) )
         {
            size_t ig = r_nmb->global_unknown_for_DOF( nn, ic, i_link ) ;
            for( size_t jl=0 ; jl<leq->nb_columns() ; ++jl )
            {
               size_t mm = leq->column_node( jl ) ;
               for( size_t j=0 ; j<c_cmps.size() ; ++j )
               {
                  size_t jc = c_cmps( j ) ;
                  if( c_lnk->DOF_is_unknown( mm, jc ) &&
                      leq->matrix_item_is_set( il, jl, i, j ) )
                  {
                     double aij = leq->matrix_item( il, jl, i, j ) ;
                     size_t jg = 
                         c_nmb->global_unknown_for_DOF( mm, jc, j_link ) ;
                     matrix->add_to_item( ig, jg, aij ) ;
                  }
               }
            }
         }
      }     
   }
   
   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_0( LA_Matrix* matrix,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* r_nmb,
                            PDE_SystemNumbering const* c_nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( r_nmb != 0 ) ;
   PEL_CHECK_PRE( c_nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == r_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == c_nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( r_nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( c_nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       r_nmb->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       c_nmb->link()->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   
   assemble_in_matrix_0( matrix, leq, r_nmb, 0, c_nmb, 0 ) ;
   
   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_0( LA_Matrix* matrix,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* nmb,
                            size_t i_link, 
                            size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( i_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( j_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link( i_link )->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       nmb->link( j_link )->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   
   assemble_in_matrix_0( matrix, leq, nmb, i_link, nmb, j_link ) ;
   
   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_matrix_0( LA_Matrix* matrix,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_matrix_0" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link()->components_table().size() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                       nmb->link()->components_table().size() ) ;
   PEL_SAVEOLD( double, matrix_state, matrix->state() ) ;
   
   assemble_in_matrix_0( matrix, leq, nmb, 0, nmb, 0 ) ;
      
   PEL_CHECK_POST( matrix->state() == OLD( matrix_state ) ||
                   matrix->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_vector_1( LA_Vector* vector,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* nmb,
                            size_t i_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_vector_1" ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( vector->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( i_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link( i_link )->components_table().size() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
     
   PDE_LinkDOF2Unknown const* r_lnk = nmb->link( i_link ) ;
   
   size_t_vector const& r_cmps = r_lnk->components_table() ;
   
   for( size_t il=0 ; il<leq->nb_rows() ; ++il )
   {
      size_t nn = leq->row_node( il ) ;
      for( size_t i=0 ; i<r_cmps.size() ; ++i )
      {
         size_t ic = r_cmps( i ) ;
         if( r_lnk->DOF_is_unknown( nn, ic ) )
         {
            size_t ig = nmb->global_unknown_for_DOF( nn, ic, i_link ) ;
            vector->add_to_item( ig, leq->vector_item( il, i ) ) ;            
         }
      }     
   }

   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
void
PDE:: assemble_in_vector_1( LA_Vector* vector,
                            PDE_LocalEquation const* leq,
                            PDE_SystemNumbering const* nmb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: assemble_in_vector_1" ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   PEL_CHECK_PRE( vector->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( vector->state() != LA::NotSync_set ) ;
   PEL_CHECK_PRE( nmb->nb_links() == 1 ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() ==
                       nmb->link()->components_table().size() ) ;
   PEL_SAVEOLD( double, vector_state, vector->state() ) ;
     
   assemble_in_vector_1( vector, leq, nmb, 0 ) ;
   
   PEL_CHECK_POST( vector->state() == OLD(vector_state) ||
                   vector->state() == LA::NotSync_add ) ; 
}

//----------------------------------------------------------------------
LA_SeqMatrix*
PDE:: create_extracted_block( PEL_Object* a_owner,
                              LA_Matrix const* matrix,
                              PDE_SystemNumbering const* nmb,
                              size_t i_link,
                              size_t j_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE:: create_extracted_block" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_CHECK_PRE( matrix->is_synchronized() ) ;
   PEL_CHECK_PRE( i_link < nmb->nb_links() ) ;
   PEL_CHECK_PRE( j_link < nmb->nb_links() ) ;
   
   PDE_LinkDOF2Unknown const* r_lnk = nmb->link( i_link ) ;
   PDE_LinkDOF2Unknown const* c_lnk = nmb->link( j_link ) ;
   
   LA_PelMatrix* result = LA_PelMatrix::create( a_owner,
                                                r_lnk->unknown_vector_size(),
                                                c_lnk->unknown_vector_size() ) ;
   
   size_t_vector r_global2local( nmb->nb_global_unknowns() ) ;
   size_t_vector c_global2local( nmb->nb_global_unknowns() ) ;
   r_global2local.set( PEL::bad_index() ) ;
   c_global2local.set( PEL::bad_index() ) ;
   
   {
      PDE_DiscreteField const* r_field = r_lnk->field() ;
      size_t_vector const& r_cmps = r_lnk->components_table() ;
      for( size_t n=0 ; n<r_field->nb_nodes() ; ++n )
      {
         for( size_t i=0 ; i<r_cmps.size() ; ++i )
         {
            size_t ic = r_cmps( i ) ;
            if( r_lnk->DOF_is_unknown( n, ic ) )
            {
               size_t ig = nmb->global_unknown_for_DOF( n, ic, i_link ) ;
               size_t il = r_lnk->unknown_linked_to_DOF( n, ic ) ;
               r_global2local( ig ) = il ;
            }
         }
      }
   }
   
   {
      PDE_DiscreteField const* c_field = c_lnk->field() ;
      size_t_vector const& c_cmps = c_lnk->components_table() ;
      for( size_t n=0 ; n<c_field->nb_nodes() ; ++n )
      {
         for( size_t i=0 ; i<c_cmps.size() ; ++i )
         {
            size_t ic = c_cmps( i ) ;
            if( c_lnk->DOF_is_unknown( n, ic ) )
            {
               size_t ig = nmb->global_unknown_for_DOF( n, ic, j_link ) ;
               size_t il = c_lnk->unknown_linked_to_DOF( n, ic ) ;
               c_global2local( ig ) = il ;
            }
         }
      }
   }
   
   LA_SeqMatrix* lmat = matrix->create_local_matrix( 0 ) ;
   LA_MatrixIterator* it = lmat->create_stored_item_iterator( lmat ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      size_t i = r_global2local( it->row() ) ;
      size_t j = c_global2local( it->col() ) ;
      if( ( i != PEL::bad_index() ) && ( j != PEL::bad_index() ) )
      {
         result->set_item( i, j, it->item() ) ;
      }
   }
   lmat->destroy() ;
   
   result->synchronize() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->is_synchronized() ) ;
   PEL_CHECK_POST( result->nb_rows() == 
                   nmb->link( i_link )->unknown_vector_size() ) ;
   PEL_CHECK_POST( result->nb_cols() == 
                   nmb->link( j_link )->unknown_vector_size() ) ;
   return( result ) ;
}

