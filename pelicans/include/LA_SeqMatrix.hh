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

#ifndef LA_SEQ_MATRIX_HH
#define LA_SEQ_MATRIX_HH

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>


class LA_MatrixIterator ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;


/*
Sequential matrices belonging to PELICANS framework such that
   - storage is not necessarily required for all coefficients ;
   - all the coefficients that are not stored are zero ;
   - there exist a traversal policy that will visit every stored
     coefficient exactly once.

PUBLISHED
*/

class PEL_EXPORT LA_SeqMatrix : public LA_Matrix
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_SeqMatrix* make( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) ;

      // Create a copy of `other' (copying all stored coefficients of `other')
      // of the same type of `self'.
      virtual LA_SeqMatrix* create_copy( PEL_Object* a_owner,
                                         LA_SeqMatrix const* other ) const ;

      virtual void re_initialize(
                        size_t a_nb_rows, size_t a_nb_cols,
                        size_t a_nb_local_rows = PEL::bad_index(),
                        size_t a_nb_local_cols = PEL::bad_index() ) ;

      virtual LA_SeqMatrix* create_matrix( PEL_Object* a_owner ) const  = 0 ;
      virtual LA_SeqVector* create_vector( PEL_Object* a_owner ) const ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

      // memory allocated by `self'
      virtual size_t allocated_memory( void ) const = 0 ;

   //-- Distributed processing

      virtual PEL_DistributedPartition const* row_distribution( void ) const ;

      virtual PEL_DistributedPartition const* col_distribution( void ) const ;

   //-- Access

      virtual LA_MatrixIterator* create_stored_item_iterator(
                                             PEL_Object* a_owner ) const = 0 ;

      // IMPLEMENTATION : copy of `self'.
      virtual LA_SeqMatrix* create_local_matrix( PEL_Object* a_owner ) const ;

      virtual size_t nb_stored_items( void ) const = 0 ;

      virtual double item( size_t i, size_t j ) const = 0 ;

      virtual void extract_diag( LA_Vector* diag ) const ;

      // Extract `i'-th row of `self' in vector `row'.
      virtual void extract_row( size_t i, LA_Vector* row ) const ;

      // Extract `j'-th column of `self' in vector `col'.
      virtual void extract_col( size_t j, LA_Vector* col ) const ;

      // Assigning the sum of the coefficients absolute value for each line
      // of `self' to the corresponding element of `lump'.
      virtual void extract_lump( LA_Vector* lump ) const ;

   //-- Element change

      virtual void nullify( void ) ;

      // Make all coefficients of the `i'-th row vanish.
      virtual void nullify_row( size_t i ) ;

      virtual void set_item( size_t i, size_t j, double x ) = 0 ;

      virtual void add_to_item( size_t i, size_t j, double x ) = 0 ;

      // Assign `val' to all the stored coefficients.
      virtual void set_stored_items( double val ) = 0 ;

      virtual void scale( double alpha ) ;

   //-- BLAS level 2 : matrix-vector operators

      virtual void multiply_vec_then_add(
                           LA_Vector const* x, LA_Vector* y,
                           double alpha = 1.0, double beta = 0.0 ) const ;

      virtual void tr_multiply_vec_then_add(
                           LA_Vector const* x, LA_Vector* y,
                           double alpha = 1.0, double beta = 0.0 ) const ;

      virtual void scale_as_diag_mat_mat( LA_Vector const* lvec ) ;

      virtual void scale_as_mat_diag_mat( LA_Vector const* rvec ) ;

      virtual void add_to_diag( LA_Vector const* vec ) ;

   //-- BLAS level 3 : matrix-matrix operators

      // Reinitialize by copying all stored coefficients of `A'.
      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

      virtual void add_Mat( LA_Matrix const* A,
                            double alpha = 1.0,
                            bool same_pattern=false ) ;

      virtual void add_tMat( LA_Matrix const* A, double alpha = 1.0 ) ;

      virtual void add_Mat_Mat( LA_Matrix const* A,
                                LA_Matrix const* B,
                                double alpha = 1.0 ) ;

      // Multiply the transpose of `A' by `B' and `alpha', then add to `self':
      //   `self' <- `self' + `alpha'*transpose(`A')*`B'.
      virtual void add_tMat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha = 1.0 ) ;

      // Multiply `A' by the transpose of `B' and `alpha', then add to `self':
      //   `self' <- `self' + `alpha'*`A'*transpose(`B').
      virtual void add_Mat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha = 1.0 ) ;

      // Multiply the transpose of `A' by the transpose of `B' and `alpha',
      // then add to `self' :
      //   `self' <- `self' + `alpha'*transpose(`A')*transpose(`B').
      virtual void add_tMat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                                  double alpha = 1.0 ) ;

   //-- System factorization(130.)

      // Factorize `self' with modified ILU 0 algorithm.
      virtual void factorize_MILU0( bool modified, double piv_min ) ;

      /* Solve a linear system of matrix L.U where L is a lower
         triangular matrix with ones on the diagonal and where U is an
         upper triangular matrix. L, appart from its diagonal, and U
         are stored in `self'. */
      virtual void solve_LU( LA_SeqVector const* rhs,
                             LA_SeqVector* sol ) const ;

      enum relaxation_mode { forward, backward, symmetric } ;

      /* Approximate the solution `sol' of the linear system:
                          `self'*`sol'=`rhs'
         performing the relaxation:
           - `mode'=forward or `mode'=symmetric
               `sol'(i) = (1.-`omega')*`sol'(i)+`omega'*r(i) i=1..`::nb_rows'
           - `mode'=backward or `mode'=symmetric
               `sol'(i) = (1.-`omega')*`sol'(i)+`omega'*r(i) i=`::nb_rows'..1
         where r(i) = ( `rhs'(i)-sum(i!=j){`self'(i,j)*`sol'(i)} )/`self'(i,i)
         The convergence is ensured for `omega' in ]0,2[.
         For optimization, the vector:
                   `omega_inv_diag'(i) = `omega'/`self'(i,i)
         is pre-built. */
      virtual void relax( double omega,
                          LA_SeqMatrix::relaxation_mode mode,
                          LA_SeqVector const* omega_inv_diag,
                          LA_SeqVector const* rhs,
                          LA_SeqVector* sol ) const ;

   //-- Input - Output

      virtual void writeMM( std::string const& file ) const ;

      // Save `self' in `mod'.
      virtual void save( PEL_Module* mod ) const ;

      // Restore matrix described by `exp'.
      virtual void restore( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

      virtual ~LA_SeqMatrix( void ) ;

      LA_SeqMatrix( PEL_Object* a_owner, std::string const& a_name ) ;

      virtual void re_initialize_with_global_sizes( size_t a_nb_rows,
                                                    size_t a_nb_cols ) = 0 ;

   //-- Plug in

      // for prototype registration only
      LA_SeqMatrix( std::string const& a_name ) ;

      virtual LA_SeqMatrix* create_replica(
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const = 0 ;

   //-- System factorization

      void raise_MILU0_zero_pivot( size_t i ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool create_copy_PRE( PEL_Object* a_owner,
                                    LA_SeqMatrix const* other ) const ;
      virtual bool create_copy_POST( LA_SeqMatrix* result,
                                     PEL_Object* a_owner,
                                     LA_SeqMatrix const* other ) const ;

      virtual bool item_PRE( size_t i, size_t j ) const ;

      virtual bool extract_row_PRE( size_t i, LA_Vector* row ) const ;
      virtual bool extract_row_POST( size_t i, LA_Vector* row ) const ;

      virtual bool extract_col_PRE( size_t j, LA_Vector* col ) const ;
      virtual bool extract_col_POST( size_t j, LA_Vector* col ) const ;

      virtual bool extract_lump_PRE( LA_Vector* lump ) const ;
      virtual bool extract_lump_POST( LA_Vector* lump ) const ;

      virtual bool nullify_POST( void ) const ;

      virtual bool nullify_row_PRE( size_t i ) const ;
      virtual bool nullify_row_POST( size_t i ) const ;

      virtual bool add_tMat_Mat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                     double alpha ) const ;
      virtual bool add_tMat_Mat_POST( void ) const ;

      virtual bool add_Mat_tMat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                     double alpha ) const ;
      virtual bool add_Mat_tMat_POST( void ) const ;

      virtual bool add_tMat_tMat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                      double alpha ) const ;
      virtual bool add_tMat_tMat_POST( void ) const ;

      virtual bool set_stored_items_POST( double val ) const ;

      virtual bool set_PRE( LA_Matrix const* A, bool same_pattern ) const ;

      virtual bool add_Mat_PRE( LA_Matrix const* A, double alpha,
                                bool same_pattern ) const ;

      virtual bool factorize_MILU0_PRE(
                                    bool modified, double piv_min ) const ;
      virtual bool factorize_MILU0_POST( void ) const ;

      virtual bool solve_LU_PRE(
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const ;
      virtual bool solve_LU_POST(
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const ;

      virtual bool relax_PRE(
                       double omega,
                       LA_SeqMatrix::relaxation_mode mode,
                       LA_SeqVector const* omega_inv_diag,
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const ;
      virtual bool relax_POST(
                       double omega,
                       LA_SeqMatrix::relaxation_mode mode,
                       LA_SeqVector const* omega_inv_diag,
                       LA_SeqVector const* rhs,
                       LA_SeqVector const* sol ) const ;

      virtual bool re_initialize_with_global_sizes_PRE(
                                                    size_t a_nb_rows,
                                                    size_t a_nb_cols ) const ;
      virtual bool re_initialize_with_global_sizes_POST(
                                                    size_t a_nb_rows,
                                                    size_t a_nb_cols ) const ;

      virtual bool save_PRE( PEL_Module const* mod ) const ;
      virtual bool save_POST( void ) const ;

      virtual bool restore_PRE( PEL_ModuleExplorer const* exp ) const ;
      virtual bool restore_POST( void ) const ;

    private: //----------------------------------------------------------

      LA_SeqMatrix( void ) ;
      LA_SeqMatrix( LA_SeqMatrix const& other ) ;
      LA_SeqMatrix& operator=( LA_SeqMatrix const& other ) ;

      void set_IMP( LA_SeqMatrix const* A, bool same_pattern ) ;

      void add_Mat_IMP( LA_SeqMatrix const* A, double alpha,
                        bool same_pattern ) ;

      void add_tMat_IMP( LA_SeqMatrix const* A, double alpha ) ;

      void add_Mat_Mat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                            double alpha ) ;

      void add_tMat_Mat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                             double alpha ) ;

      void add_Mat_tMat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                             double alpha ) ;

      void add_tMat_tMat_IMP( LA_SeqMatrix const* A, LA_SeqMatrix const* B,
                              double alpha ) ;

      bool has_same_pattern( LA_Matrix const* A ) const ;

      static PEL_ObjectRegister* sparse_mat_pluggin_map( void ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Attributes

      // Distributed partition:
      mutable PEL_DistributedPartition* ROW_DIST ;
      mutable PEL_DistributedPartition* COL_DIST ;
} ;

#endif
