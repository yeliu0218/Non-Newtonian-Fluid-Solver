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

#ifndef LA_DENSE_MATRIX_HH
#define LA_DENSE_MATRIX_HH

#include <LA_SeqMatrix.hh>
#include <LA_DenseMatrixIterator.hh>

class LA_SeqVector ;
class LA_SymmetricMatrix ;
class doubleArray2D ;
class size_t_vector ;

/*
Matrices of double, without commitment to a particular structure.

Data deck instanciation :

   MODULE LA_Matrix
      concrete_name = "LA_DenseMatrix"
   END MODULE LA_Matrix

PUBLISHED
*/

class PEL_EXPORT LA_DenseMatrix : public LA_SeqMatrix
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static LA_DenseMatrix* create( PEL_Object* a_owner,
                                     size_t a_nb_rows,
                                     size_t a_nb_cols ) ;

      // Create and return an instance.
      static LA_DenseMatrix* create( PEL_Object* a_owner,
                                     doubleArray2D const& amat ) ;

      virtual LA_DenseMatrix* create_matrix( PEL_Object* a_owner ) const ;

   //-- Characteristics

      virtual size_t allocated_memory( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

   //-- Access

      virtual LA_MatrixIterator* create_stored_item_iterator(
                                          PEL_Object* a_owner ) const ;

      virtual size_t nb_stored_items( void ) const ;

      virtual size_t nb_rows( void ) const ;

      virtual size_t nb_cols( void ) const ;

      virtual double item( size_t i, size_t j ) const ;

   //-- Element change

      virtual void set_stored_items( double val ) ;

      virtual void set_item( size_t i, size_t j, double x ) ;

      virtual void add_to_item( size_t i, size_t j, double x ) ;

      // Add `alpha'*`column' to `beta'*`self' `j'-th column.
      void add_column( size_t j,
                       LA_SeqVector const* column,
                       double alpha = 1.0,
                       double beta = 0.0  ) ;

   //-- BLAS level 3 : matrix-matrix operators

      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

      virtual void add_Mat( LA_Matrix const* A,
                            double alpha = 1.0,
                            bool same_pattern=false ) ;

      virtual void add_tMat( LA_Matrix const* A, double alpha = 1.0 ) ;

      virtual void add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                double alpha = 1.0 ) ;

      virtual void add_tMat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha = 1.0 ) ;

      virtual void add_Mat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                                 double alpha = 1.0 ) ;

      virtual void add_tMat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                                  double alpha = 1.0 ) ;

      // Replaces self by its inverse and returns its determinant.
      // Reference: Philip B Bevington,
      //             "Data Reduction and Error Analysis for the Physical
      //             Sciences", McGraw-Hill, New York, 1969, pp. 300-303.
      //             F. Murtagh, ST-ECF, Garching-bei-Muenchen,
      //             January 1986.
      void invert( double& det ) ;

      // determinant of `self'
      double determinant( void ) const ;

   //-- Statistical calculus on columns vectors(140.)

      // mean value of each column of `self'
      LA_SeqVector* col_mean( PEL_Object* a_owner ) const ;

      // standard deviation of all columns of  `self'
      LA_SeqVector* col_standard_deviation( PEL_Object* a_owner ) const ;

      // variance-covariance matrix associated to `self'
      LA_SymmetricMatrix* col_variance( PEL_Object* a_owner ) const ;

      // correlation matrix associated to column vectors of `self'
      LA_SymmetricMatrix* col_correlation( PEL_Object* a_owner ) const ;

      // reduced matrix associated to `self'
      LA_DenseMatrix*  col_reduced( PEL_Object* a_owner ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_DenseMatrix( void ) ;
      LA_DenseMatrix( LA_DenseMatrix const& other ) ;
      LA_DenseMatrix& operator=( LA_DenseMatrix const& other ) ;

      LA_DenseMatrix( PEL_Object* a_owner,
                      size_t a_nb_rows, size_t a_nb_cols ) ;

      bool largelem( double& mxel,
                     size_t k,
                     size_t_vector& ik, size_t_vector& jk ) ;

      virtual void re_initialize_with_global_sizes( size_t a_nb_rows,
                                                    size_t a_nb_cols ) ;

   //-- Plug in

      LA_DenseMatrix( void ) ;

      virtual LA_DenseMatrix* create_replica(
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const ;

   //-- BLAS level 3 : matrix-matrix operators

      void add_Mat_IMP( LA_DenseMatrix const* A, double alpha ) ;

      void add_tMat_IMP( LA_DenseMatrix const* A, double alpha ) ;

      void add_Mat_Mat_IMP( LA_DenseMatrix const* A,
                            LA_DenseMatrix const* B,
                            double alpha ) ;

      void add_tMat_Mat_IMP( LA_DenseMatrix const* A,
                             LA_DenseMatrix const* B,
                             double alpha ) ;

      void add_Mat_tMat_IMP( LA_DenseMatrix const* A,
                             LA_DenseMatrix const* B,
                             double alpha ) ;

      void add_tMat_tMat_IMP( LA_DenseMatrix const* A,
                              LA_DenseMatrix const* B,
                              double alpha ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Friends

      friend class LA_DenseMatrixIterator ;

   //-- Class attributes

      static LA_DenseMatrix* PROTOTYPE ;

   //-- Attributes

      size_t NB_ROWS ;
      size_t NB_COLS ;
      size_t SIZE ;
      double** MAT ;
} ;

#endif
