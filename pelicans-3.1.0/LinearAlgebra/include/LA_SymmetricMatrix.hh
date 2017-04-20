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

#ifndef LA_SYMMETRIC_MATRIX_HH
#define LA_SYMMETRIC_MATRIX_HH

#include <doubleVector.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SymmetricMatrixIterator.hh>
#include <LA_Vector.hh>

class LA_DenseMatrix ;

/*
Symmetric square matrices of double.

Data deck instanciation :

   MODULE LA_Matrix
      concrete_name = "LA_SymmetricMatrix"
   END MODULE LA_Matrix

PUBLISHED
*/

class PEL_EXPORT LA_SymmetricMatrix : public LA_SeqMatrix
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static LA_SymmetricMatrix* create( PEL_Object* a_owner, size_t dim ) ;

      virtual LA_SymmetricMatrix* create_matrix(
                                           PEL_Object* a_owner ) const ;

   //-- Characteristics

      virtual size_t allocated_memory( void ) const ;

      virtual bool is_symmetric( void ) const ;

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

      virtual void set_item( size_t i, size_t j, double x ) ;

      virtual void add_to_item( size_t i, size_t j, double x ) ;

      virtual void set_stored_items( double val ) ;

      virtual void scale( double alpha ) ;

   //-- BLAS level 3 : matrix-matrix operators

      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

      virtual void add_Mat( LA_Matrix const* A,
                            double alpha = 1.0,
                            bool same_pattern=false ) ;

      // Matrix-Matrix product of the form
      //    M += alpha*A*transpose( A ).
      void add_Mat_tMat( LA_SeqMatrix const* A, double alpha=1.0 ) ;

      // Replaces self by its inverse and returns its determinant.
      // Reference: Philip B Bevington,
      //             "Data Reduction and Error Analysis for the Physical
      //             Sciences", McGraw-Hill, New York, 1969, pp. 300-303.
      //             F. Murtagh, ST-ECF, Garching-bei-Muenchen,
      //             January 1986.
      void invert( double& det ) ;

      // Computes the eigenvalues and the eigenvectors of self.
      void eigen_reduce( size_t matord,
                         LA_SeqVector* egval,
                         LA_DenseMatrix* egvecs,
                         bool& success ) const ;

   //-- Input - Output

      virtual void writeMM( std::string const& file ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_SymmetricMatrix( void ) ;
      LA_SymmetricMatrix( LA_SymmetricMatrix const& other ) ;
      LA_SymmetricMatrix& operator=( LA_SymmetricMatrix const& other ) ;

      LA_SymmetricMatrix( PEL_Object* a_owner, size_t dim ) ;

      size_t index( size_t i, size_t j ) const ;

      virtual void re_initialize_with_global_sizes( size_t a_nb_rows,
                                                    size_t a_nb_cols ) ;

   //-- BLAS level 3 : matrix-matrix operators

      // Reduce a real, symmetric matrix to a symmetric, tridiagonal, matrix.
      //    matord   = order of the matrix ( will always be <= its size ) ;
      //    diagTDMA = vector of dim. matord containing, on output, diagonal
      //               elts. of trid. matrix ;
      //    subdTDMA = working vector of dim. at least matord-1 to contain
      //               subdiagonal elts. ;
      //    orth = matrix of dims. containing, on output, orthogonal
      //           transformation matrix producting the reduction.
      // Normally a call to eigenv will follow the call to reduce in order to
      // produce all eigenvectors and eigenvalues of the matrix .
      //
      // Algorithm used: Martin et al., Num. Math. 11, 181-195, 1968.
      //
      // Reference: Smith et al., Matrix Eigensystem Routines - EISPACK
      // Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 1976,
      // pp. 489-494.

      // F. Murtagh, ST-ECF Garching-bei-Muenchen, January 1986.

      //  Determine eigenvalues and eigenvectors of a symmetric,
      //  tridiagonal matrix.
      // in input :
      //   -  the current object matrix is the orthogonal matrix that
      //   produced the
      //  reduction of a precedent symmetric matrix to the tridiagonal
      //  matrix whose vector of diagonal elems is diag and subdiagonal
      //  vector is subd
      //  matord = order of the matrix ;
      //  diag = vector of dim. matord containing, on output,
      //         eigenvalues ;
      //  subd = working vector of dim. at least matord-1 ;
      //  current matrix  will contain , on output, eigenvectors as
      //  columns ;
      //  returned integer : normally 0, but 1 if no convergence.
      //  Normally the call to eigenv will be preceded by a call to
      //  reduce in order to set up the tridiagonal matrix.
      //  Algorithm used: QL method of Bowdler et al., Num. Math. 11,
      //  293-306, 1968.
      //  Reference: Smith et al., Matrix Eigensystem Routines - EISPACK
      //  Guide, Lecture Notes in Computer Science 6, Springer-Verlag,
      //  1976, pp. 468-474.
      //  F. Murtagh, ST-ECF Garching-bei-Muenchen, January 1986.
      void reduce( size_t matord,
                   LA_SeqVector* diagTDMA,
                   LA_SeqVector* subdTDMA,
		   LA_DenseMatrix* orth ) const ;

      bool eigenV( size_t matord,
                   LA_SeqVector* diagTDMA,
                   LA_SeqVector* subdTDMA,
                   LA_DenseMatrix* eigv ) const ;

   //-- Plug in

      LA_SymmetricMatrix( void ) ;

      virtual LA_SymmetricMatrix* create_replica(
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool is_symmetric_POST( bool result ) const ;

   //-- Class attributes

      static LA_SymmetricMatrix const* PROTOTYPE ;

   //-- Attributes

      size_t SIZE ;
      doubleVector DIAG ;  // Store by Diagonal
                           // Golub and Van Loan, Matrix Computation,
                           // second edition, Johns Hopkins, 1993.
} ;

#endif
