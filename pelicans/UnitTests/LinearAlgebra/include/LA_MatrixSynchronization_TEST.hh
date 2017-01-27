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

#ifndef LA_MATRIX_SYNCHRONIZATION_TEST_HH
#define LA_MATRIX_SYNCHRONIZATION_TEST_HH

#include <PEL_ObjectTest.hh>

#include <string>

class PEL_Communicator ;

class LA_SeqVector ;
class LA_SeqMatrix ;
class LA_Matrix ;
class LA_DenseMatrix ;
class LA_Scatter ;
class LA_Vector ;

class PEL_EXPORT LA_MatrixSynchronization_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      LA_MatrixSynchronization_TEST( void ) ;
     ~LA_MatrixSynchronization_TEST( void ) ;
      LA_MatrixSynchronization_TEST(
                        LA_MatrixSynchronization_TEST const& other ) ;
      LA_MatrixSynchronization_TEST& operator=(
                        LA_MatrixSynchronization_TEST const& other ) ;

   //-- Elementary tests management

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void reference_matrices_initialization( void ) ;

      void test_zero_after_make( void ) ;
      void test_zero_after_create_matrix( void ) ;
      void test_zero_after_nullify_1( void ) ;
      void test_zero_after_nullify_2( void ) ;

      void test_creation_on_only_one_process( void ) ;

      void test_set_item_1( void ) ;
      void test_set_item_2( void ) ;
      void test_mix_add_and_set( void ) ;
      void test_add_Mat_1( void ) ;
      void test_add_Mat_2( void ) ;
      void test_add_tMat_1( void ) ;
      void test_add_tMat_2( void ) ;
      void test_add_Mat_Mat_1( void ) ;
      void test_add_Mat_Mat_2( void ) ;


      bool double_equality( double x, double y ) const ;

      bool vector_equality(
                 LA_SeqVector const* vec1, LA_SeqVector const* vec2 ) const ;

      LA_Scatter* create_identity_scatter( PEL_Object* a_owner,
                                           LA_Vector const* vec ) const ;

   //-- Class attributes

      static LA_MatrixSynchronization_TEST const* PROTOTYPE ;

   //-- Attributes

      PEL_ModuleExplorer* EXP_PROTO ;

      LA_Matrix* MAT_PROTO ;
      std::string NAME ;
      size_t FIRST_ROW ;
      size_t LAST_ROW ;

      PEL_Communicator const* COM ;

      size_t DIM ;

      LA_DenseMatrix* DmatLHS ;
      LA_Matrix* matLHS ;

      LA_Vector* vecSOL ;
      double sumSOL ;
      LA_Matrix* matONES ;
      LA_SeqVector* vecRHS ;
      double sumRHS ;

      LA_Matrix* matA ;
      LA_Vector* vecX ;
      LA_Scatter* scatterX ;
      LA_SeqVector* vecSX ;


} ;

#endif
