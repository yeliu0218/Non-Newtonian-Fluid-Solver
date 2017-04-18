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

#ifndef LA_MATRIX_TEST_HH
#define LA_MATRIX_TEST_HH

#include <PEL_ObjectTest.hh>

#include <string>

class LA_SeqVector ;
class LA_SeqMatrix ;
class LA_Matrix ;
class LA_Scatter ;
class LA_Vector ;

class PEL_EXPORT LA_Matrix_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   //-- Elementary tests management
      
      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;
      
   protected: //------------------------------------------------------------

      LA_Matrix_TEST( std::string const& matrixClass,
                      std::string const& registration ) ;
      virtual ~LA_Matrix_TEST( void ) ;
 
   //-- Elementary tests management
      
      bool double_equality( double const x, double const y ) const ;
      
      bool matrix_equality(
                  LA_Matrix const* mat1, LA_SeqMatrix const* mat2 ) const ;
      
      bool vector_equality(
               LA_SeqVector const* vec1, LA_SeqVector const* vec2 ) const ;

      virtual void test( PEL_ModuleExplorer const* exp,
                         LA_Matrix const* matrix,
                         std::string const& test_matrix ) ;
      
      LA_Scatter* create_identity_scatter(
                         PEL_Object* a_owner, LA_Vector const* vec ) ;

      virtual LA_Matrix* create_matrix(
                         PEL_Object* a_owner,
                         size_t n, size_t m ) const ;
      
   private: //--------------------------------------------------------------
      
      LA_Matrix_TEST( void ) ;
      LA_Matrix_TEST( LA_Matrix_TEST const& other ) ;
      LA_Matrix_TEST& operator=( LA_Matrix_TEST const& other ) ;

      void test_vector( LA_Vector const* vec ) ;

   //-- Class attributes
      
      static LA_Matrix_TEST const* PROTOTYPE ;
      
   //-- Attributes

      double MY_EPS ;
      double MY_MIN ;
      std::string TYPE ;
      LA_Matrix* MAT_PROTO ;

} ;

#endif // LA_Matrix_TEST_HH
