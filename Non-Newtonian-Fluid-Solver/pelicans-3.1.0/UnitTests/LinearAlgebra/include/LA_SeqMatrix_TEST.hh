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

#ifndef LA_SPARSE_MATRIX_TEST_HH
#define LA_SPARSE_MATRIX_TEST_HH

#include <LA_Matrix_TEST.hh>
#include <LA_SeqMatrix.hh>

class LA_SeqMatrix ;
class doubleArray2D ;

class PEL_EXPORT LA_SeqMatrix_TEST : public LA_Matrix_TEST
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------


      virtual ~LA_SeqMatrix_TEST( void ) ;

   //-- Elementary tests management
      
      LA_SeqMatrix_TEST( std::string const& matrixClass,
                            std::string const& registration ) ;
      
      LA_SeqMatrix* create_matrix(
                         PEL_Object* a_owner,
                         size_t n, size_t m ) const ;

      void fill_sparse_matrix( LA_SeqMatrix* mat,
                               doubleArray2D& values ) const ;

      virtual void test( PEL_ModuleExplorer const* exp,
                         LA_Matrix const* matrix,
                         std::string const& test_matrix  ) ;
      
      virtual void doBasicTests( void ) ;
      virtual void doBlockVectorTests( void ) ;
      virtual void doBlas2Tests( void ) ;
      virtual void doBlas3Tests( void ) ;
      virtual void doBasicSparseTests( void ) ;
      virtual void doSolveLUTest( LA_SeqMatrix const* matA ) ;
      
   private: //--------------------------------------------------------------

      LA_SeqMatrix_TEST( void ) ;
      LA_SeqMatrix_TEST( LA_SeqMatrix_TEST const& other ) ;
      LA_SeqMatrix_TEST& operator=( LA_SeqMatrix_TEST const& other ) ;

      static LA_SeqMatrix_TEST const* PROTOTYPE ;
      bool SYM ;

} ;


#endif
