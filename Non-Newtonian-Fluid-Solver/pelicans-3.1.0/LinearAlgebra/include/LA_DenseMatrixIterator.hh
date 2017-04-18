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

#ifndef LA_DenseMatrixIterator_HH
#define LA_DenseMatrixIterator_HH

#include <LA_MatrixIterator.hh>

class LA_DenseMatrix ;

/*
`LA_MatrixIterator::' objects attached to `LA_DenseMatrix::' objects.

PUBLISHED
*/

class PEL_EXPORT LA_DenseMatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance attached to `A'.
      static LA_DenseMatrixIterator* create( PEL_Object* a_owner, 
                                             LA_DenseMatrix const* A ) ;

   //-- Status

      virtual bool is_valid( void ) const ;

   //-- Cursor movement

      virtual void start_all_items( void )  ;

      virtual void start_row_items( size_t i_row ) ;
      
      virtual void go_next( void ) ;  

   //-- Access

      virtual LA_Matrix const* matrix( void ) const ;
      
      virtual size_t nb_rows( void ) const ;

      virtual size_t row( void ) const ;

      virtual size_t col( void ) const ;

      virtual double item( void ) const ;

   //-- Hidden

      virtual void set_item( double x ) const ;

      virtual void add_to_item( double x ) const ;

   protected: //------------------------------------------------------------

   private : //-------------------------------------------------------------
  
      LA_DenseMatrixIterator( PEL_Object* a_owner, 
                            LA_DenseMatrix const* A ) ;
      
      LA_DenseMatrixIterator( void ) ;
     ~LA_DenseMatrixIterator( void ) ;
      LA_DenseMatrixIterator( LA_DenseMatrixIterator const& other ) ;
      LA_DenseMatrixIterator& operator=( LA_DenseMatrixIterator const& other ) ;
      
   //-- Attributes
    
      LA_DenseMatrix const* const MAT ;

      size_t I_ROW ;
      size_t I_COL ;
      bool STAY_IN_ROW ;
      size_t ROW_DIM ;
      size_t COL_DIM ;
} ;

#endif
