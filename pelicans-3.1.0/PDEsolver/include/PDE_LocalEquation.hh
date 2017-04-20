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

#ifndef PDE_LOCAL_EQUATION_HH
#define PDE_LOCAL_EQUATION_HH

#include <PEL_Object.hh>

#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

/*
Linear equations with two levels of indices :
           a(i,ic;j,jc) x = b(i,ic)
where i,j are the major indices, and ic,jc are the subindices.

The index i (resp. j) is the row (resp. column) index. 

Each major index is connected to a node index

The meaning of the subindices is left to the user.
*/

class PEL_EXPORT PDE_LocalEquation : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static PDE_LocalEquation* create( PEL_Object* a_owner ) ;

      // Set dimensions and initial item values.
      void initialize( size_t_vector const& row_node_connectivity,
                       size_t a_nb_row_sub_indices,
                       size_t_vector const& column_node_connectivity,
                       size_t a_nb_col_sub_indices ) ;
      
      // Re-initialize `self' and set the matrix part of as the transpose of 
      // its counterpart in `other'.
      void set_as_transpose( PDE_LocalEquation const* other ) ;

   //-- Access

      // number of rows (i.e. of row major indices)
      size_t nb_rows( void ) const ;

      // number of subindices for the rows
      size_t nb_row_sub_indices( void ) const ;

      // node connected to the i'th row index
      size_t row_node( size_t i ) const ;

      // number of columns (i.e. of columns major indices)
      size_t nb_columns( void ) const ;

      // number of subindices for the columns
      size_t nb_column_sub_indices( void ) const ;

      // node connected to the `j'-th column index
      size_t column_node( size_t j ) const ;

      // Is the left hand side matrix element of row index `i', 
      // row subindex `ic', column index `j', column subindex `jc'
      // marked as set ?
      bool matrix_item_is_set(  size_t i, size_t j,
                                size_t ic=0, size_t jc=0 ) const ;
      
      // the left hand side matrix element of row index `i', 
      // row subindex `ic', column index `j', column subindex `jc' 
      double matrix_item( size_t i, size_t j,
                          size_t ic=0, size_t jc=0 ) const ;

      // the right hand side vector element of row index `i', 
      // and row subindex `ic'
      double vector_item( size_t i, size_t ic=0 ) const ;

   //-- Modification

      // Increase by `x' the left hand side matrix element of row index `i', 
      // row subindex `ic', column index `j', column subindex `jc'.
      void add_to_matrix( double x, size_t i, size_t j, size_t ic, size_t jc ) ;
      // Increase by `x' the left hand side matrix element of row index `i' 
      // and column index `j', with the two subindices equal to 0.
      void add_to_matrix( double x, size_t i, size_t j ) ;

      // Increase by `x' the right hand side vector element of row index `i', 
      // row subindex `ic'. 
      void add_to_vector( double x, size_t i, size_t ic ) ;

      // Increase by `x' the right hand side vector element of row index `i', 
      // row subindex 0.
      void add_to_vector( double x, size_t i ) ;

      // Multiply by `coeff' all the left hand side matrix elements.
      void scale_matrix( double coeff ) ;

      // Multiply by `coeff' all the right hand side vector elements.
      void scale_vector( double coeff ) ;

      // Unset the left hand side matrix element of row index `i',
      // row subindex `ic', column index `j', column subindex `jc'.
      void nullify_then_unset_matrix_item( size_t i, size_t j,
                                           size_t ic=0, size_t jc=0 ) ;

      // Nullify, and mark unset, the left hand side matrix elements whose 
      // absolute value is lower than `relative_tolerance' times the matrix 
      // maximal absolute value.
      void nullify_then_unset_small_matrix_items( 
                                    double relative_tolerance = 1.E-8 ) ;

      // Nullify, and mark as unset, all the left hand side matrix elements.
      void nullify_then_unset_all_matrix_items( void ) ;
      
      // Mark as set all the left hand side matrix elements.
      void mark_all_matrix_items_as_set( void ) ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_LocalEquation( void ) ;
     ~PDE_LocalEquation( void ) ;     
      PDE_LocalEquation( PDE_LocalEquation const& other ) ;
      PDE_LocalEquation& operator=( PDE_LocalEquation const& other ) ;
   
      PDE_LocalEquation( PEL_Object* a_owner ) ;

      size_t mat_idx( size_t i, size_t j, size_t ic, size_t jc ) const ;
      size_t mat_idx( size_t i, size_t j ) const ;
      size_t vec_idx( size_t i, size_t ic ) const ;

      // overloaded to avoid implicit conversion
      void add_to_matrix( size_t x, size_t i, size_t j, size_t ic, double jc ) ;
      void add_to_matrix( size_t x, size_t i, double j ) ;
      void add_to_vector( size_t x, size_t i, double iComp ) ;
      void add_to_vector( size_t x, double i ) ;

   //-- Attributes

      size_t nbRows ;
      size_t nbCols ;
      size_t nbComps1 ;
      size_t nbComps2 ;

      bool* A_IS_SET ;
      double* A ;
      double* b ;
      size_t SIZE_A ;
      size_t SIZE_b ;
      
      size_t_vector loc2node1 ;
      size_t_vector loc2node2 ;
} ;

#ifndef OUTLINE
   #include <PDE_LocalEquation.icc>
#endif

#endif
