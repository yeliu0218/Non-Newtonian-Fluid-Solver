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

#ifndef GE_MATRIX_HH
#define GE_MATRIX_HH

#include <PEL_Object.hh>
#include <PEL_assertions.hh>
class doubleArray2D ;


/*
Matrices of double, specialized for geometrical manipulation. 
*/

class PEL_EXPORT GE_Matrix : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static GE_Matrix* create( PEL_Object* a_owner,
                                size_t row_nb,
                                size_t col_nb );


   //-- Access

      size_t nb_rows( void ) const ;

      size_t nb_cols( void ) const ;

      double item( size_t i, size_t j ) const ;

      // determinant of `self'
      double determinant( void ) const ;

      bool determinant_is_computed( void ) const ;
      
   //-- Element change

      void nullify( void ) ;

      // Assign `x' to the coefficient at the `i'-th row and `j'-th column.
      void set_item( size_t i, size_t j, double x ) ;

      // Add `x' to the coefficient at the `i'-th row and `j'-th column.
      void add_to_item( size_t i, size_t j, double x ) ;
      
      // Compute determinant of `self'.
      void compute_determinant( void ) ;

  //-- System resolution

      void invert( doubleArray2D const& B, doubleArray2D& X ) const ;

   //-- Input - Output

      // Print matrix line by line.
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      GE_Matrix( void ) ;
     ~GE_Matrix( void ) ;
      GE_Matrix( GE_Matrix const& other ) ;
      GE_Matrix const& operator=( GE_Matrix const& other ) ;

      GE_Matrix( PEL_Object* a_owner, size_t row_nb, size_t col_nb ) ;
      
      void gauss1x1_multi( doubleArray2D const& B, doubleArray2D& X ) const ;

      void gauss2x2_multi( doubleArray2D const& B, doubleArray2D& X ) const ;

      void gauss3x3_multi( doubleArray2D const& B, doubleArray2D& X ) const ;
      
   //-- Attributes
     
      size_t const rows ;
      size_t const cols ;
      double mat[9] ;
      double det ;
      bool det_comp ;
} ;

#ifndef OUTLINE
#include <GE_Matrix.icc>
#endif

#endif
