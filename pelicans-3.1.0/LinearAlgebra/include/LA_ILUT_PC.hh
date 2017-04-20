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

#ifndef LA_ILUT_PC_H
#define LA_ILUT_PC_H

#include <LA_Preconditioner.hh>

#include <doubleVector.hh>
#include <intVector.hh>

class LA_SeqMatrix ;

/*
Precopnditioners defined as an Incomplete LU factorization with
Threshold [Saad, 1996, section 10.4], including a  possible diagonal 
compensation  strategy [Saad, 1996, p286]  (MILU preconditioner).

Algorithm : The ILUT factorization is based on the IKJ variant of
the Gaussian elimination [Saad, 1996, p272], that reads for
a matrix a_ij (i=1..n , j=1..n) :

   1 For i = 2 -> n
   2   tol_i = ( 1-norm of row i )*tolerance
   3   For k = 1 -> i-1 with |a_ik|>tol_i
   4     a_ik = a_ik / a_kk
   5     For j = k+1 -> n with (i,j) not in P
   6       a_ij = a_ij - a_ik*a_kj
   7     End
   8   End
   9 End

where P is a subset of all possible pairs (i,j), called a 
zero pattern. The diagonal compensation strategy consist in
adding up all the elements that have been dropped at the completion
of the k loop (step 7) and then substracting this sum from the
diagonal entry in U.

The zero pattern P is such that (i,j) belongs to P if and only if :
   1.  |a_ij| > tol_i
   2.  if j>i a_ij is among the p'th elements of largest absolute value
       of the i'th line of the upper triangular part of the 
       factorized matrix
   3.  if j>i a_ij is among the p'th elements of largest absolute value
       of the i'th line of the lower triangular part of the 
       factorized matrix

Reference : [Saad, 1996]
   Y. Saad, 1996
   Iterative Methods for Sparse Linear Systems
   PWS Publishing Company 
*/
                
class PEL_EXPORT LA_ILUT_PC : public LA_Preconditioner
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /*
      Create and return an instance.
      The maximum number of nonzero non diagonal elements in each line of L 
      and of U in the LU factorized matrix is `order'.
      The elements of the a row of the factorized matrix whose absolute value 
      is lower than `tolerance' times the 1-norm of that row are dropped.
      The diagonal compensation strategy is adopted (leading to a MILU 
      preconditioner) if `diagonal_compensation' is true. `matrix_prototype'
      is used to select the type of an internal matrix that is created 
      by calling `LA_SeqMatrix::create_clone' with `matrix_prototype'
      as target object.
      */
      static LA_ILUT_PC* create( PEL_Object* a_owner,
                                 size_t order ,
                                 double tolerance,
                                 bool diagonal_compensation,
                                 double smallest_nonzero_pivot,
                                 LA_SeqMatrix const* matrix_prototype ) ;

      virtual LA_ILUT_PC* create_clone( PEL_Object* a_owner ) const ;

   //-- Status

      virtual bool is_valid( void ) const ;

      virtual size_t dimension( void ) const ;

   //-- Building

      virtual void build( LA_Matrix const* mat ) ;
      
      virtual void unbuild( void ) ;
      
   //-- Linear system solution

      virtual void solve( LA_Vector const* rhs, LA_Vector* sol ) ;

      virtual bool successful_solve( void ) const ;

   //-- Input - Output

      virtual void print_more( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_ILUT_PC( void ) ;
      LA_ILUT_PC( LA_ILUT_PC const& other ) ;
      LA_ILUT_PC& operator=( LA_ILUT_PC const& other ) ;

      LA_ILUT_PC( PEL_Object* a_owner,
                  size_t order,
                  double tolerance,
                  bool diagonal_compensation,
                  double smallest_nonzero_pivot,
                  LA_SeqMatrix const* matrix_prototype ) ;

   //-- Plug in
      
      LA_ILUT_PC( void ) ;

      virtual LA_ILUT_PC* create_replica( 
                               PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp ) const ;

   //--

      double keepLUplargest( int i,
                             double norm,
                             int lowEl,
                             int upEl ) ;

   //-- Class attributes

      static LA_ILUT_PC const* PROTOTYPE ;

   //-- Attributes

      double const RULE ;
      size_t const LFIL ;
      bool const MODIFIED ;
      double const PIV_MIN ;
      
      LA_SeqMatrix* const LU ;
      intVector COL_P_LARGEST ;
      doubleVector VAL_P_LARGEST ;
      bool BUILD_OK ;
      bool SOLVE_OK ;
} ;

#endif
