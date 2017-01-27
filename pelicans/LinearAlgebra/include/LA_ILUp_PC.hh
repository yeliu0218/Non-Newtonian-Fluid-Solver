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

#ifndef LA_ILU_P_HH
#define LA_ILU_P_HH

#include <LA_Preconditioner.hh>

class LA_SeqMatrix ;


/*
Preconditioner defined as an Incomplete LU factorization that
keeps the p'th order fill-ins [Saad, 1996, p278] (ILU(p) 
preconditioner), including a  possible diagonal compensation 
strategy [Saad, 1996, p286]  (MILU preconditioner).

Algorithm : The ILU factorization is based on the IKJ variant of
the Gaussian elimination [Saad, 1996, p272], that reads for
a matrix a_ij (i=1..n , j=1..n) :

   1 For i = 2 -> n
   2   For k = 1 -> i-1 with (i,k) not in P
   3     a_ik = a_ik / a_kk
   4     For j = k+1 -> n with (i,j) not in P
   5       a_ij = a_ij - a_ik*a_kj
   6     End
   7   End
   8 End

where P is a subset of all possible pairs (i,j), called a 
zero pattern. The diagonal compensation strategy consist in
adding up all the elements that have been dropped at the completion
of the k loop (step 7) and then substracting this sum from the
diagonal entry in U.

In the particular case of the ILU(p) factorization, 
the zero pattern P is the set of pairs (i,j) such that 
                level(i,j) > p
where level(i,j) is called the level of fill of the (i,j) element
and is computed prior to the numerical factorization using the
following method :
    
   0  Initially, level(i,j) = 0          if a_ij is no nil or i=j
                 level(i,j) = +Infinity  otherwise
   1  For i = 2 -> n
   2    For k = 1 -> i-1
   3      level(i,k) = level(i,k)/level(k,k)
   4      For j = k+1 -> n 
   5         level(i,j)=min[ level(i,j), level(i,k)+level(k,j)+1 ]
   6      End
   7    End
   8 End  

The case p=0 coincides with the ILU(0) factorization

Reference : [Saad, 1996]
   Y. Saad, 1996
   Iterative Methods for Sparse Linear Systems
   PWS Publishing Company 
*/
                
class PEL_EXPORT LA_ILUp_PC : public LA_Preconditioner
{

   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      /*
      Create an return an instance.
      All fill-in elements whose level of fill does not exceed `order' are 
      kept.
      The diagonal compensation strategy is adopted (leading to a MILU 
      preconditioner) if `diagonal_compensation' is true. `matrix_prototype'
      is used to select the type of an internal matrix that is created 
      by calling `LA_SeqMatrix::create_clone' with `matrix_prototype'
      as target object.
      */
      static LA_ILUp_PC* create( PEL_Object* a_owner,
                                 size_t order,
                                 bool diagonal_compensation,
                                 double smallest_nonzero_pivot,
                                 LA_SeqMatrix const* matrix_prototype ) ;

      virtual LA_ILUp_PC* create_clone( PEL_Object* a_owner ) const ;

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

   private: //---------------------------------------------------------

     ~LA_ILUp_PC( void ) ;
      LA_ILUp_PC( LA_ILUp_PC const& other ) ;
      LA_ILUp_PC& operator=( LA_ILUp_PC const& other ) ;

      LA_ILUp_PC( PEL_Object* a_owner,
                  size_t order,
                  bool diagonal_compensation,
                  double smallest_nonzero_pivot,
                  LA_SeqMatrix const* matrix_prototype ) ;

   //-- Plug in
      
      LA_ILUp_PC( void ) ;

      virtual LA_ILUp_PC* create_replica( 
                               PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp ) const ;

   //--

      /** Performs the incomplete factorization */
      void factorizeLU( void ) ;

      /** Compute level matrix */
      void computeLevel( LA_SeqMatrix const* mat ) ;
      
   //-- Class attributes

      static LA_ILUp_PC const* PROTOTYPE ;

   //-- Attributes

      size_t const LEVEL_ORDER ;
      bool const MODIFIED ;
      double const PIV_MIN ;

      LA_SeqMatrix* const LU ;
      LA_SeqMatrix* const LEVELS ;
      bool BUILD_OK ;
      bool SOLVE_OK ;
      
} ;

#endif

        
