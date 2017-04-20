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

#ifndef LA_DIRECT_SOLVER_GAUSS_HH
#define LA_DIRECT_SOLVER_GAUSS_HH

#include <LA_Solver.hh>

#include <size_t_vector.hh>

class LA_SeqVector ;
class LA_DenseMatrix ;
class LA_SeqMatrix ;
class LA_Vector ;

/*
Linear system solvers using the Gauss LU factorization method.

PUBLISHED
*/

class PEL_EXPORT LA_GaussLU_DS : public LA_Solver
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static LA_GaussLU_DS* create( PEL_Object* a_owner,
                                    double const a_pivot_minimal_value ) ;

      virtual LA_GaussLU_DS* create_clone( PEL_Object* a_owner ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_GaussLU_DS( void ) ;
      LA_GaussLU_DS( LA_GaussLU_DS const& other ) ;
      LA_GaussLU_DS& operator=( LA_GaussLU_DS const& other ) ;

      LA_GaussLU_DS( PEL_Object* a_owner,
                     double const a_pivot_minimal_value ) ;
      LA_GaussLU_DS( PEL_Object* a_owner, LA_GaussLU_DS const* other ) ;
      
   //-- Plug in

      LA_GaussLU_DS( void ) ;
      
      virtual LA_GaussLU_DS* create_replica( 
                                   PEL_Object* a_owner,
				   PEL_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution

      virtual void set_matrix_self( LA_Matrix const* mat, bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
      virtual void unset_matrix_self( void ) ;
      
      // Do Gauss-LU factorization with row pivoting
      //    (cf Lascaux-Theodor 2nd edition Tome1 p197).
      void factorize_LU( LA_DenseMatrix* A,
                         size_t_vector& piv,
                         bool &ok ) const ;

      // Solve linear system where A and piv are Gauss-LU factorized
      // with row pivoting of some matrix
      //    (cf Lascaux-Theodor 2nd edition Tome1 p198).
      void solve_LU( LA_SeqMatrix const* LU,
                     size_t_vector const& piv,
                     LA_SeqVector const* b,
                     LA_SeqVector* x ) const ;
      
      void gauss1x1( LA_SeqVector const* b,
                     LA_SeqVector* x ) const ;
      void gauss2x2( LA_SeqVector const* b,
                     LA_SeqVector* x ) const ;
      void gauss3x3( LA_SeqVector const* b,
                     LA_SeqVector* x ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Class attributes

      static LA_GaussLU_DS const* PROTOTYPE ;
      
   //-- Attributes

      double const PIV_MIN_VAL ;
      
      LA_DenseMatrix* const MAT ;
      size_t_vector PIV ;
      double DET ;
} ;


#endif
