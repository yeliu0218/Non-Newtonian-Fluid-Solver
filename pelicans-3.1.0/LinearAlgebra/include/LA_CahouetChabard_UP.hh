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

#ifndef LA_CAHOUET_CHABARD_UP_HH
#define LA_CAHOUET_CHABARD_UP_HH

#include <LA_UzawaPreconditioner.hh>

class LA_Implementation ;
class LA_Solver ;
class LA_Matrix ;
class LA_Vector ;
class PEL_ModuleExplorer ;

/*
Preconditionner M for the linear system

       |  A   tB  |  | X |   | f |
       |          |  |   | = |   |
       |  B   -C  |  | Y |   | g |

to be solved using `LA_UzawaSolver::' objects, where M is defined by 
       M = B diag(A)^{-1} tB+C

reference : J. Cahouet, J.-P. Chabard
            Some Fast 3D Finite Element Solvers for the Generalized
            Stokes Problem
            Int. J. for Numer. Methods Fluids, Vol.8, 869-895 (1988) 
*/

class PEL_EXPORT LA_CahouetChabard_UP : public LA_UzawaPreconditioner
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Create and return an instance according to the entries of `exp'
      static LA_CahouetChabard_UP* create( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp ) ;

      virtual LA_CahouetChabard_UP* create_clone( PEL_Object* a_owner ) const ;

   //-- Status

      virtual bool is_valid( void ) const ;
      
      virtual size_t dimension( void ) const ;
      
      virtual LA_Implementation const* implementation( void ) const ;

   //-- Building

      virtual void build( LA_Matrix const* A,
                          LA_Matrix const* B,
                          LA_Matrix const* C ) ;
      
   //-- Linear system solution

      virtual void solve( LA_Vector const* rhs, LA_Vector* sol ) ;

      virtual bool successful_solve( void ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

      LA_CahouetChabard_UP( void ) ;
     ~LA_CahouetChabard_UP( void ) ; 
      LA_CahouetChabard_UP( LA_CahouetChabard_UP const& other ) ;
      LA_CahouetChabard_UP& operator=( 
                            LA_CahouetChabard_UP const& other ) ;

      LA_CahouetChabard_UP( PEL_Object* a_owner, 
                            LA_Matrix const* a_matrix_prototype,
                            LA_Solver* a_solver_of_M ) ;

      void build_prec( LA_Matrix const* A,
                       LA_Matrix const* B,
                       LA_Matrix const* C,
                       LA_Matrix* mres ) const ;
      
   //-- Attributes

      double const MIN_DIAG ;
      LA_Matrix* const MAT_M ;
      LA_Matrix* const BT ;
      LA_Solver* const SOLVER ;
      bool BUILD_OK ;
      bool SOLVE_OK ;
} ;

#endif
