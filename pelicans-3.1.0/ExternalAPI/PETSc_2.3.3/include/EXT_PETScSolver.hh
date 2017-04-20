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

#ifndef EXT_PETSc_SOLVER_HH
#define EXT_PETSc_SOLVER_HH

#include <LA_Solver.hh>

#include <EXT_PETScAPI.hh>
#include <EXT_PETScMatrix.hh>


/*
Wrappers around PETSc iterative solvers.
*/

class EXT_PETScSolver : public LA_Solver
{
   public: //-----------------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual EXT_PETScSolver* create_clone( PEL_Object* a_owner ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~EXT_PETScSolver( void ) ;
      EXT_PETScSolver( EXT_PETScSolver const& other ) ;
      EXT_PETScSolver& operator=( EXT_PETScSolver const& other ) ;

      EXT_PETScSolver( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;
      EXT_PETScSolver( PEL_Object* a_owner, EXT_PETScSolver const* other ) ;

   //-- Plug in

      EXT_PETScSolver( void ) ;

      virtual EXT_PETScSolver* create_replica( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution
      
      virtual void set_matrix_self( LA_Matrix const* mat, bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
      virtual void unset_matrix_self( void ) ;

      void build_ksp( KSP &ksp, PEL_ModuleExplorer const* A_EXP ) ;
      
      void build_pc( KSP & ksp, PEL_ModuleExplorer const* A_EXP ) ;
      
      std::string mat_ordering( PEL_ModuleExplorer const* exp,
                                std::string const& name ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Class attributes

      static EXT_PETScSolver const* PROTOTYPE ;
      
    //-- Attributes

      PEL_ModuleExplorer const* const EXP ;
      PEL_ModuleExplorer* SUBEXP ;
      KSP MY_KSP ;
      bool HAS_TO_DESTROY_KSP ;
      EXT_PETScMatrix* MATRIX ;
      bool VERB ;
} ;

#endif
