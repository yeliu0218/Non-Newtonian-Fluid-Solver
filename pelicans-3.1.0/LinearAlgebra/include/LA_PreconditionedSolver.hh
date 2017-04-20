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

#ifndef LA_PreconditionedSolver_HH
#define LA_PreconditionedSolver_HH

#include <LA_Solver.hh>

class LA_Preconditioner ;
class LA_IterativeSolver ;
class LA_SeqVector ;
class LA_SeqMatrix ;

class PEL_EXPORT LA_PreconditionedSolver : public LA_Solver
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_PreconditionedSolver* create(
                     PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

      virtual LA_PreconditionedSolver* create_clone(
                                              PEL_Object* a_owner ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;      

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~LA_PreconditionedSolver( void ) ;
      LA_PreconditionedSolver( LA_PreconditionedSolver const& other ) ;
      LA_PreconditionedSolver& operator=(
                               LA_PreconditionedSolver const& other ) ;

      LA_PreconditionedSolver( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) ;
      LA_PreconditionedSolver( PEL_Object* a_owner,
                               LA_PreconditionedSolver const* other ) ;

   //-- Plug in

      LA_PreconditionedSolver( void ) ;

      virtual LA_PreconditionedSolver* create_replica( 
                                      PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const ;


   //-- Linear system resolution
      
      virtual void unset_matrix_self( void ) ;
      
      virtual void set_matrix_self( LA_Matrix const* mat, bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
   //-- Class attributes

      static LA_PreconditionedSolver const* PROTOTYPE ;
      
   //-- Attributes

      LA_IterativeSolver* SOLVER ;
      LA_Preconditioner* PRECOND ;
      LA_Matrix const* MY_A ;
} ;

#endif
