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

#ifndef EXT_AZTEC_SOLVER_HH
#define EXT_AZTEC_SOLVER_HH

#include <LA_Solver.hh>

#include "az_aztec.h"
#undef min
#undef max

/*
Wrappers around Aztec iterative solvers.
*/

class EXT_AztecSolver : public LA_Solver
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static EXT_AztecSolver* create( PEL_Object* a_owner, 
                                     PEL_ModuleExplorer const* exp ) ;

      virtual EXT_AztecSolver* create_clone( PEL_Object* a_owner ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;      

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~EXT_AztecSolver( void ) ;
      EXT_AztecSolver( EXT_AztecSolver const& other ) ;
      EXT_AztecSolver& operator=( EXT_AztecSolver const& other ) ;

      EXT_AztecSolver( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;
      EXT_AztecSolver( PEL_Object* a_owner, EXT_AztecSolver const* other ) ;

      void initialize( void ) ;

   //-- Plug in

      EXT_AztecSolver( void ) ;

      virtual EXT_AztecSolver* create_replica( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution      

      void build_ksp( PEL_ModuleExplorer const* exp ) ;

      void build_pc( PEL_ModuleExplorer const* exp ) ;
      
      virtual void unset_matrix_self( void ) ;
     
      virtual void set_matrix_self( LA_Matrix const* mat,
                                    bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, 
                               LA_Vector* x,
                               size_t &nb_iter,
                               bool &ok ) ;

   //-- Class attributes

      static EXT_AztecSolver const* PROTOTYPE ;

   //-- Attributes

      PEL_ModuleExplorer const* const EXP ;
      
      int nrow ;
      int proc_config[AZ_PROC_SIZE] ; // Processor information 
      int options[AZ_OPTIONS_SIZE] ;  // Array used to select solver options
      double params[AZ_PARAMS_SIZE] ; // User selected solver parameters
      int* data_org ;                 // Array to specify data layout
      double status[AZ_STATUS_SIZE] ; // Information returned from AZ_solve()
      int* update ;                   // vector elements updated on this node
      int* external ;                 // vector elements needed by this node
      int* update_index ;             // ordering of update[] and external[]
      int* extern_index ;             // locally on this processor.
      int* my_Ai ;                    // Sparse matrix to be solved is stored
      double* my_Av ;                 // in these MSR arrays.
      int N_update ;                  // # of unknowns updated on this node 
} ;

#endif
