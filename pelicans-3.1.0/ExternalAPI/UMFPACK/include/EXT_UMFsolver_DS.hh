/*
 *  Copyright 1995-2013 by IRSN
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

#ifndef EXT_UMF_SOLVER_DS_HH
#define EXT_UMF_SOLVER_DS_HH

#include <LA_Solver.hh>

class LA_Matrix ;
class LA_SeqMatrix ;
class LA_Vector ;

class PEL_ModuleExplorer ;

/*
Linear system solvers using UMFpack direct solvers.

PUBLISHED
*/

class PEL_EXPORT EXT_UMFsolver_DS : public LA_Solver
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static EXT_UMFsolver_DS* create( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;

      virtual EXT_UMFsolver_DS* create_clone( PEL_Object* a_owner ) const ;
      
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      EXT_UMFsolver_DS( void ) ;
     ~EXT_UMFsolver_DS( void ) ;
      EXT_UMFsolver_DS( EXT_UMFsolver_DS const& other ) ;
      EXT_UMFsolver_DS& operator=( EXT_UMFsolver_DS const& other ) ;

      EXT_UMFsolver_DS( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp ) ;
      EXT_UMFsolver_DS( PEL_Object* a_owner,
                        EXT_UMFsolver_DS const* other ) ;

   //-- Plug in
      
      virtual EXT_UMFsolver_DS* create_replica( 
                                   PEL_Object* a_owner,
				   PEL_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution
      
      virtual void set_matrix_self( LA_Matrix const* mat,  bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
      virtual void unset_matrix_self( void ) ;

   //-- Internals

      void initialize( void ) ;

      void inverse( bool& ok ) ;

      void raise_fatal_error( std::string const& routine, int status ) const ;
      
   //-- Class attributes

      static EXT_UMFsolver_DS const* PROTOTYPE ;
      
   //-- Attributes

      PEL_ModuleExplorer const* EXP ;
      
      size_t SIZE ;
      size_t NB_NO_ZERO ;

      int* UMF_Ap ;
      int* UMF_Ai ;
      double* UMF_Ax ;
      double* UMF_b ;
      double* UMF_x ;
      
      void* UMF_Numeric ;

      double* UMF_Control ;
      double* UMF_Info ;

      bool VERBOSE ;
} ;

#endif
