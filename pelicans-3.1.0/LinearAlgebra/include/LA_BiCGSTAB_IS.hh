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

#ifndef LA_BICGSTAB_IS_HH
#define LA_BICGSTAB_IS_HH

#include <LA_IterativeSolver.hh>

class LA_Preconditioner ;
class LA_Matrix ;
class LA_Vector ;

class PEL_ModuleExplorer ;

/*
Preconditioned BiCGSTAB iterative solvers.

PUBLISHED
*/

class LA_BiCGSTAB_IS : public LA_IterativeSolver
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_BiCGSTAB_IS( void ) ;
      LA_BiCGSTAB_IS( LA_BiCGSTAB_IS const& other ) ;
      LA_BiCGSTAB_IS& operator=( LA_BiCGSTAB_IS const& other ) ;

      LA_BiCGSTAB_IS( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

      LA_BiCGSTAB_IS( PEL_Object* a_owner, LA_BiCGSTAB_IS const* other ) ;

   //-- Instance delivery and initialization

      virtual LA_BiCGSTAB_IS* create_clone( PEL_Object* a_owner ) const ;

   //-- Plug in

      LA_BiCGSTAB_IS( void ) ;

      virtual LA_BiCGSTAB_IS* create_replica( 
                                     PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) const ;

   //-- Internals
      
      void reset_internals( LA_Vector const* prototype ) ;
      
      virtual void do_solve( LA_Matrix const* A,
                             LA_Vector const* b,
                             LA_Preconditioner* prec,
                             bool zero_initial_guess,
                             LA_Vector* x ) ;      

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static LA_BiCGSTAB_IS const* PROTOTYPE ;

   //-- Attributes

      size_t SIZE ;
      LA_Vector* R ;
      LA_Vector* RP ;
      LA_Vector* V ;
      LA_Vector* T ;
      LA_Vector* S ;
      LA_Vector* P ;
} ;


#endif
