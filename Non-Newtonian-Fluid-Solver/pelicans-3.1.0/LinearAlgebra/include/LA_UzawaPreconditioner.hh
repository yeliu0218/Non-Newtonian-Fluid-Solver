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

#ifndef LA_UZAWA_PRECONDITIONER_HH
#define LA_UZAWA_PRECONDITIONER_HH

#include <PEL_Object.hh>

class LA_Implementation ;
class LA_Matrix ;
class LA_Vector ;
class PEL_ModuleExplorer ;

/*
Preconditioners associated to `LA_UzawaSolver::' objects.
*/

class PEL_EXPORT LA_UzawaPreconditioner : public PEL_Object
{

   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance according to the entries of `exp'.
      static LA_UzawaPreconditioner* create( PEL_Object* a_owner,
                                             PEL_ModuleExplorer const* exp ) ;

      virtual LA_UzawaPreconditioner* create_clone( 
                                             PEL_Object* a_owner ) const = 0 ;

   //-- Status

      // Has `self' been built ?
      virtual bool is_valid( void ) const = 0 ;

      // dimension of `self' 
      virtual size_t dimension( void ) const = 0 ;
      
      // implementation of `self' 
      virtual LA_Implementation const* implementation( void ) const = 0 ;

   //-- Building

      // Build `self' from `A' (upper left matrix), `B' (bottom left matrix)
      // and `C' (bottom right matrix).
      virtual void build( LA_Matrix const* A, 
                          LA_Matrix const* B,
                          LA_Matrix const* C ) = 0 ;

   //-- Linear system solution

      // Solve the system defined by the preconditioner as the LHS matrix
      // and by `rhs' as the RHS, and copy the solution into `sol'.
      // At completion, the success or failure of this task is given
      // by `::successful_solve'.
      virtual void solve( LA_Vector const* rhs, LA_Vector* sol ) = 0 ;

      // Did `::solve' managed to solve successfully its linear system 
      // at completion of its last call ?
      virtual bool successful_solve( void ) const = 0 ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
            
      virtual ~LA_UzawaPreconditioner( void ) ;

      LA_UzawaPreconditioner( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool create_clone_POST( LA_UzawaPreconditioner* result,
                                      PEL_Object* a_owner ) const ;
      
      virtual bool dimension_PRE( void ) const ;
      
      virtual bool implementation_PRE( void ) const ;

      virtual bool build_PRE( LA_Matrix const* A,
                              LA_Matrix const* B,
                              LA_Matrix const* C ) const ;
      virtual bool build_POST( LA_Matrix const* A,
                               LA_Matrix const* B,
                               LA_Matrix const* C ) const ;

      virtual bool solve_PRE( LA_Vector const* rhs,
                              LA_Vector const* sol ) const ;
      virtual bool solve_POST( LA_Vector const* rhs,
                               LA_Vector const* sol ) const ;

   private: //----------------------------------------------------------

      LA_UzawaPreconditioner( void ) ;
      LA_UzawaPreconditioner( LA_UzawaPreconditioner const& other ) ;
      LA_UzawaPreconditioner& operator=( 
                              LA_UzawaPreconditioner const& other ) ;

} ;

#endif

