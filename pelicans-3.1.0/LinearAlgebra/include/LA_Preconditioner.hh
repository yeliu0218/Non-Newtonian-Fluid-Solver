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

#ifndef LA_PRECONDITIONER_HH
#define LA_PRECONDITIONER_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;
class LA_Vector ;
class LA_Matrix ;

/*
Preconditioners for iterative solvers of linear systems.

The convergence rate of iterative methods for solving a linear system
can be improved by transforming it into one that is equivalent in the
sense that it has the same solution, but has more favorable spectral
properties. A preconditioner (that is an instance of a concrete
subclass of LA_Preconditioner) is a matrix M that effects such a 
transformation. Whatever the theoretical definition of a 
preconditioner may be, it is a remarkable property of many iterative
methods that it is possible to write the steps of the method in terms
of the only preconditioner operation : solve u from Mu=v .
The present class offers an interface for such a fonctionnality.

FRAMEWORK INSTANTIATION
*/


class PEL_EXPORT LA_Preconditioner : public PEL_Object
{

   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_Preconditioner* make( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) ;

      virtual LA_Preconditioner* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Status

      virtual bool is_valid( void ) const = 0 ;

      virtual size_t dimension( void ) const = 0 ;

   //-- Building

      // Build `self' from `mat'.
      virtual void build( LA_Matrix const* mat ) = 0 ;

      // Release preconditioned matrix.
      virtual void unbuild( void ) = 0 ;
      
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
      
      virtual void print_more( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
            
      virtual ~LA_Preconditioner( void ) ;

      LA_Preconditioner( PEL_Object* a_owner ) ;

   //-- Plug in

      LA_Preconditioner( std::string const& class_name ) ;

      virtual LA_Preconditioner* create_replica( 
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool create_clone_POST( LA_Preconditioner* result,
                                      PEL_Object* a_owner ) const ;
      
      virtual bool dimension_PRE( void ) const ;

      virtual bool build_PRE( LA_Matrix const* mat ) const ;
      virtual bool build_POST( LA_Matrix const* mat ) const ;
      
      virtual bool unbuild_PRE( void ) const ;
      virtual bool unbuild_POST( void ) const ;

      virtual bool solve_PRE(
                     LA_Vector const* rhs, LA_Vector const* sol ) const ;
      virtual bool solve_POST(
                     LA_Vector const* rhs, LA_Vector const* sol ) const ;

      virtual bool create_replica_PRE( PEL_Object const* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST( LA_Preconditioner const* result,
                                        PEL_Object const* a_owner,
                                        PEL_ModuleExplorer const* exp ) const ;
      
   private: //----------------------------------------------------------

      LA_Preconditioner( void ) ;
      LA_Preconditioner( LA_Preconditioner const& other ) ;
      LA_Preconditioner& operator=( LA_Preconditioner const& other ) ;
      
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      bool const IS_PROTO ;
} ;

#endif

