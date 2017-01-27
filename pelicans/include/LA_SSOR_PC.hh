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

#ifndef LA_SSOR_PC_HH
#define LA_SSOR_PC_HH

#include <LA_Preconditioner.hh>

class LA_SOR_PC ;

/*
Symmetric Successive Over Relaxation (SSOR) sequential precontitioners.

For a matrix `A' = `L' + `D' + `U'

       `D': the diagonal of `A'     ie `D'(i,j) = 0 for i!=j
       `U': the upper part of `A'   ie `U'(i,j) = 0 for i<=j
       `L': the lower part of `A'   ie `L'(i,j) = 0 for i>=j

the solution `x' of the linear system:

                   `A' * `x' = `b'
                   
is approximated performing a fixed number of SOR "forward" then "backward"
relaxations (cf. `LA_SOR_PC::' preconditioner).

If this number of successive relaxations is one, the SSOR preconditioner
`M' related to a matrix `A' is:

     `M' = 1./`omega'/( 2.-`omega' )( `D'+`omega'*`L' ) * `D'^1 * ( `D'+`omega'*`U' )

where `omega' is a free parameter in ]0,2[.

Remark:

    For a symmetrical matrix `A' ( ie `U' = `L'^t ):

      `M' = `T' * `T'^t

    with
    
      `T' = ( `D'+`omega'*`L' ) * `D'^0.5 / ( `omega' * (2.-`omega') )^0.5

PUBLISHED
*/

class PEL_EXPORT LA_SSOR_PC : public LA_Preconditioner
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static LA_SSOR_PC* create( PEL_Object* a_owner,
                                 double omega,
                                 size_t nb_iters,
                                 double smallest_inverted_item ) ;
            
      virtual LA_SSOR_PC* create_clone( PEL_Object* a_owner ) const ;

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

     ~LA_SSOR_PC( void ) ;
      LA_SSOR_PC( LA_SSOR_PC const& other ) ;
      LA_SSOR_PC& operator=( LA_SSOR_PC const& other ) ;

      LA_SSOR_PC( PEL_Object* a_owner,
                  double omega,
                  size_t nb_iters,
                  double smallest_inverted_item ) ;

      LA_SSOR_PC( PEL_Object* a_owner, LA_SSOR_PC const* other ) ;
                  

   //-- Plug in
      
      LA_SSOR_PC( void ) ;

      virtual LA_SSOR_PC* create_replica( 
                               PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static LA_SSOR_PC const* PROTOTYPE ;

   //-- Attributes

      LA_SOR_PC* const SOR_PC ;
} ;

#endif

        
