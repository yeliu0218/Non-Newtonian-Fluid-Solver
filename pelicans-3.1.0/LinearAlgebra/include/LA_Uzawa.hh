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

#ifndef LA_UZAWA_HH
#define LA_UZAWA_HH

#include <LA_TwoBlocksMethod.hh>

class LA_Solver ;

/*
Uzawa iterative solvers for block linear systems of the form

       |  A   tB  |  | U |   | F |
       |          |  |   | = |   |
       |  B    C  |  | P |   | G |

where the sequence of approximations x(k) and y(k) are given as follows :

   for (k=0) until convergence, do
      solve    A.U(k+1) + tB.P(k) + r tB.S.( B.U(k+1) + C.P(k) - G ) = F
      compute  P(k+1) = P(k) + rho S.( B.U(k+1) + C.P(k+1) - G )
   enddo

PUBLISHED
*/

class PEL_EXPORT LA_Uzawa : public LA_TwoBlocksMethod
{
   public: //-----------------------------------------------------------
            
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_Uzawa( void ) ;
      LA_Uzawa( LA_Uzawa const& other ) ;
      LA_Uzawa& operator=( LA_Uzawa const& other ) ; 
      
      LA_Uzawa( PEL_Object* a_owner,
                PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in
      
      LA_Uzawa( void ) ;
      
      virtual LA_Uzawa* create_replica(
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;
      
   //-- System profile
         
      virtual void set_matrix_prototype_sub( LA_Matrix const* mat ) ;
         
      virtual void re_initialize_internals_sub( size_t nv_glob, 
                                                size_t np_glob,
                                                size_t nv_loc, 
                                                size_t np_loc,
                                                size_t& nv_loc_final, 
                                                size_t& np_loc_final ) ;

   //-- Auxiliary items
            
      virtual bool S_is_required( void ) const ;
         
      virtual void set_S( LA_Vector* a_S ) ;
         
   //-- System setting
            
      virtual void set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                                   LA_Vector* a_F, LA_Vector* a_G,
                                   LA_Matrix* a_C ) ;
      
      virtual void unset_system_sub( void ) ;

   //-- Estimation
            
      virtual void estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                          bool has_init_P, LA_Vector* P ) ;
         
      virtual void estimate_unknowns_sub_1( bool has_init_U, LA_Vector* U, 
                                            bool has_init_P, LA_Vector* P ) ;

      virtual void estimate_unknowns_sub_2( bool has_init_U, LA_Vector* U, 
                                            bool has_init_P, LA_Vector* P ) ;
      
   //-- Internals
      
      void augment_system( void )  ;
      
      void print_errors( size_t n, double err1, double err2, size_t nb_it ) ;      

      void print_end_errors( size_t n, double err1, double err2 ) ;

      void print_errors_2( size_t n,
                           double err1, double err2, 
                           size_t nb_it_1, size_t nb_it_2 ) ;      

      void print_end_errors_2( size_t n, double err1, double err2 ) ;

   //-- Class attributes
      
      static LA_Uzawa const* PROTOTYPE ;
      
   //-- Attributes
      
      // allocated by the client
      LA_Matrix* A ;
      LA_Matrix* B ;
      LA_Matrix* C ;
      LA_Vector* F ;
      LA_Vector* G ;
      LA_Vector* S ;
      
      // allocated internally
      LA_Matrix* BtS ;
      LA_Matrix* ImrC ;
      LA_Vector* F0 ;
      LA_Vector* PRES ;
      LA_Vector* P0 ;
      LA_Vector* DU ;
      LA_Vector* H ;
      
      LA_Solver* SOLVER_U ;
      LA_Solver* SOLVER_P ;
      
      double RR ;
      double RHO ;
      double TOL_VELO ;
      double TOL_DIV ;
      size_t MAXITS ;
} ; 

#endif
