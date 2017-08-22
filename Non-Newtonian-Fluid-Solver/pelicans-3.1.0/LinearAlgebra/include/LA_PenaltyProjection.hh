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

#ifndef LA_PENALTY_PROJECTION_HH
#define LA_PENALTY_PROJECTION_HH

#include <LA_TwoBlocksMethod.hh>

class PEL_Timer ;

class LA_Solver ;

/*
PUBLISHED
*/

class PEL_EXPORT LA_PenaltyProjection : public LA_TwoBlocksMethod
{
   public: //-----------------------------------------------------------
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_PenaltyProjection( void ) ;
      LA_PenaltyProjection( LA_PenaltyProjection const& other ) ;
      LA_PenaltyProjection& operator=( LA_PenaltyProjection const& other ) ;
      
      LA_PenaltyProjection( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in
      
      LA_PenaltyProjection( void ) ;
      
      virtual LA_PenaltyProjection* create_replica(
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
         
      virtual bool dtinv_is_required( void ) const ;

      virtual bool S_is_required( void ) const ;

      virtual void set_S( LA_Vector* a_S ) ;

      virtual bool L_is_required( void ) const ;

      virtual void set_L( LA_Matrix* a_L ) ;

      virtual bool K_is_required( void ) const ;
      
      virtual void set_K( LA_Vector* a_K ) ;
      
      virtual bool MV_is_required( void ) const ;

      virtual void set_MV( LA_Matrix* a_MV ) ;
         
   //-- System setting
         
      virtual void set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                                   LA_Vector* a_F, LA_Vector* a_G,
                                   LA_Matrix* a_C ) ;
      
      virtual void unset_system_sub( void ) ;

   //-- Estimation
            
      virtual void estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                          bool has_init_P, LA_Vector* P ) ;
      
      bool inner_solve( LA_Solver* solver, 
                        LA_Vector* rhs, 
                        LA_Vector* sol,
                        bool sol_is_init ) ;

   //-- Input - Output

      virtual void print_times( size_t indent_width ) const ;

   //-- Internals
      
      void augment_system( double augmentation_parameter )  ;
      
      void start_substep( std::string const& title ) const ;
      
   //-- Class attributes
      
      static LA_PenaltyProjection const* PROTOTYPE ;
      
   //-- Attributes
      
      // allocated by the client
      LA_Matrix* A ;
      LA_Matrix* B ;
      LA_Matrix* L ;
      LA_Matrix* M ;
      LA_Vector* F ;
      LA_Vector* G ;
      LA_Vector* S ;
      LA_Vector* K ;
      
      // allocated internally
      LA_Matrix* BtS ;
      LA_Vector* P0 ;
      
      LA_Solver* SOLVER_A ;
      LA_Solver* SOLVER_L ;
      LA_Solver* SOLVER_M ;
      
      double RR ;
      bool PROJP ;
      
      PEL_Timer* TIMER_RR ;
} ; 

#endif
