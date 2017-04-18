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

#ifndef ML_NavierStokesSystem_HH
#define ML_NavierStokesSystem_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_Timer ;
class size_t_vector ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_SystemNumbering ;

class FE_TimeIterator ;

/** \class ML_NavierStokesSystem
    \brief
 *  Adjusted version of 
 *  pelicans-3.1.0/ExamplesOfApplication/FSI/src/AP_FSINavierStokesSystem
 * 
 * Assembling and algebraic system of the discrete Navier-Stokes problem and
   solving of the discrete mixed system
   \f$
     \left( \begin{array}{cc} A & B^T \\ B & 0\end{array}\right) 
     \left( \begin{array}{c} U \\ P \end{array}\right) = 
     \left( \begin{array}{c} F \\ G \end{array}\right)
    \f$ 
 * 
 * 
 * Unknowns: U(0),P(1),T(2)
 * 
PUBLISHED
*/

class ML_NavierStokesSystem : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static ML_NavierStokesSystem* create( 
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_LinkDOF2Unknown* uu_link,
                                        PDE_LinkDOF2Unknown* pp_link,
                                        PDE_LinkDOF2Unknown* tt_link,
                                        PDE_LinkDOF2Unknown* gamma_link ) ;

      void re_initialize( void ) ;

      bool is_initialized( void ) const ;
      
      void set_calculate_LHS( bool );

   //-- Access

      PDE_SystemNumbering const* system_numbering_U( void ) const ;

      PDE_SystemNumbering const* system_numbering_P( void ) const ;
      
      PDE_SystemNumbering const* system_numbering_T( void ) const ;
      
      PDE_SystemNumbering const* system_numbering_GAMMA( void ) const ;   

      LA_SeqVector const* unknown_vector_U( void ) const ;

      LA_SeqVector const* unknown_vector_P( void ) const ;
      
      LA_SeqVector const* unknown_vector_T( void ) const ;

   //-- Element Change

      void nullify_LHS_and_RHS( void ) ;

      void nullify_for_new_internal_iteration( void ) ;

      void nullify_A_F_Inner( void ) ;

      void add_explicit_terms( void ) ;
      void add_Outer_terms( void ) ;
      void add_Inner_terms( void ) ;
      void add_Flowrate_terms( double inner_rhs_scale ) ;
      void add_BCstress_terms( void ) ;
      void add_Visco_terms( void ) ;
      
      void set_leading_BDF_over_dt( double value ) ;

      void assemble_A_F( PDE_LocalEquation const* leq ) ;
      
      void assemble_A_F_explicit( PDE_LocalEquation const* leq ) ;
      void assemble_A_F_Outer( PDE_LocalEquation const* leq ) ;
      void assemble_A_F_Inner( PDE_LocalEquation const* leq ) ;
      void assemble_F_Flowrate( PDE_LocalEquation const* leq ) ;
      void assemble_F_Stress( PDE_LocalEquation const* leq ) ;

      void assemble_B_G( PDE_LocalEquation const* leq ) ;

      void assemble_A_LagrangeMultDivMatrix( PDE_LocalEquation const* leq );

      bool MPl_is_required( void ) const ;

      void assemble_MPl( PDE_LocalEquation const* leq ) ;

      bool L_is_required( void ) const ;

      void assemble_L( PDE_LocalEquation const* leq ) ;

      bool MV_is_required( void ) const ;

      void assemble_MV( PDE_LocalEquation const* leq ) ;
      
      void set_TminusRGAMMA( double  Lagrange_parameter );

   //-- Solution

      void set_initial_guess_U( PDE_DiscreteField const* uu,
                                size_t level ) ;

      bool initial_guess_U_is_set( void ) const ;

      void set_initial_guess_P( PDE_DiscreteField const* pp,
                                size_t level ) ;

      bool initial_guess_P_is_set( void ) const ;

      void estimate_unknowns( void ) ;

      bool unknowns_are_solution( void ) const ;

   //-- Input - Output

      void set_indent( std::string const& indent ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      ML_NavierStokesSystem( void ) ;
     ~ML_NavierStokesSystem( void ) ;
      ML_NavierStokesSystem( ML_NavierStokesSystem const& other ) ;
      ML_NavierStokesSystem& operator=( 
                                ML_NavierStokesSystem const& other ) ;

      ML_NavierStokesSystem( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp,
                                PDE_LinkDOF2Unknown* uu_link,
                                PDE_LinkDOF2Unknown* pp_link,
                                PDE_LinkDOF2Unknown* tt_link,
                                PDE_LinkDOF2Unknown* gamma_link  ) ;

   //-- Internals

      enum MethodType
      { 
         AL, 
         YOS, 
         PP, 
         invalid 
      } ;

      void augment_system( double augmentation_parameter, bool set_LHS )  ;

      void inverseP( void ) ;

      void estimate_unknowns_AL( void ) ;

      void estimate_unknowns_YOS( void ) ;

      void estimate_unknowns_PP( void ) ;

      void re_initialize_matrices_and_vectors( size_t nv_glob,
                                               size_t np_glob,
                                               size_t nt_glob,
                                               size_t nv_loc,
                                               size_t np_loc,
                                               size_t nt_loc ) ;
      
   //-- Input - Output

      void print_errors_AL( std::string const& indent, size_t n,
                            double err1, double err2, size_t nb_it ) ;

      void print_end_errors_AL( std::string const& indent, size_t n,
                                double err1, double err2 ) ;

   //-- Attributes

      MethodType METH ;

      bool HAS_INIT_U ;
      bool HAS_INIT_P ;
      LA_SeqVector* U_LA ;
      LA_SeqVector* P_LA ;

      double BoverDT ;
      double RR ;
      bool CONVERGED ;

      double TOL_VELO ;
      double TOL_DIV ;

      size_t VERBOSE ;
      std::string INDENT ;

      bool INIT ;
      
      PDE_SystemNumbering* NMB_U ;
      PDE_SystemNumbering* NMB_P ;
      PDE_SystemNumbering* NMB_T ;
      PDE_SystemNumbering* NMB_GAMMA ;
      
      //??? remettre dans l'ordre
      LA_Vector* P ;
      LA_Matrix* A ;
      LA_Matrix* A_explicit ;
      LA_Matrix* M ;
      LA_Vector* F0 ;
      LA_Vector* F ;
      LA_Vector* F_explicit ;
      LA_Matrix* BtMpInv ;
      LA_Matrix* B ;
      LA_Vector* G ;
      LA_Vector* Mpl ;
      LA_Matrix* L ;
      LA_Vector* U ;
      LA_Vector* DU ;
      LA_Vector* P0 ;
      LA_Vector* PRES ;
      
      LA_Vector* T ;
      LA_Matrix* A_TtimesGAMMADOTv ;
      LA_Vector* TminusRGAMMA ;
      
      LA_Matrix* A_Outer ;
      LA_Matrix* A_Inner ;
      LA_Vector* F_Outer ; 
      LA_Vector* F_Inner ;
      LA_Vector* F_Flowrate ;    
      LA_Vector* F_Stress ;     
      LA_Vector* F_Visco ;
      
      LA_SeqVector* U_LOC ;
      LA_SeqVector* P_LOC ;
      LA_SeqVector* T_LOC ;
      
      LA_Solver* SOLVER_A ;
      LA_Solver* SOLVER_L ;
      LA_Solver* SOLVER_Mv ;
      
      bool MP_INVERTED ;
      bool LHS_flag ;
} ; 

#endif
