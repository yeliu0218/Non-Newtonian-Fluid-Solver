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

#ifndef CFV_NAVIER_STOKES_SYSTEM_HH
#define CFV_NAVIER_STOKES_SYSTEM_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_Timer ;
class size_t_vector ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_TwoBlocksMethod ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_SystemNumbering ;

class FE_TimeIterator ;

/*
PUBLISHED
*/

class AP_NavierStokes1System : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_NavierStokes1System* create( 
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_LinkDOF2Unknown* uu_link,
                                        PDE_LinkDOF2Unknown* pp_link ) ;

      void re_initialize( void ) ;

      bool is_initialized( void ) const ;

   //-- Access

      PDE_LinkDOF2Unknown const* linkDOF2Unknown_U( void ) const ;

      PDE_LinkDOF2Unknown const* linkDOF2Unknown_P( void ) const ;

      LA_SeqVector const* unknown_vector_U( void ) const ;

      LA_SeqVector const* unknown_vector_P( void ) const ;

   //-- Element Change

      size_t global_unknown_for_DOF_of_U( size_t n,
                                          size_t ic ) const ;
      
      size_t global_unknown_for_DOF_of_P( size_t n ) const ;
      
      void set_leading_BDF_over_dt( double value ) ;

      void add_to_A_item( size_t i_row, size_t j_col, double xx ) ;

      void add_to_F_item( size_t i_row, double xx ) ;

      void add_to_B_item( size_t i_row, size_t j_col, double xx ) ;

      void add_to_G_item( size_t i_row, double xx ) ;

      bool MPl_is_required( void ) const ;

      void add_to_MPl_item( size_t i_row, double xx ) ;

      void add_to_C_item( size_t i_row, size_t j_col, double xx ) ;

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

      AP_NavierStokes1System( void ) ;
     ~AP_NavierStokes1System( void ) ;
      AP_NavierStokes1System( AP_NavierStokes1System const& other ) ;
      AP_NavierStokes1System& operator=( AP_NavierStokes1System const& other ) ;

      AP_NavierStokes1System( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp,
                              PDE_LinkDOF2Unknown* uu_link,
                              PDE_LinkDOF2Unknown* pp_link  ) ;

   //-- Internals

      enum MethodType{ AL, MONO, invalid } ;

      void inverseP( void ) ;

      void estimate_unknowns_AL( void ) ;

      void estimate_unknowns_MONO( void ) ;

      void augment_system( double augmentation_parameter )  ;

   //-- Input - Output

      void print_errors_AL( std::string const& indent, size_t n,
                            double err1, double err2, 
                            size_t nit_1, size_t nit_2 ) ;

      void print_end_errors_AL( std::string const& indent, size_t n,
                                double err1, double err2 ) ;

   //-- Attributes

      MethodType METH ;

      bool HAS_INIT_U ;
      bool HAS_INIT_P ;

      double BoverDT ;
      double RR ;
      double DESC ;
      bool CONVERGED ;

      double TOL_VELO ;
      double TOL_DIV ;

      size_t VERBOSE ;
      std::string INDENT ;
      
      bool INIT ;
      
      PDE_SystemNumbering* NMB;
      PDE_SystemNumbering* NMB_U ;
      PDE_SystemNumbering* NMB_P ;
      size_t const idx_U ;
      size_t const idx_P ;
      
      LA_Vector* P ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Matrix* B ;
      LA_Vector* G ;
      LA_Vector* S ;
      LA_Matrix* C ;
      LA_Vector* U ;
      
      LA_SeqVector* U_LOC ;
      LA_SeqVector* P_LOC ;
      
      LA_Solver* SOLVER ;
      
      LA_TwoBlocksMethod* SOLVER_TB ;
} ; 

#endif
