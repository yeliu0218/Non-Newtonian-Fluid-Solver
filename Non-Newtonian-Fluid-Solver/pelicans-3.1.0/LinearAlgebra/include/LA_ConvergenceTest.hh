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

#ifndef LA_CONVERGENCE_TEST_HH
#define LA_CONVERGENCE_TEST_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class LA_IterativeSolver ;
class LA_Matrix ;
class LA_Preconditioner ;
class LA_Vector ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT LA_ConvergenceTest : public PEL_Object
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      static LA_ConvergenceTest* make( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;
      
      virtual LA_ConvergenceTest* create_clone( 
                                         PEL_Object* a_owner ) const = 0 ;
      
   //-- Characteristics
      
      enum NormType 
      { 
         Preconditioned,
         Unpreconditioned,
         Invalid
      } ;
      
      void set_norm_type( LA_ConvergenceTest::NormType nt ) ;
      
      NormType norm_type( void ) const ;
      
   //-- Convergence
      
      enum ConvergedReason
      {
         ConvergedRtol           =  2, 
         ConvergedAtol           =  3,
         ConvergedIts            =  4,
         ConvergedCGnegCurve     =  5,
         ConvergedCGconstrained  =  6,
         ConvergedStepLength     =  7,
         ConvergedHappyBreakdown =  8,
         
         DivergedNull            = -2,
         DivergedIts             = -3,
         DivergedDtol            = -4,
         DivergedBreakdown       = -5,
         DivergedBreakdownBiCG   = -6,
         DivergedNonSymmetric    = -7,
         DivergedIndefinitePC    = -8,
         DivergedNAN             = -9,
         DivergedIndefiniteMat   = -10,
         DivergedMisc            = -21,
         PreconditionerFailure   = -20 ,
         Undetermined            = -100,

         ConvergedIterating      =  0
      } ;
      
      void set_converged_reason( LA_ConvergenceTest::ConvergedReason cr ) ;
      
      ConvergedReason converged_reason( void ) const ;
      
      virtual void test_convergence( size_t iter,
                                     double r_norm,
                                     LA_Matrix const* A, 
                                     LA_Vector const* b,
                                     LA_Preconditioner* prec,
                                     bool zero_initial_guess ) = 0 ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      virtual void print_more( std::ostream& os, size_t indent_width ) const ;

      friend std::ostream& operator<<(
               std::ostream& out, LA_ConvergenceTest::ConvergedReason cr ) ;
         
   protected: //--------------------------------------------------------
      
   //-- Plug in
      
      virtual ~LA_ConvergenceTest( void ) ;
      
      LA_ConvergenceTest( std::string const& a_name ) ;

      LA_ConvergenceTest( PEL_Object* a_owner ) ;
      
      LA_ConvergenceTest( PEL_Object* a_owner,
                          LA_ConvergenceTest const* other ) ;
      
      virtual LA_ConvergenceTest* create_replica(
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const = 0 ;
      
      bool is_a_prototype( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      bool test_convergence_PRE( size_t iter,
                                 double r_norm,
                                 LA_Matrix const* A, 
                                 LA_Vector const* b,
                                 LA_Preconditioner* prec,
                                 bool zero_initial_guess ) const ;
      
      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      
      bool create_replica_POST( LA_ConvergenceTest const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;
      
   private: //----------------------------------------------------------
      
      LA_ConvergenceTest( void ) ;
      LA_ConvergenceTest( LA_ConvergenceTest const& other ) ;
      LA_ConvergenceTest& operator=( LA_ConvergenceTest const& other ) ;
            
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const IS_PROTO ;
      LA_IterativeSolver const* SOLVER ;
      NormType NORM_TYPE ;
      ConvergedReason REASON ;
} ;

#endif
