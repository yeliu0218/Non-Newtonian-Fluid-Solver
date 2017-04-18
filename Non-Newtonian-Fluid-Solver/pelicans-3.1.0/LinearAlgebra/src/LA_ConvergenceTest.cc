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

#include <LA_ConvergenceTest.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

#include <LA_Matrix.hh>
#include <LA_Preconditioner.hh>
#include <LA_Vector.hh>

#include <iostream>

using std::endl ;
using std::string ;

//---------------------------------------------------------------------------
LA_ConvergenceTest*
LA_ConvergenceTest:: make( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_ConvergenceTest:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   string name = exp->string_data( "concrete_name" ) ;
   LA_ConvergenceTest const* proto =
      static_cast<LA_ConvergenceTest const*>( plugins_map()->item( name ) ) ;
      
   LA_ConvergenceTest* result = proto->create_replica( a_owner, exp ) ;
      
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceTest:: ~LA_ConvergenceTest( void  )
//---------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
LA_ConvergenceTest:: LA_ConvergenceTest( std::string const& a_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , NORM_TYPE( Invalid )
   , REASON( Undetermined )
{
   PEL_LABEL( "LA_ConvergenceTest:: LA_ConvergenceTest" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceTest:: LA_ConvergenceTest( PEL_Object* a_owner )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NORM_TYPE( Preconditioned )
   , REASON( Undetermined )
{
}

//---------------------------------------------------------------------------
LA_ConvergenceTest:: LA_ConvergenceTest( PEL_Object* a_owner,
                                         LA_ConvergenceTest const* other )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NORM_TYPE( other->NORM_TYPE )
   , REASON( other->REASON )
{
   PEL_LABEL( "LA_ConvergenceTest:: LA_ConvergenceTest" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
void
LA_ConvergenceTest:: set_norm_type( LA_ConvergenceTest::NormType nt )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_ConvergenceTest:: set_norm_type" ) ;
   
   NORM_TYPE = nt ;
   
   PEL_CHECK_POST( norm_type() == nt ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceTest::NormType
LA_ConvergenceTest:: norm_type( void ) const
//---------------------------------------------------------------------------
{
   return( NORM_TYPE ) ;
}

//---------------------------------------------------------------------------
void
LA_ConvergenceTest:: set_converged_reason( 
                         LA_ConvergenceTest::ConvergedReason cr )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_ConvergenceTest:: set_converged_reason" ) ;
   
   REASON = cr ;
   
   PEL_CHECK_POST( converged_reason() == cr ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceTest::ConvergedReason
LA_ConvergenceTest:: converged_reason( void ) const
//---------------------------------------------------------------------------
{
   return( REASON ) ;
}

//----------------------------------------------------------------------
void
LA_ConvergenceTest:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_ConvergenceTest:: print" ) ;
   
   std::string s( indent_width, ' ') ;
   
   os << s << "LA_ConvergenceTest: \"" << type_name() << "\"" << endl ;
   os << s << "   norm type: "   ; 
   if( NORM_TYPE == Preconditioned )
   {
      os << "preconditioned" << endl ;
   }
   else if( NORM_TYPE == Unpreconditioned )
   {
      os << "unpreconditioned" << endl ;
   }
   print_more( os, indent_width+3 ) ;
}

//----------------------------------------------------------------------
void
LA_ConvergenceTest:: print_more( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, LA_ConvergenceTest::ConvergedReason cr )
//----------------------------------------------------------------------
{
   switch( cr )
   {
      case LA_ConvergenceTest::ConvergedRtol :
         out << "\"ConvergedRtol\":\n"
             << "      converged: norm(r) <= rtol*norm(b)" ;
         break ;
      case LA_ConvergenceTest::ConvergedAtol :
         out << "\"ConvergedAtol\":\n"
             << "      converged: norm(r) <= atol" ;
         break ;
      case LA_ConvergenceTest::ConvergedIts :
         out << "\"ConvergedIts\"" ;
         break ;
      case LA_ConvergenceTest::ConvergedCGnegCurve :
         out << "\"ConvergedCGnegCurve\"" ;
         break ;
      case LA_ConvergenceTest::ConvergedCGconstrained :
         out << "\"ConvergedCGconstrained\"" ;
         break ;
      case LA_ConvergenceTest::ConvergedStepLength :
         out << "\"ConvergedStepLength\"" ;
         break ;
      case LA_ConvergenceTest::ConvergedHappyBreakdown :
         out << "\"ConvergedHappyBreakdown:\""
             << "       a breakdown in the Krylov method was detected but the convergence is reached" ;
         break ;
         
      case LA_ConvergenceTest::DivergedNull :
         out << "\"DivergedNull\"" ;
         break ;
      case LA_ConvergenceTest::DivergedIts :
         out << "\"DivergedIts\":\n"
             << "      ran out of iterations before any convergence criteria was reached" ;
         break ;
      case LA_ConvergenceTest::DivergedDtol :
         out << "\"DivergedDtol\":\n"
             << "      diverged: norm(r) >= dtol*norm(b)" ;
         break ;
      case LA_ConvergenceTest::DivergedBreakdown :
         out << "\"DivergedBreakdown\":\n"
             << "       a breakdown in the Krylov method was detected so the method\n"
             << "       could not continue to enlarge the Krylov space" ;
         break ;
      case LA_ConvergenceTest::DivergedBreakdownBiCG :
         out << "\"DivergedBreakdownBiCG\":\n"
             << "       a breakdown in the BICG method was detected so the method\n"
             << "       could not continue to enlarge the Krylov space" ;
         break ;
      case LA_ConvergenceTest::DivergedNonSymmetric :
         out << "\"DivergedNonSymmetric\":\n"
             << "      it appears the operator or preconditioner is not symmetric\n"
             << "      and this Krylov method requires symmetry" ;
         break ;
      case LA_ConvergenceTest::DivergedIndefinitePC :
         out << "\"DivergedIndefinitePC\":\n"
             << "      it appears the preconditioner is indefinite\n"
             << "      (has both positive and negative eigenvalues)\n"
             << "      and this Krylov method requires it to be positive definite" ;
         break ;
      case LA_ConvergenceTest::DivergedNAN :
         out << "\"DivergedNAN\"" ;
         break ;
      case LA_ConvergenceTest::DivergedIndefiniteMat :
         out << "\"DivergedIndefiniteMat\":\n"
             << "      it appears the matrix is indefinite\n"
             << "      (has both positive and negative eigenvalues)\n"
             << "      and this Krylov method requires it to be positive definite" ;
         break ;
      case LA_ConvergenceTest::DivergedMisc :
         out << "\"DivergedMisc\"" ;
         break ;
      case LA_ConvergenceTest::PreconditionerFailure :
         out << "\"PreconditionerFailure\":\n"
             << "      preconditionner has failed" ;
         break ;
      case LA_ConvergenceTest::Undetermined :
         out << "\"Undetermined\"" ;
         break ;
         
      case LA_ConvergenceTest::ConvergedIterating :
         out << "\"ConvergedIterating\":\n"
             << "      solving is still running" ;
         break ;
         
      default:
         PEL_Error::object()-> raise_internal( "Unkwown ConvergedReason" ) ;
   }
   return( out ) ;
}

//---------------------------------------------------------------------------
bool
LA_ConvergenceTest:: is_a_prototype( void ) const
//---------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//---------------------------------------------------------------------------
bool
LA_ConvergenceTest:: test_convergence_PRE( size_t iter,
                                           double r_norm,
                                           LA_Matrix const* A, 
                                           LA_Vector const* b,
                                           LA_Preconditioner* prec,
                                           bool zero_initial_guess ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( r_norm >= 0.0 ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( b != 0 ) ;
   PEL_ASSERT( prec != 0 ) ;
   PEL_ASSERT( A->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == b->nb_rows() ) ; 
   PEL_ASSERT( prec->dimension() == A->nb_rows() ) ;
   PEL_ASSERT( converged_reason() == LA_ConvergenceTest::ConvergedIterating ) ;
   return( true ) ;  
}

//----------------------------------------------------------------------
bool
LA_ConvergenceTest:: create_replica_PRE( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
LA_ConvergenceTest:: create_replica_POST(
                                    LA_ConvergenceTest const* result,
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;  
}

//---------------------------------------------------------------------------
PEL_ObjectRegister*
LA_ConvergenceTest:: plugins_map( void )
//---------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "LA_ConvergenceTest descendant" ) ;
   return( result ) ;
}
