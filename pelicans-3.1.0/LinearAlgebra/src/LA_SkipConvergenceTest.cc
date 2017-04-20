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

#include <LA_SkipConvergenceTest.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <LA_Preconditioner.hh>
#include <LA_Vector.hh>

#include <iostream>

using std::endl ;

LA_SkipConvergenceTest const*
LA_SkipConvergenceTest:: PROTOTYPE = new LA_SkipConvergenceTest() ;

//---------------------------------------------------------------------------
LA_SkipConvergenceTest:: LA_SkipConvergenceTest( void )
//---------------------------------------------------------------------------
   : LA_ConvergenceTest( "LA_SkipConvergenceTest" )
{
}

//---------------------------------------------------------------------------
LA_SkipConvergenceTest*
LA_SkipConvergenceTest:: create( PEL_Object* a_owner,
                                 size_t a_nb_iterations )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SkipConvergenceTest:: create" ) ;
   
   LA_SkipConvergenceTest* result = 
                        new LA_SkipConvergenceTest( a_owner, a_nb_iterations ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SkipConvergenceTest:: LA_SkipConvergenceTest( 
                                      PEL_Object* a_owner,
                                      size_t a_nb_iterations )
//---------------------------------------------------------------------------
   : LA_ConvergenceTest( a_owner )
   , MAXITS( a_nb_iterations )
{
}
                                      
//---------------------------------------------------------------------------
LA_SkipConvergenceTest*
LA_SkipConvergenceTest:: create_replica( 
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SkipConvergenceTest:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_SkipConvergenceTest* result = 
                           new LA_SkipConvergenceTest( a_owner, exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SkipConvergenceTest:: LA_SkipConvergenceTest( 
                                      PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : LA_ConvergenceTest( a_owner )
   , MAXITS( exp->int_data( "nb_iterations" ) )
{
}
                                      
//---------------------------------------------------------------------------
LA_SkipConvergenceTest*
LA_SkipConvergenceTest:: create_clone( PEL_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SkipConvergenceTest:: create_clone" ) ;
   
   LA_SkipConvergenceTest* result = new LA_SkipConvergenceTest( a_owner, 
                                                                this ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SkipConvergenceTest:: LA_SkipConvergenceTest( 
                                      PEL_Object* a_owner,
                                      LA_SkipConvergenceTest const* other )
//---------------------------------------------------------------------------
   : LA_ConvergenceTest( a_owner, other )
   , MAXITS( other->MAXITS )
{
}
                                      
//---------------------------------------------------------------------------
LA_SkipConvergenceTest:: ~LA_SkipConvergenceTest( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
LA_SkipConvergenceTest:: test_convergence( size_t iter, 
                                           double r_norm,
                                           LA_Matrix const* A, 
                                           LA_Vector const* b,
                                           LA_Preconditioner* prec,
                                           bool zero_initial_guess )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SkipConvergenceTest:: test_convergence" ) ;
   PEL_CHECK_PRE( test_convergence_PRE( iter, r_norm, A, b, prec, 
                                        zero_initial_guess ) ) ;   
   
   PEL_ASSERT( iter <= MAXITS ) ;
   if( iter == MAXITS )
   {
      set_converged_reason( ConvergedIts ) ;
   }
}

//---------------------------------------------------------------------------
void
LA_SkipConvergenceTest:: print_more( std::ostream& os, 
                                     size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_SkipConvergenceTest:: print_more" ) ;
   
   std::string s( indent_width, ' ') ;
   
   os << s << "nb iterations: " << MAXITS << endl ;
}
