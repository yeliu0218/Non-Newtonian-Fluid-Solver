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

#include <LA_Identity_IS.hh>

#include <LA_Preconditioner.hh>
#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <sstream>

LA_Identity_IS const* LA_Identity_IS::PROTOTYPE = new LA_Identity_IS() ;

struct LA_Identity_IS_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
LA_Identity_IS:: LA_Identity_IS( void )
//----------------------------------------------------------------------
   : LA_IterativeSolver( "LA_Identity_IS" )
{
}

//----------------------------------------------------------------------
LA_Identity_IS*
LA_Identity_IS:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_IS:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_Identity_IS* result = new LA_Identity_IS( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Identity_IS:: LA_Identity_IS( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, exp )
{
}

//----------------------------------------------------------------------
LA_Identity_IS:: LA_Identity_IS( PEL_Object* a_owner, 
                                 LA_Identity_IS const* other )
//----------------------------------------------------------------------
   : LA_IterativeSolver( a_owner, other )
{
}

//----------------------------------------------------------------------
LA_Identity_IS:: ~LA_Identity_IS( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Identity_IS*
LA_Identity_IS:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_IS:: create_clone" ) ;
   return( new LA_Identity_IS( a_owner, this ) ) ;
}

//----------------------------------------------------------------------
void
LA_Identity_IS:: do_solve( LA_Matrix const* A,
                           LA_Vector const* b,
                           LA_Preconditioner* prec,
                           bool zero_initial_guess,
                           LA_Vector* x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Identity_IS:: do_solve" ) ;
   PEL_CHECK_PRE( do_solve_PRE( A, b, prec, zero_initial_guess, x ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   prec->solve( b, x ) ;
   if( prec->successful_solve() )
   {
      set_converged_reason( 0, LA_ConvergenceTest::ConvergedIts ) ;
   }
   else
   {
      set_converged_reason( 0, LA_ConvergenceTest::DivergedMisc ) ;
   }
}
