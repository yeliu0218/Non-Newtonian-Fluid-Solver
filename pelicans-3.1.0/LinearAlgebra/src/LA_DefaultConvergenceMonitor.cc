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

#include <LA_DefaultConvergenceMonitor.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Preconditioner.hh>
#include <LA_Vector.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;

LA_DefaultConvergenceMonitor const*
LA_DefaultConvergenceMonitor:: PROTOTYPE = new LA_DefaultConvergenceMonitor() ;

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor:: LA_DefaultConvergenceMonitor( void )
//---------------------------------------------------------------------------
   : LA_ConvergenceMonitor( "LA_DefaultConvergenceMonitor" )
{
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor*
LA_DefaultConvergenceMonitor:: create_replica( 
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   LA_DefaultConvergenceMonitor* result = 
                                 new LA_DefaultConvergenceMonitor( a_owner ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor*
LA_DefaultConvergenceMonitor:: create( PEL_Object* a_owner ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: create" ) ;
   
   LA_DefaultConvergenceMonitor* result = 
                        new LA_DefaultConvergenceMonitor( a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor:: LA_DefaultConvergenceMonitor( 
                                                    PEL_Object* a_owner )
//---------------------------------------------------------------------------
   : LA_ConvergenceMonitor( a_owner )
{
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor*
LA_DefaultConvergenceMonitor:: create_clone( PEL_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: create_clone" ) ;

   LA_DefaultConvergenceMonitor* result = 
             new LA_DefaultConvergenceMonitor( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor:: LA_DefaultConvergenceMonitor( 
                               PEL_Object* a_owner,
                               LA_DefaultConvergenceMonitor const* other )
//---------------------------------------------------------------------------
   : LA_ConvergenceMonitor( a_owner, other )
{
}

//---------------------------------------------------------------------------
LA_DefaultConvergenceMonitor:: ~LA_DefaultConvergenceMonitor( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
LA_DefaultConvergenceMonitor:: display_at_entry( 
                                       LA_Matrix const* A, 
                                       LA_Vector const* b,
                                       LA_Preconditioner* prec,
                                       bool zero_initial_guess, 
                                       LA_Vector const* x,
                                       LA_ConvergenceTest const* cvgt )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: display_at_entry" ) ;
   PEL_CHECK_PRE( display_at_entry_PRE( A, b, prec, zero_initial_guess, x, cvgt ) ) ;
   
   PEL::out() << "size = " << b->nb_rows() << std::endl ;
   if( zero_initial_guess )
   {
      PEL::out() << "zero initial guess" << std::endl ;
   }
}

//---------------------------------------------------------------------------
void
LA_DefaultConvergenceMonitor:: monitor( size_t iter, double r_norm ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: monitor" ) ;
   
   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;
   
   PEL::out() << setw( 5 ) << iter 
              << setprecision( 6 ) << setw( 15 ) << r_norm << endl ;
   
   PEL::out().flags( original_flags ) ;
}

//---------------------------------------------------------------------------
void
LA_DefaultConvergenceMonitor:: display_at_exit( LA_Matrix const* A, 
                                        LA_Vector const* b,
                                        LA_Preconditioner* prec,
                                        LA_Vector const* x,
                                        LA_ConvergenceTest const* cvgt ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "LA_DefaultConvergenceMonitor:: display_at_exit" ) ;
   PEL_CHECK_PRE( display_at_exit_PRE( A, b, prec, x, cvgt ) ) ;
   
   LA_Vector* aux = x->create_vector( 0 ) ; //??? a chaque fois ???
   
   A->multiply_vec_then_add( x, aux ) ;
   aux->sum( b, -1.0 ) ;
   PEL::out() << "                ||Ax-b||max = " ;
   PEL::out() << aux->max_norm() << std::endl ;
   PEL::out() << "                ||Ax-b||2   = " ;
   PEL::out() << aux->two_norm() << std::endl ;
   
   aux->destroy() ; aux = 0 ;
}
