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

#include <LA_RunSolver.hh>

#include <iostream>

#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Randomizer.hh>
#include <PEL_Timer.hh>
#include <PEL_assertions.hh>

#include <LA_SeqMatrix.hh>
#include <LA_Solver.hh>
#include <LA_SeqVector.hh>

using std::cout ;
using std::endl ;

LA_RunSolver const* LA_RunSolver::PROTOTYPE = new LA_RunSolver() ;

//----------------------------------------------------------------------------
LA_RunSolver:: LA_RunSolver( void )
//----------------------------------------------------------------------------
   : PEL_Application( "LA_RunSolver" )
{
}

//----------------------------------------------------------------------------
LA_RunSolver*
LA_RunSolver:: create_replica( PEL_Object* a_owner, 
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_RunSolver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_RunSolver* result = new LA_RunSolver( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
LA_RunSolver:: LA_RunSolver( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , EXP( exp->create_clone( this ) )
{
   if( exp->has_entry( "datafile" ) )
   {      
      PEL_Module* module =
         PEL_Module::create( this, "MAIN", exp->string_data( "datafile" ) ) ;
      EXP = PEL_ModuleExplorer::create( this,
                                        module->module( "PEL_Application" ) ) ;
   }                 
}

//----------------------------------------------------------------------------
LA_RunSolver:: ~LA_RunSolver( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
LA_RunSolver:: parse_arguments( stringVector& args ) 
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
LA_RunSolver:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "LA_RunSolver:: run" ) ;

   size_t nb = 1 ;
   
   if( EXP->has_entry( "nb_runs" ) )
   {
      nb = EXP->int_data( "nb_runs" )  ;
      PEL::out() << "*** Number of runs for each solving is : " << nb << std::endl ;
   }
   
   PEL_ModuleExplorer* ee = EXP->create_subexplorer( 0, "LA_Matrix" ) ;
   LA_Matrix* mat = LA_Matrix::make( 0, ee ) ;
   ee->destroy() ; ee = 0 ;

   std::string const& nom = EXP->string_data( "lhs_matrix" ) ;
   PEL::out() << "*** LHS matrix:" << endl << "      " << nom << endl << endl ;
   mat->readMM( nom ) ;

   LA_Vector* rhs = mat->create_vector( this ) ;
   LA_Vector* x = rhs->create_vector( this ) ;
   LA_Vector* res = rhs->create_vector( this ) ;
   PEL::out() << "*** RHS vector:" << endl << "      " ;
   if( EXP->has_entry( "rhs_vector" ) )
   {
      std::string const& vnom = EXP->string_data( "rhs_vector" ) ;   
      PEL::out() << vnom ;
      rhs->read( vnom ) ;  
   }
   else
   {
      PEL::out() << "random generation initialized with integer 1" ;
      PEL_Randomizer*r = PEL_Randomizer::create( this, 1 ) ;
      r->start() ;
      for( size_t i=0 ; i<rhs->nb_rows() ; i++ )
      {
         x->set_item( i, r->item() ) ;
         r->go_next() ;
      }
      x->synchronize() ;
      mat->multiply_vec_then_add( x, rhs ) ;
   }
   PEL::out() << endl << endl ;
   
   bool zero_init = true ;
   PEL::out() << "*** initial guess:" << endl << "      " ;
   std::string init_file ;
   if( EXP->has_entry( "initial_guess" ) )
   {
      init_file = EXP->string_data( "initial_guess" ) ;
      PEL::out() << init_file ;
      zero_init = false ;
   }
   else
   {
      PEL::out() << "not given" ;
   }
   PEL::out() << endl << endl ;

   ee = EXP->create_subexplorer( 0, "list_of_LA_Solvers" ) ;
   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      PEL_Timer* timer_set = PEL_Timer::create( 0 ) ;
      PEL_Timer* timer_solve = PEL_Timer::create( 0 ) ;
      PEL_Timer* timer_mult = PEL_Timer::create( 0 ) ;
      
      PEL_ModuleExplorer* se = ee->create_subexplorer( 0 ) ;
      PEL::out() << "*** MODULE " << se->name() << endl << endl ;

      LA_Solver* solver = LA_Solver::make( 0, se ) ;

      if( zero_init )
      {
         x->set( 0.0 ) ;
      }
      else
      {
         x->read( init_file ) ;
      }

      solver->print( PEL::out(), 3 ) ;
      PEL::out() << endl ;

      timer_set->start() ;
      for( size_t i= 0 ; i<nb ; i++ )
      {
         if( nb>1 ) PEL::out() << "    Setting step : " << i << std::endl ;
         if( i!=0 ) solver->unset_matrix() ;
         solver->set_matrix( mat ) ;
      }
      
      timer_set->stop() ;
      
      timer_solve->start() ;
      
      solver->set_initial_guess_nonzero( !zero_init ) ;
      for( size_t i= 0 ; i<nb ; i++ )
      {
         if( nb>1 ) PEL::out() << "    Solving step : " << i << std::endl ;
         solver->solve( rhs, x ) ;
      }
         
      timer_solve->stop() ;
      
      solver->unset_matrix() ;
      
      timer_mult->start() ;
      for( size_t i= 0 ; i<(nb==1 ? 1 : 200*nb) ; i++ )
      {
         mat->multiply_vec_then_add( x, res ) ;
      }
      timer_mult->stop() ;
      
      res->sum( rhs, -1.0 ) ;
      double residu = res->two_norm() / rhs->two_norm() ;
      PEL::out() <<"||Ax-b||2 / ||b||2 = " << res->two_norm() << " / "
                 << rhs->two_norm() << " = " << residu << endl ;
      PEL::out() << endl << endl ;

      PEL::out() << "Matrix setting timing : " ;
      timer_set->print( PEL::out(), 0 ) ;
      PEL::out() << std::endl ;
      
      PEL::out() << "Matrix solving timing : " ;
      timer_solve->print( PEL::out(), 0 ) ;
      PEL::out() << std::endl ;
      
      PEL::out() << "Matrix multipliing timing : " ;
      timer_mult->print( PEL::out(), 0 ) ;
      PEL::out() << std::endl ;
      
      solver->destroy() ;
      se->destroy() ;
      timer_set->destroy() ;
      timer_solve->destroy() ;
      timer_mult->destroy() ;
   }
   ee->destroy() ;

   mat->destroy() ;
}
