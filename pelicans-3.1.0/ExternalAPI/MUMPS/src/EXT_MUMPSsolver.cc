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

#include <EXT_MUMPSsolver.hh>

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>

#include <LA_DistVector.hh>
#include <LA_SeqVector.hh>

#include <LA_DistMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_SeqMatrix.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

EXT_MUMPSsolver const* EXT_MUMPSsolver::PROTOTYPE = new EXT_MUMPSsolver() ;
static const int USE_COMM_WORLD = -987654 ;
static const int JOB_INIT = -1 ;
static const int JOB_END = -2 ;
static const int JOB_ANALYSIS = 1 ;
static const int JOB_FACTORIZE = 2 ;
static const int JOB_SOLVE = 3 ;
static const int JOB_ANALYSIS_AND_FACTORIZE = 4 ;

//----------------------------------------------------------------------------
EXT_MUMPSsolver:: EXT_MUMPSsolver( void ) 
//----------------------------------------------------------------------------
   : LA_Solver( "EXT_MUMPSsolver" )
   , EXP( 0) 
   , EXTRA_ICNTL(0)
   , EXTRA_CNTL(0)
{   
}

//----------------------------------------------------------------------------
EXT_MUMPSsolver:: ~EXT_MUMPSsolver( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: ~EXT_MUMPSsolver" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   if( matrix_is_set() )
   {
      unset_matrix() ;
   }
}

//----------------------------------------------------------------------------
EXT_MUMPSsolver*
EXT_MUMPSsolver:: create_replica( PEL_Object* a_owner, 
                                  PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE(  a_owner, exp ) ) ;
   
   EXT_MUMPSsolver* result = new EXT_MUMPSsolver( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_MUMPSsolver:: EXT_MUMPSsolver( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner )
   , EXP( exp->create_clone( this ) )
   , EXTRA_ICNTL(40)
   , EXTRA_CNTL(15)
{
   PEL_LABEL( "EXT_MUMPSsolver:: EXT_MUMPSsolver" ) ;
   
   EXTRA_ICNTL.set(PEL::bad_int()) ;
   EXTRA_CNTL.set(PEL::bad_double()) ;
   
   // Output options
   EXTRA_ICNTL(4) = -1 ;
   if( EXP->has_entry( "MUMPS_verbosity" ) )
   {
      EXTRA_ICNTL(4) = EXP->int_data( "MUMPS_verbosity" )  ;
      EXP->test_data( "MUMPS_verbosity", "MUMPS_verbosity<=4 && MUMPS_verbosity>=-1" ) ;
   }
   
   if( EXTRA_ICNTL(4)<0  )
   {
      EXTRA_ICNTL(1) = EXTRA_ICNTL(2) = EXTRA_ICNTL(3) = EXTRA_ICNTL(4) = 0 ;
   }

   if( EXP->has_entry( "out_of_core" ) && EXP->bool_data( "out_of_core" ) )
   {
      EXTRA_ICNTL(22) = 1 ;
   }
   //
   for( size_t i=1 ; i<=40 ; i++ )
   {
      std::ostringstream os ;
      os << "icntl" << i ;
      if( EXP->has_entry( os.str() ) )
      {
         EXTRA_ICNTL(i) = EXP->int_data(os.str()) ;
      }
   }
   
   for( size_t i=1 ; i<=15 ; i++ )
   {
      std::ostringstream os ;
      os << "cntl" << i ;
      if( EXP->has_entry( os.str() ) )
      {
         EXTRA_CNTL(i) = EXP->double_data(os.str())  ;
      }
   }
   
}

//----------------------------------------------------------------------------
EXT_MUMPSsolver*
EXT_MUMPSsolver:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: create_clone" ) ;
   
   EXT_MUMPSsolver* result = new EXT_MUMPSsolver( a_owner, this ) ;

   PEL_CHECK( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_MUMPSsolver:: EXT_MUMPSsolver( PEL_Object* a_owner,
                                   EXT_MUMPSsolver const* other ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , EXP( other->EXP->create_clone( this ) )
   , SOLVER( other->SOLVER )
   , EXTRA_ICNTL(other->EXTRA_ICNTL)
   , EXTRA_CNTL(other->EXTRA_CNTL)
{
   PEL_LABEL( "EXT_MUMPSsolver:: EXT_MUMPSsolver" ) ;
}

//----------------------------------------------------------------------------
void
EXT_MUMPSsolver:: set_matrix_self( LA_Matrix const* mat,
                                   bool &ok, bool same_pattern ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: set_matrix_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   if( !same_pattern )
   {
      initialize() ;
   }
   
   ok  = true ;
   
   // Symmetry
   bool sym = mat->is_symmetric()   ;
   SOLVER.sym = ( sym ? 2 : 0 ) ;

   // Setting of matrix
   if( ok )
   {
      size_t n = mat->nb_rows() ;
      LA_MatrixIterator* it = mat->create_stored_item_iterator( 0 ) ;
      size_t nz = mat->nb_stored_items() ;
         
      int* irn_loc = ( same_pattern ? SOLVER.irn_loc : new int [nz] ) ;
      int* jcn_loc = ( same_pattern ? SOLVER.jcn_loc : new int [nz] ) ;
      double* a_loc = ( same_pattern ? SOLVER.a_loc : new double [nz] ) ;
      size_t k=0 ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         size_t i = it->row()+1 ;
         size_t j = it->col()+1 ;
         if( !sym || i<=j )
         {
            if( same_pattern )
            {
               PEL_CHECK( irn_loc[k] == it->row()+1 ) ;
               PEL_CHECK( jcn_loc[k] == it->col()+1 ) ;
            }
            else
            {
               irn_loc[k] = it->row()+1 ;
               jcn_loc[k] = it->col()+1 ;
            }
            a_loc[k] = it->item() ;
            k++ ;
         }
      }
      SOLVER.n = (int) n ;
      SOLVER.nz_loc = (int) k ;
      
      SOLVER.irn_loc = irn_loc ;
      SOLVER.jcn_loc = jcn_loc ;
      SOLVER.a_loc = a_loc ;
      it->destroy() ;

      PEL_Timer* timer = PEL_Timer::create( 0 ) ;
      if( !same_pattern )
      {
         // Analysis and Factorization
         timer->start() ;
         ok = do_the_job( JOB_ANALYSIS ) ;
         timer->stop() ;
      
         if( is_verbose() )
         {
            PEL::out() << "   [MUMPS] Minimum size required for in-core factorization (Mo) : "
                       << info(15) << " ( max over processes : " << infog(16) <<")"<<std::endl ;
            PEL::out() << "   [MUMPS] Minimum size required for out-of-core factorization (Mo) : "
                       << info(17) << " ( max over processes : " << infog(26) <<")"<<std::endl ;
         
            PEL::out() << "     Analysis time : " ;
            timer->print(PEL::out(), 0) ;
            PEL::out() << std::endl ;
         }
      }
      
      if( ok )
      {
         timer->reset() ;
         timer->start() ;
         ok = do_the_job( JOB_FACTORIZE ) ;
         timer->stop() ;
         if( is_verbose() )
         {
         
            PEL::out() << "     Numeric factorization time : " ;
            timer->print(PEL::out(), 0) ;
            PEL::out() << std::endl ;
         }
      }
      timer->destroy() ;
      
   }
   
   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------------
void
EXT_MUMPSsolver:: unset_matrix_self( void ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: unset_matrix_self" ) ;
   PEL_CHECK( unset_matrix_self_PRE() ) ;
   
   delete [] SOLVER.irn_loc ; SOLVER.irn_loc = 0 ;
   delete [] SOLVER.jcn_loc ; SOLVER.jcn_loc = 0 ;
   delete [] SOLVER.a_loc ; SOLVER.a_loc = 0 ;
   
   // Remove existing MUMPS data structures
   do_the_job( JOB_END ) ;
}

//----------------------------------------------------------------------------
void
EXT_MUMPSsolver:: solve_self( LA_Vector const* b, LA_Vector * x,
                              size_t &nb_iter, bool &ok ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   
   size_t n = b->nb_rows() ;
   size_t_vector id(n) ;
   size_t_vector null(0) ;
   for( size_t i=0 ; i<n ; i++ ) id(i)=i ;
   
   icntl(20) = 0 ; // rhs on host
   icntl(21) = 1 ; // sol distributed
   LA_Scatter* scatter = 0 ;
   LA_SeqVector* seq = 0 ;
   bool i_am_host = PEL_Exec::communicator()->rank()==0 ;
   
   if( i_am_host  )
   {
      scatter = b->create_scatter( 0, id, id ) ;
      seq = LA_SeqVector::create( 0, n ) ;
   }
   else
   {
      scatter = b->create_scatter( 0, null, null ) ;
      seq = LA_SeqVector::create( 0, 0 ) ;
   }
   scatter->get( b, seq ) ;
   if( i_am_host  )
   {
      SOLVER.rhs = seq->data() ;
      SOLVER.nrhs = 1 ;
   }
   size_t nbloc = info(23) ;
   
   SOLVER.sol_loc = new double [ nbloc ] ;
   SOLVER.lsol_loc = (int)nbloc ;
   SOLVER.isol_loc = new int [ nbloc ] ;
   
   // Solve
   ok = do_the_job( JOB_SOLVE ) ;
   for( size_t i=0 ; i<nbloc ; i++ )
      x->set_item( SOLVER.isol_loc[i]-1, SOLVER.sol_loc[i] ) ;
   x->synchronize() ;
   
   delete [] SOLVER.sol_loc ;
   delete [] SOLVER.isol_loc ;
   
   nb_iter = 1 ;

//    scatter->set( seq, x ) ;
    scatter->destroy() ;
    seq->destroy() ;
//     PEL::out() << "Solution : " << std::endl ;
//     x->print( PEL::out(), 0 ) ;
    
   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------------
void
EXT_MUMPSsolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string s( indent_width, ' ' ) ;
   os << s << "MUMPS solver interface"<< std::endl ;
}

//----------------------------------------------------------------------------
int&
EXT_MUMPSsolver:: infog( size_t i ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: icntl" ) ;
   PEL_CHECK_PRE( i>=1 ) ;
   PEL_CHECK_PRE( i<=40 ) ;

   return SOLVER.infog[i-1] ;
}

//----------------------------------------------------------------------------
int&
EXT_MUMPSsolver:: info( size_t i ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: icntl" ) ;
   PEL_CHECK_PRE( i>=1 ) ;
   PEL_CHECK_PRE( i<=40 ) ;

   return SOLVER.info[i-1] ;
}

//----------------------------------------------------------------------------
int&
EXT_MUMPSsolver:: icntl( size_t i ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: icntl" ) ;
   PEL_CHECK_PRE( i>=1 ) ;
   PEL_CHECK_PRE( i<=40 ) ;

   return SOLVER.icntl[i-1] ;
}

//----------------------------------------------------------------------------
double&
EXT_MUMPSsolver:: cntl( size_t i ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: cntl" ) ;
   PEL_CHECK_PRE( i>=1 ) ;
   PEL_CHECK_PRE( i<=15 ) ;

   return SOLVER.cntl[i-1] ;
}

//----------------------------------------------------------------------------
bool
EXT_MUMPSsolver:: do_the_job( int job ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: do" ) ;
   SOLVER.job = job ;
   dmumps_c(&SOLVER) ;
   bool result = infog(1)>=0 ;
   if( !result )
   {
      std::ostringstream msg ;
      msg << "Call to dmumps_c failed on job \'" << job << "\', error code is "<<
         infog(1) << std::endl ;

      switch(infog(1)) 
      {
         case -2 : msg << " -> NZ out of range : " << infog(2) << std::endl ; break ;
         case -3 : msg << " -> invalid job number : " << job << std::endl ; break ;
         case -10 : msg << " -> numerically singular matrix " << std::endl ; break ;
            
         default : break ;
      }
      PEL_Error::object()->display_info( msg.str() ) ;
   }
   return result ;
}

//----------------------------------------------------------------------------
void
EXT_MUMPSsolver:: initialize( void ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MUMPSsolver:: initialize" ) ;
      // Initialization
   SOLVER.par = 1 ; // host 0 participate to assembling
   SOLVER.comm_fortran=MPI_Comm_c2f(MPI_COMM_WORLD) ;
   PEL_ASSERT( do_the_job( JOB_INIT ) ) ;

   // Matrix input
   icntl(5) = 0 ; icntl(18) = 3 ;  // Distributed assembly

   // datadeck option
   for( size_t i=1 ; i<EXTRA_ICNTL.size() ; i++ )
      if( EXTRA_ICNTL(i)!=PEL::bad_int() )
      {
         if( is_verbose() ) PEL::out() << "     ....Setting icntl("<<i<<")="<<EXTRA_ICNTL(i)<<std::endl ;
         icntl(i) = EXTRA_ICNTL(i) ;
      }
   
   for( size_t i=1 ; i<EXTRA_CNTL.size() ; i++ )
      if( EXTRA_CNTL(i)!=PEL::bad_double() )
      {
         cntl(i) = EXTRA_CNTL(i) ;
         if( is_verbose() ) PEL::out() << "     ....Setting cntl("<<i<<")="<<EXTRA_CNTL(i)<<std::endl ;
      }
}



