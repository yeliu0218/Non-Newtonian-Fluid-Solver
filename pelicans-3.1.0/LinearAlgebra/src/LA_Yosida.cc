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

#include <LA_Yosida.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

LA_Yosida const* LA_Yosida:: PROTOTYPE = new LA_Yosida() ;

//----------------------------------------------------------------------
LA_Yosida:: LA_Yosida( void )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( "LA_Yosida" )
{
}

//----------------------------------------------------------------------
LA_Yosida*
LA_Yosida:: create_replica( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   LA_Yosida* result = new LA_Yosida( a_owner, exp ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Yosida:: LA_Yosida( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : LA_TwoBlocksMethod( a_owner, exp )
   , A( 0 )
   , B( 0 )
   , L( 0 )
   , F( 0 )
   , G( 0 )
   , F0( 0 )
   , P0( 0 )
   , SOLVER_A( 0 )
   , SOLVER_L( 0 )
{
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "solver_A" ) ;
   SOLVER_A = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   ee = exp->create_subexplorer( 0, "solver_L" ) ;
   SOLVER_L = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
LA_Yosida:: ~LA_Yosida( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_Yosida:: dtinv_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Yosida:: L_is_required( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: set_L( LA_Matrix* a_L )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: set_L" ) ;
   PEL_CHECK_PRE( set_L_PRE( a_L ) ) ;
   
   L = a_L ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: set_matrix_prototype_sub( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: set_matrix_prototype_sub" ) ;
   
   PEL_ASSERT( F0 == 0 ) ;
   F0 = mat->create_vector( this ) ;
   
   PEL_ASSERT( P0 == 0 ) ;
   P0 = mat->create_vector( this ) ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: re_initialize_internals_sub( size_t nv_glob, 
                                         size_t np_glob,
                                         size_t nv_loc, 
                                         size_t np_loc,
                                         size_t& nv_loc_final,
                                         size_t& np_loc_final  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: re_initialize_internals_sub" ) ;

   F0->re_initialize( nv_glob, nv_loc ) ;
   P0->re_initialize( np_glob, np_loc ) ;
   
   np_loc_final = P0->nb_local_rows() ;
   nv_loc_final = F0->nb_local_rows() ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                            LA_Vector* a_F, LA_Vector* a_G,
                            LA_Matrix* a_C )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: set_system_sub" ) ;
   
   check_zero_C( a_C ) ;
   
   if( L == 0 ) raise_invalid_usage() ;
   
   A = a_A ;
   B = a_B ;
   F = a_F ;
   G = a_G ;
   
   SOLVER_A->set_matrix( A ) ;
   SOLVER_L->set_matrix( L ) ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: unset_system_sub( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: unset_system_sub" ) ;
   
   SOLVER_A->unset_matrix() ;
   SOLVER_L->unset_matrix() ;
}

//----------------------------------------------------------------------
void
LA_Yosida:: estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                   bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Yosida:: estimate_unknowns_sub" ) ;
   
   if( has_init_P )
   {
      P0->set( P ) ;
   }
   else
   {
      P0->nullify() ;
   }
   
   if( verbose_level() > 0 ) PEL::out() << indent() << "   step 1 " ;
   
   // F0 = F
   F0->set( F ) ;
      
   // F -= BT.P0
   B->tr_multiply_vec_then_add( P0, F, -1.0, 1.0 ) ;
      
   // U = A-1.F
   SOLVER_A->set_initial_guess_nonzero( has_init_U ) ;
   SOLVER_A->solve( F, U ) ;
      
   if( ! SOLVER_A->solution_is_achieved( ) ) return ; // <---

   if( verbose_level() > 0 && SOLVER_A->is_iterative() )
      PEL::out() << ": " << SOLVER_A->nb_iterations_achieved() 
                         << " iterations" ;
   if( verbose_level() > 0 ) PEL::out() << endl << indent() << "   step 2 " ;

   // G = BU - G
   B->multiply_vec_then_add( U, G, 1.0, -1.0 ) ;
      
   // P = L-1.G
   SOLVER_L->set_initial_guess_nonzero( has_init_P ) ;
   SOLVER_L->solve( G, P ) ;
   if( ! SOLVER_L->solution_is_achieved() ) return ; // <---

   // P = P0 + beta*P/DT
   P->scale( dtinv() ) ;
   P->sum( P0 ) ;
         
   if( verbose_level() > 0 && SOLVER_L->is_iterative()  )
            PEL::out() << ": " << SOLVER_L->nb_iterations_achieved()
                               << " iterations" ;
   
   if( verbose_level() > 0 ) PEL::out() << endl << indent() << "   step 3 " ;
   // F = F0
   F->set( F0 ) ;
      
   // F -= BT.P
   B->tr_multiply_vec_then_add( P, F, -1.0, 1.0 ) ;
      
   // U = A-1.F
   SOLVER_A->solve( F, U ) ;
   if( ! SOLVER_A->solution_is_achieved() ) return ; // <---

   if( verbose_level() > 0 && SOLVER_A->is_iterative() )
      PEL::out() << ": " << SOLVER_A->nb_iterations_achieved()
                 << " iterations" ;
   if( verbose_level() > 0 ) PEL::out() << std::endl ;
   
   notify_success( true ) ;
}

