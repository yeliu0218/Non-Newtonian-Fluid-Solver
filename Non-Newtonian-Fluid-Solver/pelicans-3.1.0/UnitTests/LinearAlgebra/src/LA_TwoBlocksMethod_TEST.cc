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

#include <LA_TwoBlocksMethod_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_Matrix.hh>
#include <LA_TwoBlocksMethod.hh>
#include <LA_Vector.hh>

#include <iostream>

LA_TwoBlocksMethod_TEST const*
LA_TwoBlocksMethod_TEST:: PROTOTYPE = new LA_TwoBlocksMethod_TEST() ;

//-------------------------------------------------------------------------
LA_TwoBlocksMethod_TEST:: LA_TwoBlocksMethod_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "LA_TwoBlocksMethod", "LA_TwoBlocksMethod_TEST" )
{
}

//-------------------------------------------------------------------------
LA_TwoBlocksMethod_TEST:: ~LA_TwoBlocksMethod_TEST( void )
//-------------------------------------------------------------------------
{
   PROTOTYPE = 0 ;
}

//-------------------------------------------------------------------------
void
LA_TwoBlocksMethod_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod_TEST:: process_one_test" ) ;
   
   PEL_ModuleExplorer const* ee =
                        exp->create_subexplorer( 0, "LA_TwoBlocksMethod" ) ;
   LA_TwoBlocksMethod* solver = LA_TwoBlocksMethod::make( 0, ee ) ;
   ee->destroy() ; ee = 0 ;

   ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   LA_Matrix const* mat_proto = LA_Matrix::make( 0, ee ) ;
   ee->destroy() ; ee = 0 ;

   solver->set_matrix_prototype( mat_proto ) ;
   
   double R = 0. ;
   if( exp->has_entry( "augmentation_coefficient" ) )
   {
      R = exp->double_data( "augmentation_coefficient" ) ;
      exp->test_data( "augmentation_coefficient",
                       "augmentation_coefficient>0." ) ;
   }
   
   std::string const& n= exp->string_data( "test_name" ) ;
   process_test( n, solver, mat_proto, R,
                 exp->string_data( "A" ),
                 exp->string_data( "B" ),
                 exp->string_data( "C" ) ) ;
   
   solver->destroy() ; solver = 0 ;
   mat_proto->destroy() ; mat_proto = 0 ;
}

//-------------------------------------------------------------------------
void
LA_TwoBlocksMethod_TEST:: process_test( std::string const& test_name,
                                        LA_TwoBlocksMethod* solver,
                                        LA_Matrix const* mat_proto,
                                        double R,
                                        std::string const& Amtx,
                                        std::string const& Bmtx,
                                        std::string const& Cmtx )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod_TEST:: process_test" ) ;
   
   LA_Matrix* A = mat_proto->create_matrix( solver ) ;
   A->readMM( Amtx );
   LA_Matrix* B = mat_proto->create_matrix( solver ) ;
   B->readMM( Bmtx ) ;
   LA_Matrix* C = mat_proto->create_matrix( solver ) ;
   C->readMM( Cmtx ) ;
   
   size_t nv_glob = A->nb_rows() ;
   size_t nv_loc  = A->nb_local_rows() ;
   size_t np_glob = B->nb_rows() ;
   size_t np_loc  = B->nb_local_rows() ;
   
   PEL_ASSERT( C->nb_rows() == np_glob ) ;
   PEL_ASSERT( C->nb_cols() == np_glob ) ;
   PEL_ASSERT( C->nb_local_rows() == np_loc ) ;
   PEL_ASSERT( C->nb_local_cols() == np_loc ) ;
   
   solver->re_initialize_internals( nv_glob, np_glob, nv_loc, np_loc ) ;
   
   PEL_ASSERT( !solver->L_is_required() ) ;
   PEL_ASSERT( !solver->MV_is_required() ) ;

   if( solver->S_is_required() )
   {
      LA_Vector* S = mat_proto->create_vector( solver ) ;
      S->re_initialize( np_glob, np_loc ) ;
      S->set( 1.0 ) ;
      solver->set_S( S ) ;
   }

   LA_Vector* f = mat_proto->create_vector( solver ) ;
   f->re_initialize( nv_glob, nv_loc ) ;
   for( size_t i=0 ; i<nv_glob ; ++i )
   {
      double s = ((double) i+1)/((double) nv_glob+1) ;
      f->set_item( i, PEL::sin( s*PEL::pi() ) ) ;
   }
   f->synchronize() ;

   LA_Vector* g = mat_proto->create_vector( solver ) ;
   g->re_initialize( np_glob, np_loc ) ;
   
   LA_Vector* x = mat_proto->create_vector( solver ) ;
   x->re_initialize( nv_glob, nv_loc) ;
   LA_Vector* y = mat_proto->create_vector( solver ) ;
   y->re_initialize( np_glob, np_loc ) ;

   solver->set_system( A, B, f, g, C ) ;
   bool ok = solver->system_is_set() ;
   notify_one_test_result( test_name+" : system_is_set", ok ) ;

   // Solve :
   if( ok )
   {
      solver->estimate_unknowns( false, x, false, y ) ;
      ok = solver->successful_estimation() ;
      notify_one_test_result( test_name+" : successful_estimation", ok ) ;
   }
}
