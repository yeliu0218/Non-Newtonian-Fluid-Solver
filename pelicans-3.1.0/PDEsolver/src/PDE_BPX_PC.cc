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

#include <PDE_BPX_PC.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_LinkDOF2Unknown.hh>

#include <LA_Solver.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::set ;
using std::string ; using std::ostringstream ;

struct PDE_BPX_PC_ERROR
{
   static void n0( void ) ;
} ;

PDE_BPX_PC const* 
PDE_BPX_PC::PROTOTYPE = new PDE_BPX_PC() ;

//----------------------------------------------------------------------
PDE_BPX_PC:: PDE_BPX_PC( void )
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( "PDE_BPX_PC" )
{
}

//----------------------------------------------------------------------
PDE_BPX_PC*
PDE_BPX_PC:: create_replica( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BPX_PC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PDE_BPX_PC* result = new PDE_BPX_PC( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BPX_PC:: PDE_BPX_PC( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( a_owner, exp )
   , I_CYCLE( 0 )
   , C_SOLVER( false )
   , SOLVER( 0 )
   , SOLVE_OK( false )
   , VERBOSE( 0 )
{
   PEL_LABEL( "PDE_BPX_PC:: PDE_BPX_PC" ) ;

   if( exp->has_module( "coarse_solver" ) )
   {
      C_SOLVER = true ;
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "coarse_solver" ) ;
      SOLVER = LA_Solver::make( this, se ) ;
      se->destroy() ; se=0 ;
   }
   
   if( exp->has_entry( "verbose_level" ) )
      VERBOSE = exp->int_data( "verbose_level" ) ;
}

//----------------------------------------------------------------------
PDE_BPX_PC:: ~PDE_BPX_PC( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PDE_BPX_PC* 
PDE_BPX_PC:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BPX_PC:: create_clone" ) ;
   return( new PDE_BPX_PC( a_owner, this ) ) ;
}

//----------------------------------------------------------------------
PDE_BPX_PC:: PDE_BPX_PC( PEL_Object* a_owner, 
                                           PDE_BPX_PC const* other ) 
//----------------------------------------------------------------------
   : PDE_GeometricMultilevel_PC( a_owner, other )
{
   PEL_ASSERT( false ) ; //???????????????????????????
}

//----------------------------------------------------------------------
void
PDE_BPX_PC:: solve( LA_Vector const* rhs, LA_Vector* sol )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BPX_PC:: solve" ) ;
   PEL_CHECK_PRE( solve_PRE( rhs, sol ) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( rhs ) != 0 ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector*>( sol ) != 0 ) ;

   LA_SeqVector const* brhs = static_cast<LA_SeqVector const*>( rhs ) ;
   LA_SeqVector* bsol = static_cast<LA_SeqVector*>( sol ) ;

   //---Projections.

   //projection de rhs sur un niveau plus grossier 
   //                --> rhs_of_level(level_c)
   
   size_t level_c = nb_levels()-2 ;
   LA_SeqVector* bb_c = 
      static_cast<LA_SeqVector*>( rhs_of_level( level_c ) ) ;
   coarse_to_fine( level_c )->tr_multiply_vec_then_add( brhs, bb_c, 
                                                        1.0, 0.0 ) ;
   
   for( size_t level=nb_levels()-2 ; level!=0 ; --level )
   {
      //projection de rhs_of_level( level_f ) sur un niveau plus grossier
      //          --> rhs_of_level( level_c )
      size_t level_f = level ;
      level_c = level-1 ;
      LA_SeqVector* bb_f = 
         static_cast<LA_SeqVector*>( rhs_of_level( level_f ) ) ;
      bb_c =
         static_cast<LA_SeqVector*>( rhs_of_level( level_c ) ) ;
      coarse_to_fine(level_c)->tr_multiply_vec_then_add( bb_f, bb_c,
                                                         1.0, 0.0 ) ;
   }

   //--

   //--Actions des preconditionneurs diagonaux
   
   //Action du preconditionneur diagonal sur le niveau fin nb_levels()-1 
   //Les termes diagonaux sont tous non nuls.
   LA_SeqVector const* to_be_smoothed = 
      static_cast<LA_SeqVector const*>( unknowns_of_level( nb_levels()-2 ) ) ;
   double damped = 4./5. ;
   PEL_ASSERT( to_be_smoothed->nb_rows() == brhs->nb_rows() ) ;
   for( size_t i=0 ; i!=brhs->nb_rows() ; ++i )
   { 
      double xx ;
      if( to_be_smoothed->item( i ) > 0.9 )
      {
         LA_SeqMatrix const* mm = 
             static_cast<LA_SeqMatrix const*>( finest_mat() ) ;
         xx = damped * (brhs->item(i)) / (mm->item(i,i)) ;
         bsol -> set_item( i , xx ) ;
      } 
      else
      {
         bsol -> set_item( i , 0.0 ) ;
      }
   } 

   for( size_t level = 1 ; level!=nb_levels()-1 ; ++level )
   {
      //Action du preconditionneur diagonal sur le niveau level
      //Les termes diagonaux sont tous non nuls.
      LA_SeqVector* bb = 
         static_cast<LA_SeqVector*>( rhs_of_level( level ) ) ;
      to_be_smoothed = 
         static_cast<LA_SeqVector const*>( unknowns_of_level( level-1 ) ) ;
      for(size_t i=0 ; i!=bb->nb_rows() ; ++i )
      {
         //PEL_ASSERT( to_be_smoothed->nb_rows() == bb->nb_rows() ) ;
         if( to_be_smoothed->item( i ) > 0.9 )
         {
            double xx ;
            LA_SeqMatrix* mm = 
               static_cast<LA_SeqMatrix*>( mat_of_level( level ) ) ;
            xx = damped * bb->item( i ) / mm->item( i, i ) ;
            bb -> set_item( i , xx ) ;
         }        
         else
         {
            bb -> set_item( i , 0.0 ) ;
         }
      }
   }
   
//--Sur le niveau grossier
   //--Resolution
   if( !C_SOLVER )
   {
      //--Lissage
      LA_SeqVector* bb =
         static_cast<LA_SeqVector*>( rhs_of_level( 0 ) ) ;
      for(size_t i=0 ; i!=bb->nb_rows() ; ++i )
      {
         double xx ;
         LA_SeqMatrix* mm = 
                        static_cast<LA_SeqMatrix*>( mat_of_level( 0 ) ) ;
         xx = damped * bb->item( i ) / mm->item( i, i ) ;
         bb -> set_item( i ,xx ) ;
      }
   }
   else
   {
      LA_SeqVector* bb = 
          static_cast<LA_SeqVector*>( rhs_of_level( 0 ) ) ;
      LA_SeqVector* dummy = bb -> create_vector( 0 ) ;

      SOLVER->set_matrix( mat_of_level( 0 ) ) ;
      SOLVER->solve( dummy, bb ); 

      dummy->destroy() ;
   }

   //--

   //---Interpolations
   
   for( size_t level=0 ; level!=nb_levels()-2 ; ++level )
   {
      size_t level_f = level+1 ;
      level_c=level ; 
      LA_Vector* bb_f = rhs_of_level( level_f ) ;
      bb_c =  static_cast<LA_SeqVector*>( rhs_of_level( level_c ) ) ;
      coarse_to_fine(level_c)->multiply_vec_then_add( bb_c, bb_f,
                                                      1.0, 1.0 ) ;
   }

   level_c =nb_levels()-2 ;
   bb_c =  static_cast<LA_SeqVector*>(  rhs_of_level( level_c ) ) ;
   coarse_to_fine(level_c)->multiply_vec_then_add( bb_c, bsol, 
                                                   1.0, 1.0 ) ;

   //--

   SOLVE_OK = true ;
   PEL_CHECK_POST( solve_POST( brhs, sol ) ) ;
   PEL_CHECK_POST( successful_solve() ) ;  
}

   
//----------------------------------------------------------------------
bool
PDE_BPX_PC:: successful_solve( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BPX_PC:: successful_solve" ) ;
   return( SOLVE_OK ) ;
}

//----------------------------------------------------------------------
size_t
PDE_BPX_PC:: nb_cycles_performed( void ) const
//----------------------------------------------------------------------
{
   return( I_CYCLE ) ;
}

//----------------------------------------------------------------------
void
PDE_BPX_PC:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BPX_PC:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string s( indent_width, ' ') ;

   os << s << "preconditioner : \"PDE_BPX_PC\"" << std::endl ;
}
