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

#include <LA_ConvergenceMonitor.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

#include <LA_Matrix.hh>
#include <LA_Preconditioner.hh>
#include <LA_Vector.hh>

using std::string ;

//---------------------------------------------------------------------------
LA_ConvergenceMonitor*
LA_ConvergenceMonitor:: make( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   string name = exp->string_data( "concrete_name" ) ;
   LA_ConvergenceMonitor const* proto =
      static_cast<LA_ConvergenceMonitor const*>( plugins_map()->item( name ) ) ;
      
   LA_ConvergenceMonitor* result = proto->create_replica( a_owner, exp ) ;
      
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceMonitor:: ~LA_ConvergenceMonitor( void  )
//---------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
LA_ConvergenceMonitor:: LA_ConvergenceMonitor( std::string const& a_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "LA_ConvergenceMonitor:: LA_ConvergenceMonitor" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
LA_ConvergenceMonitor:: LA_ConvergenceMonitor( PEL_Object* a_owner )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
}

//---------------------------------------------------------------------------
LA_ConvergenceMonitor:: LA_ConvergenceMonitor( 
                                      PEL_Object* a_owner,
                                      LA_ConvergenceMonitor const* other)
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
   PEL_LABEL( "LA_ConvergenceMonitor:: LA_ConvergenceMonitor" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
bool
LA_ConvergenceMonitor:: is_a_prototype( void ) const
//---------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//---------------------------------------------------------------------------
bool 
LA_ConvergenceMonitor:: display_at_entry_PRE( 
                                   LA_Matrix const* A, 
                                   LA_Vector const* b,
                                   LA_Preconditioner* prec,
                                   bool zero_initial_guess, 
                                   LA_Vector const* x,
                                   LA_ConvergenceTest const* cvgt ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( b != 0 ) ;
   PEL_ASSERT( prec != 0 ) ;
   PEL_ASSERT( x != 0 ) ;
   PEL_ASSERT( A->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == b->nb_rows() ) ; 
   PEL_ASSERT( A->nb_cols() == x->nb_rows() ) ;
   PEL_ASSERT( prec->dimension() == A->nb_rows() ) ;
   PEL_ASSERT( cvgt != 0 ) ;
   return( true ) ;  
}

//---------------------------------------------------------------------------
bool 
LA_ConvergenceMonitor:: display_at_exit_PRE( 
                                   LA_Matrix const* A, 
                                   LA_Vector const* b,
                                   LA_Preconditioner* prec,
                                   LA_Vector const* x,
                                   LA_ConvergenceTest const* cvgt ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( b != 0 ) ;
   PEL_ASSERT( prec != 0 ) ;
   PEL_ASSERT( x != 0 ) ;
   PEL_ASSERT( A->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == b->nb_rows() ) ; 
   PEL_ASSERT( A->nb_cols() == x->nb_rows() ) ;
   PEL_ASSERT( prec->dimension() == A->nb_rows() ) ;
   PEL_ASSERT( cvgt != 0 ) ;
   return( true ) ;  
}

//----------------------------------------------------------------------
bool
LA_ConvergenceMonitor:: create_replica_PRE( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
LA_ConvergenceMonitor:: create_replica_POST(
                               LA_ConvergenceMonitor const* result,
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
LA_ConvergenceMonitor:: plugins_map( void )
//---------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "LA_ConvergenceMonitor descendant" ) ;
   return( result ) ;
}
