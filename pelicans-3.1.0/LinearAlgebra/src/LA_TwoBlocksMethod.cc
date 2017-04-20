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

#include <LA_TwoBlocksMethod.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Vector.hh>

#include <sstream>

using std::endl ; 
using std::ostringstream ;
using std::string ;

//----------------------------------------------------------------------
LA_TwoBlocksMethod:: LA_TwoBlocksMethod( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , DIST_STRAT( LA::InvalidDistribution )
   , IMPL( 0 )
{
   PEL_LABEL( "LA_TwoBlocksMethod:: LA_TwoBlocksMethod" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
LA_TwoBlocksMethod*
LA_TwoBlocksMethod:: make( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   LA_TwoBlocksMethod const* proto =
      static_cast<LA_TwoBlocksMethod const*>( plugins_map()->item( name ) ) ;
      
   LA_TwoBlocksMethod* result = proto->create_replica( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->distribution_strategy() == 
                   LA::InvalidDistribution ) ;
   PEL_CHECK_POST( result->implementation() == 0 ) ;
   PEL_CHECK_POST( result->nb_local_rows_U() == 0 ) ;
   PEL_CHECK_POST( result->nb_local_rows_P() == 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_TwoBlocksMethod:: LA_TwoBlocksMethod( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , DIST_STRAT( LA::InvalidDistribution )
   , IMPL( 0 )
   , NB_LOC_U( 0 )
   , NB_LOC_P( 0 )
   , IS_SET( false )
   , SUCCESS( false )
   , VERBOSE( exp->has_entry( "verbose_level" ) ?
                    exp->int_data( "verbose_level" ) : 0 )
   , BoverDT( PEL::bad_double() )
{
}

//----------------------------------------------------------------------
LA_TwoBlocksMethod:: ~LA_TwoBlocksMethod( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_matrix_prototype( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_matrix_prototype" ) ;
   PEL_CHECK_PRE( distribution_strategy() == LA::InvalidDistribution ) ;
   PEL_CHECK_PRE( !system_is_set() ) ;
   PEL_CHECK_PRE( mat != 0 ) ;
   PEL_CHECK_PRE( mat->is_resizable() ) ;
   PEL_CHECK_PRE( !mat->is_symmetric() ) ;
   
   DIST_STRAT = mat->distribution_strategy() ;
   IMPL = mat->implementation() ;
   
   set_matrix_prototype_sub( mat ) ;
   
   PEL_CHECK_POST( distribution_strategy() == mat->distribution_strategy() ) ;
   PEL_CHECK_POST( implementation() == mat->implementation() ) ;
}

//----------------------------------------------------------------------
LA::DistributionStrategy
LA_TwoBlocksMethod:: distribution_strategy( void ) const
//----------------------------------------------------------------------
{
   return( DIST_STRAT ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_TwoBlocksMethod:: implementation( void ) const
//----------------------------------------------------------------------
{
   return( IMPL ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: re_initialize_internals( size_t nv_glob, 
                                              size_t np_glob,
                                              size_t nv_loc, 
                                              size_t np_loc )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: re_initialize_internals" ) ;
   //   PEL_CHECK_PRE( IMPLIES( distribution_strategy() != LA::FromLocalSize,
   //                        nv_glob != PEL::bad_index() ) ) ;
   //   PEL_CHECK_PRE( IMPLIES( distribution_strategy() != LA::FromLocalSize,
   //                        np_glob != PEL::bad_index() ) ) ;
   //   
   //   PEL_CHECK_PRE( IMPLIES( distribution_strategy() == LA::FromLocalSize,
   //                        nv_loc != PEL::bad_index() ) ) ;
   //   PEL_CHECK_PRE( IMPLIES( distribution_strategy() == LA::FromLocalSize,
   //                        np_loc != PEL::bad_index() ) ) ;
   PEL_CHECK_PRE( distribution_strategy() != LA::InvalidDistribution ) ;
   PEL_CHECK_PRE( implementation() != 0 ) ;
   PEL_CHECK_PRE( !system_is_set() ) ;
   PEL_CHECK_PRE( nv_glob != PEL::bad_index() ) ;
   PEL_CHECK_PRE( np_glob != PEL::bad_index() ) ;
   PEL_CHECK_PRE( nv_loc  != PEL::bad_index() ) ;
   PEL_CHECK_PRE( np_loc  != PEL::bad_index() ) ;
   PEL_CHECK_PRE(
        IMPLIES( distribution_strategy() == LA::FromGlobalSize,
           PEL_Exec::communicator()->same_value_everywhere( (int) nv_glob ) ) ) ;
   PEL_CHECK_PRE(
        IMPLIES( distribution_strategy() == LA::FromGlobalSize,
           PEL_Exec::communicator()->same_value_everywhere( (int) np_glob ) ) ) ;
   
   re_initialize_internals_sub( nv_glob, np_glob, nv_loc, np_loc,
                                NB_LOC_U, NB_LOC_P ) ;
}

//----------------------------------------------------------------------
size_t
LA_TwoBlocksMethod:: nb_local_rows_U( void ) const
//----------------------------------------------------------------------
{
   return( NB_LOC_U ) ;
}

//----------------------------------------------------------------------
size_t
LA_TwoBlocksMethod:: nb_local_rows_P( void ) const
//----------------------------------------------------------------------
{
   return( NB_LOC_P ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: dtinv_is_required( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_dtinv( double value )
//----------------------------------------------------------------------
{
   BoverDT = value ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: S_is_required( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_S( LA_Vector* a_S )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_S" ) ;
   PEL_CHECK_PRE( set_S_PRE( a_S ) ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: L_is_required( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_L( LA_Matrix* a_L )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_L" ) ;
   PEL_CHECK_PRE( set_L_PRE( a_L ) ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: K_is_required( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_K( LA_Vector* a_K )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_K" ) ;
   PEL_CHECK_PRE( set_K_PRE( a_K ) ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: MV_is_required( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_MV( LA_Matrix* a_MV )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_MV" ) ;
   PEL_CHECK_PRE( set_MV_PRE( a_MV ) ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_system( LA_Matrix* a_A, LA_Matrix* a_B,
                                 LA_Vector* a_F, LA_Vector* a_G,
                                 LA_Matrix* a_C )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_system" ) ;
   PEL_CHECK_POST( nb_local_rows_U() != 0 ) ;
   PEL_CHECK_POST( nb_local_rows_P() != 0 ) ;
   PEL_CHECK_PRE( !system_is_set() ) ;
   PEL_CHECK_PRE( a_A != 0 ) ;
   PEL_CHECK_PRE( a_A->implementation() == implementation() ) ;
   PEL_CHECK_PRE( a_A->distribution_strategy() == distribution_strategy() ) ;
   PEL_CHECK_PRE( a_A->is_synchronized() ) ;
   PEL_CHECK_PRE( a_A->nb_local_rows() == nb_local_rows_U() ) ;
   PEL_CHECK_PRE( a_A->nb_local_cols() == nb_local_rows_U() ) ;
   PEL_CHECK_PRE( a_B != 0 ) ;
   PEL_CHECK_PRE( a_B->implementation() == implementation() ) ;
   PEL_CHECK_PRE( a_B->distribution_strategy() == distribution_strategy() ) ;
   PEL_CHECK_PRE( a_B->is_synchronized() ) ;
   PEL_CHECK_PRE( a_B->nb_local_rows() == nb_local_rows_P() ) ;
   PEL_CHECK_PRE( a_B->nb_local_cols() == nb_local_rows_U() ) ;
   PEL_CHECK_PRE( a_F != 0 ) ;
   PEL_CHECK_PRE( a_F->implementation() == implementation() ) ;
   PEL_CHECK_PRE( a_F->distribution_strategy() == distribution_strategy() ) ;
   PEL_CHECK_PRE( a_F->is_synchronized() ) ;
   PEL_CHECK_PRE( a_F->nb_local_rows() == nb_local_rows_U() ) ;
   PEL_CHECK_PRE( a_G != 0 ) ;
   PEL_CHECK_PRE( a_G->implementation() == implementation() ) ;
   PEL_CHECK_PRE( a_G->distribution_strategy() == distribution_strategy() ) ;
   PEL_CHECK_PRE( a_G->is_synchronized() ) ;
   PEL_CHECK_PRE( a_G->nb_local_rows() == nb_local_rows_P() ) ;
   PEL_CHECK_PRE( IMPLIES( a_C != 0 ,
                  a_C->implementation() == implementation() ) ) ;
   PEL_CHECK_PRE( IMPLIES( a_C != 0,
                  a_C->distribution_strategy() == distribution_strategy() ) ) ;
   PEL_CHECK_PRE( IMPLIES( a_C != 0,
                  a_C->is_synchronized() ) ) ;
   PEL_CHECK_PRE( IMPLIES( a_C != 0,
                  a_C->nb_local_rows() == nb_local_rows_P() ) ) ;
   PEL_CHECK_PRE( IMPLIES( a_C != 0,
                  a_C->nb_local_cols() == nb_local_rows_P() ) ) ;

   set_system_sub( a_A, a_B, a_F, a_G, a_C ) ;
   
   IS_SET = true ;
   
   PEL_CHECK_POST( system_is_set() ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: system_is_set( void ) const
//----------------------------------------------------------------------
{
   return( IS_SET ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: unset_system( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: unset_system" ) ;

   unset_system_sub() ;
   
   IS_SET = false ;
   NB_LOC_U = 0 ;
   NB_LOC_P = 0 ;
   
   PEL_CHECK_POST( nb_local_rows_U() == 0 ) ;
   PEL_CHECK_POST( nb_local_rows_P() == 0 ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: estimate_unknowns( bool has_init_U, LA_Vector* U, 
                                        bool has_init_P, LA_Vector* P )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: estimate_unknowns" ) ;
   PEL_CHECK_PRE( system_is_set() ) ;
   PEL_CHECK_PRE( U != 0 ) ;
   PEL_CHECK_PRE( U->is_synchronized() ) ;
   PEL_CHECK_PRE( U->nb_local_rows() == nb_local_rows_U() ) ;
   PEL_CHECK_PRE( P != 0 ) ;
   PEL_CHECK_PRE( P->is_synchronized() ) ;
   PEL_CHECK_PRE( P->nb_local_rows() == nb_local_rows_P() ) ;

   SUCCESS = false ;
   
   estimate_unknowns_sub( has_init_U, U, has_init_P, P ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: successful_estimation( void ) const
//----------------------------------------------------------------------
{
   return( SUCCESS ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: set_indent( std::string const& an_indent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_TwoBlocksMethod:: set_indent" ) ;

   INDENT = an_indent ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: print_times( size_t indent_width ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
double
LA_TwoBlocksMethod:: dtinv( void ) const
//----------------------------------------------------------------------
{
   return( BoverDT ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: check_zero_C( LA_Matrix const* a_C ) const
//----------------------------------------------------------------------
{
   if( a_C != 0 )
   {
      ostringstream mesg ;
      mesg << endl << "*** " << type_name() << ": invalid usage" << endl ;
      mesg << "   cases with nonzero C matrices are not handled" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: raise_invalid_usage( void )
//----------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << endl << "*** " << type_name() << ": invalid usage" << endl ;
   mesg << "   HINT: the following functions must have been called" << endl ;
   mesg << "         prior to \"set_system\":" << endl ;
   if( dtinv_is_required() )
      mesg << "         \"set_leading_BDF_over_dt\"" << endl ;
   if( S_is_required() )
      mesg << "         \"set_S\"" << endl ;
   if( L_is_required() )
      mesg << "         \"set_L\"" << endl ;
   if( K_is_required() )
      mesg << "         \"set_K\"" << endl ;
   if( MV_is_required() )
      mesg << "         \"set_MV\"" << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//----------------------------------------------------------------------
void
LA_TwoBlocksMethod:: notify_success( bool is_successful )
//----------------------------------------------------------------------
{
   SUCCESS = is_successful ;
}

//----------------------------------------------------------------------
size_t
LA_TwoBlocksMethod:: verbose_level( void ) const
//----------------------------------------------------------------------
{
   return( VERBOSE ) ;
}

//----------------------------------------------------------------------
std::string const&
LA_TwoBlocksMethod:: indent( void ) const
//----------------------------------------------------------------------
{
   return( INDENT ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: set_S_PRE( LA_Vector* a_S ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !system_is_set() ) ;
   PEL_ASSERT( S_is_required() ) ;
   PEL_ASSERT( a_S != 0 ) ;
   PEL_ASSERT( a_S->implementation() == implementation() ) ;
   PEL_ASSERT( a_S->distribution_strategy() == distribution_strategy() ) ;
   PEL_ASSERT( a_S->is_synchronized() ) ;
   PEL_ASSERT( a_S->nb_local_rows() == nb_local_rows_P() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: set_L_PRE( LA_Matrix* a_L ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !system_is_set() ) ;
   PEL_ASSERT( L_is_required() ) ;
   PEL_ASSERT( a_L != 0 ) ;
   PEL_ASSERT( a_L->implementation() == implementation() ) ;
   PEL_ASSERT( a_L->distribution_strategy() == distribution_strategy() ) ;
   PEL_ASSERT( a_L->is_synchronized() ) ;
   PEL_ASSERT( a_L->nb_local_rows() == nb_local_rows_P() ) ;
   PEL_ASSERT( a_L->nb_local_cols() == nb_local_rows_P() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: set_K_PRE( LA_Vector* a_K ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !system_is_set() ) ;
   PEL_ASSERT( K_is_required() ) ;
   PEL_ASSERT( a_K != 0 ) ;
   PEL_ASSERT( a_K->implementation() == implementation() ) ;
   PEL_ASSERT( a_K->distribution_strategy() == distribution_strategy() ) ;
   PEL_ASSERT( a_K->is_synchronized() ) ;
   PEL_ASSERT( a_K->nb_local_rows() == nb_local_rows_P() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: set_MV_PRE( LA_Matrix* a_MV ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !system_is_set() ) ;
   PEL_ASSERT( MV_is_required() ) ;
   PEL_ASSERT( a_MV != 0 ) ;
   PEL_ASSERT( a_MV->implementation() == implementation() ) ;
   PEL_ASSERT( a_MV->distribution_strategy() == distribution_strategy() ) ;
   PEL_ASSERT( a_MV->is_synchronized() ) ;
   PEL_ASSERT( a_MV->nb_local_rows() == nb_local_rows_U() ) ;
   PEL_ASSERT( a_MV->nb_local_cols() == nb_local_rows_U() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: create_replica_PRE( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_TwoBlocksMethod:: create_replica_POST( LA_TwoBlocksMethod const* result,
                                          PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->dtinv() == PEL::bad_double() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
LA_TwoBlocksMethod:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "LA_TwoBlocksMethod descendant" ) ;
   return( result ) ;
}

