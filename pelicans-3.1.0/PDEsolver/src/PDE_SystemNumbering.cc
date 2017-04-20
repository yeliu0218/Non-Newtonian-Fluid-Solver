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

#include <PDE_SystemNumbering.hh>

#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_Scatter.hh>
#include <LA_Vector.hh>

#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_CrossProcessUnknownNumbering.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;

struct PDE_SystemNumbering_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
PDE_SystemNumbering*
PDE_SystemNumbering::create( PEL_Object* a_owner,
                             PEL_Vector const* dof2unks,
                             std::string const& ordering,
                             size_t a_verbose_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering::create" ) ;
   PEL_CHECK_PRE( dof2unks != 0 ) ;
   PEL_CHECK_PRE( dof2unks->count() > 0 ) ;
   PEL_CHECK_PRE( dof2unks->count() == dof2unks->index_limit() ) ;
   PEL_CHECK_PRE(
     FORALL( ( size_t i=0 ; i<dof2unks->index_limit() ; ++i ),
       dynamic_cast< PDE_LinkDOF2Unknown* >( dof2unks->at(i) )!=0 ) ) ;
   PEL_CHECK_PRE(
     FORALL( ( size_t i=0 ; i<dof2unks->index_limit() ; ++i ),
       dof2unks->at( i )->owner() == 0  ) ) ;

   PDE_SystemNumbering* result = new PDE_SystemNumbering( a_owner,
                                                          dof2unks,
                                                          ordering,
                                                          a_verbose_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_links() == dof2unks->count() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<result->nb_links() ; ++i ),
              result->link( i ) == dof2unks->at( i ) ) ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<result->nb_links() ; ++i ),
              dof2unks->at( i )->is_under_ownership_of( result ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering:: PDE_SystemNumbering( PEL_Object* a_owner,
                                           PEL_Vector const* dof2unks,
                                           std::string const& ordering,
                                           size_t a_verbose_level )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_LINKS( dof2unks->index_limit() )
   , LINKS( PEL_Vector::create( this, dof2unks->index_limit() ) )
   , SHIFT( 0, 0 )
   , IDX_LOCS( new size_t_vector* [dof2unks->count()] )
   , IDX_GLOBS( new size_t_vector* [dof2unks->count()] )
   , OK_SCATTERS( false )
   , SCATTERS( PEL_Vector::create( this, dof2unks->index_limit() ) )
   , DISTRIBUTED( true )
   , NB_HANDLED_UNK( 0 )
   , MAT_SIZE( PEL::bad_index() )
   , VERB( a_verbose_level )
{
   PEL_LABEL( "PDE_SystemNumbering:: PDE_SystemNumbering" ) ;

   set_ordering_option( ordering  ) ;

   LINKS->copy( dof2unks ) ;

   for( size_t e=0 ; e<NB_LINKS ; ++e )
   {
      PDE_LinkDOF2Unknown* lnk =
          static_cast< PDE_LinkDOF2Unknown* >( dof2unks->at( e ) ) ;
      lnk->set_owner( LINKS ) ;
      LINKS->set_at( e, lnk ) ;
      IDX_LOCS[ e ]  = new size_t_vector( (size_t)0 ) ;
      IDX_GLOBS[ e ] = new size_t_vector( (size_t)0 ) ;

      if( e == 0 )
      {
         DISTRIBUTED = lnk->field()->is_distributed() ;
      }
      else if( lnk->field()->is_distributed() != DISTRIBUTED )
      {
         PDE_SystemNumbering_ERROR::n0() ;
      }
   }

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   NB_HANDLED_UNK.re_initialize( com->nb_ranks(), PEL::bad_index() ) ;

   if( ORDERING == SequenceOfUnknowns )
   {
      size_t nb_glob_unk = PEL::bad_index() ;
      for( size_t e = 0 ; e < nb_links() ; e ++ )
      {
         PDE_LinkDOF2Unknown const* lnk = link( e ) ;
         size_t nn = ( DISTRIBUTED ?
                          lnk->cross_process_numbering()->nb_global_unknowns() :
                          lnk->unknown_vector_size() ) ;
         if( e == 0 )
         {
            nb_glob_unk = nn ;
         }
         else
         {
            if( nn != nb_glob_unk )
            {
               std::ostringstream msg ;
               msg << endl << "*** PDE_SystemNumbering :" << endl ;
               msg << "    When using the option" << endl ;
               msg << "       \"sequence_of_the_unknowns\"" << endl ;
               msg << "    the number of global (cross-process) unknowns" ;
               msg << " have to be the same for all discrete fields" << endl  ;
               msg << "    (check reference elements, imposed DOFs, ...)" ;
               PEL_Error::object()->raise_plain( msg.str() ) ;
            }
         }
      }
   }

   reset_sizes( a_verbose_level ) ;
   reset_numbering() ;
}
//----------------------------------------------------------------------
PDE_SystemNumbering*
PDE_SystemNumbering::create( PEL_Object* a_owner,
                             PDE_LinkDOF2Unknown* dof2unk,
                             size_t a_verbose_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering::create" ) ;
   PEL_CHECK_PRE( dof2unk != 0 ) ;
   PEL_CHECK_PRE( dof2unk->owner() == 0 ) ;

   PDE_SystemNumbering* result = new PDE_SystemNumbering( a_owner,
                                                          dof2unk,
                                                          a_verbose_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_links() == 1 ) ;
   PEL_CHECK_POST( dof2unk->is_under_ownership_of( result) ) ;
   PEL_CHECK_POST( result->link() == dof2unk ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering:: PDE_SystemNumbering( PEL_Object* a_owner,
                                           PDE_LinkDOF2Unknown* dof2unk,
                                           size_t a_verbose_level )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_LINKS( 1 )
   , LINKS( PEL_Vector::create( this, 1 ) )
   , ORDERING( SequenceOfDiscreteFields )
   , SHIFT( 0, 0 )
   , IDX_LOCS( new size_t_vector* [ 1 ] )
   , IDX_GLOBS( new size_t_vector* [ 1 ] )
   , OK_SCATTERS( false )
   , SCATTERS( PEL_Vector::create( this, 1 ) )
   , DISTRIBUTED( true )
   , NB_HANDLED_UNK( 0 )
   , MAT_SIZE( PEL::bad_index() )
   , VERB( a_verbose_level )
{
   PEL_LABEL( "PDE_SystemNumbering:: PDE_SystemNumbering" ) ;

   dof2unk->set_owner( LINKS ) ;
   LINKS->set_at( 0, dof2unk ) ;

   IDX_LOCS[ 0 ]  = new size_t_vector( (size_t)0 ) ;
   IDX_GLOBS[ 0 ] = new size_t_vector( (size_t)0 ) ;

   DISTRIBUTED = dof2unk->field()->is_distributed() ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   NB_HANDLED_UNK.re_initialize( com->nb_ranks(), PEL::bad_index() ) ;

   reset_sizes( a_verbose_level ) ;
   reset_numbering() ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering*
PDE_SystemNumbering:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: create_clone" ) ;

   PDE_SystemNumbering* result = new PDE_SystemNumbering( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SystemNumbering:: PDE_SystemNumbering( PEL_Object* a_owner,
                                           PDE_SystemNumbering const* other )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_LINKS( other->NB_LINKS )
   , LINKS( PEL_Vector::create( this, other->NB_LINKS ) )
   , ORDERING( other->ORDERING )
   , SHIFT( other->SHIFT )
   , IDX_LOCS( new size_t_vector* [other->NB_LINKS] )
   , IDX_GLOBS( new size_t_vector* [other->NB_LINKS] )
   , OK_SCATTERS( false )                                    //??????????????
   , SCATTERS( PEL_Vector::create( this, other->NB_LINKS ) ) //??????????????
   , DISTRIBUTED( other->DISTRIBUTED )
   , NB_HANDLED_UNK( other->NB_HANDLED_UNK )
   , MAT_SIZE( other->MAT_SIZE )
   , VERB( other->VERB )
{
   PEL_LABEL( "PDE_SystemNumbering:: PDE_SystemNumbering" ) ;

   for( size_t i=0 ; i<NB_LINKS ; ++i )
   {
      PDE_LinkDOF2Unknown* lnk = other->link( i )->create_clone( this ) ;
      LINKS->set_at( i, lnk ) ;
      IDX_LOCS[i]  = new size_t_vector( other->idx_local( i )  ) ;
      IDX_GLOBS[i] = new size_t_vector( other->idx_global( i ) ) ;
   }
}


//----------------------------------------------------------------------
PDE_SystemNumbering:: ~PDE_SystemNumbering( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: ~PDE_SystemNumbering" ) ;

   for( size_t i=0 ; i<NB_LINKS ; ++i )
   {
      delete IDX_LOCS[i] ;
      delete IDX_GLOBS[i] ;
   }
   delete [] IDX_LOCS ;
   delete [] IDX_GLOBS ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_SystemNumbering:: reset( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: reset" ) ;

   if( VERB > 0 ) PEL::out() << "PDE_SystemNumbering::reset" << endl ;

   for( size_t i=0 ; i< nb_links() ; ++i )
   {
      PDE_LinkDOF2Unknown* lnk =
                           static_cast<PDE_LinkDOF2Unknown*>( LINKS->at( i ) ) ;
      lnk->reset() ;
   }

   reset_sizes( VERB ) ;
   reset_numbering() ;

   OK_SCATTERS = false ;

   PEL_CHECK_POST( !scatters_are_defined() ) ;
}


//----------------------------------------------------------------------
void
PDE_SystemNumbering:: reset( boolVector const& observed_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: reset" ) ;
   PEL_CHECK_PRE( nb_links() == 1 ) ;
   PEL_CHECK_PRE( observed_nodes.size() == link()->field()->nb_nodes() ) ;

   PDE_LinkDOF2Unknown* lnk =
                        static_cast<PDE_LinkDOF2Unknown*>( LINKS->at( 0 ) ) ;
   lnk->reset( observed_nodes ) ;

   reset_sizes( 0 ) ;
   reset_numbering() ;

   OK_SCATTERS = false ;

   PEL_CHECK_POST( !scatters_are_defined() ) ;
}

//----------------------------------------------------------------------
void
PDE_SystemNumbering:: reset( std::vector<boolVector> const& observed_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: reset" ) ;
   //??? preconditions

   for( size_t i=0 ; i< nb_links() ; ++i )
   {
      PDE_LinkDOF2Unknown* lnk =
                           static_cast<PDE_LinkDOF2Unknown*>( LINKS->at( i ) ) ;
      lnk->reset( observed_nodes[i] ) ;
   }

   reset_sizes( 0 ) ;
   reset_numbering() ;

   OK_SCATTERS = false ;

   PEL_CHECK_POST( !scatters_are_defined() ) ;
}

//----------------------------------------------------------------------
bool
PDE_SystemNumbering:: is_distributed( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: is_distributed" ) ;

   bool result = DISTRIBUTED ;

   PEL_CHECK_PRE(
      FORALL( ( size_t i_link=0 ; i_link<nb_links() ; ++i_link ),
              result == link( i_link )->field()->is_distributed() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: nb_links( void ) const
//----------------------------------------------------------------------
{
   return( NB_LINKS );
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
PDE_SystemNumbering:: link( size_t i_link ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: link" ) ;
   PEL_CHECK_PRE( i_link < nb_links() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_LinkDOF2Unknown const* result =
              static_cast<PDE_LinkDOF2Unknown*>( LINKS->at( i_link ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( IMPLIES( is_distributed(),
                            result->field()->is_distributed() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
PDE_SystemNumbering:: link( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: link" ) ;
   PEL_CHECK_PRE( nb_links() == 1 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_LinkDOF2Unknown const* result =
              static_cast<PDE_LinkDOF2Unknown*>( LINKS->at( 0 ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result == link( 0 ) ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( IMPLIES( is_distributed(),
                            result->field()->is_distributed() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SystemNumbering:: define_scatters( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: define_scatters" ) ;
   PEL_CHECK_PRE( vec->nb_rows() == nb_global_unknowns() ) ;

   for( size_t i=0 ; i<nb_links() ; ++i )
   {
      if( SCATTERS->at( i ) != 0 )
      {
         SCATTERS->destroy_possession( SCATTERS->at( i ) ) ;
      }
      SCATTERS->set_at( i,
         vec->create_scatter( SCATTERS, idx_global( i ), idx_local( i ) ) ) ;
   }
   OK_SCATTERS = true ;

   PEL_CHECK_POST( scatters_are_defined() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i_link=0 ; i_link<nb_links() ; ++i_link ),
         scatter( i_link )->implementation() == vec->implementation() ) ) ;
   //??? add all necessary postconditions to ensure that scatter( i_link )
   //??? can be used as expected
}

//----------------------------------------------------------------------
bool
PDE_SystemNumbering:: scatters_are_defined( void ) const
//----------------------------------------------------------------------
{
   return( OK_SCATTERS ) ;
}

//----------------------------------------------------------------------
LA_Scatter const*
PDE_SystemNumbering:: scatter( size_t i_link ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: scatter" ) ;
   PEL_CHECK_PRE( scatters_are_defined() ) ;
   PEL_CHECK_PRE( i_link < nb_links() ) ;

   LA_Scatter const* result =
                     static_cast< LA_Scatter* >( SCATTERS->at( i_link ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Scatter const*
PDE_SystemNumbering:: scatter( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: scatter" ) ;
   PEL_CHECK_PRE( scatters_are_defined() ) ;
   PEL_CHECK_PRE( nb_links() == 1 ) ;

   LA_Scatter const* result =
                     static_cast< LA_Scatter* >( SCATTERS->at( 0 ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result == scatter( 0 ) ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: nb_global_unknowns( void ) const
//----------------------------------------------------------------------
{
   return( MAT_SIZE );
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: nb_unknowns_on_current_process( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: nb_unknowns_on_current_process" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   return( NB_HANDLED_UNK( com->rank() ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: nb_unknowns_on_process( size_t rank ) const
//----------------------------------------------------------------------
{
   return( NB_HANDLED_UNK( rank ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: global_unknown_for_DOF( size_t n, size_t ic,
                                              size_t i_link ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: global_unknown_for_DOF" ) ;
   PEL_CHECK_PRE( i_link < nb_links() ) ;
   PEL_CHECK_PRE( n  < link( i_link )->nb_field_nodes() ) ;
   PEL_CHECK_PRE( ic < link( i_link )->field()->nb_components() ) ;
   PEL_CHECK_PRE( link( i_link )->DOF_is_unknown( n, ic ) ) ;

   PDE_LinkDOF2Unknown const* lnk = link( i_link ) ;
   size_t idx_ll = lnk->unknown_linked_to_DOF( n, ic ) ;

   size_t owning_rank = 0 ;
   size_t idx_gg = idx_ll ;
   if( DISTRIBUTED )
   {
      PDE_CrossProcessUnknownNumbering const* cpu =
                                       lnk->cross_process_numbering() ;
      owning_rank = cpu->rank_of_process_handling( idx_ll ) ;
      idx_gg = cpu->global_unknown_linked_to_DOF( n, ic ) ;
   }

   size_t result = global_unknown_index( owning_rank, i_link, idx_gg ) ;

   PEL_CHECK_POST( result < nb_global_unknowns() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: global_unknown_for_DOF( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: global_unknown_for_DOF" ) ;
   PEL_CHECK_PRE( nb_links() == 1 ) ;
   PEL_CHECK_PRE( n  < link()->nb_field_nodes() ) ;
   PEL_CHECK_PRE( ic < link()->field()->nb_components() ) ;
   PEL_CHECK_PRE( link()->DOF_is_unknown( n, ic ) ) ;

   size_t result = global_unknown_for_DOF( n, ic, 0 ) ;

   PEL_CHECK_POST( result < nb_global_unknowns() ) ;
   PEL_CHECK_POST( result == global_unknown_for_DOF( n, ic, 0 ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
PDE_SystemNumbering:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: print" ) ;

   os << "nb global unknowns = " << nb_global_unknowns() << endl ;
   for( size_t i_link=0 ; i_link<nb_links() ; ++i_link )
   {
      PDE_LinkDOF2Unknown const* lnk = link( i_link ) ;

      os << "------------------------------" << endl ;
      os << "field: \"" << lnk->field()->name() << "\"" << endl ;

      PDE_CrossProcessNodeNumbering const* cpn =
         ( is_distributed() ? lnk->field()->cross_process_numbering() : 0 ) ;
      PDE_CrossProcessUnknownNumbering const* cun =
         ( is_distributed() ? lnk->cross_process_numbering() : 0 ) ;
      os << "           unknown_vector_size  = "
         << lnk->unknown_vector_size() << endl ;
      if( is_distributed() )
      {
         os << "nb_unknowns_of_current_process  = "
            << cun->nb_unknowns_of_current_process() << endl ;
      }
      os << endl ;

      if( is_distributed() )
      {
         os << setw( 5 )  << "node"
            << setw( 5 )  << "proc"
            << setw( 5 )  << "ic"
            << setw( 8 )  << "unk_l"
            << setw( 8 )  << "unk_g"
            << setw( 10 ) << "sys_unk" ;
      }
      else
      {
         os << setw( 5 )  << "node"
            << setw( 5 )  << "ic"
            << setw( 8 )  << "unk"
            << setw( 10 ) << "sys_unk" ;
      }
      os << endl ;
      for( size_t n=0 ; n<lnk->field()->nb_nodes() ; ++n )
      {
         os << setw( 5 ) << n ;
         if( is_distributed() )
         {
            if( cpn->current_process_handles_node( n ) )
            {
               os << setw( 5 ) << "me" ;
            }
            else
            {
               ostringstream proc ;
               proc << "p" << cpn->rank_of_process_handling( n ) ;
               os << setw( 5 ) << proc.str() ;
            }
         }
         for( size_t iic=0 ; iic<lnk->field()->nb_components() ; ++iic )
         {
            size_t ic = lnk->components_table()( iic ) ;
            if( iic == 0 )
            {
               os << setw( 5 ) << ic ;
            }
            else
            {
               os << setw( 5 + 5 + 5 ) << ic ;
            }
            if( lnk->DOF_is_unknown( n, ic ) )
            {
               os << setw( 8 ) << lnk->unknown_linked_to_DOF( n, ic )
                  << setw( 8 ) ;
               if( is_distributed() )
               {
                  os << cun->global_unknown_linked_to_DOF( n, ic ) ;
               }
               os << setw( 10 )
                  << global_unknown_for_DOF( n, ic, i_link ) ;
            }
            else
            {
               os << "   not an unknown" ;
            }
            os << endl ;
         }
      }
      os << endl ;
   }
   os << "------------------------------" << endl ;
}

//----------------------------------------------------------------------
size_t_vector&
PDE_SystemNumbering:: idx_local( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: idx_local" ) ;
   PEL_CHECK( i<nb_links() ) ;

   size_t_vector* pt = IDX_LOCS[ i ] ;
   PEL_CHECK( pt != 0 ) ;

   return( *pt ) ;
}

//----------------------------------------------------------------------
size_t_vector&
PDE_SystemNumbering:: idx_global( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: idx_global" ) ;
   PEL_CHECK( i<nb_links() ) ;

   size_t_vector* pt = IDX_GLOBS[ i ] ;
   PEL_CHECK( pt != 0 ) ;

   return( *pt ) ;
}

//----------------------------------------------------------------------
void
PDE_SystemNumbering:: reset_sizes( size_t verbose_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: reset_sizes" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   MAT_SIZE = 0 ;
   if( DISTRIBUTED )
   {
      SHIFT.re_initialize( nb_links(), com->nb_ranks(), PEL::bad_index() ) ;
      size_t_array2D temp( nb_links()+1, com->nb_ranks(), 0 ) ;
      for( size_t r=0 ; r<com->nb_ranks() ; ++r )
      {
         temp( 0, r ) = 0 ;

         for( size_t e=0 ; e<nb_links() ; ++e )
         {
            PDE_LinkDOF2Unknown const* lnk = link( e ) ;
            PDE_CrossProcessUnknownNumbering const* cpu =
                                       lnk->cross_process_numbering() ;
            temp( e + 1, r ) = temp( e, r ) + cpu->nb_unknowns_on_process( r ) ;
            if( r == 0 )
            {
               MAT_SIZE       += cpu->nb_global_unknowns() ;
            }
         }
         NB_HANDLED_UNK( r ) = temp( nb_links(), r ) ;
         for( size_t e=0 ; e<nb_links() ; ++e )
         {
           if( r == 0 )
           {
              SHIFT( e, 0 ) = temp( e, 0 ) ;
           }
           else
           {
              SHIFT( e, r ) = SHIFT( e, r - 1 ) + temp( e, r )
                            + temp( nb_links(), r - 1 ) - temp( e + 1, r - 1 ) ;
           }
         }
      }
   }
   else
   {
      SHIFT.re_initialize( nb_links(), 1, PEL::bad_index() ) ;
      for( size_t e=0 ; e<nb_links() ; ++e )
      {
         PDE_LinkDOF2Unknown const* lnk = link( e  ) ;
         SHIFT( e, 0 ) = MAT_SIZE ;
         MAT_SIZE += lnk->unknown_vector_size() ;
      }
   }

   if( verbose_level > 0 )
   {
      if( DISTRIBUTED )
      {
         PEL::out() << "   Unknowns for process " << com->rank()
                    << ": " << nb_unknowns_on_current_process()
                    << " / " << MAT_SIZE
                    << std::endl ;
      }
      else
      {
         PEL::out() << "   Unknowns: " << MAT_SIZE << endl ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_SystemNumbering:: reset_numbering( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: reset_numbering" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t e=0 ; e<nb_links() ; ++e )
   {
      PDE_LinkDOF2Unknown const* lnk = link( e ) ;

      size_t const nb_tot = lnk->unknown_vector_size() ;
      size_t const nbn    = lnk->nb_field_nodes() ;
      size_t const nbc    = lnk->components_table().size() ;

      size_t_vector& idx_gg = idx_global( e ) ;
      size_t_vector& idx_ll = idx_local( e ) ;
      idx_gg.re_initialize( nb_tot, PEL::bad_index() ) ;
      idx_ll.re_initialize( nb_tot, PEL::bad_index() ) ;

      PDE_CrossProcessUnknownNumbering const* cpu =
            ( DISTRIBUTED  ? lnk->cross_process_numbering() : 0 ) ;

      size_t idx_vec = 0 ;
      for( size_t i=0 ; i<nbn ; ++i )
      {
         for( size_t iic=0 ; iic<nbc ; ++iic )
         {
            size_t ic = lnk->components_table()( iic ) ;
            if( lnk->DOF_is_unknown( i, ic ) )
            {
               size_t ll_unk = lnk->unknown_linked_to_DOF( i, ic ) ;
               idx_ll( idx_vec ) = ll_unk ;
               size_t gg_unk = ( DISTRIBUTED ?
                      cpu->global_unknown_linked_to_DOF( i, ic ) : ll_unk ) ;
               size_t rk = ( DISTRIBUTED ?
                      cpu->rank_of_process_handling( ll_unk ) : 0 ) ;
               idx_gg( idx_vec ) = global_unknown_index( rk, e, gg_unk ) ;
               idx_vec++ ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
size_t
PDE_SystemNumbering:: global_unknown_index( size_t rank,
                                            size_t i_link,
                                            size_t idx_unk ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: global_unknown_index" ) ;
   PEL_CHECK( i_link < nb_links() ) ;

   size_t result = PEL::bad_index() ;
   if( ORDERING == SequenceOfDiscreteFields )
   {
     result =  SHIFT( i_link, rank ) + idx_unk ;
   }
   else
   {
     result = idx_unk * nb_links() + i_link ;
   }
   return( result ) ;
}

//-----------------------------------------------------------------------
void
PDE_SystemNumbering:: set_ordering_option( std::string const& option )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering:: set_ordering_option" ) ;

   if( option == "sequence_of_the_discrete_fields" )
   {
      ORDERING = SequenceOfDiscreteFields ;
   }
   else if( option == "sequence_of_the_unknowns" )
   {
      ORDERING =  SequenceOfUnknowns ;
   }
   else
   {
      PEL_Error::object()->raise_plain( option + ": bad value " ) ;
   }
}

//internal--------------------------------------------------------------
void
PDE_SystemNumbering_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PDE_SystemNumbering error:" << endl
        << "    the discrete fields should be" << endl
        << "       * all distributed, or" << endl
        << "       * all not distributed" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

