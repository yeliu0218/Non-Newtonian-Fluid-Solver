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

#include <PDE_CrossProcessUnknownNumbering.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_assertions.hh>

#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_LinkDOF2Unknown.hh>

#include <string>
#include <sstream>
#include <iostream>

//----------------------------------------------------------------------
PDE_CrossProcessUnknownNumbering*
PDE_CrossProcessUnknownNumbering:: create( PEL_Object* a_owner,
                                           PDE_LinkDOF2Unknown const* a_link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: create" ) ;
   PEL_CHECK_PRE( a_link->field()->is_distributed() ) ;

   PDE_CrossProcessUnknownNumbering* result =
                    new PDE_CrossProcessUnknownNumbering( a_owner, a_link ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CrossProcessUnknownNumbering:: PDE_CrossProcessUnknownNumbering(
                                           PEL_Object* a_owner,
                                           PDE_LinkDOF2Unknown const* a_link )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , LINK( a_link )
   , GLOB_NODE( a_link->field()->cross_process_numbering() )
   , COMM( 0 )
   , NB_UNKNOWNS_ON_PROC( 0 )
   , NB_UNKNOWNS( PEL::bad_index() )
   , NB_HANDLED_UNK( PEL::bad_index() )
   , LOC_2_GLOB( 0 )
   , OWNER( 0 )
{
   PEL_CHECK_INV( invariant() ) ;

   COMM = GLOB_NODE->communicator() ;

   set_ordering_option( a_link->DOFs_ordering_in_unknown() ) ;

   size_t nb_ranks = COMM->nb_ranks() ;
   NB_UNKNOWNS_ON_PROC.re_initialize( nb_ranks , PEL::bad_int() ) ;
   
   reset() ;
}

//----------------------------------------------------------------------
PDE_CrossProcessUnknownNumbering:: ~PDE_CrossProcessUnknownNumbering( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessUnknownNumbering:: reset( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: reset" ) ;

   size_t const nbcs = LINK->components_table().size() ;
   size_t const nb_nodes = LINK->nb_field_nodes() ;

   // global_unknown(n,ic)= global index (cross-process) of the unknown
   // associated to the ic-th component of node of global number
   // (cross-process) n (PEL::bad_index() if not an unknown)
   size_t_array2D global_unknown( GLOB_NODE->nb_global_nodes(), nbcs ) ;
   globalize( global_unknown ) ;

   LOC_2_GLOB.re_initialize( LINK->unknown_vector_size() ) ;
   OWNER.re_initialize( LINK->unknown_vector_size() ) ;
   size_t_vector const& comp = LINK->components_table() ;
   for( size_t n=0 ; n<nb_nodes ; n++ )
   {
      size_t globn = GLOB_NODE->global_node_index( n ) ;
      for( size_t iic=0 ; iic<comp.size() ; ++iic )
      {
         size_t ic = comp(iic) ;
         size_t glob = global_unknown( globn, iic ) ;
         if( LINK->DOF_is_unknown(n,ic) )
         {
            size_t loc = LINK->unknown_linked_to_DOF( n, ic ) ;
            PEL_CHECK( glob != PEL::bad_index() ) ;
            PEL_CHECK( loc < LINK->unknown_vector_size()  ) ;
            LOC_2_GLOB( loc ) = glob ;
            OWNER( loc ) = GLOB_NODE->rank_of_process_handling( n ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown const*
PDE_CrossProcessUnknownNumbering:: link( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: link" ) ;

   PDE_LinkDOF2Unknown const* result = LINK ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->field()->is_distributed() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Communicator const*
PDE_CrossProcessUnknownNumbering:: communicator( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: communicator" ) ;

   PEL_Communicator const* result = COMM ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: nb_unknowns_of_current_process( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: nb_unknowns_of_current_process" ) ;

   size_t result = nb_unknowns_on_process( COMM->rank() ) ;

   PEL_CHECK_POST( result <= nb_global_unknowns() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: nb_unknowns_on_process( size_t rank ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: nb_unknowns_on_process" ) ;

   size_t result = (int) NB_UNKNOWNS_ON_PROC( rank ) ;

   PEL_CHECK_POST( result <= nb_global_unknowns() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: global_unknown_index( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: global_unknown_index" ) ;
   PEL_CHECK_PRE( i < link()->unknown_vector_size() ) ;

   size_t result = LOC_2_GLOB( i ) ;

   PEL_CHECK_POST( result < nb_global_unknowns() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: rank_of_process_handling( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: rank_of_process_handling" ) ;
   PEL_CHECK_PRE( i < link()->unknown_vector_size() ) ;

   size_t result = OWNER( i ) ;

   PEL_CHECK_POST( result < COMM->nb_ranks() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: nb_global_unknowns( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: nb_global_unknowns" ) ;

   size_t result = NB_UNKNOWNS ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessUnknownNumbering:: global_unknown_linked_to_DOF(
                                                  size_t n,
                                                  size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: global_unknown_linked_to_DOF" ) ;
   PEL_CHECK_PRE( n < link()->nb_field_nodes() ) ;
   PEL_CHECK_PRE( ic < link()->field()->nb_components() ) ;
   PEL_CHECK_PRE( link()->components_table().has( ic ) ) ;
   PEL_CHECK_PRE( link()->DOF_is_unknown( n, ic ) ) ;

   size_t loc = LINK->unknown_linked_to_DOF( n, ic ) ;

   size_t result = LOC_2_GLOB( loc ) ;

   PEL_CHECK_POST( result < nb_global_unknowns() ||
                   result == PEL::bad_index() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessUnknownNumbering:: print( std::ostream& os,
                                          size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string s( indent_width, ' ' ) ;
   os << s
      << "Local number of unknowns for " << LINK->field()->name()
      << " on P" << COMM->rank() << " is "
      << LINK->unknown_vector_size() << " / " <<  NB_UNKNOWNS << std::endl ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessUnknownNumbering:: globalize( size_t_array2D& unknown )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: globalize" ) ;
   PEL_CHECK( unknown.index_bound( 0 ) ==
              link()->field()->cross_process_numbering()->nb_global_nodes() ) ;
   PEL_CHECK( unknown.index_bound( 1 ) == link()->components_table().size() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const nb_ranks = COMM->nb_ranks() ;
   size_t const last = nb_ranks-1 ;
   size_t const rank = COMM->rank() ;
   
   int nb_unk_on_current_proc = PEL::bad_int() ;
   
   if( rank>0 )
   {
      COMM->receive( rank-1 , unknown ) ;
      COMM->receive( rank-1 , NB_UNKNOWNS ) ;
   }
   else
   {
      NB_UNKNOWNS = 0 ;
      unknown.set( PEL::bad_index() ) ;
   }
   size_t const nbn = LINK->nb_field_nodes() ;
   size_t_vector const& comp = LINK->components_table() ;
   
   size_t nb_unk = NB_UNKNOWNS ;
   
   if( ORDERING == sequenceOfTheComponents )
   {
      for( size_t il=0 ; il<nbn ; ++il )
      {
         size_t const ig = GLOB_NODE->global_node_index( il ) ;
         if( GLOB_NODE->current_process_handles_node( il ) )
         {
            for( size_t iic=0 ; iic<comp.size() ; ++iic )
            {
               PEL_CHECK( unknown( ig, iic ) == PEL::bad_index() ) ;
               size_t const ic = comp( iic ) ;
               if( LINK->DOF_is_unknown( il, ic ) )
               {
                  unknown( ig, iic ) = NB_UNKNOWNS ;
                  NB_UNKNOWNS++ ;
               }
            }
         }
      }
   }

   if( ORDERING == sequenceOfTheNodes )
   {
      for( size_t iic=0 ; iic<comp.size() ; ++iic )
      {
         for( size_t il=0 ; il<nbn ; ++il )
         {
            if( GLOB_NODE->current_process_handles_node( il ) )
            {
               size_t const ig = GLOB_NODE->global_node_index( il ) ;
               PEL_CHECK( unknown( ig, iic ) == PEL::bad_index() ) ;
               size_t const ic = comp( iic ) ;
               if( LINK->DOF_is_unknown( il, ic ) )
               {
                  unknown( ig, iic ) = NB_UNKNOWNS ;
                  NB_UNKNOWNS++ ;
               }
            }
         }
      }
   }

   nb_unk_on_current_proc = NB_UNKNOWNS - nb_unk ;
   
   if( rank!=last )
   {
      COMM->send( rank+1, unknown ) ;
      COMM->send( rank+1, NB_UNKNOWNS ) ;

      COMM->receive( last, unknown ) ;
      COMM->receive( last, NB_UNKNOWNS ) ;
   }
   else
   {
      for( size_t i=0 ; i<last ; i++ )
      {
         COMM->send( i, unknown ) ;
         COMM->send( i, NB_UNKNOWNS ) ;
      }
   }
   
   COMM->all_gather( nb_unk_on_current_proc, NB_UNKNOWNS_ON_PROC ) ;

   for( size_t il=0 ; il<nbn ; ++il )
   {
      size_t const ig = GLOB_NODE->global_node_index( il ) ;
      for( size_t iic=0 ; iic<comp.size() ; ++iic )
      {
         if( !LINK->DOF_is_unknown( il, comp(iic) ) &&
             unknown( ig, iic )!= PEL::bad_index() )
         {
            raise_globalize_error(
               "    DOF of field \""+LINK->field()->name()+"\"\n"
               "    seems not to be free for all processes",
               unknown ) ;
         }
         else if( LINK->DOF_is_unknown( il, comp(iic) ) &&
                  unknown( ig, iic ) == PEL::bad_index() )
         {
            raise_globalize_error(
               "    DOF of field \""+LINK->field()->name()+"\"\n"
               "    seems not to be imposed for all processes",
               unknown ) ;
         }
      }
   }
   
  PEL_CHECK_POST(
     FORALL(
        ( size_t i=0 ; i<link()->nb_field_nodes() ; i++ ),
        FORALL(
           ( size_t ic=0 ; ic<link()->components_table().size() ; ic++ ),
           EQUIVALENT(
              link()->DOF_is_unknown( i, link()->components_table()(ic) ),
              unknown( GLOB_NODE->global_node_index(i), ic ) != PEL::bad_index() )
            ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessUnknownNumbering:: raise_globalize_error(
           std::string const& mes, size_t_array2D const& unknown ) const
//----------------------------------------------------------------------
{
   std::ostringstream os ;

   size_t nb_unks = unknown.index_bound(0) ;

   if( nb_unks<200 )
   {
      size_t_vector const& comp = LINK->components_table() ;
      os << "Unknowns : " << std::endl ;
      for( size_t i=0 ; i<nb_unks ; ++i )
      {
         os << "   Node global : " << i ;
         for( size_t ic=0 ; ic<comp.size() ; ++ic )
         {
            if( comp.size() != 1 )
            {
               os << std::endl << "      comp : " << comp(ic) ;
            }
            if( unknown( i, ic ) == PEL::bad_index() )
            {
               os << " => fixed" << std::endl ;
            }
            else
            {
               os << " => global unk : " << unknown( i, ic )
                  << std::endl ;
            }
         }
      }
      os << std::endl ;
      os << "Field : " << std::endl ;
      LINK->field()->print( os, 3 ) ;
      os << std::endl ;
   }
   os << "*** PDE_CrossProcessUnknownNumbering error : " << std::endl ;
   os << mes ;
   PEL_Error::object()->raise_internal( os.str() ) ;
}

//-----------------------------------------------------------------------
void
PDE_CrossProcessUnknownNumbering:: set_ordering_option(
                                              std::string const& option )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessUnknownNumbering:: set_ordering_option" ) ;

   if( option == "sequence_of_the_components" )
   {
      ORDERING = sequenceOfTheComponents ;
   }
   else if( option == "sequence_of_the_nodes" )
   {
      ORDERING = sequenceOfTheNodes ;
   }
   else
   {
      PEL_Error::object()->raise_plain( option + ": bad value " ) ;
   }
}

