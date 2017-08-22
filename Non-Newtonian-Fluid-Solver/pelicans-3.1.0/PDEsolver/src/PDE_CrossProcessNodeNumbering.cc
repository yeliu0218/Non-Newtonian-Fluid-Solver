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

#include <PDE_CrossProcessNodeNumbering.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_VectorIterator.hh>
#include <PEL_assertions.hh>
#include <boolArray2D.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_GridFE.hh>
#include <PDE_MeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <string>
#include <sstream>
#include <iostream>

using std::endl ;

//----------------------------------------------------------------------
PDE_CrossProcessNodeNumbering:: PDE_CrossProcessNodeNumbering(
                                             PEL_Object* a_owner,
                                             PDE_DomainBuilder const* a_dom,
                                             PEL_Communicator const* a_com )
//----------------------------------------------------------------------
   : PEL_Object( a_owner)
   , DOM( a_dom )
   , DIM( a_dom->nb_space_dimensions() )
   , PT( 0 )
   , LOCAL_NODE( 0 )
   , GLOBAL_NODE( 0 )
   , NODE_OWNER( 0 )
   , FIELD( 0 )
   , COMM( a_com )

{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: PDE_CrossProcessNodeNumbering" ) ;

   PT = GE_Point::create( this, DIM ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_CrossProcessNodeNumbering:: ~PDE_CrossProcessNodeNumbering( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: ~PDE_CrossProcessNodeNumbering" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: attach_field( PDE_DiscreteField* a_field )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: attach_field" ) ;
   PEL_CHECK_PRE( field() == 0 ) ;

   FIELD = a_field ;
   globalize() ;
   synchronize_imposed_values() ;

   PEL_CHECK_POST( field() == a_field ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_CrossProcessNodeNumbering:: field( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: field" ) ;

   PDE_DiscreteField const* result = FIELD ;

   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Communicator const*
PDE_CrossProcessNodeNumbering:: communicator( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: communicator" ) ;

   PEL_Communicator const* result = COMM ;

   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessNodeNumbering:: nb_global_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: nb_global_nodes" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;

   return( NODE_OWNER.size() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessNodeNumbering:: global_node_index( size_t n ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: global_node_index" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;
   PEL_CHECK_PRE( n < field()->nb_nodes() ) ;

   size_t result = GLOBAL_NODE( n ) ;

   PEL_CHECK_POST( result < nb_global_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CrossProcessNodeNumbering:: local_node( size_t global_index ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: local_node" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;
   PEL_CHECK_PRE( global_index < nb_global_nodes() ) ;

   size_t result = LOCAL_NODE( global_index ) ;

   PEL_CHECK_POST( result == PEL::bad_index() ||
                   result < field()->nb_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_CrossProcessNodeNumbering:: current_process_handles_node( size_t n ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: current_process_handles_node" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;
   PEL_CHECK_PRE( n < field()->nb_nodes() ) ;

   return rank_of_process_handling(n) == COMM->rank() ;
 }

//----------------------------------------------------------------------
size_t
PDE_CrossProcessNodeNumbering:: rank_of_process_handling( size_t n ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: rank_of_process_handling" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;
   PEL_CHECK_PRE( n < field()->nb_nodes() ) ;

   size_t result = NODE_OWNER( GLOBAL_NODE( n ) ) ;
   PEL_CHECK_POST( result < COMM->nb_ranks() ) ;
   return result ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: synchronize_valid_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: synchronize_valid_nodes" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;

   boolVector DOFs_is_unknown( nb_global_nodes() ) ;

   size_t const nb_local_nodes = FIELD->nb_nodes() ;

   size_t const rank = COMM->rank() ;
   size_t const size = COMM->nb_ranks() ;
   size_t const last = size-1 ;

   if( rank>0 )
      COMM->receive( rank-1, DOFs_is_unknown ) ;
   else
      DOFs_is_unknown.set( false ) ;
   for( size_t il=0 ; il<nb_local_nodes ; ++il )
      DOFs_is_unknown( GLOBAL_NODE( il ) ) |= FIELD->node_is_active( il ) ;
   if( rank!=last )
   {
      COMM->send( rank+1, DOFs_is_unknown ) ;
      COMM->receive( last, DOFs_is_unknown ) ;
   }
   else
      for( size_t i=0 ; i<last ; i++ )
         COMM->send( i, DOFs_is_unknown ) ;

   for( size_t il=0 ; il<nb_local_nodes ; ++il )
      if( DOFs_is_unknown( GLOBAL_NODE( il ) ) && !FIELD->node_is_active( il ) )
         FIELD->set_node_active( il ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: synchronize_imposed_values( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: synchronize_imposed_values" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;

   size_t const nb_local_nodes = FIELD->nb_nodes() ;
   size_t const nb_comps = FIELD->nb_components() ;

   size_t const rank = COMM->rank() ;
   size_t const nb_ranks = COMM->nb_ranks() ;

   intVector DOFs_DBC(0) ;
   doubleVector DOFs_DBC_values(0) ;

   for( size_t il=0 ; il<nb_local_nodes ; ++il )
   {
      if( current_process_handles_node( il ) )
      {
         size_t const ig = GLOBAL_NODE(il) ;
         PEL_ASSERT( LOCAL_NODE(ig)==il ) ;

         for( size_t ic=0 ; ic<nb_comps ; ++ic )
         {
            if( FIELD->DOF_has_imposed_value( il ,ic ) )
            {
               size_t i_glob = ig*nb_comps +ic ;
               DOFs_DBC.append( i_glob ) ;
               DOFs_DBC_values.append( FIELD->DOF_imposed_value( il, ic )  ) ;
            }
         }
      }
   }
   if( rank!=0 )
   {
      COMM->send( 0, DOFs_DBC ) ;
      COMM->send( 0, DOFs_DBC_values ) ;
   }
   else
   {
      for( size_t i=1 ; i<nb_ranks ; i++ )
      {
         intVector DOFs_DBC_other(0) ;
         doubleVector DOFs_DBC_values_other(0) ;

         COMM->receive( i, DOFs_DBC_other ) ;
         COMM->receive( i, DOFs_DBC_values_other ) ;

         for( size_t j=0 ; j<DOFs_DBC_other.size() ; j++ )
         {
            DOFs_DBC.append( DOFs_DBC_other(j) ) ;
            DOFs_DBC_values.append( DOFs_DBC_values_other(j) ) ;
         }
      }
   }
   COMM->broadcast( DOFs_DBC, 0 ) ;
   COMM->broadcast( DOFs_DBC_values, 0 ) ;

   for( size_t i=0 ; i<DOFs_DBC.size() ; i++ )
   {
      size_t iglob = DOFs_DBC(i) ;
      size_t const ig = iglob / nb_comps ;
      size_t const ic = iglob % nb_comps ;
      size_t const iloc = LOCAL_NODE( ig ) ;

      if( iloc != PEL::bad_index() )
         FIELD->equip_DOF_with_imposed_value(
                                  iloc, DOFs_DBC_values( i ), ic ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: synchronize_new_nodes( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: synchronize_new_nodes" ) ;
   PEL_CHECK_PRE( field() != 0 ) ;

   size_t old_nb_nodes = GLOBAL_NODE.size() ;
   size_t old_nb_global_nodes = NODE_OWNER.size() ;
   size_t nb_new_nodes = FIELD->nb_nodes()-old_nb_nodes ;

   bool go_on = COMM->boolean_or( nb_new_nodes > 0 ) ;
   if( go_on )
   {
      // Recover all node coordinates:
      doubleArray2D coord_glob( DIM, FIELD->nb_nodes() ) ;
      boolVector is_halo_glob( FIELD->nb_nodes() ) ;
      recover_coordinates( coord_glob, is_halo_glob ) ;

      // Only take new nodes:
      doubleArray2D coord( DIM, nb_new_nodes ) ;
      boolVector is_halo( nb_new_nodes ) ;
      for( size_t i=0 ; i<nb_new_nodes ; ++i )
      {
         for( size_t ic=0 ; ic<DIM ; ++ic )
         {
            coord( ic, i ) = coord_glob( ic, i+old_nb_nodes ) ;
         }
         is_halo( i ) = is_halo_glob( i+old_nb_nodes ) ;
      }

      // Recover new node global numbering:
      size_t_vector global_node( 0 ) ;
      size_t_vector node_owner( 0 ) ;
      recover_global_numbering( coord, is_halo, global_node,
                                node_owner, false ) ;

      // Update numberings:
      GLOBAL_NODE.resize( old_nb_nodes+global_node.size() ) ;
      for( size_t i=0 ; i<global_node.size() ; ++i )
      {
         GLOBAL_NODE( old_nb_nodes+i ) = global_node(i) + old_nb_global_nodes ;
      }
      NODE_OWNER.resize( old_nb_global_nodes+node_owner.size() ) ;
      for( size_t i=0 ; i<node_owner.size() ; ++i )
      {
         NODE_OWNER( old_nb_global_nodes+i ) = node_owner(i) ;
      }
      set_global_to_local( GLOBAL_NODE, NODE_OWNER, LOCAL_NODE ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: globalize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: globalize" ) ;
   PEL_CHECK_INV( invariant() ) ;

   doubleArray2D coord( DIM, FIELD->nb_nodes() ) ;
   boolVector is_halo( FIELD->nb_nodes() ) ;

   recover_coordinates( coord, is_halo ) ;

   recover_global_numbering( coord, is_halo, GLOBAL_NODE, NODE_OWNER, true ) ;

   set_global_to_local( GLOBAL_NODE, NODE_OWNER, LOCAL_NODE ) ;

}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: recover_coordinates(
                       doubleArray2D& coord, boolVector& is_halo ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: recover_coordinates" ) ;
   PEL_CHECK( coord.index_bound(0) == DIM ) ;
   PEL_CHECK( coord.index_bound(1) == FIELD->nb_nodes() ) ;
   PEL_CHECK( is_halo.size() == FIELD->nb_nodes() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool is_boundary_field = DOM->is_defined_on_boundary( FIELD->name() ) ;

   PDE_GridFE const* grid = DOM->finite_element_grid() ;
   PEL_VectorIterator* mesh_it = 0 ;
   if( !is_boundary_field )
   {
      mesh_it = PEL_VectorIterator::create( 0, grid->cells() ) ;
   }
   else
   {
      mesh_it = PEL_VectorIterator::create( 0, grid->bounds() ) ;
   }

   recover_coordinates( mesh_it, coord, is_halo ) ;

   mesh_it->destroy() ; mesh_it = 0 ;
}

//??? done for all the fields ??? factorization ???
//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: recover_coordinates(
                                            PEL_VectorIterator* mesh_it,
                                            doubleArray2D& coord,
                                            boolVector& is_halo ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: recover_coordinates fem" ) ;
   PEL_CHECK( mesh_it != 0 ) ;
   PEL_CHECK( coord.index_bound(0) == DIM ) ;
   PEL_CHECK( coord.index_bound(1) == FIELD->nb_nodes() ) ;
   PEL_CHECK( is_halo.size() == FIELD->nb_nodes() ) ;

   boolVector is_ok( FIELD->nb_nodes() ) ; //??? for subsequent check
   is_ok.set( false ) ;

   is_halo.set( true ) ;

   for( mesh_it->start() ; mesh_it->is_valid() ; mesh_it->go_next() )
   {
      PDE_MeshFE const* mesh = static_cast< PDE_MeshFE* >( mesh_it->item() ) ;
      bool halo_mesh = ( mesh->color() == GE_Color::halo_color() ) ;
      size_t ee = mesh->index_of_reference_element( FIELD ) ;
      PDE_ReferenceElement const* elm = mesh->reference_element( ee ) ;
      for( size_t ln=0 ; ln<mesh->nb_basis_functions( ee ) ; ++ln )
      {
         PDE_BasisFunction* bf = mesh->basis_function( ee, ln ) ;
         if( bf != 0 )
         {
            GE_Point const* pt_ref = elm->node_location( ln ) ;
            mesh->polyhedron()->apply_mapping( pt_ref, PT ) ;
            size_t n = bf->node_of_DOF( FIELD ) ;        
            for( size_t d=0 ; d<DIM ; ++d )
            {
               coord( d, n ) = PT->coordinate( d ) ;
            }
            if( !halo_mesh ) is_halo( n ) = false ;
            is_ok( n ) = true ;
         }
      }
   }

   bool ok = true ;
   for( size_t n=0 ; n<FIELD->nb_nodes() ; ++n )
   {
      ok &= is_ok( n ) ;
      // if( !is_ok( n ) ) PEL::out() << "!!! not ok for node: " << n << endl ;
   }
   PEL_ASSERT( ok ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: recover_global_numbering(
                                               doubleArray2D& coord,
                                               boolVector const& is_halo,
                                               size_t_vector& global_node,
                                               size_t_vector& node_owner,
                                               bool verbose ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: recover_global_numbering" ) ;
   PEL_CHECK( coord.index_bound(0) == DIM ) ;
   PEL_CHECK( coord.index_bound(1) == is_halo.size() ) ;

   size_t N = coord.index_bound(1) ;
   
   global_node.re_initialize( N ) ;
   global_node.set( PEL::bad_index() ) ;

   size_t const rank = COMM->rank() ;
   size_t const size = COMM->nb_ranks() ;
   size_t const last = size-1 ;
   size_t nb_owned = 0 ;
   size_t nb_halo = 0 ;
   size_t nb_ext = 0 ;

   COMM->merge( GE_Point::coordinates_comparator(), coord, global_node ) ;

   if( rank>0 )
   {
      COMM->receive( rank-1, node_owner ) ;
   }
   else
   {
      node_owner.re_initialize( coord.index_bound( 1 ) ) ;
      node_owner.set( PEL::bad_index() ) ;
   }
   for( size_t i=0 ; i<global_node.size() ; ++i )
   {
      size_t global_idx = global_node( i ) ;
      if( !is_halo( i ) && ( node_owner( global_idx ) == PEL::bad_index() ) )
      {
         node_owner( global_idx ) = rank ;
         nb_owned++ ;
      }
      else
      {
         if( is_halo( i ) ) nb_halo++ ;
         else nb_ext++ ;
      }

   }
   if( rank != last )
   {
      COMM->send( rank+1, node_owner ) ;
      COMM->receive( last, node_owner ) ;
   }
   else
   {
      for( size_t i=0 ; i<last ; ++i )
         COMM->send( i, node_owner ) ;
   }

   std::ostringstream os ;
   for( size_t i=0 ; i<N ; ++i )
   {
      if( !( node_owner(global_node(i))<communicator()->nb_ranks() ) )
      {
         os << "PDE_CrossProcessNodeNumbering:: recover_global_numbering node loc " << i
            << " glob " << global_node(i) <<
            " isn't owned by any process" << std::endl ;
         if(is_halo( i ) ) os << " This node is in halo for current process" << std::endl ;
         os << " Its position is : " ;
         size_t a_idx = ( rank==0 ? global_node(i) : i ) ;
         for( size_t d=0 ; d<coord.index_bound(0) ; d++ )
            os << " " << coord(d,a_idx) ;
         os << std::endl ;
      }
   }
   
   if( !os.str().empty() )
      PEL_Error::object()->raise_internal( os.str() ) ;
   
   if( COMM->nb_ranks()>1 && verbose )
   {
      PEL::out() << "*** Cross-process node numbering for \"" 
                 << field()->name() << "\"" << std::endl ;
      PEL::out() << "    nb global (cross-process) nodes: " 
                 << node_owner.size() << std::endl ;
      PEL::out() << "    nb local (on-process) nodes: " 
                 << global_node.size() << std::endl ;
      PEL::out() << "    nb handled nodes: " << nb_owned << std::endl ;
      PEL::out() << "    nb nodes in pure halo regions: " << nb_halo << std::endl ;
      PEL::out() << "    nb nodes handled by other processes: " 
                 << nb_ext << std::endl ;
      PEL::out() << std::endl ;
   }

   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<global_node.size() ; ++i ),
              global_node(i)<node_owner.size() ) ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<global_node.size() ; ++i ),
              node_owner(global_node(i))<communicator()->nb_ranks() ) ) ;
}

//----------------------------------------------------------------------
void
PDE_CrossProcessNodeNumbering:: set_global_to_local(
                                    size_t_vector const& global_nodes,
                                    size_t_vector const& nodes_owner,
                                    size_t_vector& local_nodes ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CrossProcessNodeNumbering:: set_global_to_local" ) ;
   PEL_CHECK(
      FORALL( ( size_t i=0 ; i<global_nodes.size() ; i++ ),
              global_nodes(i)<nodes_owner.size() ) ) ;

   local_nodes.re_initialize( nodes_owner.size() ) ;
   local_nodes.set( PEL::bad_index() ) ;
   for( size_t i=0 ; i<global_nodes.size() ; i++ )
   {
      size_t ig = global_nodes( i ) ;
      local_nodes( ig ) = i ;
   }
}
