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

#include <PDE_LocalFEmulti.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_BFvalues.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_MeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ; using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;

size_t const max_nb_ns_multi = 300 ;

//----------------------------------------------------------------------
PDE_LocalFEmulti:: PDE_LocalFEmulti( PEL_Object* a_owner, size_t aNbSpDims )
//----------------------------------------------------------------------
   : PDE_LocalFE( a_owner, aNbSpDims )
   , NB_MESHES( 0 )
   , MESHES( PEL_Vector::create( this, 0 ) )
   , NB_NESTED_LEVELS( 0 )
   , MAX_NB_NODES( max_nb_ns_multi )
   , HAS_DISC( 0, 0 )
   , ELM_DERIS( 0, 0 )
   , ELM_index( 0, 0 )
   , NB_LOCAL_NODES( 0 )
   , GLOBAL_NODE( 0, 0 )
   , BASIS_FUNC( PEL_Vector::create( this, 0 ) )
   , NODE_LOCATIONS( PEL_Vector::create( this, 0 ) )
   , NODE_LOC_OK( 0 ) 
   , ELM_BF_index( 0, 0, 0 )
   , REFINEMENT_LEVEL( 0, 0, 0 )
   , HAS_CP( false )
   , CP_LOCATION( GE_Point::create( this, aNbSpDims ) )
   , CP_NB_MESHES( 0 )
   , CP_MESHES( 0 )
   , C_BF_VALS( PEL_Vector::create( this, 0 ) )
   , C_BF_VALS_POOL( PEL_Vector::create( this, 0 ) )
   , QR_PROVIDER( 0 )
   , QR_PROVIDER_STAT( PEL::bad_index() )
   , NB_IPs( 0 )
   , IP_LOCATIONS( PEL_Vector::create( this, 0 ) )
   , IP_LOCATIONS_POOL( PEL_Vector::create( this, 0 ) )
   , IP_WEIGHTS( 0 )
   , IP_NB_MESHES( 0 )
   , IP_MESHES( 0, 0 )
   , I_BF_VALS( PEL_Vector::create( this, 0 ) )
   , I_BF_VALS_POOL( PEL_Vector::create( this, 0 ) )
   , i_IP( PEL::bad_index() )
   , MASKED_VALUES( 0, 0, 0 )
   , MASKED_LEVEL( PEL::bad_index() )
   , ROW_iF( PEL::bad_index() )
   , COL_iF( PEL::bad_index() )
   , ROW_im( PEL::bad_index() )
   , COL_im( PEL::bad_index() )
   , ROW_ee( PEL::bad_index() )
   , COL_ee( PEL::bad_index() )
   , ROW_NODES( 0 )
   , COL_NODES( 0 )
   , ROW_Ns( 0 )
   , COL_Ns( 0 )
   , ROW_dNs( 0, 0 )
   , COL_dNs( 0, 0 )
   , ROW_d2Ns( 0, 0, 0 )
   , COL_d2Ns( 0, 0, 0 )
{
   PEL_LABEL( "PDE_LocalFEmulti:: PDE_LocalFEmulti" ) ;

   PEL_ASSERT( N   == PDE_BFvalues::N   ) ;
   PEL_ASSERT( dN  == PDE_BFvalues::dN  ) ;
   PEL_ASSERT( d2N == PDE_BFvalues::d2N ) ;
   
   PEL_CHECK( nb_space_dimensions()>=1 && nb_space_dimensions()<=3 ) ;
}

//----------------------------------------------------------------------
PDE_LocalFEmulti:: ~PDE_LocalFEmulti( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: nb_local_nodes( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: nb_local_nodes" ) ;
   PEL_CHECK_PRE( nb_local_nodes_PRE( ff ) ) ;

   size_t result = NB_LOCAL_NODES( field_local_index( ff ) ) ;

   PEL_CHECK_POST( nb_local_nodes_POST( result, ff ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEmulti:: local_node_location( PDE_DiscreteField const* ff,
                                        size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: local_node_location" ) ;
   PEL_CHECK_PRE( local_node_location_PRE( ff, i ) ) ;

   size_t iF = field_local_index( ff ) ;
   PEL_Vector* ff_nodes =
                   static_cast<PEL_Vector*>( NODE_LOCATIONS->at( iF ) ) ;
   if( !NODE_LOC_OK( iF ) )
   {
      size_t nb_nodes = NB_LOCAL_NODES( iF ) ;
      size_t old_size = ff_nodes->index_limit() ;
      if( nb_nodes > old_size )
      {
         ff_nodes->resize( nb_nodes ) ;
         for( size_t jj=old_size ; jj<nb_nodes ; ++jj )
         {
            ff_nodes->set_at( jj, GE_Point::create( NODE_LOCATIONS,
                                                    nb_space_dimensions() ) ) ;
         }
      }
      size_t_vector found_nodes( 0 ) ;
      size_t il = 0 ;
      for( size_t im=0 ; im<NB_MESHES ; ++im )
      {
         if( HAS_DISC( im, iF ) )
         {
            PDE_MeshFE* mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;
            size_t ee = ELM_index( im, iF ) ;
            PDE_ReferenceElement const* elm = mesh->reference_element( ee ) ;
            for( size_t cl=0 ; cl<NB_NESTED_LEVELS( im ) ; ++cl )
            {
               GE_Mpolyhedron const* poly = mesh->polyhedron() ;
               for( size_t ii=0 ; ii<mesh->nb_basis_functions(ee) ; ++ii )
               {
                  PDE_BasisFunction* bf = mesh->basis_function( ee, ii ) ;
                  if( ( bf != 0 ) &&  
                      ( ff->node_is_active( bf->node_of_DOF( ff ) ) ) )
                  {
                     size_t nn = bf->node_of_DOF( ff ) ;
                     if( !found_nodes.has( nn ) )
                     {
                        GE_Point const* pt_ref = elm->node_location( ii ) ;
                        GE_Point* pt = 
                           static_cast< GE_Point* >( ff_nodes->at( il ) ) ;
                           poly->apply_mapping( pt_ref, pt ) ;
                        PEL_CHECK( ELM_BF_index( il, im, iF ) == ii ) ;
                        PEL_CHECK( nn==GLOBAL_NODE( il, iF ) ) ;
                        ++il ;
                        found_nodes.append( nn ) ;
                     }
                  }
               }
               PEL_ASSERT( cl==NB_NESTED_LEVELS(im)-1 || mesh->parent()!=0 ) ;
               mesh = mesh->parent() ;
            }
         }
      }
      PEL_ASSERT( il == NB_LOCAL_NODES( iF ) ) ;
      NODE_LOC_OK( iF ) = true ;
   }

   GE_Point const* result = static_cast<GE_Point*>( ff_nodes->at( i ) ) ;

   PEL_CHECK_POST( local_node_location_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmulti:: local_node_is_in_mesh( PDE_DiscreteField const* ff,
                                          size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: local_node_is_in_mesh" ) ;
   PEL_CHECK_PRE( local_node_is_in_mesh_PRE( ff, i ) ) ;

   size_t ii = i * nb_handled_fields() + field_local_index( ff ) ;
   PDE_BasisFunction const* bf = 
                     static_cast<PDE_BasisFunction*>( BASIS_FUNC->at( ii ) ) ;

   bool result = is_in_mesh( ff, i, bf ) ;

   PEL_CHECK_POST( local_node_is_in_mesh_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: global_node( PDE_DiscreteField const* ff,
                                size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: global_node" ) ;
   PEL_CHECK_PRE( global_node_PRE( ff, i ) ) ;

   size_t result = GLOBAL_NODE( i, field_local_index( ff ) ) ;

   PEL_CHECK_POST( global_node_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: node_refinement_level( PDE_DiscreteField const* ff,
                                          size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: node_refinement_level" ) ;
   PEL_CHECK_PRE( node_refinement_level_PRE( ff, i ) ) ;

   size_t ii = i * nb_handled_fields() + field_local_index( ff ) ;
   PDE_BasisFunction const* bf = 
                     static_cast<PDE_BasisFunction*>( BASIS_FUNC->at( ii ) ) ;
                     
   size_t result = bf->refinement_level() ;

   PEL_CHECK_POST( node_refinement_level_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEmulti:: calculation_point( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: calculation_point" ) ;

   GE_Point const* result = 0 ;
   if( HAS_CP )
   {
      result = CP_LOCATION ;
   }

   PEL_CHECK_POST( IMPLIES( result!=0, result->is_under_ownership_of( this ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, polyhedron()->contains( result ) ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0,
                        result->nb_coordinates()==nb_space_dimensions() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: value_at_pt( PDE_DiscreteField const* ff,
                                size_t level,
                                size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: value_at_pt" ) ;
   PEL_CHECK_PRE( value_at_pt_PRE( ff, level, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t idx = ELM_BF_index( i, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
         PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                   bfv->N_at_pt( idx ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: gradient_at_pt( PDE_DiscreteField const* ff,
                                   size_t level,
                                   size_t a,
                                   size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: gradient_at_pt" ) ;
   PEL_CHECK_PRE( gradient_at_pt_PRE( ff, level, a, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t idx = ELM_BF_index( i, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
         PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                   bfv->dN_at_pt( idx, a ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: hessian_at_pt( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t a,
                                  size_t b,
                                  size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: hessian_at_pt" ) ;
   PEL_CHECK_PRE( hessian_at_pt_PRE( ff, level, a, b, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t idx = ELM_BF_index( i, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
         PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                   bfv->d2N_at_pt( idx, a, b ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: N_at_pt( PDE_DiscreteField const* ff,
                            size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: N_at_pt" ) ;
   PEL_CHECK_PRE( N_at_pt_PRE( ff, i ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;

   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;

   size_t idx = ELM_BF_index( i, im, iF ) ;
   if( idx != PEL::bad_index() ) 
   {
      size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
      result = bfv->N_at_pt( idx ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: dN_at_pt( PDE_DiscreteField const* ff,
                             size_t i,
                             size_t a ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: dN_at_pt" ) ;
   PEL_CHECK_PRE( dN_at_pt_PRE( ff, i, a ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;

   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;

   size_t idx = ELM_BF_index( i, im, iF ) ;
   if( idx != PEL::bad_index() )
   {
      size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
      result = bfv->dN_at_pt( idx, a ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: d2N_at_pt( PDE_DiscreteField const* ff,
                              size_t i,
                              size_t a,
                              size_t b ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: d2N_at_pt" ) ;
   PEL_CHECK_PRE( d2N_at_pt_PRE( ff, i, a, b ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;

   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_pt( iF, im, ee ) ;

   size_t idx = ELM_BF_index( i, im, iF ) ;
   if( idx != PEL::bad_index() )
   {
      size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( im, cl, ee ) ;
      result = bfv->d2N_at_pt( idx, a, b ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: set_row_and_col_fields( 
                                     PDE_DiscreteField const* row_field, 
                                     PDE_DiscreteField const* col_field )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: set_row_and_col_fields" ) ;
   PEL_CHECK_PRE( set_row_and_col_fields_PRE( row_field, col_field ) ) ;

   ROW_iF = field_local_index( row_field ) ;
   COL_iF = field_local_index( col_field ) ;

   build_connectivities() ;

   PEL_CHECK_POST( set_row_and_col_fields_POST( row_field, col_field ) ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_LocalFEmulti:: field( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: field" ) ;

   size_t iF = (sf == row ) ? ROW_iF : COL_iF ;
   
   PDE_DiscreteField const* result = 0 ;
   if( iF != PEL::bad_index() )
   {
      result = handled_field( iF ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEmulti:: row_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: row_field_node_connectivity" ) ;
   PEL_CHECK_PRE( row_field_node_connectivity_PRE() ) ;

   size_t_vector const& result = ROW_NODES ;

   PEL_CHECK_POST( row_field_node_connectivity_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEmulti:: col_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: col_field_node_connectivity" ) ;
   PEL_CHECK_PRE( col_field_node_connectivity_PRE() ) ;

   size_t_vector const& result = COL_NODES ;

   PEL_CHECK_POST( col_field_node_connectivity_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: nb_basis_functions( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: nb_basis_functions" ) ;
   PEL_CHECK_PRE( nb_basis_functions_PRE( sf ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = (sf == row) ? NB_LOCAL_NODES( ROW_iF ) : 
                                 NB_LOCAL_NODES( COL_iF ) ;

   PEL_CHECK_POST( nb_basis_functions_POST( result, sf ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: start_IP_iterator( GE_QRprovider const* qrp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: start_IP_iterator" ) ;
   PEL_CHECK_PRE( start_IP_iterator_PRE( qrp ) ) ;

   HAS_CP = false ;

   // si c'est le meme provider et qu'on n'a pas change de maille,
   // on va avoir la meme règle de quadrature, donc il est inutile
   // de refaire les calculs
   if( qrp != QR_PROVIDER || qrp->status() != QR_PROVIDER_STAT )
   {
      GE_QuadratureRule const* qr = quadrature_rule( qrp ) ;
      compute_itg_pts( qr ) ;
      if( qr->owner() == 0 ) qr->destroy() ;
      QR_PROVIDER = qrp ;
      QR_PROVIDER_STAT = qrp->status() ;
   }
   i_IP = 0 ;

   prepare_for_Ns_requests() ;

   PEL_CHECK_POST( start_IP_iterator_POST( qrp ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: go_next_IP( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: go_next_IP" ) ;
   PEL_CHECK_PRE( go_next_IP_PRE() ) ;

   ++i_IP ;

   prepare_for_Ns_requests() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmulti:: valid_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: valid_IP" ) ;
   PEL_CHECK_PRE( valid_IP_PRE() ) ;

   return( i_IP<NB_IPs ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: weight_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: weight_of_IP" ) ;
   PEL_CHECK_PRE( weight_of_IP_PRE() ) ;

   return( IP_WEIGHTS( i_IP ) ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEmulti:: coordinates_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: coordinates_of_IP" ) ;
   PEL_CHECK_PRE( coordinates_of_IP_PRE() ) ;

   GE_Point const* result = static_cast<GE_Point*>( IP_LOCATIONS->at(i_IP) ) ;

   PEL_CHECK_POST( coordinates_of_IP_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: value_at_IP( PDE_DiscreteField const* ff,
                                size_t level,
                                size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: value_at_IP" ) ;
   PEL_CHECK_PRE( value_at_IP_PRE( ff, level, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   if( level== MASKED_LEVEL )
   {
      result = MASKED_VALUES( iF, ic, i_IP ) ;
   }
   else
   {
      size_t im = PEL::bad_index() ;
      size_t ee = PEL::bad_index() ;
      locate_discretization_at_IP( iF, i_IP, im, ee ) ;
      for( size_t il=0 ; il<NB_LOCAL_NODES( iF ) ; ++il )
      {
         size_t idx = ELM_BF_index( il, im, iF ) ;
         if( idx != PEL::bad_index() )
         {
            size_t cl = REFINEMENT_LEVEL( il, im, iF ) ;
            PDE_BFvalues const* bfv = I_BFvalues( im, cl, ee, i_IP ) ;
            result += ff->DOF_value( level, GLOBAL_NODE( il, iF ), ic ) *
                      bfv->N_at_pt( idx ) ;
         }
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: gradient_at_IP( PDE_DiscreteField const* ff,
                                   size_t level,
                                   size_t a, 
                                   size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: gradient_at_IP" ) ;
   PEL_CHECK_PRE( gradient_at_IP_PRE( ff, level, a, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   
   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_IP( iF, i_IP, im, ee ) ;

   for( size_t il=0 ; il<NB_LOCAL_NODES( iF ) ; ++il )
   {
      size_t idx = ELM_BF_index( il, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( il, im, iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( im, cl, ee, i_IP ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( il, iF ), ic ) *
                   bfv->dN_at_pt( idx, a ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: hessian_at_IP( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t a, 
                                  size_t b, 
                                  size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: hessian_at_IP" ) ;
   PEL_CHECK_PRE( hessian_at_IP_PRE( ff, level, a, b, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   
   size_t im = PEL::bad_index() ;
   size_t ee = PEL::bad_index() ;
   locate_discretization_at_IP( iF, i_IP, im, ee ) ;

   for( size_t il=0 ; il<NB_LOCAL_NODES( iF ) ; ++il )
   {
      size_t idx = ELM_BF_index( il, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( il, im, iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( im, cl, ee, i_IP ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( il, iF ), ic ) *
                   bfv->d2N_at_pt( idx, a, b ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: N_at_IP( field_id sf, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: N_at_IP" ) ;
   PEL_CHECK_PRE( N_at_IP_PRE( sf, i ) ) ;

   double result = 0.0 ;

   if( sf == row )
   {
      size_t idx = ELM_BF_index( i, ROW_im, ROW_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, ROW_im, ROW_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( ROW_im, cl, ROW_ee, i_IP ) ;
         result = bfv->N_at_pt( idx ) ;
      }
   }
   else
   {
      size_t idx = ELM_BF_index( i, COL_im, COL_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, COL_im, COL_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( COL_im, cl, COL_ee, i_IP ) ;
         result = bfv->N_at_pt( idx ) ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: dN_at_IP( field_id sf, size_t i, size_t a )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: dN_at_IP" ) ;
   PEL_CHECK_PRE( dN_at_IP_PRE( sf, i, a ) ) ;

   double result = 0.0 ;

   if( sf == row )
   {
      size_t idx = ELM_BF_index( i, ROW_im, ROW_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, ROW_im, ROW_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( ROW_im, cl, ROW_ee, i_IP ) ;
         result = bfv->dN_at_pt( idx, a ) ;
      }
   }
   else
   {
      size_t idx = ELM_BF_index( i, COL_im, COL_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, COL_im, COL_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( COL_im, cl, COL_ee, i_IP ) ;
         result = bfv->dN_at_pt( idx, a ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEmulti:: d2N_at_IP( field_id sf,
                              size_t i, size_t a, size_t b )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: d2N_at_IP" ) ;
   PEL_CHECK_PRE( d2N_at_IP_PRE( sf, i, a, b ) ) ;
 
   double result = 0.0 ;

   if( sf == row )
   {
      size_t idx = ELM_BF_index( i, ROW_im, ROW_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, ROW_im, ROW_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( ROW_im, cl, ROW_ee, i_IP ) ;
         result = bfv->d2N_at_pt( idx, a, b ) ;
      }
   }
   else
   {
      size_t idx = ELM_BF_index( i, COL_im, COL_iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, COL_im, COL_iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( COL_im, cl, COL_ee, i_IP ) ;
         result = bfv->d2N_at_pt( idx, a, b ) ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PDE_LocalFEmulti:: Ns_at_IP( field_id sf )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: Ns_at_IP" ) ;
   PEL_CHECK_PRE( Ns_at_IP_PRE( sf ) ) ;
   
   doubleVector const& result = ( sf==row ? ROW_Ns : COL_Ns ) ;

   PEL_CHECK_POST( Ns_at_IP_POST( result, sf ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PDE_LocalFEmulti:: dNs_at_IP( field_id sf )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: dNs_at_IP" ) ;
   PEL_CHECK_PRE( dNs_at_IP_PRE( sf ) ) ;

   doubleArray2D const& result = ( sf==row ? ROW_dNs : COL_dNs ) ;

   PEL_CHECK_POST( dNs_at_IP_POST( result, sf ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PDE_LocalFEmulti:: d2Ns_at_IP( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: d2Ns_at_IP" ) ;
   PEL_CHECK_PRE( d2Ns_at_IP_PRE( sf ) ) ;

   doubleArray3D const& result = ( sf==row ? ROW_d2Ns : COL_d2Ns ) ;
   
   PEL_CHECK_POST( d2Ns_at_IP_POST( result, sf ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: mask_value_at_IP( PDE_DiscreteField const* ff,
                                     size_t level,
                                     double value,
                                     size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: mask_value_at_IP" ) ;
   PEL_CHECK_PRE( mask_value_at_IP_PRE( ff, level, value, ic ) ) ;
   
   if( MASKED_LEVEL != level )
   {
      if( MASKED_LEVEL != PEL::bad_index() )
      {
         PEL_Error::object()->raise_plain( "Only one level can be set" ) ;
      }
      MASKED_LEVEL = level ;
      MASKED_VALUES.set( 0 ) ;
   }
   
   MASKED_VALUES( field_local_index( ff ), ic, i_IP ) = value ;

   PEL_CHECK_POST( mask_value_at_IP_POST( ff, level, value, ic ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: print_local_discretization_of_fields(
                                               std::ostream& os, 
                                               size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: print_local_discretization_of_fields" ) ;
   PEL_CHECK_PRE( print_local_discretization_of_fields_PRE( os,
                                                            indent_width ) ) ;

   PDE_LocalFE::print_local_discretization_of_fields( os, indent_width ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;

   for( size_t im=0 ; im<NB_MESHES ; ++im )
   {
      PDE_MeshFE const* the_mesh = static_cast<PDE_MeshFE*>( MESHES->at(im) ) ;
      os << space << the_mesh->polyhedron()->name() 
         << " " << the_mesh->id_number() << endl ;
      for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
      {
         PDE_MeshFE const* mesh = the_mesh ;
         PDE_DiscreteField const* ff = handled_field( iF ) ;

         os << space << more << "\"" << ff->name() << "\" : " ;
	 if( HAS_DISC( im, iF ) )
	 {
            size_t ee = mesh->index_of_reference_element( ff ) ;
            os << mesh->reference_element( ee )->name() << endl ;
            for( size_t cl=1 ; cl<NB_NESTED_LEVELS( im ) ; ++cl )
            {
   	       mesh = mesh->parent() ;
               ee = mesh->index_of_reference_element( ff ) ;
               std::string spp( ff->name().size()+8, ' ' ) ;
               os << space << spp << mesh->reference_element( ee )->name() 
                  << " on parent " << mesh->polyhedron()->name() 
                  << " " << mesh->id_number() << endl ;
            }
	 }
         else
         {
            os << "no discretization" << endl ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: print_current_IP( std::ostream& os, 
                                     size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: print_current_IP" ) ;
   PEL_CHECK_PRE( print_current_IP_PRE( os, indent_width ) ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;

   ios_base::fmtflags original_flags = os.flags() ;
   os.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision( 10 ) ;

   os << space << "IP : " ;
   coordinates_of_IP()->print( os, 0 ) ;
   os << " weight : " << weight_of_IP() << endl ;
   for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
   {
      if( NB_LOCAL_NODES( iF ) != 0 )
      {
         PDE_DiscreteField const* ff = handled_field( iF ) ;
         os << space << more << "\"" << ff->name() 
            << "\" : discretization of " ;
         size_t im = PEL::bad_index() ;
         size_t e = PEL::bad_index() ;
         size_t i=0 ;
         for( ; i<IP_NB_MESHES( i_IP ) ; ++i )
         {
            im = IP_MESHES( i_IP, i ) ;
            e = ELM_index( im, iF ) ;
            if( e != PEL::bad_index() ) break ;
         }
         for( ++i ; i<IP_NB_MESHES( i_IP ) ; ++i )
         {
            PEL_CHECK( ELM_index( IP_MESHES(i_IP,i), iF ) == PEL::bad_index() ) ;
         }
         PDE_MeshFE const* mm = static_cast<PDE_MeshFE*>( MESHES->at( im ) ) ;
         os << mm->polyhedron()->name()
            << " " << mm->id_number() << endl ;

         if( handled_field_deri(iF) & PDE_LocalFEmulti::N )
         {
            for( size_t level = 0 ; level<ff->storage_depth() ; ++level )
            {
               for( size_t ic=0 ; ic < ff->nb_components() ; ++ic )
               {
                  os << space << more << more 
                     << "value( level=" << level << " ; ic= " << ic << " ) : ";
                  double val = value_at_IP( ff, level, ic ) ;
                  if( PEL::abs(val)<1.E-12 ) val = 0. ;
                  os << std::setw( 20 ) << val << endl ;
               }
            }
         }
      }
   }
   os << std::setprecision(p) ;
   os.flags( original_flags ) ;  
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: print_values_at_current_IP( std::ostream& os, 
                                               size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: print_values_at_current_IP_values" ) ;
   PEL_CHECK_PRE( print_values_at_current_IP_PRE( os, indent_width ) ) ;

   ios_base::fmtflags original_flags = os.flags() ;
   os.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision( 6 ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;

   os << space << "IP : " ;
   coordinates_of_IP()->print( os, 0 ) ;
   os << " weight : " << weight_of_IP() << endl ;

   if( COL_iF != ROW_iF )
   {
      os << space << "row field : \"" << field(row)->name() << "\"" << endl ;
      os << space << more << "   " ;
   }
   else
   {
      os << space << "row and col field : \"" << field(row)->name() << "\"" 
         << endl ;
      os << space << more << "   " ;
   }
   if( handled_field_deri(ROW_iF) & N ) os << setw( 15 ) << "N      " ;
   if( handled_field_deri(ROW_iF) & dN )
      for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
      {
         ostringstream mesg ;
         mesg << "dN(" << a << ")    " ;
         os << setw( 15 ) << mesg.str() ;
      }
   os << endl ;

   for( size_t i=0 ; i<nb_basis_functions( row ) ; ++i )
   {
      os << space << more << setw( 3 ) << ROW_NODES(i) ;
      if( handled_field_deri(ROW_iF) & N )
      {
         double n = N_at_IP( row, i ) ;
         if( PEL::abs(n)<1.E-12 ) n = 0. ;
         os << setw( 15 ) << n ;
      }
      if( handled_field_deri(ROW_iF) & dN )
	 for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
	 {
            double dn = dN_at_IP( row, i, a ) ;
            if( PEL::abs(dn)<1.E-12 ) dn = 0. ;
            os << setw( 15 ) << dn ;
         }
      os << endl ;
   }

   if( COL_iF != ROW_iF )
   {
      os << space << "col field : \"" << field(col)->name() << "\"" << endl ;
      os << space << more << "   " ;
      if( handled_field_deri(COL_iF) & N ) os << setw( 15 ) << "N      " ;
      if( handled_field_deri(COL_iF) & dN )
         for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
         {
            ostringstream mesg ;
            mesg << "dN(" << a << ")    " ;
            os << setw( 15 ) << mesg.str() ;
         }
      os << endl ;

      for( size_t i=0 ; i<nb_basis_functions( col ) ; ++i )
      {
         os << space << more << setw( 3 ) << COL_NODES(i) ;
         if( handled_field_deri(COL_iF) & N )
         {
            double n = N_at_IP( col, i ) ;
            if( PEL::abs(n)<1.E-12 ) n = 0. ;
            os << setw( 15 ) << n ;
         }
         if( handled_field_deri(COL_iF) & dN )
	    for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
   	    {
               double dn = dN_at_IP( col, i, a ) ;
               if( PEL::abs(dn)<1.E-12 ) dn = 0. ;
               os << setw( 15 ) << dn ;
            }
         os << endl ;
      }
   }

    os << std::setprecision(p) ;
    os.flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: clear_local_fields( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: clear_local_fields" ) ;

   HAS_DISC.set( false ) ;

   NB_MESHES = 0 ;
   MESHES->re_initialize( 0 ) ;

   if( HAS_DISC.index_bound(1)<nb_handled_fields() )
   {
      HAS_DISC.re_initialize( HAS_DISC.index_bound(0), nb_handled_fields() ) ;
   }
   HAS_DISC.set( false ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: append_local_field( size_t iF, PDE_MeshFE const* lf )
//----------------------------------------------------------------------
{  
   PEL_LABEL( "PDE_LocalFEmulti:: append_local_field" ) ;
   PEL_CHECK( iF < nb_handled_fields() ) ;
   PEL_CHECK( lf != 0 ) ;
   
   size_t im = local_mesh_index( lf ) ;
   if( im == PEL::bad_index() )
   {
      im = MESHES->count() ;
      MESHES->append( const_cast<PDE_MeshFE*>( lf ) ) ;
      NB_MESHES++ ;
      PEL_CHECK( local_mesh_index( lf ) == im ) ;
      
      if( im >= HAS_DISC.index_bound(0) )
      {
         HAS_DISC.raise_first_index_bound( im+1 ) ; //??? buckets
	 for( size_t i=0 ; i<nb_handled_fields() ; ++i )
	 {
	    HAS_DISC( im, i ) = false ;
	 }
      }
   }

   HAS_DISC( im, iF ) = true ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: terminate_local_fields( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: terminate_local_fields" ) ;

   size_t const old_nb_fields = ELM_index.index_bound(1) ;
   size_t const nb_fields = nb_handled_fields() ;

   // Set number of reference elements:
   size_t nb_ref_elms = 0 ;
   for( size_t im=0 ; im<NB_MESHES ; ++im )
   {
      PDE_MeshFE* the_mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;
      nb_ref_elms += the_mesh->nb_reference_elements() ;
   }

   // Resize inner table only if a new handle field has been declared:
   if( nb_fields > old_nb_fields )
   {
      for( size_t i=old_nb_fields ; i<nb_fields ; ++i )
      {
         NODE_LOCATIONS->append( PEL_Vector::create( NODE_LOCATIONS, 0 ) ) ;
      }
   }
   
   initialize_ELM_BF_vectors_1( NB_MESHES, I_BF_VALS, C_BF_VALS ) ;
   NB_NESTED_LEVELS.re_initialize( NB_MESHES ) ;
   
   ELM_index.re_initialize( NB_MESHES, nb_fields ) ;
   ELM_DERIS.re_initialize( NB_MESHES, nb_ref_elms ) ;
   ELM_index.set( PEL::bad_index() ) ;
   ELM_DERIS.set( -1 ) ;
  
//????????? limiter les reinitializations 
   NB_LOCAL_NODES.re_initialize( nb_fields ) ;
   GLOBAL_NODE.re_initialize( MAX_NB_NODES, nb_fields ) ;
   BASIS_FUNC->re_initialize( MAX_NB_NODES*nb_fields ) ;
   NODE_LOC_OK.re_initialize( nb_fields ) ;
   ELM_BF_index.re_initialize( MAX_NB_NODES, NB_MESHES, nb_fields ) ;
   ELM_BF_index.set( PEL::bad_index() ) ;
   REFINEMENT_LEVEL.re_initialize( MAX_NB_NODES, NB_MESHES, nb_fields ) ;

   std::vector< size_t_vector > found_nodes( nb_fields, 0 ) ;
   std::vector< size_t > il( nb_fields ) ;
   for( size_t im=0 ; im<NB_MESHES ; ++im )
   {
      PDE_MeshFE* the_mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;

      NB_NESTED_LEVELS( im ) = 0 ;
      for( size_t iF=0 ; iF<nb_fields ; ++iF )
      {
         PDE_DiscreteField const* ff = handled_field( iF ) ;

         if( HAS_DISC( im, iF ) )
         {
            PDE_MeshFE* mesh = the_mesh ;

            size_t ee = mesh->index_of_reference_element( ff ) ;

            if( ELM_DERIS( im, ee ) == -1 )
            {
               ELM_DERIS( im, ee ) = handled_field_deri( iF ) ;
            }
            else
            {
               ELM_DERIS( im, ee ) |= handled_field_deri( iF ) ;
            }
            ELM_index( im, iF ) = ee ;

            size_t cl = 0 ;
            do
            {
               for( size_t i=0 ; i<mesh->nb_basis_functions(ee) ; ++i )
               {
                  PDE_BasisFunction* bf = mesh->basis_function( ee, i ) ;
                  if( ( bf != 0 ) && 
                      ( ff->node_is_active( bf->node_of_DOF( ff ) ) ) )
                  {
                     size_t nn = bf->node_of_DOF( ff ) ;
                     if( !found_nodes[iF].has( nn ) )
                     {
                        found_nodes[iF].append( nn ) ;
                        GLOBAL_NODE( il[iF], iF ) = nn ;
                        BASIS_FUNC->set_at( il[iF]*nb_fields+iF, bf ) ;
                        ELM_BF_index( il[iF], im, iF ) = i ;
                        REFINEMENT_LEVEL( il[iF], im, iF ) = cl ;
                        if( cl+1 > NB_NESTED_LEVELS( im ) )
                                 NB_NESTED_LEVELS( im ) = cl+1 ;
                        ++il[iF] ;
                        if( il[iF] == MAX_NB_NODES ) 
                        {
                           MAX_NB_NODES += 10 ;
                           GLOBAL_NODE.raise_first_index_bound( 
                                                         MAX_NB_NODES ) ;
                           BASIS_FUNC->resize( MAX_NB_NODES*nb_fields ) ;
                           ELM_BF_index.raise_first_index_bound( 
                                                         MAX_NB_NODES ) ;
                           for( size_t j=il[iF] ; j<MAX_NB_NODES ; ++j )
                              for( size_t jm=0 ; jm<NB_MESHES ; ++jm )
                                 for( size_t jF=0 ; jF<nb_fields ; ++jF )
                                    ELM_BF_index( j, jm, jF ) = PEL::bad_index() ;
                           REFINEMENT_LEVEL.raise_first_index_bound(  
                                                         MAX_NB_NODES ) ;
                        }
                     }
                     else 
                     {
                        bool found = false ;
                        for( size_t j=0 ; !found && j<=il[iF] ; ++j )
                        {
                           if( GLOBAL_NODE( j, iF ) == nn )
                           {
                              found = true ;
                              ELM_BF_index( j, im, iF ) = i ;
                              REFINEMENT_LEVEL( j, im, iF ) = cl ;
                              if( cl+1 > NB_NESTED_LEVELS( im ) )
                                 NB_NESTED_LEVELS( im ) = cl+1 ;
                           }
                        }
                        PEL_ASSERT( found ) ;
                     }
                  }
               }
               ++cl ;
               mesh = mesh->parent() ;
            } while( mesh != 0 ) ;
         }
      }
   }

   for( size_t iF = 0 ; iF<nb_fields ; ++iF )
   {
      NB_LOCAL_NODES( iF ) = il[ iF ] ;
      NODE_LOC_OK( iF ) = false ;
   }

   QR_PROVIDER = 0 ;
   QR_PROVIDER_STAT = PEL::bad_index() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEmulti:: is_in_mesh( PDE_DiscreteField const* ff,
                               size_t i,
                               PDE_BasisFunction const* bf ) const
//----------------------------------------------------------------------
{
   PEL_Error::object()->raise_not_implemented( this, "is_in_mesh" ) ;
   return( false ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: nb_nested_meshes( PDE_MeshFE const* mesh ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: nb_nested_meshes" ) ;

   size_t im = local_mesh_index( mesh ) ;
   size_t result = ( im == PEL::bad_index() ? 0 : NB_NESTED_LEVELS( im ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: reset_calculation_point( void )
//----------------------------------------------------------------------
{
   HAS_CP = false ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: reset_row_and_col_fields( void )
//----------------------------------------------------------------------
{
   ROW_iF = PEL::bad_index() ;
   COL_iF = PEL::bad_index() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: reset_IP_iterator( void )
//----------------------------------------------------------------------
{
   i_IP = PEL::bad_index() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: set_CP( GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: set_CP" ) ;

   HAS_CP = true ;

   CP_LOCATION->set( pt ) ;

   initialize_BF_vector_2( 1, C_BF_VALS, C_BF_VALS_POOL ) ;

   CP_NB_MESHES = 0 ;
   CP_MESHES.re_initialize( NB_MESHES ) ; //surdimensionne
   CP_MESHES.set( PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: set_nb_IPs( size_t aNbItgPts )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: set_nb_IPs" ) ;

   NB_IPs = aNbItgPts ;
   IP_WEIGHTS.re_initialize( NB_IPs ) ;

   initialize_IP_vector( NB_IPs, nb_space_dimensions(),
                         IP_LOCATIONS, IP_LOCATIONS_POOL ) ;
   
   MASKED_LEVEL = PEL::bad_index() ;
   
   MASKED_VALUES.re_initialize( nb_handled_fields(), 
                                handled_fields_max_nb_comps(),
                                NB_IPs ) ;

   initialize_BF_vector_2( NB_IPs, I_BF_VALS, I_BF_VALS_POOL ) ;

   IP_NB_MESHES.re_initialize( NB_IPs ) ;
   IP_MESHES.re_initialize( NB_IPs, NB_MESHES ) ; //surdimensionne
   IP_MESHES.set( PEL::bad_index() ) ;

   i_IP = PEL::bad_index() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: append_IP( GE_Point const* pt, double weight )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: append_IP" ) ;

   PEL_CHECK( pt!=0 ) ;
   PEL_CHECK( i_IP==PEL::bad_index() || i_IP<NB_IPs ) ;

   if( i_IP == PEL::bad_index() )
   { 
      i_IP = 0 ;
   }
   else
   {
      ++i_IP ;
   }
   
   GE_Point* ip_pt = static_cast<GE_Point*>( IP_LOCATIONS->at(i_IP) ) ;
   ip_pt->set( pt ) ;
   
   IP_WEIGHTS( i_IP ) = weight ;

}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: connect_CP( PDE_MeshFE const* mesh,
                               size_t coarse_level,
                               GE_Point const* pt_ref,
                               GE_Matrix const* tr_jac,
                               doubleArray3D const* hessian )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: connect_CP" ) ;
   PEL_CHECK( mesh!=0 && pt_ref!=0 ) ;

   size_t im = local_mesh_index( mesh ) ;
   if( im != PEL::bad_index() )
   {
      if( coarse_level == 0 )
      {
         CP_MESHES( CP_NB_MESHES ) = im ;
         ++CP_NB_MESHES ;
      }
      else PEL_ASSERT( CP_MESHES( CP_NB_MESHES-1 ) == im ) ;

      for( size_t e=0 ; e<mesh->nb_reference_elements() ; ++e ) 
      {
         int deri = ELM_DERIS( im, e ) ;
         if( deri > 0 )
         {
            PDE_ReferenceElement const* elm = mesh->reference_element( e ) ;
            PDE_BFvalues* bfv = C_BFvalues( im, coarse_level, e ) ;
            bfv->re_initialize( elm, deri,  pt_ref, tr_jac, hessian ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: connect_IP( PDE_MeshFE const* mesh,
                               size_t coarse_level,
                               GE_Point const* pt_ref,
                               GE_Matrix const* tr_jac,
                               doubleArray3D const* hessian )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: connect_IP" ) ;
   PEL_CHECK( mesh!=0 && pt_ref!=0 ) ;

   size_t im = local_mesh_index( mesh ) ;
   if( im != PEL::bad_index() )
   {
      if( coarse_level == 0 )
      {
         IP_MESHES( i_IP, IP_NB_MESHES( i_IP ) ) = im ;
         ++( IP_NB_MESHES( i_IP ) ) ;
      }
      else PEL_ASSERT( IP_MESHES( i_IP, IP_NB_MESHES( i_IP )-1 ) == im ) ;

      for( size_t e=0 ; e<mesh->nb_reference_elements() ; ++e ) 
      {
         int deri = ELM_DERIS( im, e ) ;
         if( deri > 0 )
         {
            PDE_ReferenceElement const* elm = mesh->reference_element( e ) ;
            PDE_BFvalues* bfv = I_BFvalues( im, coarse_level, e, i_IP ) ;
            bfv->re_initialize( elm, deri, pt_ref, tr_jac, hessian ) ;
         }
      }
   }
}

//-----------------------------------------------------------------------
GE_QuadratureRule const*
PDE_LocalFEmulti:: quadrature_rule( GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: quadrature_rule" ) ;
   PEL_CHECK_PRE( quadrature_rule_PRE( qrp ) ) ;

   GE_QuadratureRule const* result = 
                qrp->quadrature_rule( polyhedron()->reference_polyhedron() ) ;

   PEL_CHECK_POST( quadrature_rule_POST( result, qrp ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
bool
PDE_LocalFEmulti:: quadrature_rule_PRE( GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( !valid_IP() ) ;
   PEL_ASSERT( qrp!=0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
PDE_LocalFEmulti:: quadrature_rule_POST( GE_QuadratureRule const* result,
                                         GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: locate_discretization_at_pt( size_t iF,
                                                size_t& im, size_t& ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: locate_discretization_at_pt" ) ;
   PEL_CHECK( iF < nb_handled_fields() ) ;
   
   im = PEL::bad_index() ;
   ee = PEL::bad_index() ;
   size_t i = 0 ;
   for( ; i<CP_NB_MESHES ; ++i )
   {
      im = CP_MESHES( i ) ;
      ee = ELM_index( im, iF ) ;
      if( ee != PEL::bad_index() ) break ;
   }
   for( ++i ; i<CP_NB_MESHES ; ++i )
   {
      PEL_CHECK( ELM_index( CP_MESHES(i), iF ) == PEL::bad_index() ) ;
   }

   PEL_CHECK( im < NB_MESHES ) ;
   PEL_CHECK( ee < 
    static_cast< PDE_MeshFE* >( MESHES->at( im ) )->nb_reference_elements() ) ;
}

//----------------------------------------------------------------------
PDE_BFvalues*
PDE_LocalFEmulti:: C_BFvalues( size_t im, size_t cl, size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: C_BFvalues" ) ;
   PEL_CHECK( im < NB_MESHES ) ;
   PEL_CHECK( cl < NB_NESTED_LEVELS( im ) ) ;
   
   PDE_MeshFE const* mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;
   PEL_CHECK( ee < mesh->nb_reference_elements() ) ;

   PEL_Vector const* vv = static_cast<PEL_Vector*>( C_BF_VALS->at( im ) ) ;
   size_t idx = ( cl*mesh->nb_reference_elements() + ee ) ;
   PDE_BFvalues* result = static_cast<PDE_BFvalues*>( vv->at( idx ) ) ;

   PEL_CHECK( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: locate_discretization_at_IP( size_t iF, size_t ip,
                                                size_t& im, size_t& ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: locate_discretization_at_IP" ) ;
   PEL_CHECK( iF < nb_handled_fields() ) ;
   PEL_CHECK( ip < NB_IPs ) ;

   im = PEL::bad_index() ;
   ee = PEL::bad_index() ;
   size_t i = 0 ;
   for( ; i<IP_NB_MESHES( ip ) ; ++i )
   {
      im = IP_MESHES( ip, i ) ;
      ee = ELM_index( im, iF ) ;
      if( ee != PEL::bad_index() ) break ;
   }
   for( ++i ; i<IP_NB_MESHES( ip ) ; ++i )
   {
      PEL_CHECK( ELM_index( IP_MESHES(ip,i), iF ) == PEL::bad_index() ) ;
   }

   PEL_CHECK( im < NB_MESHES ) ;
   PEL_CHECK( ee < 
    static_cast< PDE_MeshFE* >( MESHES->at( im ) )->nb_reference_elements() ) ;
}

//----------------------------------------------------------------------
PDE_BFvalues*
PDE_LocalFEmulti:: I_BFvalues( size_t im, 
                               size_t cl, size_t ee, size_t ip ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: I_BFvalues" ) ;
   PEL_CHECK( im < NB_MESHES ) ;
   PEL_CHECK( cl < NB_NESTED_LEVELS( im ) ) ;
   PEL_CHECK( ip < NB_IPs ) ;

   PDE_MeshFE const* mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;
   PEL_CHECK( ee < mesh->nb_reference_elements() ) ;

   PEL_Vector const* vv = static_cast<PEL_Vector*>( I_BF_VALS->at( im ) ) ;
   size_t idx = ( cl*mesh->nb_reference_elements() + ee )*NB_IPs + ip ;
   PDE_BFvalues* result = static_cast<PDE_BFvalues*>( vv->at( idx ) ) ;

   PEL_CHECK( result != 0 ) ;
   return result ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEmulti:: local_mesh_index( PDE_MeshFE const* mesh ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: local_mesh_index" ) ;

   size_t result = PEL::bad_index() ;
   bool found = false ;
   size_t i=0 ;
   for( ; i<NB_MESHES ; ++i )
   {
      PDE_MeshFE const* mm = static_cast<PDE_MeshFE*>( MESHES->at( i ) ) ;
      if( mm == mesh )
      {
         found = true ;
	 break ;
      }
   }
   if( found ) result = i ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: build_connectivities( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: build_connectivities" ) ;

   ROW_NODES.re_initialize( NB_LOCAL_NODES( ROW_iF ) ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( ROW_iF ) ; i++ )
   {
      ROW_NODES(i) = GLOBAL_NODE( i, ROW_iF ) ;
   }
   
   COL_NODES.re_initialize( NB_LOCAL_NODES( COL_iF ) ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( COL_iF ) ; i++ )
   {
      COL_NODES(i) = GLOBAL_NODE( i, COL_iF ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: prepare_for_Ns_requests( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: prepare_for_Ns_requests" ) ;

   if( ROW_iF != PEL::bad_index() && valid_IP() )
   {
      locate_discretization_at_IP( ROW_iF, i_IP, ROW_im, ROW_ee ) ;

      build_Ns( ROW_iF, ROW_im, ROW_ee,
                ROW_Ns, ROW_dNs, ROW_d2Ns ) ;

      locate_discretization_at_IP( COL_iF, i_IP, COL_im, COL_ee ) ;

      build_Ns( COL_iF, COL_im, COL_ee,
                COL_Ns, COL_dNs, COL_d2Ns ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEmulti:: build_Ns( size_t iF, size_t im, size_t ee,
                             doubleVector& val,
                             doubleArray2D& grad,
                             doubleArray3D& grad_grad ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: build_Ns" ) ;

   size_t dim = nb_space_dimensions() ;
   size_t deri = handled_field_deri( iF ) ;
   size_t nb_bfs = NB_LOCAL_NODES( iF ) ;
   if( deri & N )
   {
      val.re_initialize( nb_bfs ) ;
   }
   if( deri & dN )
   {
      grad.re_initialize( nb_bfs, dim ) ;
   }
   if( deri & d2N )
   {
      grad_grad.re_initialize( nb_bfs, dim, dim ) ;
   }

   for( size_t i=0 ; i<nb_bfs ; ++i )
   {
      size_t idx = ELM_BF_index( i, im, iF ) ;
      if( idx != PEL::bad_index() )
      {
         size_t cl = REFINEMENT_LEVEL( i, im, iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( im, cl, ee, i_IP ) ;
         if( deri & N )
         {
            val( i ) = bfv->N_at_pt( idx ) ;
         }
         if( deri & dN )
         {
            for( size_t a=0 ; a<dim ; a++ )
               grad( i, a ) = bfv->dN_at_pt( idx, a ) ;
         }
         if( deri & d2N )
         {
            for( size_t a=0 ; a<dim ; a++ )
               for( size_t b=0 ; b<dim ; b++ )
                  grad_grad( i, a, b ) = bfv->d2N_at_pt( idx, a, b ) ;
         }
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_LocalFEmulti:: initialize_IP_vector( size_t new_length,
                                         size_t nb_dims,
                                         PEL_Vector* vec,
                                         PEL_Vector* pool )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: initialize_IP_vector" ) ;

   if( pool->index_limit() <  new_length )
   {
      size_t old_count = pool->index_limit() ;
      pool->resize( new_length ) ;
      for( size_t i=old_count ; i<new_length ; ++i )
      {
         pool->set_at( i, GE_Point::create( pool, nb_dims ) ) ;
      }
   }
   PEL_CHECK( pool->index_limit() >= new_length) ;
   
   if( vec->count() != new_length )
   {
      vec->re_initialize( new_length ) ;
      for( size_t i=0 ; i<new_length ; ++i )
      {
         vec->set_at( i, pool->at(i) ) ;
      }
   }

   PEL_CHECK( vec->count()== new_length ) ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEmulti:: initialize_ELM_BF_vectors_1( size_t new_length,
                                                PEL_Vector* vec2,
                                                PEL_Vector* vec3 )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: initialize_ELM_BF_vectors_1" ) ;
   PEL_CHECK( vec2->count() == vec3->count() ) ;

   size_t old_length = vec2->count() ;

   if( old_length != new_length )
   {
      for( size_t i=0 ; i<old_length ; ++i )
      {
         vec2->destroy_possession( vec2->at(i) ) ;
         vec3->destroy_possession( vec3->at(i) ) ;
      }
      vec2->re_initialize( new_length ) ;
      vec3->re_initialize( new_length ) ;
      for( size_t i=0 ; i<new_length ; ++i )
      {
         vec2->set_at( i, PEL_Vector::create( vec2, 0 ) ) ;
         vec3->set_at( i, PEL_Vector::create( vec3, 0 ) ) ;
      }
   }
   else
   {
      for( size_t i=0 ; i<new_length ; ++i )
      {
         static_cast<PEL_Vector*>( vec2->at(i) )->re_initialize( 0 ) ;
         static_cast<PEL_Vector*>( vec3->at(i) )->re_initialize( 0 ) ;
      }
   }

   PEL_CHECK( vec2->count() == new_length ) ;
   PEL_CHECK( vec3->count() == new_length ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<new_length ; ++i ),
                dynamic_cast<PEL_Vector*>( vec2->at(i) ) != 0 && 
                dynamic_cast<PEL_Vector*>( vec3->at(i) ) != 0  ) ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<new_length ; ++i ),
                static_cast<PEL_Vector*>( vec2->at(i) )->count()==0 ) ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<new_length ; ++i ),
                static_cast<PEL_Vector*>( vec3->at(i) )->count()==0 ) ) ;
}
   
//---------------------------------------------------------------------------
void
PDE_LocalFEmulti:: initialize_BF_vector_2( size_t nb_pts,
                                           PEL_Vector* vec,
                                           PEL_Vector* pool )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEmulti:: initialize_BF_vector_2" ) ;
   PEL_CHECK( vec->count() == NB_MESHES ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<NB_MESHES ; ++i ),
                 dynamic_cast<PEL_Vector*>( vec->at(i) ) != 0 ) ) ;

   size_t ii = 0 ;
   for( size_t im=0 ; im<NB_MESHES ; ++im )
   {
      PDE_MeshFE const* mesh = static_cast< PDE_MeshFE* >( MESHES->at( im ) ) ;
      size_t nb_elms = mesh->nb_reference_elements() ;
      size_t last_index = ii + NB_NESTED_LEVELS(im)*nb_elms*nb_pts ;
      if( pool->index_limit() < last_index )
      {
         size_t old_count = pool->index_limit() ;
         pool->resize( last_index ) ;
         for( size_t i=old_count ; i<last_index ; ++i )
         {
            pool->set_at( i, PDE_BFvalues::create( pool ) ) ;
         }
      }
      PEL_CHECK( pool->index_limit()>=last_index ) ;
      
      PEL_Vector* vv = static_cast<PEL_Vector*>( vec->at( im ) ) ;
      if( vv->count() != NB_NESTED_LEVELS(im)*nb_elms*nb_pts )
      {
         vv->re_initialize( NB_NESTED_LEVELS(im)*nb_elms*nb_pts ) ;
         for( size_t i=0 ; i<NB_NESTED_LEVELS(im)*nb_elms*nb_pts ; ++i )
         {
            vv->set_at( i, pool->at(ii) ) ;
            ++ii ;
         }
      }
   }
   
   PEL_CHECK( FORALL( ( size_t i=0 ; i<NB_MESHES ; ++i ),
                      dynamic_cast<PEL_Vector*>( vec->at( i ) ) != 0 ) ) ;
}
