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

#include <PDE_LocalFEsingle.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_BFvalues.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_MeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ; using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;
using std::vector ;

//----------------------------------------------------------------------
PDE_LocalFEsingle:: PDE_LocalFEsingle( PEL_Object* a_owner, size_t aNbSpDims )
//----------------------------------------------------------------------
   : PDE_LocalFE( a_owner, aNbSpDims )
   , NB_NESTED_LEVELS( PEL::bad_index() )
   , MAX_NB_NODES( 300 )
   , ELM_DERIS( 0 )
   , ELM_index( 0 )
   , NB_LOCAL_NODES( 0 )
   , GLOBAL_NODE( 0, 0 )
   , BASIS_FUNC( PEL_Vector::create( this, 0 ) )
   , ELM_BF_index( 0, 0 )
   , NODE_LOCATIONS( PEL_Vector::create( this, 0 ) )
   , NODE_LOC_OK( 0 ) 
   , REFINEMENT_LEVEL( 0, 0 )
   , HAS_CP( false )
   , CP_LOCATION( GE_Point::create( this, aNbSpDims ) )
   , C_BF_VALS( PEL_Vector::create( this, 0 ) )
   , QR_PROVIDER( 0 )
   , QR_PROVIDER_STAT( PEL::bad_index() )
   , QR( 0 )
   , NB_IPs( 0 )
   , IP_LOCATION( GE_Point::create( this, aNbSpDims ) )
   , IP_WEIGHTS( 0 )
   , I_BF_VALS( PEL_Vector::create( this, 0 ) )
   , i_IP( PEL::bad_index() )
   , MASKED_VALUES( 0, 0, 0 )
   , MASKED_LEVEL( PEL::bad_index() )
   , SF_2_iF( 2 )
   , ROW_NODES( 0 )
   , COL_NODES( 0 )
   , ROW_bfv( 0 )
   , COL_bfv( 0 )
   , ROW_Ns( 0 )
   , COL_Ns( 0 )
   , ROW_dNs( 0, 0 )
   , COL_dNs( 0, 0 )
   , ROW_d2Ns( 0, 0, 0 )
   , COL_d2Ns( 0, 0, 0 )
   , PT( GE_Point::create( this, aNbSpDims ) )
{
   PEL_LABEL( "PDE_LocalFEsingle:: PDE_LocalFEsingle" ) ;

   PEL_ASSERT( node == 0  ) ;
   PEL_ASSERT( N    == PDE_BFvalues::N   ) ;
   PEL_ASSERT( dN   == PDE_BFvalues::dN  ) ;
   PEL_ASSERT( d2N  == PDE_BFvalues::d2N ) ;
   
   PEL_CHECK( aNbSpDims>=1 && aNbSpDims<=3 ) ;

   reset_row_and_col_fields() ;

   reset_IP_iterator() ;
}

//----------------------------------------------------------------------
PDE_LocalFEsingle:: ~PDE_LocalFEsingle( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEsingle:: nb_local_nodes( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: nb_local_nodes" ) ;
   PEL_CHECK_PRE( nb_local_nodes_PRE( ff ) ) ;

   size_t result = NB_LOCAL_NODES( field_local_index( ff ) ) ;

   PEL_CHECK_POST( nb_local_nodes_POST( result, ff ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEsingle:: local_node_location( PDE_DiscreteField const* ff,
                                         size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: local_node_location" ) ;
   PEL_CHECK_PRE( local_node_location_PRE( ff, i ) ) ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   PDE_ReferenceElement const* elm = THE_MESH->reference_element( ee ) ;
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
      PDE_MeshFE* mesh = THE_MESH ;
      size_t il = 0 ;
      for( size_t cl=0 ; cl<NB_NESTED_LEVELS ; ++cl )
      {
         GE_Mpolyhedron const* poly = mesh->polyhedron() ;
//         bool found = false ;
         for( size_t ii=0 ; ii<mesh->nb_basis_functions(ee) ; ++ii )
         {
            PDE_BasisFunction* bf = mesh->basis_function( ee, ii ) ;
            if( ( bf != 0 ) && ( ff->node_is_active( bf->node_of_DOF( ff ))))
            {
//               found = true ;
               GE_Point const* pt_ref = elm->node_location( ii ) ;
               GE_Point* pt = static_cast< GE_Point* >( ff_nodes->at( il ) ) ;
               poly->apply_mapping( pt_ref, pt ) ;
               ++il ;
            }
         }
         PEL_ASSERT( cl==NB_NESTED_LEVELS-1 || mesh->parent()!=0 ) ;
         mesh = mesh->parent() ;
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
PDE_LocalFEsingle:: local_node_is_in_mesh( PDE_DiscreteField const* ff,
                                           size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: local_node_is_in_mesh" ) ;
   PEL_CHECK_PRE( local_node_is_in_mesh_PRE( ff, i ) ) ;

   size_t ii = i * nb_handled_fields() + field_local_index( ff ) ;
   PDE_BasisFunction const* bf = 
                     static_cast<PDE_BasisFunction*>( BASIS_FUNC->at( ii ) ) ;

   bool result = is_in_mesh( bf ) ;

   PEL_CHECK_POST( local_node_is_in_mesh_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEsingle:: global_node( PDE_DiscreteField const* ff,
                                 size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: global_node" ) ;
   PEL_CHECK_PRE( global_node_PRE( ff, i ) ) ;

   size_t result = GLOBAL_NODE( i, field_local_index( ff ) ) ;

   PEL_CHECK_POST( global_node_POST( result, ff, i ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEsingle:: node_refinement_level( PDE_DiscreteField const* ff,
                                           size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: node_refinement_level" ) ;
   PEL_CHECK_PRE( node_refinement_level_PRE( ff, i ) ) ;

   size_t ii = i * nb_handled_fields() + field_local_index( ff ) ;
   PDE_BasisFunction const* bf = 
                     static_cast<PDE_BasisFunction*>( BASIS_FUNC->at( ii ) ) ;
   size_t result = bf->refinement_level() ;

   PEL_CHECK_POST( node_refinement_level_POST( result, ff, i ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEsingle:: calculation_point( void ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: calculation_point" ) ;

   GE_Point const* result = 0 ;
   if( HAS_CP )
   {
      result = CP_LOCATION ;
   }

   PEL_CHECK_POST( calculation_point_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: value_at_pt( PDE_DiscreteField const* ff,
                                 size_t level,
                                 size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: value_at_pt" ) ;
   PEL_CHECK_PRE( value_at_pt_PRE( ff, level, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
      result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                bfv->N_at_pt( ELM_BF_index( i, iF ) ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: gradient_at_pt( PDE_DiscreteField const* ff,
                                    size_t level,
                                    size_t a,
                                    size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: gradient_at_pt" ) ;
   PEL_CHECK_PRE( gradient_at_pt_PRE( ff, level, a, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
      result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                bfv->dN_at_pt( ELM_BF_index( i, iF ), a ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: hessian_at_pt( PDE_DiscreteField const* ff,
                                   size_t level,
                                   size_t a,
                                   size_t b,
                                   size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: hessian_at_pt" ) ;
   PEL_CHECK_PRE( hessian_at_pt_PRE( ff, level, a, b, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
      result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                bfv->d2N_at_pt( ELM_BF_index( i, iF ), a, b ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: N_at_pt( PDE_DiscreteField const* ff,
                             size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: N_at_pt" ) ;
   PEL_CHECK_PRE( N_at_pt_PRE( ff, i ) ) ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
   double result = bfv->N_at_pt( il ) ;

   return( result );
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: dN_at_pt( PDE_DiscreteField const* ff,
                              size_t i, 
                              size_t a ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: dN_at_pt" ) ;
   PEL_CHECK_PRE( dN_at_pt_PRE( ff, i, a ) ) ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
   double result = bfv->dN_at_pt( il, a ) ;

   return( result );
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: d2N_at_pt( PDE_DiscreteField const* ff,
                               size_t i, 
                               size_t a,
                               size_t b ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: d2N_at_pt" ) ;
   PEL_CHECK_PRE( d2N_at_pt_PRE( ff, i, a, b ) ) ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = C_BFvalues( cl, ee ) ;
   double result = bfv->d2N_at_pt( il, a, b ) ;
   
   return( result );
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: set_row_and_col_fields( 
                                      PDE_DiscreteField const* row_field, 
                                      PDE_DiscreteField const* col_field )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: set_row_and_col_fields" ) ;
   PEL_CHECK_PRE( set_row_and_col_fields_PRE( row_field, col_field ) ) ;

   SF_2_iF( row ) = field_local_index( row_field ) ;
   SF_2_iF( col ) = field_local_index( col_field ) ;

   build_connectivities() ;

   PEL_CHECK_POST( set_row_and_col_fields_POST( row_field, col_field ) ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_LocalFEsingle:: field( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: field" ) ;

   PDE_DiscreteField const* result = 0 ;
   size_t iF = SF_2_iF( sf ) ;
   
   if( iF != PEL::bad_index() )
   {
      result = handled_field( iF ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEsingle:: row_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: row_field_node_connectivity" ) ;
   PEL_CHECK_PRE( row_field_node_connectivity_PRE() ) ;

   size_t_vector const& result = ROW_NODES ;

   PEL_CHECK_POST( row_field_node_connectivity_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LocalFEsingle:: col_field_node_connectivity( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: col_field_node_connectivity" ) ;
   PEL_CHECK_PRE( col_field_node_connectivity_PRE() ) ;

   size_t_vector const& result = COL_NODES ;

   PEL_CHECK_POST( col_field_node_connectivity_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEsingle:: nb_basis_functions( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: nb_basis_functions" ) ;
   PEL_CHECK_PRE( nb_basis_functions_PRE( sf ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t result = NB_LOCAL_NODES( SF_2_iF( sf ) ) ;

   PEL_CHECK_POST( nb_basis_functions_POST( result, sf ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: start_IP_iterator( GE_QRprovider const* qrp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: start_IP_iterator" ) ;
   PEL_CHECK_PRE( start_IP_iterator_PRE( qrp ) ) ;

   // si c'est le meme provider et qu'on n'a pas change de maille,
   // on va avoir la meme règle de quadrature, donc il est inutile
   // de refaire les calculs
   if( qrp != QR_PROVIDER || qrp->status() != QR_PROVIDER_STAT )
   {
      if( QR != 0 && QR->owner() == this ) destroy_possession(QR) ;
      QR = quadrature_rule( qrp ) ;
      if( QR->owner() == 0 )
      {
         const_cast<GE_QuadratureRule*>( QR )->set_owner(this) ;
      }
      compute_itg_pts( QR ) ;
      QR_PROVIDER = qrp ;
      QR_PROVIDER_STAT = qrp->status() ;
   }
   i_IP = 0 ;

   prepare_for_Ns_requests() ;

   PEL_CHECK_POST( start_IP_iterator_POST( qrp ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: go_next_IP( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: go_next_IP" ) ;
   PEL_CHECK_PRE( go_next_IP_PRE() ) ;

   ++i_IP ;

   prepare_for_Ns_requests() ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFEsingle:: valid_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: valid_IP" ) ;
   PEL_CHECK_PRE( valid_IP_PRE() ) ;

   return( i_IP<NB_IPs ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: weight_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: weight_of_IP" ) ;
   PEL_CHECK_PRE( weight_of_IP_PRE() ) ;

   return( IP_WEIGHTS( i_IP ) ) ;
}

//----------------------------------------------------------------------
GE_Point const*
PDE_LocalFEsingle:: coordinates_of_IP( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: coordinates_of_IP" ) ;
   PEL_CHECK_PRE( coordinates_of_IP_PRE() ) ;

   polyhedron()->apply_mapping( QR->point( i_IP ), IP_LOCATION ) ;
   GE_Point* result = IP_LOCATION ;
   
   PEL_CHECK_POST( coordinates_of_IP_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: value_at_IP( PDE_DiscreteField const* ff,
                                 size_t level,
                                 size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: value_at_IP" ) ;
   PEL_CHECK_PRE( value_at_IP_PRE( ff, level, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   if( level == MASKED_LEVEL )
   {
      result = MASKED_VALUES( iF, ic, i_IP ) ;
   }
   else
   {
      size_t ee = ELM_index( iF ) ;
      for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
      {
         size_t cl = REFINEMENT_LEVEL( i, iF ) ;
         PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;
         result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                   bfv->N_at_pt( ELM_BF_index( i, iF ) ) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: gradient_at_IP( PDE_DiscreteField const* ff,
                                    size_t level,
                                    size_t a, 
                                    size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: gradient_at_IP" ) ;
   PEL_CHECK_PRE( gradient_at_IP_PRE( ff, level, a, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;
      result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                bfv->dN_at_pt( ELM_BF_index( i, iF ), a ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: hessian_at_IP( PDE_DiscreteField const* ff,
                                   size_t level,
                                   size_t a, 
                                   size_t b, 
                                   size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: hessian_at_IP" ) ;
   PEL_CHECK_PRE( hessian_at_IP_PRE( ff, level, a, b, ic ) ) ;

   double result = 0.0 ;

   size_t iF = field_local_index( ff ) ;
   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( iF ) ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;
      result += ff->DOF_value( level, GLOBAL_NODE( i, iF ), ic ) *
                bfv->d2N_at_pt( ELM_BF_index( i, iF ), a, b ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: N_at_IP( field_id sf, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: N_at_IP" ) ;
   PEL_CHECK_PRE( N_at_IP_PRE( sf, i ) ) ;
   
   size_t iF = SF_2_iF( sf ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;

   double result = bfv->N_at_pt( il ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: dN_at_IP( field_id sf, size_t i, size_t a )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: dN_at_IP" ) ;
   PEL_CHECK_PRE( dN_at_IP_PRE( sf, i, a ) ) ;

   size_t iF = SF_2_iF( sf ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;

   double result = bfv->dN_at_pt( il, a ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_LocalFEsingle:: d2N_at_IP( field_id sf,
                               size_t i, size_t a, size_t b )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: d2N_at_IP" ) ;
   PEL_CHECK_PRE( d2N_at_IP_PRE( sf, i, a, b ) ) ;
 
   size_t iF = SF_2_iF( sf ) ;
   size_t ee = ELM_index( iF ) ;
   size_t cl = REFINEMENT_LEVEL( i, iF ) ;
   size_t il = ELM_BF_index( i, iF ) ;
   PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;

   double result = bfv->d2N_at_pt( il, a, b ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PDE_LocalFEsingle:: Ns_at_IP( field_id sf )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: Ns_at_IP" ) ;
   PEL_CHECK_PRE( Ns_at_IP_PRE( sf ) ) ;

   doubleVector const* result = 0 ;
   if( NB_NESTED_LEVELS == 1 )
   {
      result = ( sf==row ? &( ROW_bfv->Ns_at_pt() ) 
                         : &( COL_bfv->Ns_at_pt() ) ) ;
   }
   else
   {
      result = ( sf==row ? &ROW_Ns : &COL_Ns ) ;
   }

   PEL_CHECK_POST( Ns_at_IP_POST( *result, sf ) ) ;
   return( *result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PDE_LocalFEsingle:: dNs_at_IP( field_id sf )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: dNs_at_IP" ) ;
   PEL_CHECK_PRE( dNs_at_IP_PRE( sf ) ) ;

   doubleArray2D const* result = 0 ;
   if( NB_NESTED_LEVELS == 1 )
   {
      result = ( sf==row ? &( ROW_bfv->dNs_at_pt() ) 
                         : &( COL_bfv->dNs_at_pt() ) ) ;
   }
   else
   {
      result = ( sf==row ? &ROW_dNs : &COL_dNs ) ;
   }

   PEL_CHECK_POST( dNs_at_IP_POST( *result, sf ) ) ;
   return( *result ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PDE_LocalFEsingle:: d2Ns_at_IP( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: d2Ns_at_IP" ) ;
   PEL_CHECK_PRE( d2Ns_at_IP_PRE( sf ) ) ;

   doubleArray3D const* result = 0 ;

   if( NB_NESTED_LEVELS == 1 )
   {
      result = ( sf==row ? &( ROW_bfv->d2Ns_at_pt() ) 
                         : &( COL_bfv->d2Ns_at_pt() ) ) ;
   }
   else
   {
      result = ( sf==row ? &ROW_d2Ns : &COL_d2Ns ) ;
   }

   PEL_CHECK_POST( d2Ns_at_IP_POST( *result, sf ) ) ;
   return( *result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: mask_value_at_IP( PDE_DiscreteField const* ff,
                                      size_t level,
                                      double value,
                                      size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: mask_value_at_IP" ) ;
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
PDE_LocalFEsingle:: print_local_discretization_of_fields(
                                               std::ostream& os, 
                                               size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: print_local_discretization_of_fields" ) ;
   PEL_CHECK_PRE( print_local_discretization_of_fields_PRE( os,
                                                            indent_width ) ) ;

   PDE_LocalFE::print_local_discretization_of_fields( os, indent_width ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;

//???????????? !!!!!!!!!!! ??????????????????????
   PDE_MeshFE const* the_mesh = 0 ;
   for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
   {
      PDE_DiscreteField const* ff = handled_field( iF ) ;

      PDE_MeshFE const* lf = THE_MESH ;
      if( the_mesh == 0 )
      {
         the_mesh = lf ;
         os << space << the_mesh->polyhedron()->name() 
            << " " << the_mesh->id_number() << endl ;
      }
      else
      {
         PEL_ASSERT( lf == the_mesh ) ;
      }
      os << space << more << "\"" << ff->name() << "\" : " ;

      if( lf->has_discretization( ff ) )
      {
         size_t ee = lf->index_of_reference_element( ff ) ;
         os << lf->reference_element( ee )->name() << endl ;
         for( size_t cl=1 ; cl<NB_NESTED_LEVELS ; ++cl )
         {
            lf = lf->parent() ;
            PEL_ASSERT( lf != 0 ) ;
            ee = lf->index_of_reference_element( ff ) ;
            std::string spp( ff->name().size()+8, ' ' ) ;
            os << space << spp << lf->reference_element( ee )->name() 
               << " on parent " << lf->polyhedron()->name() 
               << " " << lf->id_number() << endl ;
         }
      }
      else
      {
         os << " no discretization" << endl ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: print_current_IP( std::ostream& os, 
                                      size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: print_current_IP" ) ;
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
      PDE_DiscreteField const* ff = handled_field( iF ) ;
      os << space << more << "\"" << ff->name() 
         << "\" : discretization of " ;
//      PDE_MeshFE const* lf = LOCAL_FIELDS[iF][0] ;
      PDE_MeshFE const* lf = THE_MESH ;
      os << lf->polyhedron()->name()
         << " " << lf->id_number() << endl ;
      
      if( handled_field_deri(iF) & PDE_LocalFE::N )
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

   os << std::setprecision( p ) ;
   os.flags( original_flags ) ;
}

//????????? comme PDE_LocalFEmulti.... remonter dans PDE_LocalFE
//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: print_values_at_current_IP( std::ostream& os, 
                                                size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: print_values_at_current_IP" ) ;
   PEL_CHECK_PRE( print_values_at_current_IP_PRE( os, indent_width ) ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;

   ios_base::fmtflags original_flags = os.flags() ;
   os.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision( 6 ) ;
   
   os << space << "IP : " ;
   coordinates_of_IP()->print( os, 0 ) ;
   os << " weight : " << weight_of_IP() << endl ;

   if( SF_2_iF(0) != SF_2_iF(1) )
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
   if( handled_field_deri(SF_2_iF(0)) & N ) os << setw( 15 ) << "N      " ;
   if( handled_field_deri(SF_2_iF(0)) & dN )
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
      if( handled_field_deri(SF_2_iF(0)) & N )
      {
         double n = N_at_IP( row, i ) ;
         if( PEL::abs(n)<1.E-12 ) n = 0. ;
         os  << setw( 15 ) << n ;
      }
      if( handled_field_deri(SF_2_iF(0)) & dN )
	 for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
	 {
            double dn = dN_at_IP( row, i, a ) ;
            if( PEL::abs(dn)<1.E-12 ) dn = 0. ;
            os << setw( 15 ) << dn ;
         }
      os << endl ;
   }

   if( SF_2_iF(0) != SF_2_iF(1) )
   {
      os << space << "col field : \"" << field(col)->name() << "\"" << endl ;
      os << space << more << "   " ;
      if( handled_field_deri(SF_2_iF(1)) & N ) os << setw( 15 ) << "N      " ;
      if( handled_field_deri(SF_2_iF(1)) & dN )
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
         if( handled_field_deri(SF_2_iF(1)) & N )
         {
            double n = N_at_IP( col, i ) ;
            if( PEL::abs(n)<1.E-12 ) n = 0. ;
            os <<  setw( 15 ) << n ;
         }
         if( handled_field_deri(SF_2_iF(1)) & dN )
            for( size_t a=0 ; a<polyhedron()->nb_space_dimensions() ; ++a )
            {
               double dn = dN_at_IP( col, i, a ) ;
               if( PEL::abs(dn)<1.E-12 ) dn = 0. ;
               os << setw( 15 ) << dn ;
            }
         os << endl ;
      }
   }

   os << std::setprecision( p ) ;
   os.flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: set_single_mesh( PDE_MeshFE* a_mesh  )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: set_single_mesh" ) ;
   PEL_CHECK( a_mesh != 0 ) ;

   THE_MESH = a_mesh ;
   PDE_MeshFE* mesh = THE_MESH ;

   size_t const nb_fields = nb_handled_fields() ;
   size_t const nb_ref_elms = mesh->nb_reference_elements() ;
  
   // Resize inner table only if a new handle field has been declared:
   size_t const old_nb_fields = ELM_index.size() ;
   if( nb_fields > old_nb_fields )
   {
      ELM_index.re_initialize( nb_fields ) ;
      NB_LOCAL_NODES.re_initialize( nb_fields ) ;
      GLOBAL_NODE.re_initialize( MAX_NB_NODES, nb_fields ) ;
      BASIS_FUNC->re_initialize( MAX_NB_NODES*nb_fields ) ;
      ELM_BF_index.re_initialize( MAX_NB_NODES, nb_fields ) ;
      REFINEMENT_LEVEL.re_initialize( MAX_NB_NODES, nb_fields ) ;
      for( size_t i=old_nb_fields ; i<nb_fields ; ++i )
      {
         NODE_LOCATIONS->append( PEL_Vector::create( NODE_LOCATIONS, 0 ) ) ;
      }
      NODE_LOC_OK.re_initialize( nb_fields) ;
   }

   // Resize new number of reference elements:
   size_t const old_nb_ref_elms = ELM_DERIS.size() ;
   if( nb_ref_elms > old_nb_ref_elms )
   {
      ELM_DERIS.re_initialize( nb_ref_elms ) ;
   }
  
   // Set inner tables:
   ELM_DERIS.set( -1 ) ;
   NB_NESTED_LEVELS = 0 ;
   NB_LOCAL_NODES.set( 0 ) ;
   NODE_LOC_OK.set( false ) ;
   for( size_t iF=0 ; iF<nb_fields ; ++iF )
   {
      mesh = a_mesh ;
      PDE_DiscreteField const* ff = handled_field( iF ) ;

      size_t const ee = mesh->index_of_reference_element( ff ) ;
      if( ee != PEL::bad_index() )
      {
         if( ELM_DERIS( ee ) == -1 )
         {
            ELM_DERIS( ee ) = handled_field_deri( iF ) ;
         }
         else
         {
            ELM_DERIS( ee ) |= handled_field_deri( iF ) ;
         }
         ELM_index( iF ) = ee ;

         size_t cl = 0 ;
         size_t il = 0 ;
      
         do
         {
            size_t const nb_bfs = mesh->nb_basis_functions(ee) ;
            for( size_t ii=0 ; ii<nb_bfs ; ++ii )
            {
               PDE_BasisFunction* bf = mesh->basis_function( ee, ii ) ;
               if( ( bf != 0 ) && ( ff->node_is_active( bf->node_of_DOF( ff ))))
               {
                  GLOBAL_NODE( il, iF ) = bf->node_of_DOF( ff ) ;
                  BASIS_FUNC->set_at( il*nb_fields+iF, bf ) ;
                  ELM_BF_index( il, iF ) =  ii ;
                  REFINEMENT_LEVEL( il, iF ) = cl ;
                  if( cl+1 > NB_NESTED_LEVELS ) NB_NESTED_LEVELS = cl+1 ;
                  ++il ; 
                  if( il == MAX_NB_NODES ) 
                  {
                     MAX_NB_NODES += 10 ;
                     GLOBAL_NODE.raise_first_index_bound( MAX_NB_NODES ) ;
                     BASIS_FUNC->resize( MAX_NB_NODES*nb_fields ) ;
                     ELM_BF_index.raise_first_index_bound( MAX_NB_NODES ) ;
                     REFINEMENT_LEVEL.raise_first_index_bound( MAX_NB_NODES ) ;
                  }
               }
            }
            ++cl ;
            mesh = mesh->parent() ;
         } while( mesh != 0 ) ;

         NB_LOCAL_NODES( iF ) = il ;
      }
   }

   QR_PROVIDER = 0 ;
   QR_PROVIDER_STAT = PEL::bad_index() ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFEsingle:: nb_nested_meshes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: nb_nested_meshes" ) ;

   return( NB_NESTED_LEVELS ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: reset_row_and_col_fields( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: reset_row_and_col_fields" ) ;

   SF_2_iF.set( PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: reset_calculation_point( void )
//----------------------------------------------------------------------
{
   HAS_CP = false ;
}

//-----------------------------------------------------------------------
void
PDE_LocalFEsingle:: set_CP( GE_Point const* pt )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: set_CP" ) ;

   HAS_CP = true ;

   CP_LOCATION->set( pt ) ;

   init_C_BFvalues( NB_NESTED_LEVELS, THE_MESH->nb_reference_elements() ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: connect_CP( size_t coarse_level,
                                GE_Point const* pt_ref,
                                GE_Matrix const* tr_jac,
                                doubleArray3D const* hessian )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: connect_CP" ) ;

   for( size_t e = 0 ; e<THE_MESH->nb_reference_elements() ; ++e )
   {
      int deri = ELM_DERIS( e ) ;
      if( deri > 0 )
      {
         PDE_ReferenceElement const* elm = THE_MESH->reference_element( e ) ;
         PDE_BFvalues* bfv = C_BFvalues( coarse_level, e ) ;
         bfv->re_initialize( elm, deri, pt_ref, tr_jac, hessian ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: reset_IP_iterator( void )
//----------------------------------------------------------------------
{
   i_IP = PEL::bad_index() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: set_nb_IPs( size_t a_nb_ips )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: set_nb_IPs" ) ;
   PEL_CHECK( a_nb_ips > 0 ) ;
   PEL_CHECK( THE_MESH != 0 ) ;

   NB_IPs = a_nb_ips ;
   IP_WEIGHTS.re_initialize( NB_IPs ) ;
   
   MASKED_LEVEL = PEL::bad_index() ;
   
   MASKED_VALUES.re_initialize( nb_handled_fields(), 
                                handled_fields_max_nb_comps(),
                                NB_IPs ) ;
 
   init_I_BFvalues( NB_NESTED_LEVELS, 
                    THE_MESH->nb_reference_elements(), 
                    NB_IPs ) ;

   i_IP = PEL::bad_index() ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: append_IP( double weight )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: append_IP" ) ;
   PEL_CHECK( i_IP==PEL::bad_index() || i_IP<NB_IPs ) ;

   if( i_IP == PEL::bad_index() )
   { 
      i_IP = 0 ;
   }
   else
   {
      ++i_IP ;
   }
   IP_WEIGHTS( i_IP ) = weight ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: connect_IP( size_t coarse_level,
                                GE_Point const* pt_ref,
                                GE_Matrix const* tr_jac,
                                doubleArray3D const* hessian )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: connect_IP" ) ;

   for( size_t e = 0 ; e<THE_MESH->nb_reference_elements() ; e++ )
   {
      int deri = ELM_DERIS( e ) ;
      if( deri > 0 )
      {
         PDE_ReferenceElement const* elm = THE_MESH->reference_element( e ) ;
         PDE_BFvalues* bfv = I_BFvalues( coarse_level, e, i_IP ) ;
         //??? tout est recalcule, meme pour les fonctions de base inactives
         bfv->re_initialize( elm, deri, pt_ref, tr_jac, hessian ) ;
      }
   }
}

//-----------------------------------------------------------------------
GE_QuadratureRule const*
PDE_LocalFEsingle:: quadrature_rule( GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: quadrature_rule" ) ;
   PEL_CHECK_PRE( quadrature_rule_PRE( qrp ) ) ;
   
   GE_QuadratureRule const* result = 
                qrp->quadrature_rule( polyhedron()->reference_polyhedron() ) ;
   
   PEL_CHECK_POST( quadrature_rule_POST( result, qrp ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
bool
PDE_LocalFEsingle:: quadrature_rule_PRE( GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( !valid_IP() ) ;
   PEL_ASSERT( qrp!=0 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
PDE_LocalFEsingle:: quadrature_rule_POST( GE_QuadratureRule const* result,
                                          GE_QRprovider const* qrp ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: init_C_BFvalues( size_t nb_cl, size_t nb_ee )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: init_C_BFvalues" ) ;

   size_t old_nb_bfvs = C_BF_VALS->index_limit() ;
   size_t nb_bfvs = nb_cl * nb_ee ;
   if( nb_bfvs > old_nb_bfvs )
   {
      C_BF_VALS->resize( nb_bfvs ) ;
      for( size_t i=old_nb_bfvs ; i<nb_bfvs ; ++i )
      {
         C_BF_VALS->set_at( i, PDE_BFvalues::create( C_BF_VALS ) ) ;
      }
   }
}

//----------------------------------------------------------------------
PDE_BFvalues*
PDE_LocalFEsingle:: C_BFvalues( size_t cl, size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: C_BFvalues" ) ;
   PEL_CHECK( cl < NB_NESTED_LEVELS ) ;
   PEL_CHECK( ee < THE_MESH->nb_reference_elements() ) ;

   size_t idx = cl*THE_MESH->nb_reference_elements() + ee  ;   
   PDE_BFvalues* result = static_cast<PDE_BFvalues*>( C_BF_VALS->at( idx ) ) ;

   PEL_CHECK( result != 0 ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: init_I_BFvalues( size_t nb_cl, size_t nb_ee, size_t nb_ip )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: init_I_BFvalues" ) ;

   size_t old_nb_bfvs = I_BF_VALS->index_limit() ;
   size_t nb_bfvs = nb_cl * nb_ee * nb_ip ;
   if( nb_bfvs > old_nb_bfvs )
   {
      I_BF_VALS->resize( nb_bfvs ) ;
      for( size_t i=old_nb_bfvs ; i<nb_bfvs ; ++i )
      {
         I_BF_VALS->set_at( i, PDE_BFvalues::create( I_BF_VALS ) ) ;
      }
   }
}

//----------------------------------------------------------------------
PDE_BFvalues* 
PDE_LocalFEsingle:: I_BFvalues( size_t cl, size_t ee, size_t ip ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: I_BFvalues" ) ;
   PEL_CHECK( cl < NB_NESTED_LEVELS ) ;   
   PEL_CHECK( ee < THE_MESH->nb_reference_elements() ) ;
   PEL_CHECK( ip < NB_IPs ) ;

   size_t idx = ( cl*THE_MESH->nb_reference_elements() + ee )*NB_IPs + ip ;   
   PDE_BFvalues* result = static_cast<PDE_BFvalues*>( I_BF_VALS->at( idx ) ) ;

   PEL_CHECK( result != 0 ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: build_connectivities( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: build_connectivities" ) ;

   ROW_NODES.re_initialize( NB_LOCAL_NODES( SF_2_iF( 0 ) ) ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( SF_2_iF( 0 ) ) ; i++ )
   {
      ROW_NODES(i) = GLOBAL_NODE( i, SF_2_iF( 0 ) ) ;
   }
   
   COL_NODES.re_initialize( NB_LOCAL_NODES( SF_2_iF( 1 ) ) ) ;
   for( size_t i=0 ; i<NB_LOCAL_NODES( SF_2_iF( 1 ) ) ; i++ )
   {
      COL_NODES(i) = GLOBAL_NODE( i, SF_2_iF( 1 ) ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: prepare_for_Ns_requests( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: prepare_for_Ns_requests" ) ;

   if( SF_2_iF( row ) != PEL::bad_index() && valid_IP() )
   {
      if( NB_NESTED_LEVELS == 1 )
      {
         ROW_bfv = I_BFvalues( 0, ELM_index( SF_2_iF(row) ), i_IP ) ;
         COL_bfv = I_BFvalues( 0, ELM_index( SF_2_iF(col) ), i_IP ) ;
      }
      else
      {
         build_Ns( SF_2_iF(row), ROW_Ns, ROW_dNs, ROW_d2Ns ) ;
         build_Ns( SF_2_iF(col), COL_Ns, COL_dNs, COL_d2Ns ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEsingle:: build_Ns( size_t iF, 
                              doubleVector& val, 
                              doubleArray2D& grad,
                              doubleArray3D& grad_grad )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEsingle:: build_Ns" ) ;

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

   size_t ee = ELM_index( iF ) ;
   for( size_t i=0 ; i<nb_bfs ; ++i )
   {
      size_t cl = REFINEMENT_LEVEL( i, iF ) ;
      size_t il = ELM_BF_index( i, iF ) ;
      PDE_BFvalues const* bfv = I_BFvalues( cl, ee, i_IP ) ;
      if( deri & N )
      {
         val( i ) = bfv->N_at_pt( il ) ;
      }
      if( deri & dN )
      {
         for( size_t a=0 ; a<dim ; ++a )
            grad( i, a ) = bfv->dN_at_pt( il, a ) ;
      }
      if( deri & d2N )
      {
         for( size_t a=0 ; a<dim ; ++a )
            for( size_t b=0 ; b<dim ; ++b )
               grad_grad( i, a, b ) = bfv->d2N_at_pt( il, a, b ) ;
      }
   }
}

