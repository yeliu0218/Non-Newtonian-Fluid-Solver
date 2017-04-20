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

#include <PDE_FluxCellIndicator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>
#include <string>

using std::string ;
using std::cout ; using std::endl ;

PDE_FluxCellIndicator const*
PDE_FluxCellIndicator:: PROTOTYPE = new PDE_FluxCellIndicator() ;

//---------------------------------------------------------------------------
PDE_FluxCellIndicator:: PDE_FluxCellIndicator( void )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( "PDE_FluxCellIndicator" )
   , CELL_ERRORS( 0, 0 )
{
}

//---------------------------------------------------------------------------
PDE_FluxCellIndicator*
PDE_FluxCellIndicator:: create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FluxCellIndicator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp, a_verbose_level ) ) ;

   PDE_FluxCellIndicator* result = 
       new PDE_FluxCellIndicator( a_owner, dom, exp, a_verbose_level ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, 
                                   dom, exp, a_verbose_level ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_FluxCellIndicator:: PDE_FluxCellIndicator( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp,
                                           size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( a_owner, a_verbose_level )
   , NB_REFS( 10000 )
   , PT( GE_Point::create( this, dom->nb_space_dimensions() ) )
   , ICALL( 0 )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data( "field" ) ) )
   , L_UU( exp->int_data( "level_of_field" ) )
   , BCs( dom->set_of_boundary_conditions() )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP( GE_QRprovider::object( 
                           exp->string_data( "quadrature_rule_provider" ) ) )
   , NB_DIMS( dom->nb_space_dimensions() )
   , CELL_ERRORS( 0, 0 )
   , MAX_ERR( exp->double_data( "maximum_error" ) )
   , MIN_ERR( -1.E+10 )
{
   sFE->include_color( GE_Color::halo_color() ) ;
   cFE->include_color( GE_Color::halo_color() ) ;
   bFE->include_color( GE_Color::halo_color() ) ;

   sFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;

   CELL_ERRORS.re_initialize( cFE->nb_meshes(), UU->nb_components() ) ;

   if( exp->has_entry( "minimum_error" ) ) 
         MIN_ERR = exp->double_data( "minimum_error" ) ;

   if( exp->has_entry( "max_nb_steps" ) )
         NB_REFS = exp->int_data( "max_nb_steps" ) ;
}

//---------------------------------------------------------------------------
PDE_FluxCellIndicator:: ~PDE_FluxCellIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_FluxCellIndicator:: reset( void )
//---------------------------------------------------------------------------
{
   ICALL = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_FluxCellIndicator:: build( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FluxCellIndicator:: build" ) ;
   
   ++ICALL ;
   CELL_ERRORS.set( 0.0 ) ;
   
   size_t nbc = UU->nb_components() ;

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Vector const* normal = ( sFE->is_periodic() ? 
                                  sFE->normal( 0 ) : sFE->normal() ) ;
      size_t id_0 = fe_K->mesh_id() ;
      if( id_0 >= CELL_ERRORS.index_bound( 0 ) ) 
      {
         CELL_ERRORS.raise_first_index_bound( id_0 + 1 ) ;
      }
      size_t id_1 = fe_L->mesh_id() ;
      if( id_1 >= CELL_ERRORS.index_bound( 0 ) )
      {
         CELL_ERRORS.raise_first_index_bound( id_1 + 1) ;
      }
      
      sFE->start_IP_iterator( QRP ) ;
      doubleVector jump( nbc ) ;
      for( ; sFE->valid_IP() ; sFE->go_next_IP() )
      {
         double ww = ( sFE->is_periodic() ? 
                       sFE->weight_of_IP( 0 ) : sFE->weight_of_IP() ) ;
         for( size_t ic=0 ; ic < UU->nb_components() ; ++ic )
         {
            double gg = 0 ;
            for( size_t d=0 ; d<NB_DIMS ; ++d )
            {
               gg += ( fe_K->gradient_at_IP( UU, L_UU, d, ic ) -
                       fe_L->gradient_at_IP( UU, L_UU, d, ic ) ) 
                  * normal->component( d ) ;
            }
            jump( ic ) += 0.25 * gg*gg * ww ;
         }
      }
      
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         CELL_ERRORS( id_0, ic ) += jump( ic ) ;
         CELL_ERRORS( id_1, ic ) += jump( ic ) ;
      }
   }
}

//---------------------------------------------------------------------------
bool
PDE_FluxCellIndicator:: to_be_refined( double bf_indicator,
                                       GE_Mpolyhedron const* poly,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FluxCellIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( bf_indicator, poly, elm, local_node ) ) ;

   bool result = false ;
   if( ICALL <= NB_REFS )
   {
      if( bf_indicator > MAX_ERR ) result = true ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PDE_FluxCellIndicator:: to_be_unrefined( double bf_indicator,
                                         GE_Mpolyhedron const* poly,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FluxCellIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( bf_indicator, poly, elm, local_node ) );

   bool result = false ;
   if( MIN_ERR > 0.0 )
   {
      if( ICALL <= NB_REFS )
      {
         if( bf_indicator < MIN_ERR ) result = true ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
double
PDE_FluxCellIndicator:: cell_indicator( size_t cell_id ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FluxCellIndicator:: cell_indicator" ) ;

   //??? these computation may be stored (lazy evaluation) if 
   //??? this function is to be called many times for the same cell_id 
   //??? between two calls of build
   
   double max_cell_err = CELL_ERRORS( cell_id, 0 ) ;
   for( size_t ic=1 ; ic<UU->nb_components() ; ++ic )
   {
      double xx = CELL_ERRORS( cell_id, ic ) ;
      if( xx > max_cell_err ) max_cell_err = xx ;
   }
   
   cFE->go_i_th( cell_id ) ;
   PEL_ASSERT( cFE->is_valid() ) ;

   double hh = cFE->polyhedron()->inter_vertices_maximum_distance() ;
   double result = PEL::sqrt( hh/24.0 * max_cell_err ) ;
   return( result ) ;
}
