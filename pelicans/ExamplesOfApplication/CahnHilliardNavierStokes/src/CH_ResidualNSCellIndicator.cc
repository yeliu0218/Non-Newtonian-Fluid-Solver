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

#include <CH_ResidualNSCellIndicator.hh>

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

#include <FE_SetOfParameters.hh>
#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>
#include <iostream>
#include <string>
#include <math.h>

using std::string ;
using std::cout ; using std::endl ;

CH_ResidualNSCellIndicator const*
CH_ResidualNSCellIndicator:: PROTOTYPE = new CH_ResidualNSCellIndicator() ;

//---------------------------------------------------------------------------
CH_ResidualNSCellIndicator:: CH_ResidualNSCellIndicator( void )
//---------------------------------------------------------------------------
   : FE_AdaptationIndicator( "CH_ResidualNSCellIndicator" )
   , CELL_ERRORS( 0 )
{
}

//---------------------------------------------------------------------------
CH_ResidualNSCellIndicator*
CH_ResidualNSCellIndicator:: create_replica(
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_ResidualNSCellIndicator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp, a_verbose_level ) ) ;

   CH_ResidualNSCellIndicator* result =
                    new CH_ResidualNSCellIndicator( a_owner, dom, prms,
                                                    exp, a_verbose_level ) ;

   PEL_CHECK( create_replica_POST( result, a_owner,
                                   dom, prms, exp, a_verbose_level ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
CH_ResidualNSCellIndicator:: CH_ResidualNSCellIndicator(
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer const* exp,
                                           size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : FE_AdaptationIndicator( a_owner, a_verbose_level )
   , NB_REFS( 10000 )
   , PT( GE_Point::create( this, dom->nb_space_dimensions() ) )
   , ICALL( 0 )
   , UU( dom->set_of_discrete_fields()->item(
                     exp->string_data( "velocity" ) ) )
   , PP( dom->set_of_discrete_fields()->item(
                     exp->string_data( "pressure" ) ) )
   , L_UU( exp->int_data( "level_of_velocity" ) )
   , L_PP( exp->int_data( "level_of_pressure" ) )
   , DENS( prms->item( exp->string_data( "param_density" ) ) )
   , MU( prms->item( exp->string_data( "param_viscous" ) ) )
   , RHS( prms->item( exp->string_data( "param_source" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP( GE_QRprovider::object(
                           exp->string_data( "quadrature_rule_provider" ) ) )
   , NB_DIMS( dom->nb_space_dimensions() )
   , CELL_ERRORS( 0 )
   , MAX_ERR( exp->double_data( "maximum_error" ) )
   , MIN_ERR( -1.E+10 )
   , BUILD_OK( false )
{
   sFE->include_color( GE_Color::halo_color() ) ;
   sFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
   sFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   DENS->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
   MU->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;

   bFE->include_color( GE_Color::halo_color() ) ;
   bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   DENS->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
   MU->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;

   cFE->include_color( GE_Color::halo_color() ) ;
   cFE->require_field_calculation( PP, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::d2N ) ;
   DENS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   MU->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   RHS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   CELL_ERRORS.re_initialize( cFE->nb_meshes() ) ;
   BUILD_OK = false ;

   if( exp->has_entry( "minimum_error" ) )
         MIN_ERR = exp->double_data( "minimum_error" ) ;

   if( exp->has_entry( "max_nb_steps" ) )
         NB_REFS = exp->int_data( "max_nb_steps" ) ;
}

//---------------------------------------------------------------------------
CH_ResidualNSCellIndicator:: ~CH_ResidualNSCellIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
CH_ResidualNSCellIndicator:: reset_self( void )
//---------------------------------------------------------------------------
{
   ICALL = 0 ;
}

//---------------------------------------------------------------------------
void
CH_ResidualNSCellIndicator:: build_self( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_ResidualNSCellIndicator:: build_self" ) ;

   ++ICALL ;
   CELL_ERRORS.set( 0.0 ) ;

   size_t nbc = UU->nb_components() ;
   PEL_ASSERT( nbc == cFE->nb_space_dimensions() ) ; //cf ugradu et div

   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Vector const* normal = sFE->normal() ;
      size_t id_0 = fe_K->mesh_id() ;
      if( id_0 >= CELL_ERRORS.size() )
      {
         CELL_ERRORS.resize( id_0 + 1 ) ;
      }
      size_t id_1 = fe_L->mesh_id() ;
      if( id_1 >= CELL_ERRORS.size() )
      {
         CELL_ERRORS.resize( id_1 + 1) ;
      }

      sFE->start_IP_iterator( QRP ) ;
      double jump = 0.0 ;
      for( ; sFE->valid_IP() ; sFE->go_next_IP() )
      {
         for( size_t ic=0 ; ic < UU->nb_components() ; ++ic )
         {
            double gg = 0 ;
            for( size_t d=0 ; d<NB_DIMS ; ++d )
            {
               gg += ( MU->cell_value_at_IP( t_it, fe_K )
                     * fe_K->gradient_at_IP( UU, L_UU, d, ic )
                     - MU->cell_value_at_IP( t_it, fe_L )
                     * fe_L->gradient_at_IP( UU, L_UU, d, ic ) )
                     * normal->component( d ) ;
            }
            gg -= ( fe_K->value_at_IP( PP, L_PP, 0 )
                  - fe_L->value_at_IP( PP, L_PP, 0 ) )
                  * normal->component( ic ) ;
            jump += 0.25 * PEL::sqr( gg ) ;
         }
      }

      double hh_0 = fe_K->polyhedron()->inter_vertices_maximum_distance() ;
      CELL_ERRORS( id_0 ) += hh_0 * jump  ;
      double hh_1 = fe_L->polyhedron()->inter_vertices_maximum_distance() ;
      CELL_ERRORS( id_1 ) += hh_1* jump  ;
   }

   //???????????????Boucle sur les bords neumann

         for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
         {
            size_t id = cFE->mesh_id() ;
            PEL_ASSERT( id < CELL_ERRORS.size() ) ;

            CELL_ERRORS( id ) = PEL::sqrt( CELL_ERRORS( id ) ) ;

            cFE->start_IP_iterator( QRP ) ;
            double res = 0.0 ;
            double div = 0.0 ;
            for( ; cFE->valid_IP() ; cFE->go_next_IP() )
            {
               double div_ip = 0.0 ;
               for( size_t ic=0 ; ic<nbc ; ++ic )
               {
                  double lapl = 0.0 ; //????? Assembler div( 2 eta DU )
                  double ugradu = 0.0 ;
                  for( size_t a=0 ; a<cFE->nb_space_dimensions() ; ++a )
                  {
                     lapl += cFE->hessian_at_IP( UU, L_UU, a, a, ic ) ;
                     ugradu += cFE->value_at_IP( UU, L_UU, a )
                     * cFE->gradient_at_IP( UU, L_UU, a, ic )  ;
                  }
                  res += PEL::sqr( RHS->cell_value_at_IP( t_it, cFE, ic )
                                 + MU->cell_value_at_IP( t_it, cFE ) *  lapl
                                 - DENS ->cell_value_at_IP( t_it, cFE ) * ugradu
                                 - cFE->gradient_at_IP( PP, L_PP, ic, 0 ) ) ;
                  div_ip +=  cFE->gradient_at_IP( UU, L_UU, ic, ic )  ;
               }
               div += PEL::sqr( div_ip ) ;
            }
            double hh = cFE->polyhedron()->inter_vertices_maximum_distance() ;
            CELL_ERRORS( id ) += hh* ( PEL::sqrt( res ) + PEL::sqrt( div ) ) ;
            //hh en 3d aussi  ???
         }
         CELL_ERRORS_MAX = CELL_ERRORS( 0 ) ;
         for( size_t i = 1 ; i < CELL_ERRORS.size() ; ++i)
         {
            if( CELL_ERRORS( i ) > CELL_ERRORS_MAX )
               CELL_ERRORS_MAX = CELL_ERRORS( i ) ;
         }

         if( !t_it->is_started() ||
             !(t_it->iteration_number() > t_it->initial_iteration_number()) )
         {
            PEL::out() << " (CH_ResidualCellIndicator not applied) " ;
            CELL_ERRORS.set( 0.0 ) ;
            CELL_ERRORS_MAX = 0.0 ;
         }
         BUILD_OK = true ;
}

//---------------------------------------------------------------------------
bool
CH_ResidualNSCellIndicator:: to_be_refined( double bf_indicator,
                                       GE_Mpolyhedron const* poly,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_ResidualNSCellIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( bf_indicator, poly, elm, local_node ) ) ;

   bool result = false ;
   if( ICALL <= NB_REFS )
   {
      if( bf_indicator/CELL_ERRORS_MAX > MAX_ERR ) result = true ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
CH_ResidualNSCellIndicator:: to_be_unrefined( double bf_indicator,
                                         GE_Mpolyhedron const* poly,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_ResidualNSCellIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( bf_indicator, poly, elm, local_node ) );

   bool result = false ;
   if( MIN_ERR > 0.0 )
   {
      if( ICALL <= NB_REFS )
      {
         if( bf_indicator/CELL_ERRORS_MAX < MIN_ERR ) result = true ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
double
CH_ResidualNSCellIndicator:: cell_indicator( size_t cell_id ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_ResidualNSCellIndicator:: cell_indicator" ) ;
   PEL_ASSERT( BUILD_OK ) ;
   double result = CELL_ERRORS( cell_id );
   return( result ) ;
}
