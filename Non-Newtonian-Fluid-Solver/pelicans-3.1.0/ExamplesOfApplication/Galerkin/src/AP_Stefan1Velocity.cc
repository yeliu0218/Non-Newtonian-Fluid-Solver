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

#include <AP_Stefan1Velocity.hh>

#include <PEL_ModuleExplorer.hh>

#include <GE_Vector.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
using namespace std ;

AP_Stefan1Velocity const* 
AP_Stefan1Velocity:: PROTOTYPE = new AP_Stefan1Velocity() ;

//-------------------------------------------------------------------------
AP_Stefan1Velocity:: AP_Stefan1Velocity( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "AP_Stefan1Velocity" )
{
}

//-------------------------------------------------------------------------
AP_Stefan1Velocity*
AP_Stefan1Velocity:: create_replica( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "AP_Stefan1Velocity:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   AP_Stefan1Velocity* result = new AP_Stefan1Velocity( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
AP_Stefan1Velocity:: AP_Stefan1Velocity( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , TT( dom->set_of_discrete_fields()->item( 
                                         exp->string_data( "temperature" ) ) )
   , L_TT( exp->int_data( "level_of_temperature" ) )
   , EXP( exp->create_clone( this ) )
   , KAPPA( 0 )
   , S_DENS( exp->double_data( "solid_density" ) )
   , LHEAT( exp->double_data( "h_in_minus_h_out" ) )
   , FLUX_EXT( 0 )
   , NB_DIMS( dom->nb_space_dimensions() )
   , TMEL( 0 )
   , CHECK_BOUNDING( false )
   , TT_MIN( 0.0 )
   , TT_MAX( 0.0 )
{
   if( exp->has_module( "check_bounding" ) )
   {
      PEL_ModuleExplorer const* se = 
                             exp->create_subexplorer( 0, "check_bounding" ) ;
      TT_MIN = se->double_data( "value_min" ) ;
      TT_MAX = se->double_data( "value_max" ) ;
      CHECK_BOUNDING = true ;
      se->destroy() ;
   }
}

//-------------------------------------------------------------------------
AP_Stefan1Velocity:: ~AP_Stefan1Velocity( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
AP_Stefan1Velocity:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   size_t result = NB_DIMS ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
AP_Stefan1Velocity:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "AP_Stefan1Velocity:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

   KAPPA = prms->item( EXP->string_data( "conductivity" ) ) ;
   FLUX_EXT = prms->item( EXP->string_data( "flux_from_out" ) ) ;
}

//-------------------------------------------------------------------------
double
AP_Stefan1Velocity:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "AP_Stefan1Velocity:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;

   if( CHECK_BOUNDING )
   {
      double VAL = fe->value_at_IP( TT, L_TT ) ;
      if( (VAL <= TT_MIN) || (VAL >= TT_MAX) )
      {
         return( result ) ;
	 // -------------
      }
   }

   GE_Vector const* nn = fe->outward_normal() ;
   double flux_in = 0.0 ;
   for( size_t d=0 ; d<NB_DIMS ; ++d )
   {
      flux_in -= fe->gradient_at_IP( TT, L_TT, d ) * nn->component( d ) ;
   }
   flux_in *= KAPPA->bound_value_at_IP( t_it, fe ) ;

   result = ( FLUX_EXT->bound_value_at_IP( t_it, fe ) + flux_in ) 
            / S_DENS / LHEAT ;

   return( result ) ;
   // -------------
}

//-------------------------------------------------------------------------
void
AP_Stefan1Velocity:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   if( CHECK_BOUNDING ) fe->require_field_calculation( TT, PDE_LocalFE::N )  ;
   fe->require_field_calculation( TT, PDE_LocalFE::dN ) ;
   KAPPA->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   FLUX_EXT->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
}

