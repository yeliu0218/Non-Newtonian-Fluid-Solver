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

#include <FE_SurfaceForce.hh>

#include <iostream>
#include <cmath>

#include <PEL.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <sstream>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

struct FE_SurfaceForce_ERROR {
   static void n0( std::string const& a_name ) ;
   static void n1( std::string const& a_name ) ;
} ;

std::map< std::string, FE_SurfaceForce* > FE_SurfaceForce::OBJS ;

FE_SurfaceForce const* FE_SurfaceForce::PROTOTYPE = new FE_SurfaceForce() ;

//---------------------------------------------------------------------------
FE_SurfaceForce:: FE_SurfaceForce( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_SurfaceForce" )
   , FORCE( 0 )
   , FORCE_P( 0 )
   , FORCE_VISC( 0 )
{
}

//---------------------------------------------------------------------------
FE_SurfaceForce*
FE_SurfaceForce:: create_replica( PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_SurfaceForce* result = new FE_SurfaceForce( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_SurfaceForce:: FE_SurfaceForce( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , NAME( exp->has_entry( "name") ? exp->string_data( "name" ) : "unknown" )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data("velocity") ) )
   , L_UU( exp->int_data( "level_of_velocity" )  )
   , PP( dom->set_of_discrete_fields()->item( exp->string_data("pressure") ) )
   , L_PP( exp->int_data( "level_of_pressure" ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , QRP( GE_QRprovider::object( "GE_QRprovider_5" ) )
   , MU( prms->item( exp->string_data( "param_visc" ) ) )
   , SURFCOL( GE_Color::object( exp->string_data( "surface_color" ) ) )
   , FORCE( dom->nb_space_dimensions() )
   , FORCE_P( dom->nb_space_dimensions() )
   , FORCE_VISC( dom->nb_space_dimensions() )
{
   PEL_LABEL( "FE_SurfaceForce:: FE_SurfaceForce" ) ;

   if( OBJS.count( NAME ) != 0 ) FE_SurfaceForce_ERROR::n0( NAME ) ;
   OBJS[ NAME ] = this ;
   
   check_param_nb_components( MU, "param_visc", 1 ) ;
   MU->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;

   bFE->require_field_calculation( UU, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( PP, PDE_LocalFE::N ) ;
}

//---------------------------------------------------------------------------
FE_SurfaceForce:: ~FE_SurfaceForce( void )
//---------------------------------------------------------------------------
{
}

//--------------------------------------------------------------------------
FE_SurfaceForce*
FE_SurfaceForce:: object( std::string const& a_name )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: object" ) ;
   std::map<std::string,FE_SurfaceForce*>::const_iterator it = 
                                                         OBJS.find( a_name ) ;
   if( it == OBJS.end() ) FE_SurfaceForce_ERROR::n1( a_name ) ;
   
   FE_SurfaceForce* result = (*it).second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const&
FE_SurfaceForce:: name( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: name" ) ;

   return( NAME ) ;
}

//---------------------------------------------------------------------------
void
FE_SurfaceForce:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
FE_SurfaceForce:: save_other_than_time_and_fields( 
                                            FE_TimeIterator const* t_it, 
                                            PDE_ResultSaver* rs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   start_total_timer( "FE_SurfaceForce:: save_other_than_time_and_fields" ) ;

   compute_force( t_it ) ;

   rs->save_variable( FORCE( 0 ), "XFX" ) ;
   rs->save_variable( FORCE( 1 ), "XFY" ) ;
   if( bFE->nb_space_dimensions() == 3 )
   {
      rs->save_variable( FORCE( 2 ), "XFZ" ) ;
   }

   if( verbose_level() >= 2 )
   {
      cout << indent() << "   f_p_x = " << FORCE_P( 0 ) << endl ;
      cout << indent() << "   f_p_y = " << FORCE_P( 1 ) << endl ;
      if( bFE->nb_space_dimensions() == 3 )
      {
         cout << indent() << "   f_p_z = " << FORCE_P( 2 ) << endl ;
      }
      cout << indent() << "   f_visc_x = " << FORCE_VISC( 0 ) << endl ;
      cout << indent() << "   f_visc_y = " << FORCE_VISC( 1 ) << endl ;
      if( bFE->nb_space_dimensions() == 3 )
      {
         cout << indent() << "   f_visc_z = " << FORCE_VISC( 2 ) << endl ;
      }
      cout << indent() << "   fx = " << FORCE( 0 ) << endl ;
      cout << indent() << "   fy = " << FORCE( 1 ) << endl ;
      if( bFE->nb_space_dimensions() == 3 )
      {
         cout << indent() << "   fz = " << FORCE( 2 ) << endl ;
      }
   }

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
double
FE_SurfaceForce:: force( size_t d ) const
//---------------------------------------------------------------------------
{
   return( FORCE( d ) ) ;
}

//---------------------------------------------------------------------------
void
FE_SurfaceForce:: compute_force( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SurfaceForce:: compute_force" ) ;

   size_t nb_dims = bFE->nb_space_dimensions() ;
   FORCE.re_initialize( nb_dims ) ;
   FORCE_P.re_initialize( nb_dims ) ;
   FORCE_VISC.re_initialize( nb_dims ) ;

   bFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      if( SURFCOL->is_overlapping( bFE->color() ) )
      {
         GE_Vector const* nor = bFE->outward_normal() ;
         bFE->start_IP_iterator( QRP ) ;
         for( ; bFE->valid_IP() ; bFE->go_next_IP() )
         { 
            double weight = bFE->weight_of_IP() ;
            if( FE::geometry() == FE::axisymmetrical )
            {
               weight *= bFE->coordinates_of_IP()->coordinate( 0 ) ;
               weight *= 2.0 * PEL::pi() ;
            }
            double mu = MU->bound_value_at_IP( t_it, bFE ) ;
            for( size_t d=0 ; d < nb_dims ; d++ )
            {
               FORCE( d ) += weight * bFE->value_at_IP( PP, L_PP ) 
		                              * nor->component( d ) ;
               FORCE_P( d ) += weight * bFE->value_at_IP( PP, L_PP ) 
		                                * nor->component( d ) ;
               for( size_t i=0 ; i < nb_dims ; i++ )
               {                    
                  FORCE( d ) += weight * -nor->component( i ) * mu
                              *( bFE->gradient_at_IP( UU, L_UU, d, i ) 
                               + bFE->gradient_at_IP( UU, L_UU, i, d ) ) ;
                  FORCE_VISC( d ) += weight * -nor->component( i ) * mu
			                     *( bFE->gradient_at_IP( UU, L_UU, d, i ) 
			                      + bFE->gradient_at_IP( UU, L_UU, i, d ) ) ;
               }
            }
         }
      }
   }
}

//internal-------------------------------------------------------------------
void
FE_SurfaceForce_ERROR:: n0( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_SurfaceForce :" << endl ;
   mesg << "   attempt to create two instances with the" << endl ;
   mesg << "   same name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_SurfaceForce_ERROR:: n1( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_SurfaceForce :" << endl ;
   mesg << "   there is no recorded instance of name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

