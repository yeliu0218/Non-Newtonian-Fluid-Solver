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

#include <CH_CapillaryForce.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::ostringstream ;
using std::string ;

CH_CapillaryForce const* 
CH_CapillaryForce:: PROTOTYPE = new CH_CapillaryForce() ;

//-------------------------------------------------------------------------
CH_CapillaryForce:: CH_CapillaryForce( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "CH_CapillaryForce" )
{
}

//-------------------------------------------------------------------------
CH_CapillaryForce*
CH_CapillaryForce:: create_replica( PEL_Object* a_owner,
				       PDE_DomainAndFields const* dom,
				       PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CapillaryForce:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   CH_CapillaryForce* result =
                             new CH_CapillaryForce( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
CH_CapillaryForce:: CH_CapillaryForce( PEL_Object* a_owner,
					     PDE_DomainAndFields const* dom,
					     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , NBCS( dom->nb_space_dimensions() )
   , S1( exp->double_data( "coef_sigma_1" ) )
   , S2( exp->double_data( "coef_sigma_2" ) )
   , S3( exp->double_data( "coef_sigma_3" ) )
{
   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;
   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
   se->start_module_iterator() ;
   for( ; se->is_valid_module() ; se->go_next_module() )
   {
      PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
      PDE_DiscreteField const* ff = 
                               dfs->item( sse->string_data( "phase_field" ) ) ;
      CCs.push_back( ff ) ;
      ff = dfs->item( sse->string_data( "generalized_potential" ) ) ;
      MMs.push_back( ff ) ;
      L_CCs.push_back( sse->int_data( "level_of_fields" ) ) ;
      sse->destroy() ;
   }
   se->destroy() ;

}

//-------------------------------------------------------------------------
CH_CapillaryForce:: ~CH_CapillaryForce( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
CH_CapillaryForce:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
CH_CapillaryForce:: cell_value_at_IP( FE_TimeIterator const* t_it,
					 PDE_LocalFEcell const* fe,
					 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CapillaryForce:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ; 

   double m =  fe->value_at_IP( MMs[0], L_CCs[0] ) ;
   double n =  fe->value_at_IP( MMs[1], L_CCs[1] ) ;
   double p = - ( m*S3/S1 + n*S3/S2 ) ;

   double grad_c = fe->gradient_at_IP( CCs[0], L_CCs[0], ic ) ;
   double grad_d = fe->gradient_at_IP( CCs[1], L_CCs[1], ic ) ;
   double grad_f = - ( grad_c + grad_d ) ;
   
   result = m*grad_c + n*grad_d + p*grad_f ;

   return( result ) ;
}
//-------------------------------------------------------------------------
double
CH_CapillaryForce:: cell_value_at_pt( FE_TimeIterator const* t_it,
					 PDE_LocalFEcell const* fe,
					 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_CapillaryForce:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ; 

   double m =  fe->value_at_pt( MMs[0], L_CCs[0] ) ;
   double n =  fe->value_at_pt( MMs[1], L_CCs[1] ) ;
   double p = - ( m*S3/S1 + n*S3/S2 ) ;

   double grad_c = fe->gradient_at_pt( CCs[0], L_CCs[0], ic ) ;
   double grad_d = fe->gradient_at_pt( CCs[1], L_CCs[1], ic ) ;
   double grad_f = - ( grad_c + grad_d ) ;
   
   result = m*grad_c + n*grad_d + p*grad_f ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
CH_CapillaryForce:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   for( size_t i=0 ; i<CCs.size() ; ++i ) 
   {
      fe->require_field_calculation( CCs[i], PDE_LocalFE::dN )  ;   
      fe->require_field_calculation( MMs[i], PDE_LocalFE::N )  ;   
   }
}
