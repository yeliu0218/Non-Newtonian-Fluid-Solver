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

#include <FE_NormalVelocityBC.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <FE_TimeIterator.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;
using std::ostringstream ;
using std::map ; using std::vector ;

FE_NormalVelocityBC const* 
FE_NormalVelocityBC:: PROTOTYPE = new FE_NormalVelocityBC() ;

//-------------------------------------------------------------------------
FE_NormalVelocityBC:: FE_NormalVelocityBC( void )
//-------------------------------------------------------------------------
   : FE_OneBCbuilder( "FE_NormalVelocityBC" )
{
}

//---------------------------------------------------------------------------
FE_NormalVelocityBC*
FE_NormalVelocityBC:: create_replica( PEL_Object* a_owner,
     			              PDE_DiscreteField const* ff,
			              GE_Color const* color ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_NormalVelocityBC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, ff, color ) ) ;

   FE_NormalVelocityBC* result = 
                              new FE_NormalVelocityBC( a_owner, ff, color ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, ff, color ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_NormalVelocityBC:: FE_NormalVelocityBC( PEL_Object* a_owner,
				           PDE_DiscreteField const* ff,
                                           GE_Color const* color )
//---------------------------------------------------------------------------
   : FE_OneBCbuilder( a_owner, ff, color )
{
}

//---------------------------------------------------------------------------
FE_NormalVelocityBC:: ~FE_NormalVelocityBC( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_NormalVelocityBC:: read_boundary_condition( size_t idx,
                                               PDE_DomainAndFields const* dom,
					       FE_SetOfParameters const* prms,
					       PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_NormalVelocityBC:: read_boundary_condition" ) ;
   PEL_CHECK_PRE( read_boundary_condition_PRE( idx, dom, prms, exp ) ) ;

   PEL_ASSERT( idx == VN.size()   &&
               idx == PCOEF.size() ) ;

   VN.push_back( prms->item( exp->string_data( "param_normal_velocity" ) ) ) ;
   PCOEF.push_back( exp->double_data( "penalization_coefficient" ) ) ;
}
   
//---------------------------------------------------------------------------
void
FE_NormalVelocityBC:: transfer_calculation_requirements( 
                                                 PDE_LocalFEbound* fe ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_NormalVelocityBC:: transfer_calculation_requirements" ) ;
   PEL_CHECK_PRE( transfer_calculation_requirements_PRE( fe ) ) ;

   vector< FE_Parameter* >::const_iterator it ;
   for( it = VN.begin() ; it != VN.end() ; ++it )
   {
      (*it)->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   }
   fe->require_field_calculation( field(), PDE_LocalFE::N ) ;
}

//---------------------------------------------------------------------------
void
FE_NormalVelocityBC:: build( PDE_LocalEquation* leq, 
                             PDE_LocalFEbound* fe, 
   	          	     FE_TimeIterator const* t_it,
                             GE_QRprovider const* qrp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_NormalVelocityBC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( leq, fe, t_it, qrp ) ) ;

   size_t i = index_of_color( fe->color() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   static doubleVector cimp( nb_dims ) ;

   for( fe->start_IP_iterator( qrp ) ; fe->valid_IP() ; fe->go_next_IP() )
   {
      double vn = VN[i]->bound_value_at_IP( t_it, fe ) ;

      for( size_t d=0 ; d<nb_dims ; ++d )
            cimp( d ) = PCOEF[i] * vn * fe->outward_normal()->component( d ) ;

      FE::add_row_col_S( leq, fe, PCOEF[i] ) ;
      FE::add_row( leq, fe, cimp ) ;
   }
}

