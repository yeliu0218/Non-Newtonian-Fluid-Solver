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

#include <FE_RobinBC.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>

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

FE_RobinBC const* FE_RobinBC:: PROTOTYPE = new FE_RobinBC() ;

//-------------------------------------------------------------------------
FE_RobinBC:: FE_RobinBC( void )
//-------------------------------------------------------------------------
   : FE_OneBCbuilder( "FE_RobinBC" )
{
}

//---------------------------------------------------------------------------
FE_RobinBC*
FE_RobinBC:: create_replica( PEL_Object* a_owner, 
  	                     PDE_DiscreteField const* ff,
                             GE_Color const* color ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RobinBC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, ff, color ) ) ;

   FE_RobinBC* result = new FE_RobinBC( a_owner, ff, color ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, ff, color ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_RobinBC:: FE_RobinBC( PEL_Object* a_owner,
	                 PDE_DiscreteField const* ff,
                         GE_Color const* color )
//---------------------------------------------------------------------------
   : FE_OneBCbuilder( a_owner, ff, color )
{
}

//---------------------------------------------------------------------------
FE_RobinBC:: ~FE_RobinBC( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_RobinBC:: read_boundary_condition( size_t idx,
                                      PDE_DomainAndFields const* dom,
			              FE_SetOfParameters const* prms,
				      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RobinBC:: read_boundary_condition" ) ;
   PEL_CHECK_PRE( read_boundary_condition_PRE( idx, dom, prms, exp ) ) ;

   PEL_ASSERT( idx == HH.size()   && 
               idx == UU_OUT.size() ) ;

   HH.push_back( prms->item( exp->string_data( "prefactor" ) ) ) ;
   UU_OUT.push_back( prms->item( exp->string_data( "field_out_value" ) ) ) ;
}

//---------------------------------------------------------------------------
void
FE_RobinBC:: transfer_calculation_requirements( PDE_LocalFEbound* fe ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RobinBC:: transfer_calculation_requirements" ) ;
   PEL_CHECK_PRE( transfer_calculation_requirements_PRE( fe ) ) ;

   vector< FE_Parameter* >::const_iterator it ;
   for( it = HH.begin() ; it != HH.end() ; ++it )
   {
      (*it)->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   }
   for( it = UU_OUT.begin() ; it != UU_OUT.end() ; ++it )
   {
      (*it)->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
   }
   fe->require_field_calculation( field(), PDE_LocalFE::N ) ;
}

//---------------------------------------------------------------------------
void
FE_RobinBC:: build( PDE_LocalEquation* leq, 
                    PDE_LocalFEbound* fe, 
	            FE_TimeIterator const* t_it,
                    GE_QRprovider const* qrp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RobinBC:: build" ) ;
   PEL_CHECK_PRE( build_PRE( leq, fe, t_it, qrp ) ) ;

   size_t i = index_of_color( fe->color() ) ;

   for( fe->start_IP_iterator( qrp ) ; fe->valid_IP() ; fe->go_next_IP() )
   {
      double hh = HH[i]->bound_value_at_IP( t_it, fe ) ;
      double uu_out = UU_OUT[i]->bound_value_at_IP( t_it, fe ) ;
      FE::add_row_col_S( leq, fe, hh ) ;
      FE::add_row( leq, fe, hh*uu_out ) ;
   }
}

