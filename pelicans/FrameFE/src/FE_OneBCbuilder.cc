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

#include <FE_OneBCbuilder.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfBCs.hh>
#include <FE_TimeIterator.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;
using std::ostringstream ;
using std::map ;

//---------------------------------------------------------------------------
FE_OneBCbuilder*
FE_OneBCbuilder:: make( PEL_Object* a_owner,
                        std::string const& type,
			PDE_DiscreteField const* ff,
			GE_Color const* color,
                        PDE_DomainAndFields const* dom,
			FE_SetOfParameters const* prms,
                        PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneBCbuilder:: make" ) ;
   PEL_CHECK_PRE( !type.empty() ) ;

   FE_OneBCbuilder const* proto =
      static_cast<FE_OneBCbuilder const*>( plugins_map()->item( type ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   FE_OneBCbuilder* result = proto->create_replica( a_owner, ff, color ) ;
   result->read_boundary_condition( 0, dom, prms, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == ff ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_OneBCbuilder:: FE_OneBCbuilder( std::string const& a_type )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "FE_OneBCbuilder:: FE_OneBCbuilder" ) ;

   plugins_map()->register_item( a_type, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
FE_OneBCbuilder:: FE_OneBCbuilder( PEL_Object* a_owner,
				   PDE_DiscreteField const* ff,
                                   GE_Color const* color )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , FF( ff )
   , IDX( 0 )
{
   COLOR_2_BCindex[ color ] = 0 ;

   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
FE_OneBCbuilder:: ~FE_OneBCbuilder( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_OneBCbuilder:: extend( GE_Color const* color,
			  PDE_DomainAndFields const* dom,
                          FE_SetOfParameters const* prms,
                          PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneBCbuilder:: extend" ) ;

   map<GE_Color const*,size_t>::const_iterator it ;
   it = COLOR_2_BCindex.find( color ) ;
   if( it == COLOR_2_BCindex.end() )
   {
      ++IDX ;
      COLOR_2_BCindex[ color ] = IDX ;
      read_boundary_condition( IDX, dom, prms, exp ) ;
   }
}

//---------------------------------------------------------------------------
PDE_DiscreteField const*
FE_OneBCbuilder:: field( void ) const
//---------------------------------------------------------------------------
{
   PDE_DiscreteField const* result = FF ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneBCbuilder:: is_a_prototype( void ) const
//---------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//---------------------------------------------------------------------------
size_t
FE_OneBCbuilder:: index_of_color( GE_Color const* color ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneBCbuilder:: index_of_color" ) ;

   map<GE_Color const*,size_t>::const_iterator it = 
                                        COLOR_2_BCindex.find( color ) ;
   PEL_ASSERT( it != COLOR_2_BCindex.end() ) ;
   size_t result = it->second  ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool 
FE_OneBCbuilder:: read_boundary_condition_PRE( 
                                   size_t idx,
  		                   PDE_DomainAndFields const* dom,
				   FE_SetOfParameters const* prms,
				   PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool 
FE_OneBCbuilder:: transfer_calculation_requirements_PRE( 
                                                PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool 
FE_OneBCbuilder:: build_PRE( PDE_LocalEquation* leq,
 		             PDE_LocalFEbound* fe,
		             FE_TimeIterator const* t_it,
                             GE_QRprovider const* qrp ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->field( PDE_LocalFE::row ) == field() ) ;
   PEL_ASSERT( fe->field( PDE_LocalFE::col ) == field() ) ;
   PEL_ASSERT( !fe->valid_IP() ) ;
   PEL_ASSERT( leq != 0 ) ;
   PEL_ASSERT( leq->nb_rows() == 
               fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_ASSERT( leq->nb_columns() == 
               fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_ASSERT( leq->nb_row_sub_indices() == field()->nb_components() ) ;
   PEL_ASSERT( leq->nb_column_sub_indices() == field()->nb_components() ) ;
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( qrp != 0 ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_OneBCbuilder:: create_replica_PRE( PEL_Object* a_owner,
	 	   		      PDE_DiscreteField const* ff,
  				      GE_Color const* color ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( color != 0 ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_OneBCbuilder:: create_replica_POST( FE_OneBCbuilder const* result,
                                       PEL_Object* a_owner,
				       PDE_DiscreteField const* ff,
  				       GE_Color const* color ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   PEL_ASSERT( result->field() == ff ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
FE_OneBCbuilder:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "FE_OneBCbuilder descendant" ) ;
   return( result ) ;
}
