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

#include <FE_LocalBCsBuilder.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_OneBCbuilder.hh>
#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;
using std::ostringstream ;
using std::map ;

struct FE_LocalBCsBuilder_ERROR
{
   static void n0( std::string const& field_name,
                   std::string const& bc_type ) ;
} ;

//---------------------------------------------------------------------------
FE_LocalBCsBuilder*
FE_LocalBCsBuilder:: create( PEL_Object* a_owner,
			     PDE_DomainAndFields const* dom,
                             std::string const& field_name,
			     stringVector const& bc_types,
			     FE_SetOfParameters const* prms )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_LocalBCsBuilder:: create" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;
   PEL_CHECK_PRE( bc_types.size() != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;

   FE_LocalBCsBuilder* result = new FE_LocalBCsBuilder( a_owner, 
                                                        dom,
                                                        field_name,
                                                        bc_types, 
                                                        prms ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->current_BC_type_is_ok() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_LocalBCsBuilder:: FE_LocalBCsBuilder( PEL_Object* a_owner,
					 PDE_DomainAndFields const* dom,
                                         std::string const& field_name,
					 stringVector const& bc_types,
				         FE_SetOfParameters const* prms )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FF( dom->set_of_discrete_fields()->item( field_name ) )
{
   PEL_LABEL( "FE_LocalBCsBuilder:: FE_LocalBCsBuilder" ) ;

   boolVector ok_types( bc_types.size() ) ;

   map< string, FE_OneBCbuilder* >::const_iterator it ;
   PDE_SetOfBCs const* bcs = dom->set_of_boundary_conditions() ;

   PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Color const* color = fe->color() ;
      if( bcs->has_BC( fe->color(), FF ) )
      {
	 PEL_ModuleExplorer const* ee = bcs->BC_explorer( color, FF ) ;
         string const& type = ee->string_data( "type" ) ;
         if( bc_types.has( type ) ) 
	 {
            it = BUILDERS.find( type ) ;
            if( it == BUILDERS.end() )
	    {
               ok_types( bc_types.index_of( type ) ) = true ;
               BUILDERS[ type ] = FE_OneBCbuilder::make( this, type, 
                                                         FF, color, 
                                                         dom, prms, ee ) ;
	    }
	    else
	    {
               it->second->extend( color, dom, prms, ee ) ;
	    }
	 }
      }
   }
   fe->destroy() ;

   for( size_t i=0 ; i<bc_types.size() ; ++i )
   {
      if( !ok_types(i) ) 
         FE_LocalBCsBuilder_ERROR::n0( field_name, bc_types(i) ) ;
   }
   IT_BC_to_build = BUILDERS.end() ;
}

//---------------------------------------------------------------------------
FE_LocalBCsBuilder:: ~FE_LocalBCsBuilder( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
PDE_DiscreteField const*
FE_LocalBCsBuilder:: field( void ) const
//---------------------------------------------------------------------------
{
   return( FF ) ;
}

//---------------------------------------------------------------------------
void
FE_LocalBCsBuilder:: transfer_calculation_requirements( 
                                                 PDE_LocalFEbound* fe ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_LocalBCsBuilder:: transfer_calculation_requirements" ) ;
   PEL_CHECK_PRE( !current_BC_type_is_ok() ) ;

   map< string, FE_OneBCbuilder* >::const_iterator it ;
   for( it = BUILDERS.begin() ; it != BUILDERS.end() ; ++it )
   {
      it->second->transfer_calculation_requirements( fe ) ;
   }

   PEL_CHECK_POST( !current_BC_type_is_ok() ) ;
}

//---------------------------------------------------------------------------
void
FE_LocalBCsBuilder:: set_current_BC_type( std::string const& bc_type )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_LocalBCsBuilder:: set_current_BC_type" ) ;
   PEL_CHECK_PRE( !current_BC_type_is_ok() ) ;

   IT_BC_to_build = BUILDERS.find( bc_type ) ;
}

//---------------------------------------------------------------------------
bool
FE_LocalBCsBuilder:: current_BC_type_is_ok( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_LocalBCsBuilder:: current_BC_type_is_ok" ) ;

   bool result = ( IT_BC_to_build != BUILDERS.end() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE_LocalBCsBuilder:: build_current_BC( PDE_LocalEquation* leq,
                                       PDE_LocalFEbound* fe,
	            	               FE_TimeIterator const* t_it,
                                       GE_QRprovider const* qrp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_LocalBCsBuilder:: build_current_BC" ) ;
   PEL_CHECK_PRE( current_BC_type_is_ok() ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == field() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) == field() ) ;
   PEL_CHECK_PRE( !fe->valid_IP() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == field()->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == field()->nb_components() ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   PEL_CHECK_PRE( qrp != 0 ) ;

   IT_BC_to_build->second->build( leq, fe, t_it, qrp ) ;

   IT_BC_to_build = BUILDERS.end() ;

   PEL_CHECK_POST( !current_BC_type_is_ok() ) ;
}

//internal-------------------------------------------------------------------
void
FE_LocalBCsBuilder_ERROR:: n0( std::string const& field_name,
                               std::string const& bc_type )
//internal-------------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_LocalBCsBuilder:" << endl << endl ;
   msg << "no boundary_condition of type \"" << bc_type << "\"" << endl ;
   msg << "has been found for field \"" << field_name << "\"" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
