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

#include <FE_GravityParameter.hh>

#include <PEL_ModuleExplorer.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

FE_GravityParameter const* 
FE_GravityParameter:: PROTOTYPE = new FE_GravityParameter() ;

//-------------------------------------------------------------------------
FE_GravityParameter:: FE_GravityParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_GravityParameter" )
   , GRAV( 0 )
{
}

//-------------------------------------------------------------------------
FE_GravityParameter*
FE_GravityParameter:: create_replica( PEL_Object* a_owner,
        		              PDE_DomainAndFields const* dom,
                                      PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_GravityParameter* result = new FE_GravityParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_GravityParameter:: FE_GravityParameter( PEL_Object* a_owner,
        		                   PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , GRAV( exp->doubleVector_data( "gravity" ) )
   , RHOS( exp->double_data( "rho_shift" ) )
   , RHO( 0 )
   , EXP( exp->create_clone( this ) ) 
{
}

//-------------------------------------------------------------------------
FE_GravityParameter:: ~FE_GravityParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_GravityParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   size_t result = GRAV.size() ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_GravityParameter:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

   RHO = prms->item( EXP->string_data( "rho" ) ) ;
}

//-------------------------------------------------------------------------
double
FE_GravityParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = ( RHO->cell_value_at_pt( t_it, fe ) - RHOS ) * GRAV(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_GravityParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = ( RHO->cell_value_at_IP( t_it, fe ) - RHOS ) * GRAV(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_GravityParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   RHO->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//-------------------------------------------------------------------------
void
FE_GravityParameter:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GravityParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;
   RHO->transfer_side_calculation_requirements(  fe, FE_Parameter::Val ) ;
}
