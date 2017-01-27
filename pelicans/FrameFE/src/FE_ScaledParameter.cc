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

#include <FE_ScaledParameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <FE_SetOfParameters.hh>

FE_ScaledParameter const* 
FE_ScaledParameter:: PROTOTYPE = new FE_ScaledParameter() ;

//-------------------------------------------------------------------------
FE_ScaledParameter:: FE_ScaledParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_ScaledParameter" )
   , PRM_NAME()
   , PRM( 0 )
   , COEF( PEL::bad_double() )
{
}

//-------------------------------------------------------------------------
FE_ScaledParameter*
FE_ScaledParameter:: create_replica( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_ScaledParameter* result = new FE_ScaledParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_ScaledParameter:: FE_ScaledParameter( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , PRM_NAME( exp->string_data( "parameter" ) )
   , PRM( 0 )
   , COEF( exp->double_data( "coefficient" ) )
{
}

//-------------------------------------------------------------------------
FE_ScaledParameter:: ~FE_ScaledParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

   PRM = prms->item( PRM_NAME ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: reset( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: reset" ) ;
   PEL_CHECK_PRE( reset_PRE( exp ) ) ;
   
   PRM->reset( exp ) ;
}

//-------------------------------------------------------------------------
size_t
FE_ScaledParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( PRM->nb_components() ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_value" ) ;
   PEL_CHECK_PRE( cell_value_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->cell_value( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_gradient( FE_TimeIterator const* t_it,
                                    PDE_LocalFEcell const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_gradient" ) ;
   PEL_CHECK_PRE( cell_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->cell_gradient( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->cell_value_at_pt( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->cell_gradient_at_pt( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->cell_value_at_IP( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->cell_gradient_at_IP( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_value( FE_TimeIterator const* t_it,
                                  PDE_LocalFEbound const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_value" ) ;
   PEL_CHECK_PRE( bound_value_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->bound_value( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_gradient( FE_TimeIterator const* t_it,
                                     PDE_LocalFEbound const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_gradient" ) ;
   PEL_CHECK_PRE( bound_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->bound_gradient( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->bound_value_at_pt( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEbound const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->bound_gradient_at_pt( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->bound_value_at_IP( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->bound_gradient_at_IP( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_value( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_value" ) ;
   PEL_CHECK_PRE( side_value_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->side_value( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_gradient( FE_TimeIterator const* t_it,
                                    PDE_CursorFEside const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_gradient" ) ;
   PEL_CHECK_PRE( side_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->side_gradient( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->side_value_at_pt( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_gradient_at_pt" ) ;
   PEL_CHECK_PRE( side_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->side_gradient_at_pt( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = COEF * PRM->side_value_at_IP( t_it, fe, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ScaledParameter:: side_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: side_gradient_at_IP" ) ;
   PEL_CHECK_PRE( side_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = COEF * PRM->side_gradient_at_IP( t_it, fe, a, ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   PRM->transfer_cell_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_gradient_on_cells" ) ;
   PEL_CHECK( prepare_for_gradient_on_cells_PRE( fe ) ) ;
   PRM->transfer_cell_calculation_requirements( fe, FE_Parameter::Grad ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;
   PRM->transfer_bound_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_gradient_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_gradient_on_bounds" ) ;
   PEL_CHECK( prepare_for_gradient_on_bounds_PRE( fe ) ) ;
   PRM->transfer_bound_calculation_requirements( fe, FE_Parameter::Grad ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;
   PRM->transfer_side_calculation_requirements( fe, FE_Parameter::Val ) ;
}

//-------------------------------------------------------------------------
void
FE_ScaledParameter:: prepare_for_gradient_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ScaledParameter:: prepare_for_gradient_on_sides" ) ;
   PEL_CHECK( prepare_for_gradient_on_sides_PRE( fe ) ) ;
   PRM->transfer_side_calculation_requirements( fe, FE_Parameter::Grad ) ;
}
