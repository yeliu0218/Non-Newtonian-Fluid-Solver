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

#include <FE_AdvVelocityParameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <sstream>

FE_AdvVelocityParameter const* 
FE_AdvVelocityParameter:: PROTOTYPE_1 = new FE_AdvVelocityParameter(
                              order1, "FE_AdvVelocityParameter#order1" ) ;
FE_AdvVelocityParameter const* 
FE_AdvVelocityParameter:: PROTOTYPE_2 = new FE_AdvVelocityParameter(
                              order2, "FE_AdvVelocityParameter#order2" ) ;
FE_AdvVelocityParameter const* 
FE_AdvVelocityParameter:: PROTOTYPE_U_DEF = new FE_AdvVelocityParameter(
                  user_defined, "FE_AdvVelocityParameter#user_defined" ) ;

struct FE_AdvVelocityParameter_ERROR
{
   static void n0( std::string const& velo_name,
                   size_t level0, size_t level1 ) ;
   static void n1( std::string const& velo_name ) ;
} ;

//-------------------------------------------------------------------------
FE_AdvVelocityParameter:: FE_AdvVelocityParameter(
             FE_AdvVelocityLaw const id, std::string const& concrete_name )
//-------------------------------------------------------------------------
   : FE_Parameter( concrete_name )
   , ADV_VELO_LAW( id )
   , VELOCITY( 0 )
   , LEVELS( 0 )
   , COEFFS( 0 )
{
}

//-------------------------------------------------------------------------
FE_AdvVelocityParameter*
FE_AdvVelocityParameter:: create_replica(
                               PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   size_t_vector levels(0) ;
   doubleVector coeffs(0) ;
   
   if( ADV_VELO_LAW == order1 )
   {
      coeffs.re_initialize(1) ;
      coeffs(0) = 1. ;
      levels.re_initialize(1) ;
      levels(0) = (size_t) exp->int_data( "velocity_level" ) ;
   }
   else if( ADV_VELO_LAW == order2 )
   {
      coeffs.re_initialize(2) ;
      coeffs(0) = -1. ;
      coeffs(1) =  2. ;
      levels.re_initialize(2) ;
      levels(0) = (size_t) exp->int_data( "previous_velocity_level" ) ;
      levels(1) = (size_t) exp->int_data( "initial_velocity_level" ) ;
   }
   else
   {
      coeffs = exp->doubleVector_data( "coefficients_table" ) ;
      intVector const& l = exp->intVector_data( "velocity_levels_table" ) ;
      if( coeffs.size() != l.size() )
      {
         PEL_Error::object()->raise_plain(
            "Incompatible sizes for the FE_AdvVelocityParameter tables :\n"
            "   - \"coefficients_table\"\n"
            "   - \"velocity_levels_table\"" ) ;
      }
      levels.re_initialize( l.size() ) ;
      for( size_t i=0 ; i<l.size() ; ++i )
      {
         levels(i) = (size_t) l(i) ;
      }
   }
   
   FE_AdvVelocityParameter* result =
      new FE_AdvVelocityParameter( a_owner,
                                   ADV_VELO_LAW,
                                   levels, coeffs,
                                   dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_AdvVelocityParameter:: FE_AdvVelocityParameter(
                               PEL_Object* a_owner,
                               FE_AdvVelocityLaw const id,
                               size_t_vector const& levels,
                               doubleVector const& coeffs,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , ADV_VELO_LAW( id )
   , VELOCITY(
      dom->set_of_discrete_fields()->item( exp->string_data( "velocity_name" ) ) )
   , LEVELS( levels )
   , COEFFS( coeffs )
{   
   size_t max_level = 0 ;
   for( size_t i=0 ; i<levels.size() ; ++i )
   {
      max_level = PEL::max( max_level, levels(i) ) ;
   }
   
   if( VELOCITY->storage_depth() <= max_level )
   {
      FE_AdvVelocityParameter_ERROR::n0( VELOCITY->name(),
                                         VELOCITY->storage_depth(),
                                         max_level+1 ) ;
   }
   if( VELOCITY->nb_components() != dom->nb_space_dimensions() )
   {
      FE_AdvVelocityParameter_ERROR::n1( VELOCITY->name() ) ;
   }
}

//-------------------------------------------------------------------------
FE_AdvVelocityParameter:: ~FE_AdvVelocityParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_AdvVelocityParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( VELOCITY->nb_components() ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                            PDE_LocalFEcell const* fe,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->value_at_pt( VELOCITY, LEVELS(i), ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                               PDE_LocalFEcell const* fe,
                                               size_t a,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a,  ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->gradient_at_pt( VELOCITY, LEVELS(i), a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                            PDE_LocalFEcell const* fe,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->value_at_IP( VELOCITY, LEVELS(i), ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                               PDE_LocalFEcell const* fe,
                                               size_t a,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a,  ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->gradient_at_IP( VELOCITY, LEVELS(i), a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                             PDE_LocalFEbound const* fe,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->value_at_pt( VELOCITY, LEVELS(i), ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                                PDE_LocalFEbound const* fe,
                                                size_t a,
                                                size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->gradient_at_pt( VELOCITY, LEVELS(i), a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                             PDE_LocalFEbound const* fe,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->value_at_IP( VELOCITY, LEVELS(i), ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_AdvVelocityParameter:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                                PDE_LocalFEbound const* fe,
                                                size_t a,
                                                size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = 0. ;
   for( size_t i=0 ; i<LEVELS.size() ; ++i )
   {
      result += COEFFS(i)*fe->gradient_at_IP( VELOCITY, LEVELS(i), a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_value_on_cells(
                                                PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::N ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_gradient_on_cells(
                                                PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_gradient_on_cells" ) ;
   PEL_CHECK( prepare_for_gradient_on_cells_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::dN ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_value_on_bounds(
                                               PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::N ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_gradient_on_bounds(
                                               PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_gradient_on_bounds" ) ;
   PEL_CHECK( prepare_for_gradient_on_bounds_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::dN ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_value_on_sides(
                                                PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::N ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_AdvVelocityParameter:: prepare_for_gradient_on_sides(
                                                PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdvVelocityParameter:: prepare_for_gradient_on_sides" ) ;
   PEL_CHECK( prepare_for_gradient_on_sides_PRE( fe ) ) ;
   
   if( VELOCITY != 0 )
   {
      fe->require_field_calculation( VELOCITY, PDE_LocalFE::dN ) ;
   }
}

//internal--------------------------------------------------------------
void 
FE_AdvVelocityParameter_ERROR:: n0( std::string const& velo_name,
                                    size_t level0, size_t level1 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** FE_AdvVelocityParameter error\n" ;
   mesg << "    level of the discrete field \""
        << velo_name << "\"\n" ;
   mesg << "    is " << level0 << " and at least " << level1 << " is expected." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_AdvVelocityParameter_ERROR:: n1( std::string const& velo_name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** FE_AdvVelocityParameter error\n" ;
   mesg << "    bad number of components of the discrete field \""
        << velo_name << "\"\n" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
