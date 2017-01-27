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

#include <FE_DivergenceParameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <sstream>

FE_DivergenceParameter const* 
FE_DivergenceParameter:: PROTOTYPE = new FE_DivergenceParameter() ;

struct FE_DivergenceParameter_ERROR
{
   static void n0( std::string const& f_name,
                   size_t current, size_t expected ) ;
} ;

//-------------------------------------------------------------------------
FE_DivergenceParameter:: FE_DivergenceParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_DivergenceParameter" )
   , FIELD( 0 )
   , FIELD_LEVEL( PEL::bad_index() )
   , COEF( PEL::bad_double() )
{
}

//-------------------------------------------------------------------------
FE_DivergenceParameter*
FE_DivergenceParameter:: create_replica( 
                             PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;
   
   FE_DivergenceParameter* result = new FE_DivergenceParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_DivergenceParameter:: FE_DivergenceParameter( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , FIELD( dom->set_of_discrete_fields()->item( 
                                           exp->string_data( "field_name" ) ) )
   , FIELD_LEVEL( exp->int_data( "field_level" ) )
   , COEF( exp->has_entry( "coefficient" ) ? 
                exp->double_data( "coefficient" ) : 1.0 )
{
   if( FIELD_LEVEL>=FIELD->storage_depth() )
   {
      FE_DivergenceParameter_ERROR::n0( exp->string_data( "name" ),
                                   FIELD->storage_depth(), FIELD_LEVEL ) ;
   }
}

//-------------------------------------------------------------------------
FE_DivergenceParameter:: ~FE_DivergenceParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
double
FE_DivergenceParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                      PDE_LocalFEcell const* fe,
                                      size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;
 
 double result = PEL::bad_double();
   if ( ic==0 )
   {
   result = COEF *(fe->gradient_at_pt( FIELD, FIELD_LEVEL, 0, 0 )+ fe->gradient_at_pt( FIELD, FIELD_LEVEL, 1, 2 ));
    }
    if ( ic==1 )
    {
    result = COEF *(fe->gradient_at_pt( FIELD, FIELD_LEVEL, 1, 1 )+ fe->gradient_at_pt( FIELD, FIELD_LEVEL, 0, 2 ));
     }
     if ( ic==2 )
    {
    result = 0.0 ;
     }

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_DivergenceParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                      PDE_LocalFEcell const* fe,
                                      size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

    double result = PEL::bad_double();
   if ( ic==0 )
   {
   result =COEF *( fe->gradient_at_IP( FIELD, FIELD_LEVEL, 0, 0 )+ fe->gradient_at_IP( FIELD, FIELD_LEVEL, 1, 2 ));
    }
    if ( ic==1 )
    {
    result =COEF *( fe->gradient_at_IP( FIELD, FIELD_LEVEL, 1, 1 )+ fe->gradient_at_IP( FIELD, FIELD_LEVEL, 0, 2 ));
     }
     if ( ic==2 )
    {
    result = 0.0 ;
     }
   return( result ) ;
}



//-------------------------------------------------------------------------
double
FE_DivergenceParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEbound const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = PEL::bad_double();
   if ( ic==0 )
   {
   result =COEF *(fe->gradient_at_pt( FIELD, FIELD_LEVEL, 0, 0 )+ fe->gradient_at_pt( FIELD, FIELD_LEVEL, 1, 2 ));
    }
    if ( ic==1 )
    {
    result =COEF *( fe->gradient_at_pt( FIELD, FIELD_LEVEL, 1, 1 )+ fe->gradient_at_pt( FIELD, FIELD_LEVEL, 0, 2 ));
     }
     if ( ic==2 )
    {
    result = 0.0 ;
     }
   return( result ) ;
}


//-------------------------------------------------------------------------
double
FE_DivergenceParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEbound const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = PEL::bad_double();
   if ( ic==0 )
   {
   result =COEF *(fe->gradient_at_IP( FIELD, FIELD_LEVEL, 0, 0 )+ fe->gradient_at_IP( FIELD, FIELD_LEVEL, 1, 2 ));
    }
    if ( ic==1 )
    {
    result = COEF *(fe->gradient_at_IP( FIELD, FIELD_LEVEL, 1, 1 )+ fe->gradient_at_IP( FIELD, FIELD_LEVEL, 0, 2 ));
     }
     if ( ic==2 )
    {
    result = 0.0 ;
     }
   return( result ) ;
}


//-------------------------------------------------------------------------
size_t
FE_DivergenceParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: nb_components" ) ;
   
   size_t result = FIELD->nb_components() ;
   
   PEL_CHECK_POST( result == discrete_field()->nb_components() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DiscreteField const*
FE_DivergenceParameter:: discrete_field( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: discrete_field" ) ;
   
   PDE_DiscreteField const* result = FIELD ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
FE_DivergenceParameter:: field_level( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: field_level" ) ;
   
   size_t result = FIELD_LEVEL ;
   
   PEL_CHECK_POST( result < discrete_field()->storage_depth() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_DivergenceParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;

  fe->require_field_calculation( FIELD, PDE_LocalFE::dN )  ;
}

//-------------------------------------------------------------------------
void
FE_DivergenceParameter:: prepare_for_value_on_bounds(
                                               PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

fe->require_field_calculation( FIELD, PDE_LocalFE::dN )  ;
}


//-------------------------------------------------------------------------
void
FE_DivergenceParameter:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DivergenceParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;

fe->require_field_calculation( FIELD, PDE_LocalFE::dN )  ;
}


//internal--------------------------------------------------------------
void 
FE_DivergenceParameter_ERROR:: n0( std::string const& f_name,
                              size_t current, size_t expected )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** FE_DivergenceParameter error\n" ;
   mesg << "    bad storage depth for the field \""
        << f_name << "\"" << std::endl ;
   mesg << "       expected : " << expected << std::endl ;
   mesg << "       current  : " << current << std::endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
