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

#include <FE_FieldCompositionParameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

FE_FieldCompositionParameter const* 
FE_FieldCompositionParameter:: PROTOTYPE = new FE_FieldCompositionParameter() ;

struct FE_FieldCompositionParameter_ERROR
{
   static void n0( std::string const& param_name,
                   std::string const& law_name ) ;
   static void n1( std::string const& param_name,
                   std::string const& field_name ) ;
} ;

//-------------------------------------------------------------------------
FE_FieldCompositionParameter*
FE_FieldCompositionParameter:: create( PEL_Object* a_owner,
                                       std::string const& a_name,
                                       PDE_DomainAndFields const* a_dom,
                                       PDE_FieldComposition* a_compo,
                                       size_t a_field_level )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_dom != 0 ) ;
   PEL_CHECK_PRE( a_compo != 0 ) ;

   FE_FieldCompositionParameter* result = 
      new FE_FieldCompositionParameter( a_owner,
                                        a_name, a_dom, a_compo, a_field_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_FieldCompositionParameter:: FE_FieldCompositionParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_FieldCompositionParameter" )
   , LAW( 0 )
   , FIELDS_LEVEL( PEL::bad_index() )
{
}

//-------------------------------------------------------------------------
FE_FieldCompositionParameter*
FE_FieldCompositionParameter:: create_replica(
                                      PEL_Object* a_owner,
        		              PDE_DomainAndFields const* dom,
                                      PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   // Name :
   std::string const& a_name = exp->string_data( "name" ) ;

   // `PDE_FieldComposition::' object :
   PDE_SetOfFieldCompositions const* compos =
                                          dom->set_of_field_compositions() ;
   if( !compos->has( exp->string_data( "field_composition_name" ) ) )
   {
      FE_FieldCompositionParameter_ERROR::n0(
         a_name,
         exp->string_data( "field_composition_name" ) ) ;
   }
   PDE_FieldComposition* a_compo =
               compos->item( exp->string_data( "field_composition_name" ) ) ;

   // Fields level :
   size_t a_field_level = exp->int_data( "fields_level" ) ;
   
   FE_FieldCompositionParameter* result =
      new FE_FieldCompositionParameter( a_owner,
                                        a_name, dom, a_compo, a_field_level ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_FieldCompositionParameter:: FE_FieldCompositionParameter(
                               PEL_Object* a_owner,
                               std::string const& a_name,
                               PDE_DomainAndFields const* a_dom,
                               PDE_FieldComposition* a_compo,
                               size_t a_field_level )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, a_name )
   , LAW( a_compo )
   , FIELDS_LEVEL( a_field_level )
{
   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      if( field->storage_depth()<=FIELDS_LEVEL )
      {
         FE_FieldCompositionParameter_ERROR::n1( name(), field->name() ) ;
      }
   }
}

//-------------------------------------------------------------------------
FE_FieldCompositionParameter:: ~FE_FieldCompositionParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_FieldCompositionParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( LAW->nb_components() ) ;
}

//-------------------------------------------------------------------------
double
FE_FieldCompositionParameter:: cell_value_at_IP(
                                               FE_TimeIterator const* t_it,
                                               PDE_LocalFEcell const* fe,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      for( size_t i=0 ; i<field->nb_components() ; ++i )
      {
         LAW->set_variable_value(
            field, i, fe->value_at_IP( field, FIELDS_LEVEL, i ) ) ;
      }
   }
   LAW->compute() ;
   
   double result = LAW->value( ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_FieldCompositionParameter:: cell_value_at_pt(
                                               FE_TimeIterator const* t_it,
                                               PDE_LocalFEcell const* fe,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      for( size_t i=0 ; i<field->nb_components() ; ++i )
      {
         LAW->set_variable_value(
            field, i, fe->value_at_pt( field, FIELDS_LEVEL, i ) ) ;
      }
   }
   LAW->compute() ;
   
   double result = LAW->value( ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_FieldCompositionParameter:: bound_value_at_IP(
                                               FE_TimeIterator const* t_it,
                                               PDE_LocalFEbound const* fe,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      for( size_t i=0 ; i<field->nb_components() ; ++i )
      {
         LAW->set_variable_value(
            field, i, fe->value_at_IP( field, FIELDS_LEVEL, i ) ) ;
      }
   }
   LAW->compute() ;
   
   double result = LAW->value( ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_FieldCompositionParameter:: bound_value_at_pt(
                                               FE_TimeIterator const* t_it,
                                               PDE_LocalFEbound const* fe,
                                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      for( size_t i=0 ; i<field->nb_components() ; ++i )
      {
         LAW->set_variable_value(
            field, i, fe->value_at_pt( field, FIELDS_LEVEL, i ) ) ;
      }
   }
   LAW->compute() ;
   
   double result = LAW->value( ic ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_FieldCompositionParameter:: prepare_for_value_on_cells(
                                                PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   
   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      fe->require_field_calculation( field, PDE_LocalFE::N ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_FieldCompositionParameter:: prepare_for_value_on_bounds(
                                               PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;
   
   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      fe->require_field_calculation( field, PDE_LocalFE::N ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_FieldCompositionParameter:: prepare_for_value_on_sides(
                                                PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldCompositionParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;
   
   LAW->start_variable_iterator() ;
   for( ; LAW->valid_variable() ; LAW->go_next_variable() )
   {
      PDE_DiscreteField const* field = LAW->variable() ;
      fe->require_field_calculation( field, PDE_LocalFE::N ) ;
   }
}

//internal--------------------------------------------------------------
void 
FE_FieldCompositionParameter_ERROR:: n0( std::string const& param_name,
                                         std::string const& law_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** FE_FieldCompositionParameter of name : \""+param_name +"\"\n" ;
   mesg += "    unexpected field composition \""+law_name+"\"" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldCompositionParameter_ERROR:: n1( std::string const& param_name,
                                         std::string const& field_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** FE_FieldCompositionParameter of name : \""+param_name +"\"\n" ;
   mesg += "    incompatible storage depth for the field \""+field_name+"\"" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}
