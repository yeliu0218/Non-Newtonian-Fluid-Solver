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

#include <PDE_FieldCompositionExp.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Data.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <iostream>

struct PDE_FieldCompositionExp_ERROR {
   static void n0( std::string const& name, std::string const& var_name ) ;
   static void n1( std::string const& name, std::string const& law_name ) ;
   static void n2( std::string const& name ) ;
} ;

//----------------------------------------------------------------------
PDE_FieldCompositionExp*
PDE_FieldCompositionExp:: create( PEL_Object* a_owner,
                                  std::string const& a_name,
                                  PEL_Data* expression,
                                  PDE_SetOfDiscreteFields const* dfs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldCompositionExp:: create" ) ;
   PEL_CHECK_PRE( expression != 0 ) ;
   PEL_CHECK_PRE( expression->owner()==0 ) ;
   PEL_CHECK_PRE( dfs != 0 ) ;
   
   PDE_FieldCompositionExp* result = new PDE_FieldCompositionExp( a_owner, 
                                                                  a_name, 
                                                                  expression,
                                                                  dfs ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( expression->is_under_ownership_of( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_FieldComposition*
PDE_FieldCompositionExp:: create_replica(
                             PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const
//----------------------------------------------------------------------
{
    PEL_LABEL( "PDE_FieldCompositionExp:: create_replica" ) ;
    PEL_CHECK( create_replica_PRE( a_owner, exp, nb_sp_dims, dfs ) ) ;
    
    PDE_FieldComposition* result = 0 ;
    PEL_Error::object()->raise_internal(
                "Class \"PDE_FieldCompositionExp\" is not pluggable" ) ;
    
    PEL_CHECK( create_replica_POST( result, a_owner, exp, nb_sp_dims, dfs ) ) ;
    return( result ) ;
}

//----------------------------------------------------------------------
PDE_FieldCompositionExp:: PDE_FieldCompositionExp(
                                       PEL_Object* a_owner,
                                       std::string const& a_name,
                                       PEL_Data* expression,
                                       PDE_SetOfDiscreteFields const* dfs )
//----------------------------------------------------------------------
   : PDE_FieldComposition( a_owner, a_name )
   , IS_CSTE( false )
   , EXPR( expression )
   , CONTEXT( PEL_ContextSimple::create( this ) )
   , FIELDS( PEL_Vector::create( this, 0 ) )
   , VAL( 1 )
   , INNER_LAW_NAMES( 0 )
   , INNER_LAWS( PEL_Vector::create( this, 0 ) )
{
   expression->set_owner( this ) ;
   if( expression->data_type() != PEL_Data::DoubleVector )
   {
      PDE_FieldCompositionExp_ERROR:: n2( name() ) ;
   }
   PEL_List* lst = PEL_List::create( 0 ) ;
   EXPR->declare( lst ) ;
   if( lst->count() == 0 )
   {
      VAL = EXPR->to_double_vector() ;
      IS_CSTE = true ;
      compute() ;
   }
   else
   {
      for( size_t i=0 ; i<lst->count() ; i++ )
      {
         PEL_Variable const* var = static_cast<PEL_Variable*>( lst->at( i ) ) ;
         std::string var_name = var->name().substr( 3, var->name().length()-3 ) ;
         if( var->name().substr(0,3) != "DV_" )
         {
            PDE_FieldCompositionExp_ERROR:: n0( name(), var->name() ) ;
         }
         if( dfs->has( var_name ) )
         {
            PDE_DiscreteField const* field = dfs->item( var_name ) ;
            add_one_variable( field ) ;
            FIELDS->append( const_cast<PDE_DiscreteField*>( field ) ) ;
            doubleVector v( field->nb_components() ) ;
            v.set( PEL::bad_double() ) ;
            PEL_DoubleVector* val = PEL_DoubleVector::create( CONTEXT, v ) ;
            CONTEXT->extend( var, val ) ;
         }
         else
         {
            INNER_LAW_NAMES.append( var_name ) ;
         }
      }
   }
   lst->destroy() ; lst = 0 ;
}

//----------------------------------------------------------------------
PDE_FieldCompositionExp:: ~PDE_FieldCompositionExp( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldCompositionExp:: ~PDE_FieldCompositionExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void 
PDE_FieldCompositionExp:: do_the_links( PDE_SetOfFieldCompositions const* fcs )
//----------------------------------------------------------------------
{
    PEL_LABEL( "PDE_FieldCompositionExp:: do_the_links" ) ;
    PEL_CHECK_PRE( do_the_links_PRE( fcs )  ) ;

    if( !IS_CSTE )
    {
       for( size_t j=0 ; j<INNER_LAW_NAMES.size() ; j++ )
       {
          std::string var_name = INNER_LAW_NAMES(j) ;
          if( !fcs->has( var_name ) )
          {
             PDE_FieldCompositionExp_ERROR::n1( name(), var_name ) ;
          }
          PDE_FieldComposition* compo = fcs->item( var_name ) ;
          INNER_LAWS->append( compo ) ;
          add_one_composition( compo ) ;
          doubleVector v( compo->nb_components() ) ;
          v.set( PEL::bad_double() ) ;
          PEL_DoubleVector* val = PEL_DoubleVector::create( CONTEXT, v ) ;
       
          CONTEXT->extend( PEL_Variable::object( "DV_"+var_name ), val ) ;
       }
    }
  
    PEL_CHECK_POST( do_the_links_POST( fcs )  ) ;
}

//----------------------------------------------------------------------
size_t
PDE_FieldCompositionExp:: nb_components( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldCompositionExp:: nb_components" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( VAL.size() ) ;
}

//----------------------------------------------------------------------
void 
PDE_FieldCompositionExp:: compute_self( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldCompositionExp:: compute_self" ) ;

   if( !IS_CSTE )
   {
      for( size_t i=0 ; i<FIELDS->index_limit() ; i++ )
      {
         PDE_DiscreteField const* f = 
            static_cast<PDE_DiscreteField const*>( FIELDS->at(i) ) ;
         doubleVector v( f->nb_components() ) ;
         for( size_t j=0 ; j<f->nb_components() ; j++ )
         {
            v( j ) = variable_value( f, j ) ;
         }
         PEL_Variable const* var = PEL_Variable::object( "DV_"+f->name() ) ;
         PEL_DoubleVector* vec =
            static_cast<PEL_DoubleVector*>( CONTEXT->value( var ) ) ;
         vec->set( v ) ;
      }
      for( size_t i=0 ; i<INNER_LAWS->index_limit() ; i++ )
      {
         PDE_FieldComposition const* l = 
            static_cast<PDE_FieldComposition*>( INNER_LAWS->at( i ) ) ;
         doubleVector v( l->nb_components() ) ;
         for( size_t j=0 ; j<l->nb_components() ; j++ )
         {
            v( j ) = l->value( j ) ;
         }
         PEL_Variable const* var =
                          PEL_Variable::object( "DV_"+INNER_LAW_NAMES(i) ) ;
         PEL_DoubleVector* vec =
            static_cast<PEL_DoubleVector*>( CONTEXT->value( var ) ) ;
         vec->set( v ) ;
      }
      VAL = EXPR->to_double_vector( CONTEXT ) ;
   }
}

//----------------------------------------------------------------------
double 
PDE_FieldCompositionExp:: value( size_t ic ) const
//----------------------------------------------------------------------
{
    PEL_LABEL( "PDE_FieldCompositionExp:: value" ) ;
    PEL_CHECK_PRE( value_PRE( ic ) ) ;
    
    return( VAL( ic ) ) ;
}

//internal--------------------------------------------------------------
void
PDE_FieldCompositionExp_ERROR:: n0( std::string const& name,
                                    std::string const& var_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** field composition of name : \""+name+"\"\n" ;
   mesg += "    invalid variable name : \"$"+var_name+"\"\n" ;
   mesg += "    Rem : discrete field or field composition variable\n" ;
   mesg += "          should start with \"$DV_\"" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void
PDE_FieldCompositionExp_ERROR:: n1( std::string const& name,
                                    std::string const& law_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** field composition of name : \""+name+"\"\n" ;
   mesg += "    variable $DV_"+law_name+" can not be evaluated :\n" ;
   mesg += "       - unknown discrete field of name \"" + law_name + "\"\n" ;
   mesg += "       - unknown field composition of name \"" + law_name + "\"" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void
PDE_FieldCompositionExp_ERROR:: n2( std::string const& name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** field composition of name : \""+name+"\"\n" ;
   mesg += "    a doubleVector expression is expected" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}
