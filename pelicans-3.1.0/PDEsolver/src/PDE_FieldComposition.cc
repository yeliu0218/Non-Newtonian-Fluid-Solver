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

#include <PDE_FieldComposition.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>

#include <iostream>

//-------------------------------------------------------------------------
PDE_FieldComposition*
PDE_FieldComposition:: make( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dimensions,
                             PDE_SetOfDiscreteFields const* dfs )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( nb_sp_dimensions == 0 || nb_sp_dimensions == 1 ||
                  nb_sp_dimensions == 2 || nb_sp_dimensions == 3 ) ;
   PEL_CHECK_PRE( dfs != 0 ) ;

   std::string const& n = exp->string_data( "concrete_name" ) ;
   PDE_FieldComposition const* proto =
      static_cast<PDE_FieldComposition const*>(
                                    plugins_map()->item( n ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PDE_FieldComposition* result =
            proto->create_replica( a_owner, exp, nb_sp_dimensions, dfs ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_FieldComposition:: PDE_FieldComposition(
                                    std::string const& a_concrete_name )
//----------------------------------------------------------------------
  : PEL_Object( plugins_map() )
  , NAME( a_concrete_name )
  , PROTO( true )
  , VARS( 0 )
  , NB_VARS( 0 )
  , VAR_IDX( PEL::bad_index() )
  , GLOB_2_iVAR( 0 )
  , VAR_VALUES( 0, 0 )
  , LAWS( 0 )
  , COMPUTED( false )
  , VALUES_PROVIDER( 0 )
{
   PEL_LABEL( "PDE_FieldComposition:: PDE_FieldComposition(a_concrete_name)" ) ;
   
   plugins_map()->register_item( a_concrete_name, this ) ;

   PEL_CHECK_POST( nb_variables() == 0 ) ;
   PEL_CHECK_POST( name() == a_concrete_name ) ;
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PDE_FieldComposition:: PDE_FieldComposition( PEL_Object* a_owner,
                                             std::string const& a_name )
//----------------------------------------------------------------------
  : PEL_Object( a_owner )
  , NAME( a_name )
  , PROTO( false )
  , VARS( PEL_Vector::create( this, 0 ) )
  , NB_VARS( 0 )
  , VAR_IDX( PEL::bad_index() )
  , GLOB_2_iVAR( 0 )
  , VAR_VALUES( 0, 0 )
  , LAWS( PEL_Vector::create( this, 0 ) )
  , COMPUTED( false )
  , VALUES_PROVIDER( this )
{
   PEL_LABEL( "PDE_FieldComposition:: PDE_FieldComposition(a_owner,a_name)" ) ;
   
   PEL_CHECK_POST( nb_variables() == 0 ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PDE_FieldComposition:: ~PDE_FieldComposition( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: ~PDE_FieldComposition" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_FieldComposition:: do_the_links( PDE_SetOfFieldCompositions const* fcs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( fcs ) ) ;

   PEL_CHECK_POST( do_the_links_POST( fcs ) ) ;
}

//----------------------------------------------------------------------
bool 
PDE_FieldComposition:: has_consistent_internal_dependencies( 
                                PDE_FieldComposition const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: has_consistent_internal_dependencies" ) ;
   bool result = true ;

   for( size_t i=0 ; i<NB_VARS ; i++ )
   {
      PDE_DiscreteField const* f = 
                      static_cast<PDE_DiscreteField const*>( VARS->at( i ) ) ;
      result = result && other->is_a_variable( f ) ;
   }
   
   for( size_t i=0 ; i<LAWS->count() ; ++i )
   {
      PDE_FieldComposition const* fc = 
                     static_cast<PDE_FieldComposition const*>( LAWS->at(i) ) ;
      result = result && fc->has_consistent_internal_dependencies( other ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void 
PDE_FieldComposition:: complete_internal_dependencies( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: complete_internal_dependencies" ) ;
   
   for( size_t i=0 ; i<LAWS->count() ; ++i )
   {
      PDE_FieldComposition* fc = 
                          static_cast<PDE_FieldComposition*>( LAWS->at(i) ) ;
      
      fc->complete_internal_dependencies() ;
      fc->start_variable_iterator() ;      
      for( ; fc->valid_variable() ; fc->go_next_variable() )
      {
         PDE_DiscreteField const* f = fc->variable() ;
         if( !is_a_variable( f ) )
         {
            add_one_variable( f ) ;
         }
      }
   }
   PEL_CHECK_POST( has_consistent_internal_dependencies( this ) ) ;
}

//----------------------------------------------------------------------
std::string const&
PDE_FieldComposition:: name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( NAME ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: is_constant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: is_constant" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( nb_variables() == 0 ) ;
   for( size_t i=0 ; result && i<LAWS->index_limit() ; ++i )
   {
      PDE_FieldComposition const* c =
         static_cast<PDE_FieldComposition const*>( LAWS->at(i) ) ;
      result &= c->is_constant() ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_FieldComposition:: nb_variables( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: nb_variables" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( NB_VARS ) ;
}

//----------------------------------------------------------------------
bool 
PDE_FieldComposition:: is_a_variable( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: is_a_variable" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( VARS->item( ff ) != 0 ) ;
   return( result );
}

//----------------------------------------------------------------------
void
PDE_FieldComposition:: start_variable_iterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: start_variable_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   VAR_IDX = 0 ;
}

//----------------------------------------------------------------------
void
PDE_FieldComposition:: go_next_variable( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: go_next_variable" ) ;
   PEL_CHECK_PRE( valid_variable() ) ;
   PEL_CHECK_INV( invariant() ) ;

   ++VAR_IDX ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: valid_variable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: valid_variable" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( VAR_IDX < NB_VARS ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_FieldComposition:: variable( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: variable" ) ;
   PEL_CHECK_PRE( valid_variable() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_DiscreteField const* result =
               static_cast<PDE_DiscreteField const*>( VARS->at( VAR_IDX ) ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( is_a_variable( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: variable_value_is_set( PDE_DiscreteField const* ff,
                                              size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: variable_value_is_set" ) ;
   PEL_CHECK_PRE( is_a_variable( ff ) ) ;
   PEL_CHECK_PRE( ic < ff->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t idx = VALUES_PROVIDER->GLOB_2_iVAR( ff->id_number() ) ;
   bool result = VALUES_PROVIDER->VAR_VALUES( idx, ic ) != PEL::bad_double() ;
   
   return result ;
}

//----------------------------------------------------------------------
void
PDE_FieldComposition:: set_variable_value( PDE_DiscreteField const* ff,
                                           size_t ic,
                                           double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: set_variable_value" ) ;
   PEL_CHECK_PRE( is_a_variable( ff ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   size_t idx = GLOB_2_iVAR( ff->id_number() ) ;
   VAR_VALUES( idx, ic ) = x ;
   COMPUTED = false ;
   
   PEL_CHECK_POST( !is_computed() ) ;
   PEL_CHECK_POST( variable_value_is_set( ff, ic ) ) ;  
}

//----------------------------------------------------------------------
double
PDE_FieldComposition:: variable_value( PDE_DiscreteField const* ff,
                                       size_t ic ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: variable_value" ) ;
   PEL_CHECK_PRE( ff!=0 ) ;
   PEL_CHECK_PRE( is_a_variable( ff ) ) ;
   PEL_CHECK_PRE( variable_value_is_set( ff, ic ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_ASSERT( VALUES_PROVIDER!=0 ) ;
   
   size_t idx = VALUES_PROVIDER->GLOB_2_iVAR( ff->id_number() ) ;
   double result = VALUES_PROVIDER->VAR_VALUES( idx , ic ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
void 
PDE_FieldComposition:: compute( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: compute" ) ;
   PEL_CHECK_PRE( has_consistent_internal_dependencies( this ) ) ;
   PEL_CHECK_PRE( 
 FORALL( ( start_variable_iterator() ; valid_variable() ; go_next_variable() ),
    FORALL( ( size_t ic=0 ; ic<variable()->nb_components() ; ic++ ),
       variable_value_is_set( variable(), ic ) ) ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   compute_from_provider( this ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_computed() ) ;
}

//----------------------------------------------------------------------
void 
PDE_FieldComposition:: compute_from_provider( 
                                     PDE_FieldComposition const* provider )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: compute_from_provider" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   VALUES_PROVIDER = provider ;
   
   for( size_t i=0 ; i<LAWS->count() ; ++i )
   {
      PDE_FieldComposition* fc =
                  static_cast<PDE_FieldComposition*>( LAWS->at(i) ) ;
      PEL_CHECK( dynamic_cast<PDE_FieldComposition*>( LAWS->at(i) ) !=0 ) ;
      
      fc->compute_from_provider( VALUES_PROVIDER ) ;
   }
   compute_self() ;
   COMPUTED = true ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_computed() ) ;
}

//----------------------------------------------------------------------
bool 
PDE_FieldComposition:: is_computed( void ) const
//----------------------------------------------------------------------
{
    PEL_LABEL( "PDE_FieldComposition:: is_computed" ) ;

    return COMPUTED ;
}

//----------------------------------------------------------------------
void 
PDE_FieldComposition:: add_one_composition( PDE_FieldComposition* fc ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: add_one_composition" ) ;
   PEL_CHECK_PRE( fc != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( fc->depends_on( this ) )
   {
      PEL_Error::object()->raise_plain(
         "*** \""+NAME+"\" field composition error:\n"
         "    circular definition (recursive calls)" ) ;
   }
   LAWS->append( fc ) ;
   fc->start_variable_iterator() ;
   for( ; fc->valid_variable() ; fc->go_next_variable() )
   {
      PDE_DiscreteField const* f = fc->variable() ;
      if( !is_a_variable( f ) )
      {
         add_one_variable( f ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( depends_on( fc ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( fc->start_variable_iterator() ; fc->valid_variable() ;
                fc->go_next_variable() ),
              is_a_variable( fc->variable() ) ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: depends_on( PDE_FieldComposition const* fc ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: depends_on" ) ;
   PEL_CHECK_PRE( fc != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = ( fc == this ) ;
   for( size_t i=0 ; !result && i<LAWS->count() ; ++i )
   {
      PDE_FieldComposition const* other = 
                     static_cast<PDE_FieldComposition const*>( LAWS->at(i) ) ;
      result = other->depends_on( fc ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FieldComposition:: add_one_variable( PDE_DiscreteField const* df )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition:: add_one_variable" ) ;
   PEL_CHECK_PRE( df != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( !is_a_variable( df ) )
   {
      VARS->append( const_cast<PDE_DiscreteField*>( df ) ) ;
      size_t glob = df->id_number() ;
      size_t idx = VARS->index_of( df ) ;
      if( GLOB_2_iVAR.size()<=glob )
      {
         GLOB_2_iVAR.resize( glob+1 ) ;
      }
      GLOB_2_iVAR( glob ) = idx ;
      
      size_t nb_comp = PEL::max( df->nb_components(), VAR_VALUES.index_bound(1) ) ;
      VAR_VALUES.re_initialize( idx+1, nb_comp ) ;
      VAR_VALUES.set( PEL::bad_double() ) ;
      NB_VARS = VARS->count() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_a_variable( df ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( PROTO ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: do_the_links_PRE(
                           PDE_SetOfFieldCompositions const* fcs ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( fcs != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: do_the_links_POST(
                           PDE_SetOfFieldCompositions const* fcs ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: value_PRE( size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( ic<nb_components() ) ;
   PEL_ASSERT( is_computed() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: create_replica_PRE(
                             PEL_Object const* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( exp != 0 ) ;
   PEL_ASSERT( nb_sp_dims == 0 || nb_sp_dims == 1 ||
               nb_sp_dims == 2 || nb_sp_dims == 3 ) ;
   PEL_ASSERT( dfs != 0 ) ;
   PEL_ASSERT( is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: create_replica_POST(
                             PDE_FieldComposition const* result,
                             PEL_Object const* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_FieldComposition:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( EQUIVALENT( !is_a_prototype(),
                           VARS!=0 && LAWS!=0 ) ) ;
   PEL_ASSERT( IMPLIES( !is_a_prototype(),
                        VARS->count() == NB_VARS ) ) ;
   PEL_ASSERT( IMPLIES( !is_a_prototype(),
                        NB_VARS == VAR_VALUES.index_bound(0) ) ) ;
   PEL_ASSERT( IMPLIES( !is_a_prototype(),
                        NB_VARS == VARS->index_limit() ) ) ;
   PEL_ASSERT( IMPLIES( !is_a_prototype(),
                        LAWS->count() == LAWS->index_limit() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PDE_FieldComposition:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "PDE_FieldComposition descendant" ) ;
   return( result ) ;
}
