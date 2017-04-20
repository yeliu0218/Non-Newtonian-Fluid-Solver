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

#include <PEL_ContextPair.hh>

#include <PEL_assertions.hh>
#include <PEL_Vector.hh>
#include <PEL_Variable.hh>

//----------------------------------------------------------------------
PEL_ContextPair*
PEL_ContextPair:: create( PEL_Object* a_owner,
                          PEL_Context const* first,
                          PEL_Context const* second ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: create" ) ;
   PEL_CHECK_PRE( EQUIVALENT( first == 0, second == 0 ) ) ;

   PEL_ContextPair* result = new PEL_ContextPair( a_owner, first, second  ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( IMPLIES( first == 0 || second == 0,
                            result->nb_variables() == 0 ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; first != 0 && i<first->nb_variables() ; ++i ),
         result->has_variable( first->variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; second != 0 && i<second->nb_variables() ; ++i ),
         result->has_variable( second->variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         first->has_variable( result->variable( i ) ) ||
         second->has_variable( result->variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         IMPLIES(
            second->has_variable( result->variable( i ) ),
            result->value( result->variable( i ) ) == second->value( result->variable( i ) ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         IMPLIES(
            !second->has_variable( result->variable( i ) ),
            first->has_variable( result->variable( i ) ) &&
            result->value( result->variable( i ) ) == first->value( result->variable( i ) ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ContextPair:: PEL_ContextPair( PEL_Object* a_owner,
                                   PEL_Context const* first,
                                   PEL_Context const* second ) 
//----------------------------------------------------------------------
   : PEL_Context( a_owner )
   , CT1( 0 )
   , CT2( 0 )
   , VALUES( PEL_Vector::create( this, PEL_Variable::nb_objects() ) )
{
   PEL_CHECK( EQUIVALENT( first != 0, second != 0 ) ) ;
   if( first != 0 || second != 0 )
   {
      re_initialize( first, second ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PEL_ContextPair:: re_initialize( PEL_Context const* first,
                                 PEL_Context const* second )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: re_initialize" ) ;
   PEL_CHECK_PRE( first != 0 && second != 0 ) ;
   PEL_CHECK_PRE( first != second ) ;
   PEL_CHECK_PRE( first != this ) ;
   PEL_CHECK_PRE( this != second ) ;
   
   if( CT1!=first || CT2!=second )
   {
      if( (CT1!=0) && (CT1!=first) )
      {
         CT1->detach_observer( this ) ;
      }
      if( (CT2!=0) && (CT2!=second) )
      {
         CT2->detach_observer( this ) ;
      }
      CT1 = first ;
      CT2 = second ;
      if( CT1!=0 )
      {
         CT1->attach_observer( this ) ;
      }
      if( CT2!=0 )
      {
         CT2->attach_observer( this ) ;
      }
      update() ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<first->nb_variables() ; ++i ),
         has_variable( first->variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<second->nb_variables() ; ++i ),
         has_variable( second->variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         first->has_variable( variable( i ) ) ||
         second->has_variable( variable( i ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         IMPLIES(
            second->has_variable( variable( i ) ),
            value( variable( i ) ) == second->value( variable( i ) ) ) ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         IMPLIES(
            !second->has_variable( variable( i ) ),
            first->has_variable( variable( i ) ) &&
            value( variable( i ) ) == first->value( variable( i ) ) ) ) ) ;
}

//----------------------------------------------------------------------
PEL_ContextPair*
PEL_ContextPair:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: create_clone" ) ;

   PEL_ContextPair* result = new PEL_ContextPair( a_owner, CT1, CT2  ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_variables() == nb_variables() ) ;
   PEL_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         result->has_variable( variable( i ) ) &&
         result->value( variable( i ) ) == value( variable( i ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ContextPair:: ~PEL_ContextPair( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: ~PEL_ContextPair" ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   if( CT1!=0 )
   {
      CT1->detach_observer( this ) ;
   }
   if( CT2!=0 )
   {
      CT2->detach_observer( this ) ;
   }
//   PEL_Context::update_dependencies(this) ;
   
   notify_observers_of_my_destruction() ;
}

//----------------------------------------------------------------------
size_t
PEL_ContextPair:: nb_variables( void ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: nb_variables" ) ;
   
   size_t result = VALUES->count() ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Variable const* 
PEL_ContextPair:: variable( size_t i ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: variable" ) ;
   PEL_CHECK_PRE( variable_PRE(i) ) ;
   
   PEL_Variable const* result = 0 ;
   for( size_t j=0 ; j<VALUES->index_limit() ; j++ )
   {
      if( VALUES->at(j)!=0 )
      {
         if( i==0 )
         {
            result = PEL_Variable::object(j) ;
            break ;
         }
         i-- ;
      }
   }

   PEL_CHECK_POST( variable_POST( result ) ) ;
   return( result ) ;
      
}

//----------------------------------------------------------------------
bool
PEL_ContextPair:: has_variable( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: has_variable" ) ;
   PEL_CHECK_PRE( has_variable_PRE( var ) ) ;
   
   size_t idx = var->id_number() ;
   bool result = ( VALUES->index_limit()>idx ) &&
                 ( VALUES->at( idx )!=0 ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data* 
PEL_ContextPair:: value( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: value" ) ;
   PEL_CHECK_PRE( value_PRE( var ) ) ;
   
   PEL_Data* result = 
                   static_cast<PEL_Data*>( VALUES->at( var->id_number() ) ) ;

   PEL_CHECK_POST( value_POST( result, var ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ContextPair:: update( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ContextPair:: update" ) ;
   
   VALUES->resize(PEL_Variable::nb_objects()) ;
   if( CT1!=0 )
   {
      for( size_t i=0 ; i<CT1->nb_variables() ; ++i )
      {
	 PEL_Variable const* var = CT1->variable( i ) ;
	 size_t idx = var->id_number() ;
	 VALUES->set_at( idx, CT1->value(var) ) ;
      }
   }
   if( CT2!=0 )
   {
      for( size_t i=0 ; i<CT2->nb_variables() ; ++i )
      {
	 PEL_Variable const* var = CT2->variable( i ) ;
	 size_t idx = var->id_number() ;
	 VALUES->set_at( idx, CT2->value( var ) ) ;
      }
   }
   update_observers() ;
}

//----------------------------------------------------------------------
void
PEL_ContextPair:: update_for_destruction_of( PEL_Context const* subject )
//----------------------------------------------------------------------
{
   if( subject == CT1 )
   {
      CT1 = 0 ;
   }
   else
   {
      PEL_ASSERT( subject == CT2 ) ;
      CT2 = 0 ;
   }
   VALUES->clear() ;
}

//----------------------------------------------------------------------
bool
PEL_ContextPair:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( VALUES != 0 ) ;
   return( true ) ;
}




