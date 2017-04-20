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

#include <PEL_Context.hh>

#include <PEL_assertions.hh>
#include <PEL_List.hh>
#include <PEL_Variable.hh>

#include <iostream>
#include <string>

//----------------------------------------------------------------------
PEL_Context:: PEL_Context( PEL_Object* a_owner ) 
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , OBSERVERS( PEL_List::create( this ) )
{
}

//----------------------------------------------------------------------
PEL_Context:: ~PEL_Context( void ) 
//----------------------------------------------------------------------
{
   OBSERVERS = 0 ;
}

//----------------------------------------------------------------------
void
PEL_Context:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Context:: print" ) ;
   print( os, indent_width, (PEL_Context const*) 0 ) ;
}

//----------------------------------------------------------------------
void
PEL_Context:: print( std::ostream& os, size_t indent_width,
                     PEL_Context const* ctx ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Context:: print" ) ;
   
   std::string bl( indent_width, ' ' ) ;
   for( size_t i=0 ; i<nb_variables() ; ++i  )
   {
      PEL_Variable const* var = variable(i) ;
      PEL_Data const* dat = value( var ) ;
      os << bl <<  "$" << var->name() ;
      if( ctx == 0 || !ctx->has_variable( var ) )
      {
         os << " = " ;
      }
      else
      {
         os << " == " ;
      }
      dat->print( os, 0 ) ;
      os << std::endl ;
   }   
}

//----------------------------------------------------------------------
void
PEL_Context:: attach_observer( PEL_Context* observer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Context:: attach_observer" ) ;
   PEL_ASSERT( observer != this ) ;
   OBSERVERS->extend( observer ) ;
}

//----------------------------------------------------------------------
void
PEL_Context:: detach_observer( PEL_Context* observer ) const
//----------------------------------------------------------------------
{
   OBSERVERS->remove( observer ) ;
}


//----------------------------------------------------------------------
void
PEL_Context:: update_observers( void ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<OBSERVERS->index_limit() ; i++ )
   {
      PEL_Context* ct = static_cast<PEL_Context*>( OBSERVERS->at(i) ) ;
      ct->update() ;
   }
}

//----------------------------------------------------------------------
void
PEL_Context:: notify_observers_of_my_destruction( void ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<OBSERVERS->index_limit() ; i++ )
   {
      PEL_Context* ct = static_cast<PEL_Context* >(OBSERVERS->at(i)) ;
      ct->update_for_destruction_of( this ) ;
   }
}

//----------------------------------------------------------------------
bool
PEL_Context:: variable_PRE( size_t i ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_variables() ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Context:: variable_POST( PEL_Variable const* result ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( has_variable( result ) ) ;
   PEL_ASSERT( result!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Context:: has_variable_PRE( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( var!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Context:: value_PRE( PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( var!=0 ) ;
   PEL_ASSERT( has_variable( var ) ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
PEL_Context:: value_POST( PEL_Data* result,
                          PEL_Variable const* var ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( var->data_type()==result->data_type() ) ;
   return true ;
}

