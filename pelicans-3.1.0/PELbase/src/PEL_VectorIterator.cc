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

#include <PEL_VectorIterator.hh>

#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

//----------------------------------------------------------------------
PEL_VectorIterator*
PEL_VectorIterator:: create( PEL_Object* a_owner,
                             PEL_Vector const* vector )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: create" ) ;
   PEL_CHECK_PRE( vector!=0 ) ;
   
   PEL_VectorIterator* result =  new PEL_VectorIterator( a_owner, vector ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( EQUIVALENT( vector->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_VectorIterator:: PEL_VectorIterator( PEL_Object* a_owner,
                                         PEL_Vector const* vector )
//----------------------------------------------------------------------
   : PEL_Iterator( a_owner, vector ), 
     vec( vector ),
     valid( false ),
     current_i( static_cast<size_t>(~0) ),
     current_item( 0 )
{
   PEL_LABEL( "PEL_VectorIterator:: PEL_VectorIterator" ) ;
   start() ;
}



//----------------------------------------------------------------------
PEL_VectorIterator:: ~PEL_VectorIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: ~PEL_VectorIterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PEL_VectorIterator:: go_i_th( size_t i )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: go_i_th" ) ;
   PEL_CHECK_INV( invariant() ) ;

   current_i = i ;
   current_item = vec->at( current_i ) ;
   if( current_item != 0 )
   { 
      valid = true ;
   }
   else
   {
      valid = false ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( is_valid(), index_of_item()==i ) ) ;
}

//----------------------------------------------------------------------
void
PEL_VectorIterator:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   valid = false ;
   current_item = 0 ;

   for( current_i=0 ; current_i<vec->index_limit() ; current_i++ )
   {
      current_item = vec->at( current_i ) ;
      if( current_item != 0 )
      {
         valid = true ;
         break ;
      }
   }
   record_container_state_id() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
bool
PEL_VectorIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: is_valid" ) ;
   PEL_CHECK_PRE( is_valid_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( valid ) ;
}



//----------------------------------------------------------------------
void
PEL_VectorIterator:: go_next( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: go_next" ) ;
   PEL_CHECK_PRE( go_next_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   valid = false ;
   current_i++ ;
   current_item = 0 ;
   for( ; current_i<vec->index_limit() ; current_i++ )
   {
      current_item = vec->at( current_i ) ;
      if( current_item != 0 )
      {
         valid = true ;
         break ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PEL_VectorIterator:: index_of_item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: index_of_item" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( current_i ) ;
}

//----------------------------------------------------------------------
PEL_Object*
PEL_VectorIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VectorIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK_POST( item_POST( current_item ) ) ;

   return( current_item ) ;
}



//----------------------------------------------------------------------
bool
PEL_VectorIterator:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Iterator::invariant() ) ;

   PEL_ASSERT( IMPLIES( !container_has_been_modified() && is_valid(),
                        ( current_item != 0 &&
                          vec->at( current_i ) == current_item ) ) ) ;

   return( true ) ;
}



