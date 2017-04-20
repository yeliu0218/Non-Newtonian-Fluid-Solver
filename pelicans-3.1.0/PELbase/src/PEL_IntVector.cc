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

#include <PEL_IntVector.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL_Container.hh>
#include <PEL_Iterator.hh>
#include <intVector.hh>

//----------------------------------------------------------------------------
PEL_IntVector*
PEL_IntVector:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_IntVector* result = 0 ;
   result = new PEL_IntVector( a_owner, dv ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_IntVector*
PEL_IntVector:: create( PEL_Object* a_owner,
                        intVector const& aIntVector )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: create(intVector)" ) ;
   PEL_IntVector* result = new PEL_IntVector( a_owner, aIntVector ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->to_int_vector() == aIntVector ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_IntVector:: PEL_IntVector( PEL_Object* a_owner,
                                 intVector const& aIntVector )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , dv( aIntVector )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_IntVector:: PEL_IntVector( PEL_Object* a_owner,
                               PEL_Container const* aDataList,
                               bool& error )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , dv( aDataList->count() )
{
   PEL_Iterator* list_iterator = aDataList->create_iterator( 0 ) ;
   size_t cpt=0 ;
   for( list_iterator->start() ;
        list_iterator->is_valid() && !error ;
        list_iterator->go_next() )
   {
      PEL_Data const* val =
         dynamic_cast<PEL_Data const*>( list_iterator->item() ) ;
      error |= ( val==0  ||
                 val->data_type()!=Int ||
                 !val->value_can_be_evaluated( 0 ) ) ;
      if( !error ) dv(cpt++) = val->to_int() ;
   }
   
   list_iterator->destroy() ; list_iterator=0 ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_IntVector*
PEL_IntVector:: create( PEL_Object* a_owner,
                        PEL_Container const* aDataList ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: create(PEL_Container)" ) ;
   PEL_CHECK_PRE( aDataList!=0 ) ;
   PEL_CHECK_PRE( aDataList->count()!=0 ) ;

   bool error = false ;
   PEL_IntVector* result = new PEL_IntVector( 0, aDataList, error ) ;
   if( error )
   {
      result->destroy() ;
      result = 0 ;
   }
   else if( a_owner != 0 )
   {
      result->set_owner( a_owner ) ;
   }

   PEL_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->to_int_vector().size() == aDataList->count() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_IntVector:: ~PEL_IntVector( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
PEL_IntVector:: set( intVector const& other )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: set" ) ;
   dv=other ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<other.size() ; i++ ),
                           to_int_vector()(i)==other(i) ) ) ;
}

//----------------------------------------------------------------------------
PEL_Data::Type
PEL_IntVector:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( IntVector ) ;
}

//----------------------------------------------------------------------------
intVector const& 
PEL_IntVector:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE(ct) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( dv ) ;
}

//----------------------------------------------------------------------------
void
PEL_IntVector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntVector:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   PEL_CHECK_INV( invariant() ) ;
   os << space << "< " ;
   for( size_t i=0 ; i<dv.size() ; i++ )
   {
       os << dv(i) << " " ;
   }
   os << " > " ;
}

//-------------------------------------------------------------------------
bool
PEL_IntVector:: invariant( void ) const 
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Data::invariant() ) ;
   return true ;
}
