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

#include <PEL_ListIdentity.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL_ListIterator.hh>

using std::ostream ;
using std::endl ;



//----------------------------------------------------------------------
PEL_ListIdentity*
PEL_ListIdentity:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ListIdentity:: create" ) ;
   PEL_ListIdentity* result = new PEL_ListIdentity( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->index_limit() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
PEL_ListIdentity:: PEL_ListIdentity( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_List( a_owner )
{
   PEL_LABEL( "PEL_ListIdentity:: PEL_ListIdentity" ) ;
}



//-------------------------------------------------------------------------
PEL_ListIdentity:: ~PEL_ListIdentity( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ListIdentity:: ~PEL_ListIdentity" ) ;
}



//-------------------------------------------------------------------------
PEL_ListIdentity*
PEL_ListIdentity:: create_clone( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ListIdentity:: create_clone" ) ;
   PEL_ListIdentity* result = create( a_owner ) ;
   
   PEL_Iterator* it = create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      result->append( it->item() ) ;
   }
   it->destroy() ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
bool
PEL_ListIdentity:: matching_items( PEL_Object const* object1,
                                   PEL_Object const* object2 ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ListIdentity:: matching_items" ) ;
   PEL_CHECK_PRE( matching_items_PRE( object1, object2 ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( object1 == object2 ) ;
}



