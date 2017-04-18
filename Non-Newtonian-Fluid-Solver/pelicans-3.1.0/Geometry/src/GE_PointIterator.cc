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

#include <GE_PointIterator.hh>

#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_assertions.hh>
#include <GE_Point.hh>



//----------------------------------------------------------------------
GE_PointIterator*
GE_PointIterator:: create( PEL_Object* a_owner, PEL_List const* a_list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: create" ) ;
   PEL_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<a_list->count() ; ++i ),
              dynamic_cast<GE_Point const*>( a_list->at(i) )!=0 ) ) ;
   
   GE_PointIterator* result =  new GE_PointIterator( a_owner, a_list ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( EQUIVALENT( a_list->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
GE_PointIterator:: GE_PointIterator( PEL_Object* a_owner,
                                     PEL_List const* aList )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     list_it( PEL_ListIterator::create( this, aList ) )
{
   PEL_LABEL( "GE_PointIterator:: GE_PointIterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
GE_PointIterator:: ~GE_PointIterator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: ~GE_PointIterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
bool
GE_PointIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: is_valid" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( list_it->is_valid() ) ;
}



//----------------------------------------------------------------------
void
GE_PointIterator:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: start" ) ;
   PEL_CHECK_INV( invariant() ) ;

   list_it->start() ;
   
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
void
GE_PointIterator:: go_next( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: go_next" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ;

   list_it->go_next() ;

   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
GE_Point const*
GE_PointIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointIterator:: item" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_INV( invariant() ) ; 

   GE_Point const* result =
                       static_cast<GE_Point const*>( list_it->item() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 ) ;

   return( result ) ;
}



//---------------------------------------------------------------------
bool
GE_PointIterator:: invariant( void ) const
//---------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;

   PEL_ASSERT( list_it!=0 && list_it->owner()==this ) ;

   return( true ) ;
}
