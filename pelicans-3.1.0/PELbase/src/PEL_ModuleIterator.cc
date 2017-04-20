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

#include <PEL_ModuleIterator.hh>

#include <PEL_assertions.hh>
#include <PEL_List.hh>
#include <PEL_ListItem.hh>
#include <PEL_Module.hh>



//----------------------------------------------------------------------
PEL_ModuleIterator*
PEL_ModuleIterator:: create( PEL_Object* a_owner, PEL_List const* a_list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleIterator:: create" ) ;
   PEL_CHECK_PRE( a_list != 0 ) ;

   PEL_ModuleIterator* result = new PEL_ModuleIterator( a_owner, a_list ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_ModuleIterator:: PEL_ModuleIterator( PEL_Object* a_owner,
                                         PEL_List const* a_list )
//----------------------------------------------------------------------
   : PEL_ListIterator( a_owner, a_list )
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
PEL_ModuleIterator:: ~PEL_ModuleIterator( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
PEL_Module*
PEL_ModuleIterator:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleIterator:: item" ) ;
   PEL_CHECK_PRE( item_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ; 

   PEL_Object* obj = PEL_ListIterator::item() ;
   
   PEL_Module* resu = dynamic_cast<PEL_Module*>( obj ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( item_POST( resu ) ) ;

   return( resu ) ;
}



//---------------------------------------------------------------------
bool
PEL_ModuleIterator:: invariant( void ) const
//---------------------------------------------------------------------
{
   PEL_ASSERT( PEL_ListIterator::invariant() ) ;

   return( true ) ;
}
