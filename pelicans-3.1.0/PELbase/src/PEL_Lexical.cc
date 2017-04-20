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

#include <PEL_Lexical.hh>

#include <PEL_assertions.hh>
#include <PEL_Module.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_String.hh>

#include <iostream>

//----------------------------------------------------------------------
PEL_Lexical* PEL_Lexical::the_common_owner = 0 ;
//----------------------------------------------------------------------



//----------------------------------------------------------------------
PEL_Lexical*
PEL_Lexical:: common_owner( void )
//----------------------------------------------------------------------
{
   if( the_common_owner==0 )
   {
      the_common_owner = new PEL_Lexical() ;
   }
   return the_common_owner ;
}



//----------------------------------------------------------------------
void
PEL_Lexical:: remove_all_lexical( void )
//----------------------------------------------------------------------
{
   if( the_common_owner!=0 )
   {
      the_common_owner->destroy() ;
      the_common_owner=0 ;
   }
}



//----------------------------------------------------------------------
PEL_Lexical*
PEL_Lexical:: create( PEL_Module* aModule )
//----------------------------------------------------------------------
{
   return new PEL_Lexical( common_owner(), aModule ) ;
}



//----------------------------------------------------------------------
PEL_Lexical*
PEL_Lexical:: create( PEL_Data* aSimpleValue )
//----------------------------------------------------------------------
{
   return new PEL_Lexical( common_owner(), aSimpleValue ) ;
}



//----------------------------------------------------------------------
PEL_Lexical*
PEL_Lexical:: create( PEL_List* aSimpleValue )
//----------------------------------------------------------------------
{
   return new PEL_Lexical( common_owner(), aSimpleValue ) ;
}



//----------------------------------------------------------------------
PEL_Lexical:: PEL_Lexical( void )
//----------------------------------------------------------------------
      :   PEL_Object( 0 ),
          myModule( 0 ),
          myData( 0 ),
          myList( 0 )
{
}



//----------------------------------------------------------------------
PEL_Lexical:: PEL_Lexical( PEL_Object* anOwner,
                           PEL_Module* aModule )
//----------------------------------------------------------------------
      :   PEL_Object( anOwner ),
          myModule( aModule ),
          myData( 0 ),
          myList( 0 )
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
PEL_Lexical:: PEL_Lexical( PEL_Object* anOwner,
                           PEL_Data* aSimpleValue )
//----------------------------------------------------------------------
      :   PEL_Object( anOwner ),
          myModule( 0 ),
          myData( aSimpleValue ),
          myList( 0 )
{
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( aSimpleValue->owner()==0 ) ;
   aSimpleValue->set_owner( this ) ;
}




//----------------------------------------------------------------------
PEL_Lexical:: PEL_Lexical( PEL_Object* anOwner,
                           PEL_List* aSimpleValue )
//----------------------------------------------------------------------
      :   PEL_Object( anOwner ),
          myModule( 0 ),
          myData( 0 ),
          myList( aSimpleValue )
{
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( aSimpleValue->owner()==0 ) ;
   aSimpleValue->set_owner( this ) ;
}




//----------------------------------------------------------------------
PEL_Lexical:: ~PEL_Lexical( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
bool
PEL_Lexical:: is_list( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: is_list" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return myList!=0 ;
}



//----------------------------------------------------------------------
bool
PEL_Lexical:: is_module( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: is_module" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return myModule!=0 ;
}



//----------------------------------------------------------------------
bool
PEL_Lexical:: is_data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: is_data" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return myData!=0 ;
}



//-------------------------------------------------------------------------
PEL_List *
PEL_Lexical:: to_list( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: to_module" ) ;
   PEL_CHECK_PRE( is_list() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return myList ;
}



//-------------------------------------------------------------------------
PEL_Module *
PEL_Lexical:: to_module( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: to_module" ) ;
   PEL_CHECK_PRE( is_module() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return myModule ;
}



//-------------------------------------------------------------------------
PEL_Data *
PEL_Lexical:: to_data( void )  
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Lexical:: to_data" ) ;
   PEL_CHECK_PRE( is_data() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Data *resu = myData ;
   return resu ;
}



//----------------------------------------------------------------------
bool
PEL_Lexical:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( this==the_common_owner || myModule!=0 || myData!=0 || myList!=0 ) ;
   PEL_ASSERT( myModule==0 || myModule->owner()!=0 ) ;
   return( true ) ;
}





