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

#include <PEL_Root.hh>

#include <PEL_assertions.hh>
#include <PEL_Vector.hh>

#include <iostream>

//-------------------------------------------------------------------------
PEL_Object*
PEL_Root:: object( size_t i ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Root:: object" ) ;

   PEL_Vector* objs = objects() ;
   if( objs->index_limit()<=i )
   {
      objs->resize( i+1 ) ;
   }
   if( objs->at(i) == 0 )
   {
      objs->set_at( i, new PEL_Root( 0 ) ) ;
   }
   
   PEL_Object* result = objs->at(i) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_Root:: PEL_Root( PEL_Object* a_owner ) 
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//-------------------------------------------------------------------------
PEL_Root:: ~PEL_Root( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_Root:: cleanup( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Root:: cleanup" ) ;
   
   PEL_Vector* objs = objects() ;
   objs->destroy_items_and_remove_section( 0, objs->index_limit() ) ;
   objects( true ) ;
}

//-------------------------------------------------------------------------
PEL_Vector*
PEL_Root:: objects( bool destroy )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Root:: objects" ) ;

   static PEL_Vector* result = 0 ;
   if( result == 0 )
   {
      result = PEL_Vector::create( 0, 0 ) ;
   }
   if( destroy )
   {
      result->destroy() ;
      result = 0 ;
   }
   
   PEL_CHECK_POST( IMPLIES( destroy, result == 0 ) ) ;
   PEL_CHECK_POST( IMPLIES( !destroy, result != 0 && result->owner() == 0 ) ) ;
   return( result ) ;
}


