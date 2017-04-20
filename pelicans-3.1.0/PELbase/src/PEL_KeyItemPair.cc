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

#include <PEL_KeyItemPair.hh>

#include <PEL_assertions.hh>

//-------------------------------------------------------------------------
PEL_KeyItemPair*
PEL_KeyItemPair:: create( PEL_Object* a_owner,
                          PEL_Object* a_key,
                          PEL_Object* a_item  ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: create" ) ;
   return( new PEL_KeyItemPair( a_owner, a_key, a_item ) ) ;
}



//-------------------------------------------------------------------------
PEL_KeyItemPair:: PEL_KeyItemPair( PEL_Object* a_owner,
                                   PEL_Object* a_key,
                                   PEL_Object* a_item ) 
//-------------------------------------------------------------------------
   : PEL_Object( a_owner ), 
     the_key( a_key ),
     the_item( a_item )
{
   PEL_LABEL( "PEL_KeyItemPair:: PEL_KeyItemPair" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_KeyItemPair:: ~PEL_KeyItemPair( void ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: ~PEL_KeyItemPair" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_KeyItemPair:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: hash_code" ) ;
   return( key()->hash_code() ) ;  
}


//-------------------------------------------------------------------------
bool
PEL_KeyItemPair:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( key() != 0 ) ;

   return( true ) ;
} 


//-------------------------------------------------------------------------
PEL_Object*
PEL_KeyItemPair:: item( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: item" ) ;
   return( the_item ) ;
} 


//-------------------------------------------------------------------------
void
PEL_KeyItemPair:: set_item( PEL_Object* object )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: set_item" ) ;
   the_item = object ;
} 


//-------------------------------------------------------------------------
PEL_Object*
PEL_KeyItemPair:: key( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: key" ) ;
   return( the_key ) ;
}



//-------------------------------------------------------------------------
bool
PEL_KeyItemPair:: is_equal( PEL_Object const* other ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_KeyItemPair const* other_as_pair =
                           dynamic_cast<PEL_KeyItemPair const*>( other ) ;
   bool resu ;
   if( other_as_pair != 0 )
   {
      resu = key()->is_equal( other_as_pair->key() ) ;
   }
   else
   {
      resu = key()->is_equal( other ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_equal_POST( resu, other ) ) ;

   return( resu ) ;
}



//-------------------------------------------------------------------------
bool
PEL_KeyItemPair:: comparable( PEL_Object const* other ) const 
//-------------------------------------------------------------------------
{
   return PEL_Object::comparable( other ) ||
      key()->comparable( other ) ;
}



//-------------------------------------------------------------------------
int
PEL_KeyItemPair:: three_way_comparison( PEL_Object const* other ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeyItemPair:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_KeyItemPair const* other_as_pair =
                             dynamic_cast<PEL_KeyItemPair const*>( other ) ;
   int resu ;
   if( other_as_pair != 0 )
   {
      resu = key()->three_way_comparison( other_as_pair->key() ) ;
   }
   else
   {
      resu = key()->three_way_comparison( other ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( three_way_comparison_POST( resu, other ) ) ;

   return( resu );   
}



