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

#include <PEL_KeywordDataPair.hh>

#include <doubleVector.hh>
#include <intVector.hh>

#include <PEL_Error.hh>
#include <PEL_KeyItemPair.hh>
#include <PEL_Root.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_Data.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>
#include <PEL_String.hh>

#include <iostream>

//----------------------------------------------------------------------
PEL_KeywordDataPair*
PEL_KeywordDataPair:: create( PEL_Object* a_owner,
                              PEL_String const* a_keyword, 
                              PEL_Data const* a_data )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: create" ) ;
   return( new PEL_KeywordDataPair( a_owner, a_keyword, a_data ) ) ;
} 



//----------------------------------------------------------------------
PEL_KeywordDataPair:: PEL_KeywordDataPair( PEL_Object* a_owner,
                                           PEL_String const* a_keyword, 
                                           PEL_Data const* a_data )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     pair( 0 )
{
   pair = PEL_KeyItemPair::create( this,
                                   const_cast<PEL_String*>(a_keyword), 
                                   const_cast<PEL_Data*>(a_data) ) ;

   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
std::string const&
PEL_KeywordDataPair:: keyword( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: keyword" ) ;
   return( static_cast<PEL_String*>( pair->key() )->to_string() ) ;
}



//----------------------------------------------------------------------
PEL_Data const* 
PEL_KeywordDataPair:: data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: data" ) ;
   return( static_cast<PEL_Data const*>( pair->item() ) ) ;
}



//----------------------------------------------------------------------
void
PEL_KeywordDataPair:: replace_data( PEL_Data const* a_data )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: replace_data" ) ;
   PEL_CHECK_INV( invariant() ) ;

   pair->set_item( const_cast<PEL_Data*>(a_data) ) ;
}



//----------------------------------------------------------------------
PEL_KeywordDataPair:: ~PEL_KeywordDataPair( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
size_t
PEL_KeywordDataPair:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: hash_code" ) ;
   return( pair->hash_code() ) ;  
}


//----------------------------------------------------------------------
bool
PEL_KeywordDataPair:: comparable( PEL_Object const* other ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: comparable" ) ;
   return PEL_Object::comparable( other ) ||
      pair->comparable( other ) ;
}


//----------------------------------------------------------------------
bool
PEL_KeywordDataPair:: is_equal( PEL_Object const* other ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_KeywordDataPair const* other_as_assignment =
                           dynamic_cast<PEL_KeywordDataPair const*>( other ) ;
   bool result ;
   if( other_as_assignment != 0 )
   {
      result = pair->is_equal( other_as_assignment->pair ) ;
   }
   else
   {
      result = pair->is_equal( other ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;

   return( result ) ;
}



//-------------------------------------------------------------------------
int
PEL_KeywordDataPair:: three_way_comparison( PEL_Object const* other ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: three_way_comparison" ) ;
   PEL_Error::object()->raise_not_implemented( this, "three_way_comparison" ) ;
   return( 1 ) ;
}



//----------------------------------------------------------------------
void
PEL_KeywordDataPair:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_KeywordDataPair:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << keyword() ;
   os << " = " ;
   bool is_str = data()->data_type()==PEL_Data::String ;
   
   if( is_str )
   {
      os << "\"" ;
   }
   data()->print( os, 0 ) ;
   if( is_str )
   {
      os << "\"" ;
   }
}
