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

#include <PEL_String.hh>

#include <iostream>
#include <numeric>

#include <PEL_assertions.hh>

//----------------------------------------------------------------------------
PEL_String*
PEL_String:: create( PEL_Object* a_owner, std::string const& a )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: create" ) ;
   return( new PEL_String( a_owner, a ) ) ;
}



//----------------------------------------------------------------------------
PEL_String:: PEL_String( PEL_Object* a_owner, std::string const& a )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner ),
     str( a )
{
   PEL_LABEL( "PEL_String:: PEL_String" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
PEL_String:: ~PEL_String( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: ~PEL_String" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
PEL_Data::Type
PEL_String:: data_type( void) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: data_type" ) ;
   return( PEL_Data::String ) ;
}



//----------------------------------------------------------------------------
PEL_String*
PEL_String:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_String* resu = new PEL_String( a_owner, str ) ;
   PEL_CHECK_POST( resu!=0 ) ;
   return resu ;
}



//----------------------------------------------------------------------------
std::string const&
PEL_String:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE(ct) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( str ) ;
}



//----------------------------------------------------------------------------
bool
PEL_String:: is_equal( PEL_Object const* other ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   
   bool result = three_way_comparison( other ) == 0 ;
   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   
   return( result ) ;
}



//----------------------------------------------------------------------------
int
PEL_String:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   PEL_String const* tString = static_cast<PEL_String const*>( other ) ;

   int result = str.compare( tString->str ) ;
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return result ;
}



//-------------------------------------------------------------------------
size_t
PEL_String:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: hash_code" ) ;
   size_t resu = std::accumulate( str.begin(), str.end(), 0 ) ;
   return( resu ) ;
}



//----------------------------------------------------------------------
void
PEL_String:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_String:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   if(str.find('"')<str.length())
   {
      std::string dup(str) ;
      size_t idx ;
      while( (idx=dup.find('"'))<dup.length()) dup.replace( idx, 1, "'" ) ;
      os << space << "\"" << dup << "\"" ;
   }
   else
   {
      os << space << "\"" << str << "\"" ;
   }
   
}



//----------------------------------------------------------------------------
void
PEL_String:: set( std::string const& val )
//----------------------------------------------------------------------------
{
   str = val ;
}



