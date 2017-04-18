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

#include <GE_Color.hh>

#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::string ;
using std::ostream ; using std::endl ;
using std::ostringstream ;

int GE_Color::NEXT_ID = 0 ;

stringVector GE_Color::NAMES = stringVector( 0)  ;
intArray2D GE_Color::CONNECTIVITY = intArray2D( 0, 0 ) ;

struct GE_Color_ERROR
{
   static void n0( string const& name ) ;
   static void n1( string const& name ) ;
   static void n2( string const& composite, string const& component ) ;
} ;

//-------------------------------------------------------------------------
GE_Color:: GE_Color( PEL_Object* a_owner, std::string const& a_name )
//------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , MY_NAME( a_name )
   , ID( NEXT_ID++ )
   , IS_COMPOSITE( false )
   , COMPOSING_COLORS( 0 )
{
   COMPOSING_COLORS = PEL_ListIdentity::create( this ) ;
   NAMES.append( MY_NAME ) ;
}

//------------------------------------------------------------------------
GE_Color:: ~GE_Color( void )
//------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
GE_Color const*
GE_Color:: object( std::string const& a_name )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: object" ) ;

   GE_Color* result = static_cast<GE_Color*>( colors()->item( a_name ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_Color:: null_color( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: null_color" ) ;
   
   static GE_Color const* result = 0 ;
   if( result==0 )
   {
      PEL_ASSERT( !exist( "null color" ) ) ;
      extend( "null color" ) ;
      result = object( "null color" ) ;
   }
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( !result->is_composite() ) ;
   PEL_CHECK_POST( result->name()=="null color" ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_Color:: halo_color( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: halo_color" ) ;
   
   static GE_Color const* result = 0 ;
   if( result==0 )
   {
      PEL_ASSERT( !exist( "halo color" ) ) ;
      extend( "halo color" ) ;
      result = object( "halo color" ) ;
   }
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( !result->is_composite() ) ;
   PEL_CHECK_POST( result->name()=="halo color" ) ;
   return result ;
}

//------------------------------------------------------------------------
bool
GE_Color:: exist( std::string const& a_name )
//------------------------------------------------------------------------
{
    PEL_LABEL( "GE_Color:: exist" ) ;
    PEL_CHECK_PRE( a_name.length() > 0 ) ;
    
    bool result = colors()->has( a_name ) ;
    
    PEL_CHECK_POST( IMPLIES( result, object(a_name)->name()==a_name ) ) ;
    return( result ) ;
}

//-------------------------------------------------------------------------
void
GE_Color:: extend( std::string const& a_name )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: extend" ) ;
   PEL_CHECK_PRE( a_name.length() > 0 ) ;
   
   if( !colors()->has( a_name ) )
   {
      GE_Color* color = new GE_Color( colors(), a_name ) ; 
      color->IS_COMPOSITE = false ;
      colors()->register_item( a_name, color ) ;
   }
   else
   {
      GE_Color const* color =
         static_cast<GE_Color const*>( colors()->item( a_name ) ) ;
      if( color->is_composite() ) GE_Color_ERROR::n0( a_name ) ;
   }
   
   PEL_CHECK_POST( !object(a_name)->is_composite() ) ;
}

//-------------------------------------------------------------------------
void
GE_Color:: extend( std::string const& a_name,
                   stringVector const& a_name_list )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: extend" ) ;
   PEL_CHECK_PRE( a_name.length() > 0 ) ;
   PEL_CHECK_PRE( a_name_list.size() > 0 ) ;
   
   if( colors()->has( a_name ) )
   {
      GE_Color const* color =
         static_cast<GE_Color const*>( colors()->item( a_name ) ) ;
      bool same = ( a_name_list.size() == color->COMPOSING_COLORS->count() ) ;
      for( size_t i=0 ; i<a_name_list.size() ; ++i )
      {
         same = same && color->has( a_name_list(i) ) ;
      }
      if( !same ) GE_Color_ERROR::n1( a_name ) ;
   }
   else
   {
      GE_Color* color = new GE_Color( colors(), a_name ) ;
      color->IS_COMPOSITE = true ;
      for( size_t i=0 ; i<a_name_list.size() ; ++i )
      {
         GE_Color const* cmp = object( a_name_list( i ) ) ;
         if( cmp->is_composite() ) GE_Color_ERROR::n2( a_name, 
                                                       a_name_list( i ) ) ;
         color->COMPOSING_COLORS->append( const_cast<GE_Color*>( cmp ) ) ;
      }
      colors()->register_item( a_name, color ) ;
   }
   
   PEL_CHECK_POST( object(a_name)->is_composite() ) ;
   PEL_CHECK_POST( FORALL( (size_t i=0 ; i<a_name_list.size() ; ++i),
                           object(a_name)->has( a_name_list(i) ) ) ) ;
}

//------------------------------------------------------------------------
stringVector const&
GE_Color:: color_table( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: color_table" ) ;
   return( NAMES ) ;
}

//------------------------------------------------------------------------
intArray2D const&
GE_Color:: color_table_connectivity( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: color_table_connectivity" ) ;
   size_t n = NAMES.size() ;

   if( CONNECTIVITY.index_bound(0) != n )
   {
      CONNECTIVITY.re_initialize(n,n) ;
      for( size_t i=0 ; i<n; i++ )
      {
         GE_Color const* ci = GE_Color::object( NAMES(i) ) ;
         if( ci->is_composite() ) 
         {
            for( size_t j=0 ; j<n; j++ )
            {
               GE_Color const* cj = GE_Color::object( NAMES(j) ) ;
               if( !cj->is_composite() && ci->is_matching( cj ) )
               {
                  CONNECTIVITY(i,j) = 1 ;
               }
            }
         }
         else
            CONNECTIVITY(i,i) = 1 ;
      }
   }
   return( CONNECTIVITY ) ;
}

//------------------------------------------------------------------------
int
GE_Color:: identifier( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: identifier" ) ;

   int result = ID ;

   PEL_CHECK_POST( result >= 0 ) ;
   PEL_CHECK_POST( color_table()(result)==name() ) ;   
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_Color:: name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: name" ) ;
   return( MY_NAME ) ;
}

//------------------------------------------------------------------------
bool
GE_Color:: is_composite( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: is_composite" ) ;
   return( IS_COMPOSITE ) ;
}

//------------------------------------------------------------------------
bool
GE_Color:: has( std::string const& a_name ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: has" ) ;
   PEL_CHECK_PRE( is_composite() ) ;

   bool result = false ;
   if( colors()->has( a_name ) )
   {
      GE_Color const* col =
               static_cast<GE_Color const*>( colors()->item( a_name ) ) ;
      result = COMPOSING_COLORS->has( col ) ;
   }

   return( result ) ;
}

//------------------------------------------------------------------------
bool
GE_Color:: is_overlapping( GE_Color const* other ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: is_overlapping" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   
   bool result = false ;
   if( is_composite() && other->is_composite() )
   {
      PEL_Iterator* it = COMPOSING_COLORS->create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() && !result ; it->go_next() )
      {
         result = result || other->is_matching(
            static_cast<GE_Color const*>( it->item() ) ) ;
         
      }
      it->destroy() ; it=0 ;
   }
   else
   {
      result = is_matching( other ) ;
   }
   return( result ) ;
}

//------------------------------------------------------------------------
bool
GE_Color:: is_matching( GE_Color const* other ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: is_matching" ) ;
   PEL_CHECK_PRE( other != 0 ) ;
   PEL_CHECK_PRE( !is_composite() || !other->is_composite() ) ;
   
   bool result = ( this == other )
      || ( is_composite() && has( other->name() ) )
      || ( other->is_composite() && other->has( name() ) ) ;
              
   PEL_CHECK_POST( result == ( 
     ( name() == other->name() ) ||
     (  is_composite() && !other->is_composite() && has(other->name()) ) ||
     ( !is_composite() &&  other->is_composite() && other->has(name()) ) 
                             ) );

   return( result ) ;
}

//-------------------------------------------------------------------------
void
GE_Color:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << name() ;
}

//-------------------------------------------------------------------------
PEL_ObjectRegister*
GE_Color:: colors( void )
//-------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
             PEL_ObjectRegister::create( PEL_Root::object(), "GE_Color" ) ;
   return( result ) ;
}

//internal-----------------------------------------------------------------
void
GE_Color_ERROR:: n0( string const& name )
//internal-----------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Attempt to define a composite color and a leaf color" << endl ;
   mesg << "with the same : \"" << name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-----------------------------------------------------------------
void
GE_Color_ERROR:: n1( string const& name )
//internal-----------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "The definition of the composite color" << endl ;
   mesg << "   of name \"" << name << "\"" << endl ;
   mesg << "is inconsistent with a previous definition." << endl ;
   mesg << "(There might be an existing leaf color with" << endl ;
   mesg << "the same name)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-----------------------------------------------------------------
void
GE_Color_ERROR:: n2( string const& composite, string const& component )
//internal-----------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Attempt to use the composite color" << endl ;
   mesg << "   of name \"" << component << "\"" << endl ;
   mesg << "in the definition of the composite color" << endl ;
   mesg << "   of name \"" << composite << "\"" << endl ;
   mesg << "Composite colors cannot be defined from " << endl << endl ;
   mesg << "other composite colors" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
