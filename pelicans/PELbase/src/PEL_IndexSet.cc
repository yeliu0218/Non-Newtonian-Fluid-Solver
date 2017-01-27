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

#include <PEL_IndexSet.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <iostream>
#include <string>

//-------------------------------------------------------------------------
PEL_IndexSet*
PEL_IndexSet:: create( PEL_Object* a_owner ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: create" ) ;
   
   PEL_IndexSet* result = new PEL_IndexSet( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->elements().size() == 0 ) ;
   PEL_CHECK_POST( result->id() == PEL::bad_index() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_IndexSet:: PEL_IndexSet( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ID( PEL::bad_index() )
   , SET( 0 )
{
}

//-------------------------------------------------------------------------
PEL_IndexSet*
PEL_IndexSet:: create( PEL_Object* a_owner,
                       size_t_vector const& vec,
                       size_t a_id ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: create" ) ;
   
   PEL_IndexSet* result = new PEL_IndexSet( a_owner, vec, a_id ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->elements().size() == vec.size() ) ;
   PEL_CHECK_POST( FORALL( (size_t i=0;i<vec.size()-1;i++),
                           result->elements()(i)<=result->elements()(i+1) )) ;
   PEL_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           vec.has(result->elements()(i) )) ) ;
   PEL_CHECK_POST( result->id() == a_id ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_IndexSet:: PEL_IndexSet( PEL_Object* a_owner,
                             size_t_vector const& vec,
                             size_t a_id )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ID( a_id )
   , SET( vec )
{
   PEL_LABEL( "PEL_IndexSet:: PEL_IndexSet" ) ;
   SET.sort_increasingly() ;
}

//-------------------------------------------------------------------------
PEL_IndexSet:: ~PEL_IndexSet( void )
//-------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_IndexSet:: re_initialize( size_t_vector const& vec, 
                              size_t a_id )
//----------------------------------------------------------------------
{
   SET = vec ;
   SET.sort_increasingly() ;
   ID = a_id ;
}

//----------------------------------------------------------------------
size_t
PEL_IndexSet:: id( void ) const
//----------------------------------------------------------------------
{
   return( ID ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PEL_IndexSet:: elements( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: elements" ) ;
   
   return( SET ) ;
}

//----------------------------------------------------------------------
bool
PEL_IndexSet:: is_equal( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   
   bool result = ( three_way_comparison( other ) == 0 ) ;
   
   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
PEL_IndexSet:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   PEL_CHECK( dynamic_cast<PEL_IndexSet const*>( other ) != 0 ) ;
   PEL_IndexSet const* oo = static_cast<PEL_IndexSet const*>( other ) ;
                  
   size_t const ss = SET.size() ;
   size_t const os = oo->SET.size() ;

   int result = 0 ;
   if( ss<os )
   {
      result = -1 ;
   }
   else if( ss>os )
   {
      result = 1 ;
   }
   else if( ss>0 )
   {
      for( size_t i=0 ; i<ss ; i++ )
      {
         size_t si = SET(i) ;
         size_t osi = oo->SET(i) ;
         
         if( si < osi )
         {
            result = -1 ;
            break ;
         }
         else if( si > osi )
         {
            result = 1 ;
            break ;
         }
      }
      
   }
   
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PEL_IndexSet:: hash_code( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: hash_code" ) ;
   size_t result = 0 ;
   for( size_t i=0 ; i<SET.size() ; i++ )
   {
      result += SET(i) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_IndexSet:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IndexSet:: print" ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "index: " << ID << ",  " << SET ;
}
