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

#include <DOC_Type.hh>
#include <DOC_Typedef.hh>
#include <DOC_Class.hh>

#include <map>

//--------------------------------------------------------------------
static std::map< std::string, DOC_Type* > DUMMY_MAP ;
//--------------------------------------------------------------------

//--------------------------------------------------------------------
DOC_Type*
DOC_Type::to_type( void ) 
//--------------------------------------------------------------------
{
   return this ;
}


//--------------------------------------------------------------------
DOC_Type*
DOC_Type::create( PEL_Object* a_owner,
              string const& a_name ) 
//--------------------------------------------------------------------
{
   return new DOC_Type( a_owner, a_name ) ;
}



//--------------------------------------------------------------------
DOC_Type::DOC_Type( PEL_Object* a_owner,
            std::string const& a_name ) 
//--------------------------------------------------------------------
   : DOC_Symbol( a_owner ),
     identificateur( a_name ),
     estPointeur( false ), estReference( false ), is_constant( false ),
     type_complet( a_name )
{
}



//--------------------------------------------------------------------
DOC_Type::~DOC_Type( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
void 
DOC_Type::set_pointer( void ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::set_pointer" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   estPointeur = true ;
}



//--------------------------------------------------------------------
void 
DOC_Type::set_reference( void ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::set_reference" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   estReference = true ;
}



//--------------------------------------------------------------------
void 
DOC_Type::set_constant( void ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::set_constant" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   is_constant = true ;
}



//--------------------------------------------------------------------
string 
DOC_Type::modifiers( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::modifiers" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   string modif ;
   
   if( is_constant ) 
   {
      modif += " const" ;
   }
   if( estReference ) 
   {
      modif += "&" ;
   }
   if( estPointeur ) 
   {
      modif += "*" ;
   }
   return modif ;
}



//--------------------------------------------------------------------
string const&
DOC_Type::name( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::name" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   return identificateur ;
}



//--------------------------------------------------------------------
string const&
DOC_Type::full_type_name( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::full_type_name" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   type_complet = identificateur + modifiers() ;
   return type_complet ;
}



//--------------------------------------------------------------------
string 
DOC_Type::type_reference( DOC_Writer const& sullizer ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::type_reference" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   string ret = full_type_name() ;
   
   DOC_Class const* cl = DOC_Class::search( identificateur ) ;
   if( cl!=0 && sullizer.source_reference() )
   {
      ret = sullizer.reference( cl, "", identificateur ) + modifiers() ;
   }
   // Search des typedef
   DOC_ClassItem const* elem ;   
   if( ( elem = DOC_Class::search_and_find( identificateur ) ) != 0 )
   {
      ret = sullizer.reference( elem->def_class(),
                                elem->name(),
                                identificateur ) + modifiers() ;
   }
   
   return ret ;
}



//--------------------------------------------------------------------
void 
DOC_Type::display( ostream& out ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::display" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   out << full_type_name() ;
}



//--------------------------------------------------------------------
bool 
DOC_Type::is_equal( PEL_Object const* other ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Type::is_equal" ) ;
   PEL_CHECK_PRE( const_cast<DOC_Type*>(this)->to_type()!=0 ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   DOC_Type const* t = static_cast<DOC_Type const*>( other ) ;
   bool result = ( full_type_name() == t->full_type_name() ) ;
   
   if( DOC_Tools::message() && identificateur!=t->identificateur )
   {
      string s1 = identificateur ;
      string s2 = t->identificateur ;
      
      if( ( s1.find( s2 )< s1.length() || s2.find( s1 )< s2.length() ) &&
          ( s1.find( "::" )< s1.length() || s2.find( "::" )< s2.length() ) )
      {
         DOC_Tools::warning( "Comparing same type from different namespaces : " + s1
                  + " " + s2 ) ;
      }
   }
   
   return result ;
}
