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

#include <DOC_Attribute.hh>

#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Tools.hh>
#include <DOC_Writer.hh>
#include <DOC_Type.hh>
#include <PEL_assertions.hh>

//--------------------------------------------------------------------
DOC_Attribute*
DOC_Attribute::create( std::string const& a_name,
                  DOC_Type* typ,
                  Protection laProtection,
                  DOC_Category const* laDOC_Category,
                  DOC_Text const* leComment ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Attribute::create" ) ;
   
   return new DOC_Attribute( a_name, typ, laProtection,
                        laDOC_Category, leComment ) ;
}



//--------------------------------------------------------------------
DOC_Attribute::DOC_Attribute( std::string const& a_name,
		    DOC_Type* typ,
		    Protection laProtection,
                    DOC_Category const* laDOC_Category,
		    DOC_Text const* leComment ) 
//--------------------------------------------------------------------
   : DOC_ClassItem( laProtection, laDOC_Category, leComment ),
     myName( a_name ),
     myDOC_Type( typ ),
     estStatique( false ),
     init( 0 )
{
   PEL_CHECK( typ!=0 ) ;
}



//--------------------------------------------------------------------
DOC_Attribute::~DOC_Attribute( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
void
DOC_Attribute::set_static( void ) 
//--------------------------------------------------------------------
{
   estStatique = true ;
}



//--------------------------------------------------------------------
void
DOC_Attribute::initialize( DOC_Symbol const* s ) 
//--------------------------------------------------------------------
{
   init = s ;
}



//--------------------------------------------------------------------
bool
DOC_Attribute::is_static( void ) const
//--------------------------------------------------------------------
{
   return estStatique ;
}



//--------------------------------------------------------------------
string const&
DOC_Attribute::name( void ) const
//--------------------------------------------------------------------
{
   return myName ;
}



//--------------------------------------------------------------------
string 
DOC_Attribute::prototype( DOC_Writer& sullizer ) const
//--------------------------------------------------------------------
{
   string ret =  myDOC_Type->type_reference( sullizer ) + ' ' + myName ;
   if( estStatique )
   {
      ret = "static " + ret ;
   }
   if( init!=0 ) ret += " = " + init->text() ;
   
   return( ret ) ;
}



//--------------------------------------------------------------------
string
DOC_Attribute::signature( void ) const
//--------------------------------------------------------------------
{
   string ret =  myDOC_Type->full_type_name() + ' ' + myName ;
   if( estStatique )
   {
      ret = "static " + ret ;
   }
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_Attribute::display( std::ostream& out ) const
//--------------------------------------------------------------------
{
   if( !comment().empty() )
   {
      out << "      " << comment() << std::endl ;
   }
   if( estStatique )
   {
      out << " CLASS " ;
   }
   switch( protection() )
   {
      case Public : out << " PUBLIC " ;
	 break ;
      case Protected : out << " PROTECTED " ;
	 break ;
      case Private : out << " PRIVATE " ;
	 break ;
   }

   out  << signature() ;
}



//--------------------------------------------------------------------
void
DOC_Attribute:: display_info( std::ostream& os, size_t indent_width ) const
//--------------------------------------------------------------------
{
   string w( ' ',indent_width ) ;
   os << w ;
   display(os) ;
}



//--------------------------------------------------------------------
DOC_Attribute*
DOC_Attribute::attribute( void )
//--------------------------------------------------------------------
{
   return this ;
}



//--------------------------------------------------------------------
DOC_Type const* 
DOC_Attribute::type( void ) const 
//--------------------------------------------------------------------
{
   return myDOC_Type ;
}



//--------------------------------------------------------------------
bool
DOC_Attribute::is_equal( PEL_Object const* other ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Attribute::is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   DOC_Attribute const* att = static_cast<DOC_Attribute const*>( other ) ;
   
   return name()==att->name() ;
}
