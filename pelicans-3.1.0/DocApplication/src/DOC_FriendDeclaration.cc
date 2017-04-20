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

#include <DOC_FriendDeclaration.hh>

#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Tools.hh>
#include <DOC_Writer.hh>
#include <DOC_Type.hh>
#include <PEL_assertions.hh>

//--------------------------------------------------------------------
DOC_FriendDeclaration*
DOC_FriendDeclaration::create( std::string const& a_name,
                  Protection laProtection,
                  DOC_Category const* laDOC_Category,
                  DOC_Text const* leComment ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_FriendDeclaration::create" ) ;
   
   return new DOC_FriendDeclaration( a_name, laProtection,
                        laDOC_Category, leComment ) ;
}



//--------------------------------------------------------------------
DOC_FriendDeclaration::DOC_FriendDeclaration( std::string const& a_name,
		    Protection laProtection,
                    DOC_Category const* laDOC_Category,
		    DOC_Text const* leComment ) 
//--------------------------------------------------------------------
   : DOC_ClassItem( laProtection, laDOC_Category, leComment ),
     myName( a_name )
{
}



//--------------------------------------------------------------------
DOC_FriendDeclaration::~DOC_FriendDeclaration( void ) 
//--------------------------------------------------------------------
{
}




//--------------------------------------------------------------------
string const&
DOC_FriendDeclaration::name( void ) const
//--------------------------------------------------------------------
{
   return myName ;
}



//--------------------------------------------------------------------
string 
DOC_FriendDeclaration::prototype( DOC_Writer& sullizer ) const
//--------------------------------------------------------------------
{
   string ret =  myName ;   
   return( ret ) ;
}



//--------------------------------------------------------------------
string
DOC_FriendDeclaration::signature( void ) const
//--------------------------------------------------------------------
{
   string ret = myName ;
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_FriendDeclaration::display( std::ostream& out ) const
//--------------------------------------------------------------------
{
   if( !comment().empty() )
   {
      out << "      " << comment() << std::endl ;
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
DOC_FriendDeclaration:: display_info( std::ostream& os, size_t indent_width ) const
//--------------------------------------------------------------------
{
   string w( ' ',indent_width ) ;
   os << w ;
   display(os) ;
}





//--------------------------------------------------------------------
bool
DOC_FriendDeclaration::is_equal( PEL_Object const* other ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_FriendDeclaration::is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   DOC_FriendDeclaration const* att = static_cast<DOC_FriendDeclaration const*>( other ) ;
   
   return name()==att->name() ;
}
