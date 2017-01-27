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

#include <DOC_ClassItem.hh>

#include <PEL.hh>
#include <DOC_Category.hh>
#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Tools.hh>
#include <DOC_Writer.hh>
#include <DOC_Type.hh>

size_t DOC_ClassItem::global_index = 0 ;

//--------------------------------------------------------------------
DOC_ClassItem*
DOC_ClassItem::to_class_item( void ) 
//--------------------------------------------------------------------
{
   return this ;
}

//--------------------------------------------------------------------
bool
DOC_ClassItem::declare( std::string const& a__name ) const
//--------------------------------------------------------------------
{
   return name()==a__name ;
}



//--------------------------------------------------------------------
DOC_ClassItem::DOC_ClassItem( Protection maprotection,
                      DOC_Category const* macategory,
                      DOC_Text const* mycomment ) 
//--------------------------------------------------------------------
      : DOC_Symbol( 0 ), 
     niveauProtection( maprotection ),
     maDOC_Category( macategory ),
     maDOC_Class( 0 ),
     my_index( global_index++ )
{
   static std::string aliased = "@SIGNET " ;
   
   if( mycomment!=0 )
   {
      myComment = DOC_Tools::text_comment( mycomment->text() ) ;
      size_t idx = myComment.find( aliased ) ;
      if( idx < myComment.length() ) 
      {
         size_t end = myComment.find( "\n", idx+ aliased.length() ) ;
         if( end >= myComment.length() ) end = myComment.length()-1 ;
         PEL_ASSERT( end > idx + aliased.length() ) ;
         
         alias = myComment.substr( idx+aliased.length(), end+1 ) ;
         myComment.replace( idx, end+aliased.length()+1, "" ) ;
      }
   }
   if( macategory==0 )
   {
      maDOC_Category = DOC_Category::create( "" ) ;
   }

}



//--------------------------------------------------------------------
DOC_ClassItem::~DOC_ClassItem( void ) 
//--------------------------------------------------------------------
{
}

//--------------------------------------------------------------------
size_t
DOC_ClassItem:: index( void ) const
//--------------------------------------------------------------------
{
   return( my_index ) ;
}



//--------------------------------------------------------------------
DOC_Category const*
DOC_ClassItem:: category( void ) const
//--------------------------------------------------------------------
{
   return maDOC_Category ;
}



//--------------------------------------------------------------------
string
DOC_ClassItem:: category_displaye( void ) const
//--------------------------------------------------------------------
{
   string result = maDOC_Category->name() ;
   if( protection()==Protected ) 
   {
      result += " (protected)" ;
   }
   return( result ) ;
}



//--------------------------------------------------------------------
bool
DOC_ClassItem::has_to_document( DOC_Class const* cl ) const
//--------------------------------------------------------------------
{
   Protection p = ( cl->is_plug_point() ?
                    Protected : Public ) ;
   bool result = DOC_Tools::private_doc() ||
      ( visible( p ) && category()->name()!="Hidden" ) ;
   
   return result ;
}



//--------------------------------------------------------------------
DOC_Method*
DOC_ClassItem::method( void )
//--------------------------------------------------------------------
{
   return 0 ;
}




//--------------------------------------------------------------------
DOC_Attribute*
DOC_ClassItem::attribute( void )
//--------------------------------------------------------------------
{
   return 0 ;
}



//--------------------------------------------------------------------
DOC_ClassItem::Protection
DOC_ClassItem::protection( void ) const
//--------------------------------------------------------------------
{
   return niveauProtection ;
}



//--------------------------------------------------------------------
void
DOC_ClassItem::attach( DOC_Class* classe ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_ClassItem::attach" ) ;
   PEL_CHECK_PRE( classe!=0 ) ;
   
   maDOC_Class = classe ;
   
}



//--------------------------------------------------------------------
DOC_Class const* 
DOC_ClassItem::def_class( void ) const 
//--------------------------------------------------------------------
{
   return maDOC_Class ;
}



//--------------------------------------------------------------------
bool 
DOC_ClassItem::visible( Protection level ) const 
//--------------------------------------------------------------------
{
   bool ret = niveauProtection<=level ;
   return ret ;
}



//--------------------------------------------------------------------
bool 
DOC_ClassItem::is_inheritable( void ) const 
//--------------------------------------------------------------------
{
   bool ret = niveauProtection!=Private ;
   return ret ;
}



//--------------------------------------------------------------------
string 
DOC_ClassItem::comment( void ) const 
//--------------------------------------------------------------------
{
   return myComment ;
}



//--------------------------------------------------------------------
void
DOC_ClassItem::inherit_category( DOC_ClassItem const* other) 
//--------------------------------------------------------------------
{
   maDOC_Category = other->category() ;
}



//--------------------------------------------------------------------
int
DOC_ClassItem:: three_way_comparison( PEL_Object const* other ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_ClassItem:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   DOC_ClassItem const* e2 = static_cast<DOC_ClassItem const*>( other ) ;
   
   int result = 0 ;
   if( protection() != e2->protection() )
   {
      result = ( protection() < e2->protection() ? -1 : 1 ) ;
   }
   else
   {
      result = category()->three_way_comparison( e2->category() ) ;
      if( result==0 )
      {
         int res = index() - e2->index() ;
         if( res<0 )
            result=-1;
         else if( res>0 )
            result=1;
      }
   }   
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return result ;
}


//--------------------------------------------------------------------
void
DOC_ClassItem:: set_bookmark( std::string const& a_name )
//--------------------------------------------------------------------
{
   alias = a_name ;
}


//--------------------------------------------------------------------
std::string const&
DOC_ClassItem:: bookmark( void ) const
//--------------------------------------------------------------------
{
   return( alias ) ;
}



