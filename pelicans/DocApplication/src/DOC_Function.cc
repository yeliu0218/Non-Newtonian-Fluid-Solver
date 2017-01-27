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

#include <DOC_Function.hh>

#include <DOC_Class.hh>
#include <DOC_Method.hh>
#include <DOC_Type.hh>

//--------------------------------------------------------------------
int DOC_Function::nb_blancs = 0 ;
//--------------------------------------------------------------------



//--------------------------------------------------------------------
DOC_Function*
DOC_Function:: to_function( void )
//--------------------------------------------------------------------
{
   return this ;
}



//--------------------------------------------------------------------
DOC_Function*
DOC_Function::create( PEL_Object* a_owner,
                  std::string const& theDOC_Text ) 
//--------------------------------------------------------------------
{
   return new DOC_Function( a_owner, theDOC_Text ) ;
}



//--------------------------------------------------------------------
DOC_Function::DOC_Function( PEL_Object* a_owner,
                    std::string const& theDOC_Text ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        myDOC_Symbol( theDOC_Text ),
        mesDOC_Arguments( 0 ),
        infixe( false ),
        symb( 0 ),
        unaire( false ),
        identificateur( true )
{
}



//--------------------------------------------------------------------
DOC_Function*
DOC_Function::create( PEL_Object* a_owner,
                  std::string const& theDOC_Text,
                  DOC_Symbol * arg ) 
//--------------------------------------------------------------------
{
   return new DOC_Function( a_owner, theDOC_Text, arg ) ;
}



//--------------------------------------------------------------------
DOC_Function::DOC_Function( PEL_Object* a_owner,
                    std::string const& theDOC_Text,
                    DOC_Symbol * arg ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        myDOC_Symbol( theDOC_Text ),
        infixe( false ),
        symb( 0 ),
        unaire( true ),
        identificateur( false ),
        mesDOC_Arguments( PEL_List::create( this ) )
{
   mesDOC_Arguments->append( arg ) ;   
}



//--------------------------------------------------------------------
DOC_Function*
DOC_Function::create( PEL_Object* a_owner,
                  std::string const& theDOC_Text,
                  DOC_Symbol * arg1,
                  DOC_Symbol * arg2,
                  bool est_infixe ) 
//--------------------------------------------------------------------
{
   return new DOC_Function( a_owner, theDOC_Text, arg1, arg2, est_infixe ) ;
}



//--------------------------------------------------------------------
DOC_Function::DOC_Function( PEL_Object* a_owner,
                    std::string const& theDOC_Text,
                    DOC_Symbol * arg1,
                    DOC_Symbol * arg2,
                    bool est_infixe ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        myDOC_Symbol( theDOC_Text ),
        mesDOC_Arguments( PEL_List::create(this) ),
        infixe( est_infixe ),
        symb( 0 ),
        unaire( false ),
        identificateur( false )
{
   mesDOC_Arguments->append( arg1 ) ;
   mesDOC_Arguments->append( arg2 ) ;
}



//--------------------------------------------------------------------
DOC_Function*
DOC_Function::create( PEL_Object* a_owner,
                  std::string const& theDOC_Text,
                  DOC_Symbol * arg1,
                  DOC_Symbol * arg2,
                  DOC_Symbol * arg3 ) 
//--------------------------------------------------------------------
{
   return new DOC_Function( a_owner, theDOC_Text, arg1, arg2, arg3 ) ;
}



//--------------------------------------------------------------------
DOC_Function::DOC_Function( PEL_Object* a_owner,
                    std::string const& theDOC_Text,
                    DOC_Symbol * arg1,
                    DOC_Symbol * arg2,
                    DOC_Symbol * arg3 ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        myDOC_Symbol( theDOC_Text ),
        mesDOC_Arguments( PEL_List::create( this ) ),
        infixe( false ),
        symb( 0 ),
        unaire( false ),
        identificateur( false )
{
   mesDOC_Arguments->append( arg1 ) ;
   mesDOC_Arguments->append( arg2 ) ;
   mesDOC_Arguments->append( arg3 ) ;
}



//--------------------------------------------------------------------
DOC_Function*
DOC_Function::create( PEL_Object* a_owner,
                  DOC_Function * fonc,
                  PEL_List * lst ) 
//--------------------------------------------------------------------
{
   return new DOC_Function( a_owner, fonc, lst ) ;
}



//--------------------------------------------------------------------
DOC_Function::DOC_Function( PEL_Object* a_owner,
                    DOC_Function * fonc,
                    PEL_List* lst ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        mesDOC_Arguments( lst ),
        infixe( false ),
        symb( fonc ),
        unaire( false ),
        identificateur( false )
{
   myDOC_Symbol = fonc->text() ;
}



//--------------------------------------------------------------------
DOC_Function::~DOC_Function( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
std::string
DOC_Function::text( void ) const
//--------------------------------------------------------------------
{
   return referenced_string( 0, 0, 0 ) ;
}



//--------------------------------------------------------------------
std::string
DOC_Function:: referenced_string( DOC_Class const* classe,
                              DOC_Method const* method,
                              DOC_Writer const* red ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Function:: referenced_string" ) ;
   
   PEL_CHECK_PRE( red!=0 || ( classe==0 && method==0 ) ) ;
   
   string ret ;
   string symbol_fonc = myDOC_Symbol ;
   
   DOC_ClassItem const* elem = 0 ;
   
   
   if( classe!=0 && ( elem = classe->find( myDOC_Symbol ) ) != 0 )
   {
      DOC_Class const* classe_src = classe ;
      if( elem->def_class()->is_upper_class() )
         classe_src = elem->def_class() ;
      symbol_fonc  = red->reference( classe_src,
                                      elem->name(),
                                      myDOC_Symbol ) ;
   }
   else if( red!=0 && ( elem = DOC_Class::search_and_find( myDOC_Symbol ) ) != 0 )
   {
      symbol_fonc  = red->reference( elem->def_class(),
                                      elem->name(),
                                      myDOC_Symbol ) ;
   }
   if( elem==0 )
   {
      if( symb!=0 )
      {
         symbol_fonc = symb->referenced_string( classe,method,red ) ;
      }
      else
      {
         if( method!=0 && method->has_argument( myDOC_Symbol ) )
         {
            symbol_fonc = red->underlined_text( 1, myDOC_Symbol ) ;
         }
      }
   }
   
   if( identificateur )
   {
      
      /* DOC_Argument de la method */
      if( elem==0 && method!=0 && method->has_argument( myDOC_Symbol ) )
      {
         ret = red->underlined_text( 1, symbol_fonc ) ;
      }
      else
      {
         ret = symbol_fonc ;
      }
   }
   else
   {
      if( infixe )
      {
         DOC_Symbol const* arg1 = arg(0) ;
         DOC_Symbol const* arg2 = arg(mesDOC_Arguments->index_limit()-1);
         string bl = " " ;
         if( myDOC_Symbol=="->" || myDOC_Symbol=="." )
         {
            bl = "" ;
         }
         ret += arg1->referenced_string( classe, method, red ) ;
         ret += bl + myDOC_Symbol + bl ;
         ret += arg2->referenced_string( classe, method, red ) ;
      }
      else if( unaire )
      {
         DOC_Symbol const* aarg = arg(0);
         if( myDOC_Symbol.empty() )
         {
            ret =  "(" + aarg->referenced_string( classe, method, red ) + ")" ;
         }
         else
         {
            ret = myDOC_Symbol + aarg->referenced_string( classe, method, red ) ;
         }
      }
      else if( myDOC_Symbol=="IMPLIES" || myDOC_Symbol=="EQUIVALENT" )
      {
         PEL_ASSERT( mesDOC_Arguments!= 0 && mesDOC_Arguments->count()==2 ) ;
         DOC_Symbol const* arg1 = arg(0);
         DOC_Symbol const* arg2 = arg(mesDOC_Arguments->count()-1);
         ret = arg1->referenced_string( classe, method, red ) ;
         if( myDOC_Symbol=="IMPLIES" )
         {
            ret += " ==&gt; " ;
         }
         else if( myDOC_Symbol=="EQUIVALENT" )
         {
            ret += " &lt;==&gt; " ;
         }
         ret += arg2->referenced_string( classe, method, red ) ;
      }     
      else if( myDOC_Symbol.find( "_cast" )< myDOC_Symbol.length() )
      {
         PEL_ASSERT( mesDOC_Arguments!= 0 && mesDOC_Arguments->count()==2 ) ;
         DOC_Type const* arg1 = arg(0)->to_type() ;
         DOC_Symbol const* arg2 = arg(mesDOC_Arguments->count()-1) ;
         
         ret = myDOC_Symbol + "< " ;
         ret += arg1->type_reference( *red ) ;
         ret += " >(" ;
         ret += arg2->referenced_string( classe, method, red ) ;
         ret += ")" ;
      }     
      else if( myDOC_Symbol=="FORALL" || myDOC_Symbol=="EXISTS" )
      {
         PEL_ASSERT( mesDOC_Arguments!= 0 && mesDOC_Arguments->count()==2 ) ;
         DOC_Symbol const* arg1 = arg(0) ;
         DOC_Symbol const* arg2 = arg(mesDOC_Arguments->count()-1) ;
         ret = myDOC_Symbol + arg1->referenced_string( classe, method, red ) ;
         ret += "\n" + string( 2*(nb_blancs+1), ' ' ) ;
         nb_blancs ++ ;
         ret += arg2->referenced_string( classe, method, red ) ;
         nb_blancs -- ;
      }     
      else if( mesDOC_Arguments!= 0 && mesDOC_Arguments->count()>0 )
      {
         ret = symbol_fonc + "( " ;
         if( mesDOC_Arguments!= 0 )
         {
            bool prem = true ;
            PEL_Iterator* it = mesDOC_Arguments->create_iterator(0) ;
            
            for( it->start() ; it->is_valid() ; it->go_next() )
            {
               if( !prem ) ret += ", " ;
               prem = false ;
               DOC_Symbol*s = static_cast<DOC_Symbol*>( it->item() ) ;
               PEL_CHECK( dynamic_cast<DOC_Symbol*>( s ) !=0 ) ;
               ret += s->referenced_string( classe, method, red ) ;
            }
            it->destroy() ;
         }
         ret += " )" ;
      }
      else
      {
         ret = symbol_fonc + "()" ;
      }
      
   }
   return ret ;
}



//--------------------------------------------------------------------
std::string
DOC_Function::name( void ) const
//--------------------------------------------------------------------
{
   return myDOC_Symbol ;
}
 


//--------------------------------------------------------------------
DOC_ClassItem const*
DOC_Function::is_class_item( DOC_Class const* a_class ) const
//--------------------------------------------------------------------
{
   string symbol ;
   if( !identificateur )
   {
      symbol = myDOC_Symbol ;
      if( symb!=0 )
      {
         symbol = symb->myDOC_Symbol ;
      }
   }
   else
   {
      if( myDOC_Symbol.substr( 0, 4 ) == "old_" )
      {
         symbol = myDOC_Symbol.substr( 4, myDOC_Symbol.length()-4 ) ;
      }
   }
   DOC_ClassItem * el = a_class->find( symbol ) ;
   
   return el ;
}



 
//--------------------------------------------------------------------
DOC_Symbol*
DOC_Function::arg( size_t n ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Function::arg" ) ;
   PEL_CHECK( mesDOC_Arguments!=0 && (n<mesDOC_Arguments->count() ) ) ;
   DOC_Symbol* result =
      static_cast<DOC_Symbol* >( mesDOC_Arguments->at(n) ) ;
   PEL_CHECK_POST( dynamic_cast<DOC_Symbol* >( result )!=0 ) ;
   return result ;
}


