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

#include <DOC_Writer.hh>

#include <algorithm>
#include <sstream>
#include <functional>
#include <vector>

#include <DOC_WriterHtml.hh>
#include <DOC_WriterJavascript.hh>
#include <DOC_WriterLaTex.hh>

#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <DOC_Function.hh>
#include <DOC_Method.hh>

//#include <unistd.h>
#include <errno.h>
#include <PEL_System.hh>
#include <PEL_ModuleExplorer.hh>

using::std::endl ;

//--------------------------------------------------------------------
DOC_Writer::Format DOC_Writer::format_de_sortie = DocumentationHtml ;
bool DOC_Writer::interactif = false ;
//--------------------------------------------------------------------


//--------------------------------------------------------------------
void
DOC_Writer::parse_explorer( PEL_ModuleExplorer const* exp ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::parse_explorer" ) ;
   
   if( exp->has_entry( "interactif" ) )
   {
      interactif = exp->bool_data( "interactif" ) ;
   }
   std::string const& format = exp->string_data( "format" ) ;
   if( format=="html" )
   {
      format_de_sortie = DOC_Writer::DocumentationHtml ;
   }
   else if( format=="js" )
   {
      format_de_sortie = DOC_Writer::DocumentationJS ;
   }
   else if( format=="latex" )
   {
      format_de_sortie = DOC_Writer::DocumentationLaTex ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "Unkown output format : "+format ) ;
   } 
}



//--------------------------------------------------------------------
DOC_Writer*
DOC_Writer::create( PEL_Object* a_owner, string const& name ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::create" ) ;
   
   DOC_Writer* result = 0 ;
   string chaine = name ;
   switch( format_de_sortie )
   {
      case DocumentationHtml :
         result =  DOC_WriterHtml::create( a_owner, chaine ) ;
         break ;
      case DocumentationJS :
         result =  DOC_WriterJavascript::create( a_owner, chaine ) ;
         break ;
      case DocumentationLaTex :
         result =  DOC_WriterLaTex::create( a_owner, chaine ) ;
         break ;
      default :
         Erreur( "Bad configuration" ) ;
         break ;
   }
   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   return result ;
}



//--------------------------------------------------------------------
DOC_Writer::DOC_Writer( PEL_Object* a_owner, string const& chaine ) 
//--------------------------------------------------------------------
      : PEL_Object( a_owner ),
        leStream( chaine.c_str() )
{
   if( !leStream )
   {
      PEL_Error::object()->raise_plain( "Unable to open "+chaine ) ;
   }
}



//--------------------------------------------------------------------
DOC_Writer::~DOC_Writer( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
std::ostream&
DOC_Writer::out( void ) const
//--------------------------------------------------------------------
{
   return const_cast<std::ofstream&>( leStream ) ;
}



//--------------------------------------------------------------------
string
DOC_Writer:: itoa( int i ) const
//--------------------------------------------------------------------
{
   std::ostringstream os ;
   os << i ;
   
   return( os.str() ) ;
}





//--------------------------------------------------------------------
void
DOC_Writer::process_documentation_package( DOC_Package const* pack,
                              size_t level ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::process_documentation_package" ) ;
   PEL_CHECK_PRE( pack!=0 ) ;
   
   doc_class_from( pack, level+1 ) ;
   PEL_Iterator* it = pack->sub_packages()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Package* p = static_cast<DOC_Package*>( it->item() ) ;
      if( p->own_non_external_classes() )
      {
         this->process_documentation_package( p, level+1 ) ;
      }

   }
   it->destroy() ;
}




//--------------------------------------------------------------------
void
DOC_Writer::doc_class_from( DOC_Package const* pack,
                                         int level )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::doc_class_from" ) ;
   PEL_CHECK_PRE( pack!=0 ) ;

   bool fini = false ;
   PEL_List* liste_classes = DOC_Class::list()->create_clone(this) ;
   while( !fini )
   {
      DOC_Class const* candidate = 0 ;
      PEL_Iterator* it = liste_classes->create_iterator(0) ;
      for( it->start() ; it->is_valid() && candidate==0 ; it->go_next() )

      {
         DOC_Class* cl=static_cast<DOC_Class*>( it->item() ) ;
         if( cl->package() == pack )
         {
            candidate = cl ;
         }
      }
      it->destroy() ;
      fini = candidate==0 ;
      if( !fini )
      {
         while( candidate->mother() !=0 &&
                !candidate->mother()->is_upper_class() &&
                candidate->mother()->package() == pack )
            candidate = candidate->mother() ;
         
         disp_classes_tree( candidate, 0, liste_classes ) ;
      }
   }
}

//--------------------------------------------------------------------
void
DOC_Writer::disp_classes_tree( DOC_Class const* classe,
                                  int niveau,
                                  PEL_List* liste_classes )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::disp_classes_tree" ) ;
   PEL_CHECK_PRE( classe!=0 ) ;
   PEL_CHECK_PRE( liste_classes!=0 ) ;
   PEL_CHECK_PRE( liste_classes->has( classe ) ) ;
   
   size_t idx = liste_classes->index_of( classe ) ;
   PEL_ASSERT( idx<liste_classes->index_limit() ) ;
   liste_classes->remove_at( idx ) ;
   
   if( !classe->is_upper_class() )
   {
      bool fini = false ;
      while( !fini )
      {
         fini = true ;
         PEL_Iterator* it=liste_classes->create_iterator(0) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Class* cl = static_cast<DOC_Class*>( it->item() ) ;
            if( cl->mother() == classe &&
                cl->package()==classe->package() )
            {
               disp_classes_tree( cl, niveau+1, liste_classes ) ;
               fini = false ;
               break ;
            }
         }
         it->destroy() ;
      }
   }
   
} 
   



//--------------------------------------------------------------------
void
DOC_Writer:: copy( string const& source,
                    string const& target ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::copy" ) ;
   
   std::string src = PEL_System::getenv( "PELICANSHOME" ) ;
   if( src.empty() )
   {
      Erreur( "PELICANSHOME environment variable not properly defined !" ) ;
   }
   else
   {
      src += PEL_System::path_name_separator() ;
      src += source ;      
      if( ! PEL_System::copy( src, target ) )
      {
         if( errno==EEXIST )
         {
         }
         else
         {
            Erreur( "Can't execute a copy of " << src << " in " << target ) ;
            Erreur(  "Unknown error refers to errno : " << errno << endl ) ;
         }
         
      }
   }
   
   
}
