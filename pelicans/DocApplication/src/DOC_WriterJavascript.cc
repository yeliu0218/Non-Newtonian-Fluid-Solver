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

#include <DOC_WriterJavascript.hh>
#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <algorithm>
#include <DOC_ClassItem.hh>
#include <DOC_Function.hh>
#include <DOC_Method.hh>

#include <PEL_System.hh>

using std::endl ;

//--------------------------------------------------------------------
stringVector DOC_WriterJavascript::menu = stringVector(0) ;
//--------------------------------------------------------------------



//--------------------------------------------------------------------
DOC_WriterJavascript*
DOC_WriterJavascript:: create( PEL_Object* a_owner, std::string const& name ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterJavascript:: create" ) ;
   
   return new DOC_WriterJavascript( a_owner, name ) ;
}



//--------------------------------------------------------------------
DOC_WriterJavascript::DOC_WriterJavascript( PEL_Object* a_owner,
                                    std::string const& name ) 
//--------------------------------------------------------------------
      : DOC_Writer( a_owner, name + ".html" ),
        file( name )
{
   out() << "<html>" << endl ;
   out() << "<head>" << endl ;
   out() << "<script LANGUAGE=\"JavaScript\" src=\"./js/peldocs.js\"></script>" << endl ;
   out() << "<script LANGUAGE=\"JavaScript\">" << endl ;
   out() << "<!--" << endl ;
   out() << "var parent = 0;" << endl ;
   out() << "var object = 0;" << endl ;
   
}



//--------------------------------------------------------------------
DOC_WriterJavascript::~DOC_WriterJavascript( void ) 
//--------------------------------------------------------------------
{
   out() << "// end -->" << endl ;
   out() << "</script>" << endl ;
   out() << "<HEAD> <title>" << file << "</title></HEAD>" << endl ;
   out() << "<BODY TEXT=\"#000000\" BGCOLOR=\"#FFFFFF\" LINK=\"#0000EE\" VLINK=\"#551A8B\" ALINK=\"#FF0000\">" << endl ;
   out() << "<script LANGUAGE=\"JavaScript\">" << endl ;
   out() << " <!--" << endl ;
   out() << "drawAll(\"class " << file << "\")" << endl ;
   out() << "// end -->" << endl ;
   out() << "</script>" << endl ;
   out() << "</BODY>" << endl ;
   out() << "</HTML>" << endl ;

}



//--------------------------------------------------------------------
string
DOC_WriterJavascript::underlined_text( size_t marqueur,
                               std::string const& text_initial ) const 
//--------------------------------------------------------------------
{
   string ret ;
   
   ret = "<span class=\"" ;
   switch( marqueur )
   {
      case 0 : ret+="uid0" ;
         break ;
      case 1 : ret+="uid1" ;
         break ;
      case 2 : ret+="uid2" ;
         break ;
      case 10 : ret+="uid10" ;
         break ;
      default : ret+="uidd" ;
         break ;
   }
   ret += "\">" + text_initial + "</span></B>" ;
//    if( marqueur == 10 )
//    {
//       ret += "<img SRC=\"http://sand/~chailan/Images/hot-spot.gif\" height=30 width=30>" ;
//    }
      
   return ret ;
}






//--------------------------------------------------------------------
string 
DOC_WriterJavascript::reference( DOC_Class const* cl,
                                 std::string const& label,
                                 std::string const& display ) const
//--------------------------------------------------------------------
{
   string ret = display ;
   if( !cl->is_external() )
   {
      ret = "<A HREF=\"" + cl->name() + ".html#" + label + "\">" + display + "</A>" ;
   }
   return ret ;
}



// //--------------------------------------------------------------------
// void
// DOC_WriterJavascript::ecrit_code( string const& com )
// //--------------------------------------------------------------------
// {
//    out() << "<TT><B>" ;
//    for( string::const_iterator iter = com.begin() ;
//         iter != com.end() ;
//         ++iter )
//    {
//       if( *iter=='\n' )
//       {
//          out() << "<BR>\n" ;
//       }
//       else if( *iter==' ' && iter+1 != com.end() && *(iter+1)==' ' )
//       {
//          out() << "&nbsp; " ;
//          ++iter ;
//       }
//       else
//       {
//           out() << *iter ;
//       }
//    }
//    out() << "</B></TT>" << endl ;
// }



//--------------------------------------------------------------------
void
DOC_WriterJavascript::ecrit_comment( std::string const& com )
//--------------------------------------------------------------------
{
   string deb = "object.add(\'" ;
   string fin = " \')\n" ;
   string mot ;
   
   string::const_iterator iter = com.begin() ;
   for(  ;
         iter != com.end() ;
         ++iter )
   {
      if( *iter!='\n' && *iter!=' ' ) break ;
   }
   
   string::const_iterator iter_fin = com.end()-1 ;
   for(  ;
         iter_fin > iter ;
         --iter_fin )
   {
      if( *iter_fin!='\n' && *iter_fin!=' ' ) break ;
   }
   bool first = true ;
   
   if( iter < iter_fin )
   {
      out() << "object = parent.add(new comment())" << endl ;
      out() << deb ;
      for( ;
           iter <= iter_fin ;
           ++iter )
      {
         if( *iter=='.' && first )
         {
            out() << '.' << fin << deb ;
            first = false ;
         }
         else if( *iter=='\n' )
         {
            if( !first )
            {
               out() << fin << deb ;
            }
            else
            {
               out() << "\\n" ;
            }
            
         }
         else if( *iter=='\'' )
         {
            out() << " " ;
         }
         else
         {
            out() << *iter ;
         }
      }
      out() << fin ;
   }
}

   
//--------------------------------------------------------------------
void
DOC_WriterJavascript::ecrit( std::string const& com )
//--------------------------------------------------------------------
{
   for( string::const_iterator iter = com.begin() ;
        iter != com.end() ;
        ++iter )
   {
      if( *iter=='\n' )
      {
      }
      else if( *iter=='\'' )
      {
         out() << "\\'" ;
      }
      else
      {
         out() << *iter ;
      }
   }
}

   
//--------------------------------------------------------------------
void
DOC_WriterJavascript::finalize( void )
//--------------------------------------------------------------------
{
   std::string src = "tools" ;
   src += PEL_System::path_name_separator() ;
   src += "peldoc" ;
   src += PEL_System::path_name_separator() ;
   src += "js" ;
   std::string dest = "." ;
   dest += PEL_System::path_name_separator() ;
   dest += "js" ;
   copy( src, dest ) ; 
//   copy( "tools/peldoc/js",
//           "./js" ) ;
   
   process_documentation_package( DOC_Package::master() ) ;
}



//--------------------------------------------------------------------
bool
DOC_WriterJavascript::source_reference( void ) const
//--------------------------------------------------------------------
{
   return true ;
}



//--------------------------------------------------------------------
void
DOC_WriterJavascript::process_documentation_package( DOC_Package const* pack,
                                   size_t level ) 
//--------------------------------------------------------------------
{   
   if( level==0 )
   {
      out() << "var level = new Array()" << endl ;
      out() << "level[0] = db" << endl ;
   }
   else
   {
      out() << "level["<<level<<"] = level["<<level-1
            <<"].add(new Entry(\"\",\"" << pack->comment() << "\".fontcolor('black')))" << endl ;  
   }
   package_level = level ;
   DOC_Writer::process_documentation_package( pack, level ) ;
   
}



//--------------------------------------------------------------------
void
DOC_WriterJavascript::disp_classes_tree( DOC_Class const* classe,
                                     int niveau,
                                     PEL_List* liste_classes )
//--------------------------------------------------------------------
{
   size_t level = niveau+package_level+1 ;
   string color = "blue" ;
   if( classe->is_plug_point() ) color = "red" ;
   
   out() << "level["<<level<<"] = level["<<level-1
         <<"].add(new Entry(\"\",\"" << classe->name()
         << "\".fontcolor('" << color << "')"
         << ".link('" << classe->name() << ".html') ))" << endl ;     
   
   DOC_Writer::disp_classes_tree( classe,
                                     niveau,
                                     liste_classes ) ;
} 




//--------------------------------------------------------------------
void
DOC_WriterJavascript::process_documentation( DOC_Class const* classe )
//--------------------------------------------------------------------
{
   if( classe->mother()!=0 ) 
   {
      out() << " parent = db.add(new DOC_Text( 'Ininheritd from " 
            << reference( classe->mother(),
                          "",
                          classe->mother()->name() )
            << "'))" << endl ;
   }
   if( classe->is_virtual() &&
       classe->inherited_classes()->count()>0 )
   {
      out() << " parent = db.add(new DOC_Text('Derived classe(s) : " ;
      PEL_List const* inherited_classes = classe->inherited_classes() ;
      
      string heritage ;
      PEL_Iterator* it = inherited_classes->create_iterator(0) ;
      bool prem=true ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Class*cl = static_cast<DOC_Class*>( it->item() ) ;
         if( !prem ) out() << " , " ; else prem=false ;
         out() << reference( cl,
                             "",
                             cl->name() ) ;
      }
      out() << "'))" <<  endl ;
      it->destroy() ;
   }
   out() << "parent = db.add(new Description())" << endl ;
   ecrit_comment(
      DOC_Method::format_comment( *this,
                                    1,
                                    classe->comment(),
                                    '`',
                                    '\'',
                                    classe ) )  ;
   PEL_List const* elements = classe->owned_components() ;
   
   if( elements->count()>0 )
   {
      string category = "no init cat" ;
      bool premierDomaine = true ;
      
      string heritage ;
      PEL_Iterator* it = elements->create_iterator(0) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_ClassItem*elem = static_cast<DOC_ClassItem*>( it->item() ) ;
         if( classe->has_to_document( elem ) )
         {
            bool nouvelleDOC_Category = elem->category_displaye() != category ;
            if( nouvelleDOC_Category )
            {
               premierDomaine = false ;
               category = elem->category_displaye() ;
               out() << "parent = db.add(new Section('" ;
               if( category.empty() )
               {
                  out() << "Miscellaneous" ;
               }
               else
               {
                  out() << DOC_Method::format_comment( *this,
                                                         1,
                                                         category,
                                                         '`',
                                                         '\'',
                                                         elem->def_class() ) ;
               }
               out() << "'))" << endl ;
            }
            process_documentation( classe, elem ) ;            
         }
      }
      it->destroy() ;
   }  
}



//--------------------------------------------------------------------
void
DOC_WriterJavascript::process_documentation( DOC_Class const* classeAProcess_Documentationr,
                           DOC_ClassItem * elem )  
//--------------------------------------------------------------------
{
   string proto = elem->prototype( *this ) ;
   for( string::iterator iter = proto.begin() ;
        iter!=proto.end() ;
        ++iter )
      if( *iter=='\n' ) *iter = ' ' ;
   
   
   out() << "parent = db.add(new Entry('"
         << elem->name() << "','" 
         << proto << "'))" << endl ;
   ecrit_comment(
      DOC_Method::format_comment( *this,
                                    1,
                                    elem->comment(),
                                    '`',
                                    '\'',
                                    elem->def_class() ) ) ;
   DOC_Method const* method = elem->method() ;
   if( method != 0 )
   {
      PEL_List const* preConditions =  method->pre_conditions() ;
      PEL_List const* postConditions =  method->post_conditions() ;
      if( preConditions->count()>0 )
      {
         string title = "Preconditions" ;
         if( method->is_precondition() ) title+= " implementation for method " +
                                         method->short_name();
         out() << "object = parent.add(new List(\""
               << title << "\"))" << endl ;
         
         PEL_Iterator* it = preConditions->create_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
            out() << "object.add('" ;
            ecrit( f->referenced_string(
               method->def_class(),
               method,
               this ) ) ;
            out()   << "')" << endl ;
         }
         it->destroy() ;
      }
      if( postConditions->count()>0 )
      {
         string title = "Postconditions" ;
         if( method->is_postcondition() ) title+= " implementation for method " +
                                          method->short_name();
         out() << "object = parent.add(new List(\""
               << title << "\"))" << endl ;
         PEL_Iterator* it = postConditions->create_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
            out() << "object.add('" ;
            ecrit( f->referenced_string(
               method->def_class(),
               method,
               this ) ) ;    
            out() << "')" << endl ;
         }
         it->destroy() ;
      }
      DOC_Method const* surdefinition = method->overiden_method() ;
   
      if( classeAProcess_Documentationr!=method->def_class() )
      {
         out() << "object = parent.add(new DOC_Text('"
               << " Extension from "
               << reference( method->def_class(),
                             "",
                             method->def_class()->name() )
               << "'))" << endl ;
      } else if( surdefinition!=0 )
      {
         out() << "object = parent.add(new DOC_Text('"
               << " Overridden from "
               << reference( surdefinition->def_class(),
                             "",
                             surdefinition->def_class()->name() )
               << "'))" << endl ;
      }
   }
}

