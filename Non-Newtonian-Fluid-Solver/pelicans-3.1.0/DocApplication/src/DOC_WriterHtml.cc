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

#include <DOC_WriterHtml.hh>
#ifdef OUTLINE
#define inline
#include <DOC_WriterHtml.icc>
#undef inline
#endif

#include <PEL.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <algorithm>
#include <DOC_ClassItem.hh>
#include <DOC_Function.hh>
#include <DOC_Method.hh>
#include <DOC_WebView.hh>
#include <fstream>

using std::endl ;

//--------------------------------------------------------------------
stringVector DOC_WriterHtml::menu = stringVector(0) ;
//--------------------------------------------------------------------



//--------------------------------------------------------------------
DOC_WriterHtml*
DOC_WriterHtml::create( PEL_Object* a_owner, std::string const& str )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::create" ) ;
   
   return new DOC_WriterHtml( a_owner, str ) ;
}



//--------------------------------------------------------------------
DOC_WriterHtml::DOC_WriterHtml( PEL_Object* a_owner,
std::string const& str )
//--------------------------------------------------------------------
: DOC_Writer( a_owner, str + ".html" )
, file( str ) {
   PEL_LABEL( "DOC_WriterHtml::DOC_WriterHtml" ) ;
   std::string strtitle = file ;
   if( file == "index" ) strtitle = DOC_Package::application() + " class reference" ;
   head( out(), strtitle, file !="index" ) ; /** <HTML /> */
   
   /** <MODIF> - Create frameset */
   if( file == "index" ) {
      DOC_Package const* pack = DOC_Package::master() ;
      PEL_Iterator* it = pack->sub_packages()->create_iterator(this) ;
      string def = "MASTER-tree" ;
      if( it->is_valid() ) {
         DOC_Package const* first = dynamic_cast<DOC_Package*>( it->item() ) ;
         PEL_CHECK( dynamic_cast<DOC_Package const*>( first )!=0 ) ;
         def = first->name() ;
      }
      /** <HTML> */
      out() << "<frameset cols=\"20%,*\" frameborder=\"1\">" << endl ;
      out() << "  <frameset rows=\"30%,*\">" << endl ;
      out() << "    <frame src=\"overview-frame.html\" name=\"LibListFrame\">" << endl ;
      out() << "    <frame src=\"" << def << ".html\" name=\"LibFrame\">" << endl ;
      out() << "  </frameset>" << endl ;
      out() << "  <frame src=\"" + DOC_Package::master()->name() + "-tree.html\" name=\"classFrame\">" << endl ;
      out() << "  <noframes>" << endl ;
      out() << "  <body>" << endl ;
      out() << "    <p>Sorry, a frames-capable browser is required to access the PELICANS documentation.</p>" << endl ;
      out() << "    <a href=\"" << def << ".html\">.: No frames :.</a>" << endl ;
      out() << "  </body>" << endl ;
      out() << "  </noframes>" << endl ;
      out() << "</frameset>" << endl ;
      /** </HTML> */
   }
   /** </MODIF> - Create frameset */
   package_level = 0 ;
   class_level = 0 ;
}



//--------------------------------------------------------------------
DOC_WriterHtml::~DOC_WriterHtml( void )
//--------------------------------------------------------------------
{
   if( file != "index" ) {
      out() << "\n</body>" << endl ; /** <HTML /> */
   }
   out() << "</html>" << endl ; /** </HTML /> */
   
}



//--------------------------------------------------------------------
string
DOC_WriterHtml::underlined_text( size_t marqueur,
std::string const& text_initial ) const
//--------------------------------------------------------------------
{
   string ret ;
   /** <MODIF> - Highlighting styles definition */
   /** <HTML> */
   ret = "<span class=\"" ;
   switch( marqueur ) {
      /** style =  */
      case 0 : ret+="uid0" ;
      break ;
      /** style =  */
      case 1 : ret+="uid1" ;
      break ;
      /** style =  */
      case 2 : ret+="uid2" ;
      break ;
      /** style =  */
      case 10 : ret+="uid10" ;
      break ;
      /** style =  */
      case 20 : ret+="uid20" ;
      break ;
      /** style = default */
      default : ret+="uidd" ;
      break ;
   }
   ret += "\">" + text_initial + "</span>" ;
   /** </HTML> */
   /** </MODIF> - Highlighting styles definition */
   return ret ;
}

//--------------------------------------------------------------------
void
DOC_WriterHtml::hruler( void )
//--------------------------------------------------------------------
{
   out() << "<hr />" << endl ;
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::start_category( std::string const& name,
std::string const& protection )
//--------------------------------------------------------------------
{
   /** <MODIF> - category table heading */
   /** <HTML> */
   out() << "\n<!-- ========== NEW CATEGORY ========== -->" << endl ;
   out() << "<table class=category cellpadding=2 cellspacing=0 rules=rows>" << endl ;
   out() << "  <tr class=heading><th class=name>" << name << "</th><th class=protection>" << protection << "</th>" << "</tr>" << endl ;
   /** </HTML> */
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::title( std::string const& tit, int niveau )
//--------------------------------------------------------------------
{
   if( tit.length() > 0 ) {
      out() << "<h1 id=maintitle>" << tit << "</h1>" << endl;
   }
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::comment( std::string const& com )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::comment" ) ;
   
   if( com.length()>0 ) {
      // Searching for html syntax
      std::string to_output = expand( com, "http://" ) ;
      to_output = expand( to_output, "mailto:" ) ;
      out() << "<pre>" << to_output << "</pre>" << endl ;
   }
   else {
      out() << "<br />" << endl ;
   }
}


//--------------------------------------------------------------------
std::string
DOC_WriterHtml::expand( std::string const& com, std::string const& html_cmd )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml:: expand" ) ;
   std::string result ;
   size_t idx_deb = 0 ;
   size_t idx ;
   size_t l = com.length() ;
   
   while( ( idx = com.find( html_cmd, idx_deb ) ) < l ) {
      result += com.substr( idx_deb, idx-idx_deb ) ;
      
      size_t idx_fin = idx+1 ;
      for( ; idx_fin<l ; idx_fin++ ) {
         if( com[idx_fin]==' ' || com[idx_fin]=='\n' )
            break ;
      }
      std::string link = com.substr( idx+html_cmd.length(), idx_fin-(idx+html_cmd.length()) ) ;
      
      result += "<a href=\"" + html_cmd + link + "\">" + link + "</a>" ; /** <HTML /> */
      idx_deb = idx_fin ;
   }
   result += com.substr( idx_deb, com.length()-idx_deb ) ;
   
   return result ;
   
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::list( std::string const& id )
//--------------------------------------------------------------------
{
   if( !id.empty() ) {
      out() << id << endl ;
      out() << "<ul>" << endl ;
   }
   else {
      out() << "</ul>" << endl ;
   }
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::new_element( void )
//--------------------------------------------------------------------
{
   write( "<li>" ) ;
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::write_label( std::string const& label )
//--------------------------------------------------------------------
{
   out() << "<a name=\"" << label << "\"></a>" ;
}



//--------------------------------------------------------------------
string
DOC_WriterHtml::reference( DOC_Class const* cl,
std::string const& label,
std::string const& display ) const
//--------------------------------------------------------------------
{
   string ret = display ;
   string ext ;
   if( !cl->is_external() ) {
      ret = "<a href=\"" + cl->name() + ".html#" + label
      + "\" target=\"classFrame\">" + display + "</a>" ;
   }
   else if( ! ( ext = DOC_Tools::external( cl->name() + ".html" ) ).empty() ) {
      ret = "<a href=\"" + ext + "#" + label + "\" target=\"_blank\">"
      + underlined_text( 2, display ) + "</a>" ;
   }
   
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::write_code( std::string const& com )
//--------------------------------------------------------------------
{
   out() << "<code>" ;
   for( string::const_iterator iter = com.begin() ;
   iter != com.end() ;
   ++iter ) {
      /** ecrire de manière propre le traitement des caractères spéciaux */
      if( *iter=='\n' ) {
         out() << "<br />\n" ;
      }
      else if( *iter==' ' && iter+1 != com.end() && *(iter+1)==' ' ) {
         out() << "&nbsp;" ;
         ++iter ;
      }
      else {
         out() << *iter ;
      }
   }
   out() << "</code>" ;
}

//--------------------------------------------------------------------
void
DOC_WriterHtml::write_prototype( std::string const& com )
//--------------------------------------------------------------------
{
   out() << "<table class=prototype>\n  <tr><td><code>" ;
   for( string::const_iterator iter = com.begin() ;
   iter != com.end() ;
   ++iter ) {
      if( *iter=='\n' ) {
         out() << "</tr>\n  <tr><td /><td><code>" ;
      }
      else if( *iter == 40 ) {
         out() << *iter << "</code></td><td><code>" ;
         ++iter ;
      }
      else if( *iter==' ' && iter+1 != com.end() && *(iter+1)==' ' ) {
         ++iter ;
      }
      /** réutiliser le traitement de caracteres speciaux de write_code */
      else {
         out() << *iter ;
      }
   }
   out() << "</code></td></tr>\n</table>" << endl ;
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::write( std::string const& com )
//--------------------------------------------------------------------
{
   
   for( string::const_iterator iter = com.begin() ;
   iter != com.end() ;
   ++iter ) {
      if( *iter=='\n' ) {
         out() << "<br />\n" ;
      }
      else {
         out() << *iter ;
      }
   }
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::head( std::ostream& ovf, std::string const& title, bool set_body )
//--------------------------------------------------------------------
{
   std::string nm = title ;
   size_t idx=nm.find_last_of( "/" ) ;
   if( idx<nm.length() )
      nm = nm.substr( idx+1, nm.length()-idx-1 ) ;
   /** <HTML> */
   ovf << "<!DOCTYPE html" << endl ;
   ovf << "     PUBLIC \"-//W3C//DTD XHTML 1.0 Frameset//EN\"" << endl ;
   ovf << "     \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-frameset.dtd\">" << endl ;
   ovf << "<html>" << endl ;
   ovf << "<head>" << endl ;
   ovf << "  <title>" << nm << "</title>" << endl ;
   ovf << "  <link rel =\"stylesheet\" type=\"text/css\" href=\"stylesheet.css\" title=\"Style\">" << endl ;
   ovf << "</head>" << endl ;
   if( set_body ) {
      ovf << "<script>" << endl ;
      ovf << "  function asd()" << endl ;
      ovf << "  {" << endl ;
      ovf << "    parent.document.title=\"" <<
      nm << "_" << DOC_Package::application() << " \";" << endl ;
      ovf << "  }" << endl ;
      ovf << "</script> " << endl ;
      ovf << "\n<body onload=\"asd();\">" << endl ;
   }
   /** <HTML> */
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::finalize( void )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::finalize" ) ;
   
   std::string src = "doc" ;
   src += PEL_System::path_name_separator() ; 
   src += "share" ;
   src += PEL_System::path_name_separator() ; 
   src += "stylesheet.css" ;
   copy( src, "." ) ;
//   copy( "doc/share/stylesheet.css", "." ) ;
   src = "doc" ;
   src += PEL_System::path_name_separator() ; 
   src += "share" ;
   src += PEL_System::path_name_separator() ;
   src += "img" ;
   src += PEL_System::path_name_separator() ;
   src += "arrow.png" ;
   copy( src, "." ) ;
//   copy( "doc/share/img/arrow.png", "." ) ;
   
   // process_documentation_package( DOC_Package::master() ) ;
   
   DOC_WriterHtml* red = new DOC_WriterHtml( 0, DOC_Package::master()->name()+"-tree" ) ;
   red->start_package_doc( DOC_Package::master(), 0 ) ;
   red->title( "Hierarchy for "+DOC_Package::application()+" libraries" , 0 ) ;
   red->process_documentation_package( DOC_Package::master() ) ;
   red->destroy() ;
   
   //??????????????? numero de version ???????????????????????,
   std::ofstream ovf( "overview-frame.html" ) ;
   head( ovf, DOC_Package::application() + " API Specification : Overview", false  ) ;
   ovf << "<body>" << endl ;
   ovf << "<h3>"<<DOC_Package::application() <<"</h3>" << endl ;
   ovf << "<h4><a href=\""+DOC_Package::master()->name()+"-tree.html\" target=\"classFrame\">All Classes</a></h4>" << endl ;
   ovf << "<h4>Libraries</h4>" << endl ;
   ovf << "<div id=\"liblist\">" << endl ;
   ovf << "<ul>" << endl ;
   
   PEL_List const* classes = DOC_Class::list() ;
   DOC_Package const* pack = DOC_Package::master() ;
   PEL_Iterator* it = pack->sub_packages()->create_iterator(this) ;
   for( it->start() ; it->is_valid() ; it->go_next() ) {
      DOC_Package* p = static_cast<DOC_Package*>( it->item() ) ;
      PEL_CHECK( dynamic_cast<DOC_Package*>(p)!=0 ) ;
      
      if( ! p->own_non_external_classes() ) continue ;
      
      string nn = p->name() ;
      
      ovf << "<li><a href=\"" << nn << ".html\""
      << " target=\"LibFrame\">" << p->comment()  << "</a></li>" << endl ;
      
      red =  DOC_WriterHtml::create( 0, nn ) ;
      red->out() << "<h3><a href=\"" << p->name() << "-tree.html\" target=\"classFrame\">" << p->comment() << "</a></h3>" << endl ;
      red->out() << "<h4>Classes</h4>" << endl ;
      red->out() << "<div id=\"classlist\">" << endl ;
      
      PEL_Iterator* cl_it = classes->create_iterator( this) ;
      red->out() << "<ul>" << endl ;
      for( cl_it->start() ; cl_it->is_valid() ; cl_it->go_next() ) {
         DOC_Class* cl = static_cast<DOC_Class* >( cl_it->item() ) ;
         if( p->own_class( cl  ) && !cl->is_external() ) {
            string display = cl->name() ;
            if( cl->is_plug_point() )
               display = underlined_text(10, cl->name()) ;
            red->out() << "<li>" << reference( cl, "", display ) << "</li>" << endl ;
         }
      }
      red->out() << "</ul>\n</div>" << endl ;
      red->destroy() ;
      
      nn = p->name()+"-tree" ;
      red = DOC_WriterHtml::create( 0, nn ) ;
      red->start_package_doc( p, 0 ) ;
      red->title( "Hierarchy for Library "+p->comment(), 0 ) ;
      red->process_documentation_package( p, 0 ) ;
      red->destroy() ;
   }
   
   ovf << "</ul>\n</div>" << endl ;
   ovf << "</body>\n</html>" << endl ;
}



//--------------------------------------------------------------------
bool
DOC_WriterHtml::source_reference( void ) const
//--------------------------------------------------------------------
{
   return true ;
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::start_package_doc( DOC_Package const* pack,
size_t level )
//--------------------------------------------------------------------
{
   stringVector items(1), refs(1) ;
   items(0) = "Tree" ;
   refs(0) = pack->name()+"-tree.html" ;
   navbar( out(), items, refs, 0 ) ;
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::process_documentation_package( DOC_Package const* pack,
size_t level )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::process_documentation_package" ) ;
   
   string indent = string( 2*level, ' ' ) ;
   /** ??? -> moche ??? */
   if(class_level == 1) out() << indent << "</ul>" << endl ;
   if( level > 0 ) {
      class_level = 0 ;
      out() << indent << "\n<!-- ========== NEW PACKAGE ========== -->" << endl ;
      out() << indent << "<h2 class=name>" << pack->comment() << "</h2>" << endl ;
   }
   
   out() << indent << "<ul class=packagetree>" << endl ;
   package_level = level ;
   DOC_Writer::process_documentation_package( pack, level ) ;
   out() << indent << "</ul>" << endl ;
}


//--------------------------------------------------------------------
void
DOC_WriterHtml::disp_classes_tree( DOC_Class const* classe,
int niveau,
PEL_List* list_classes )
//--------------------------------------------------------------------
{
   
   if( !classe->is_external() ) {
      string indent = string( 2*(niveau+package_level), ' ' ) ;
      string style ;
      /** select list element style depending on location in the class tree */
      if ( niveau > 0 ) {
         if ( classe->inherited_classes()->count() == 0 ) {
            style = "leaf" ;
         }
         else {
            style = "node" ;
         }
      }
      else {
         if ( classe->inherited_classes()->count() == 0 ) {
            style = "nochild" ;
         }
         else {
            style = "root" ;
         }
      }
      if ( class_level+1 == (size_t) niveau ) out() << indent << "<ul class=branch>" << endl ;
      /** write matching end-tags for open lists */
      if ( class_level > (size_t) niveau ) {
         for (size_t i = class_level - (size_t) niveau; i > 0; i--) {
            out() << indent << "  </ul>" << endl ;
         }
      }
      string display = classe->name() ;
      if( classe->is_plug_point() ) display = underlined_text(10, classe->name()) ;
      string comm ;
      string comm_court = classe->short_description() ;
      if( !comm_court.empty() ) {
         int st = 30 - classe->name().length() ;
         if( st<0 ) st=2 ;
         comm = string( st, '.' ) + comm_court ;
      }
      /** write parent class name if applicable */
      if ( (niveau == 0) && 
           (classe->mother() != NULL) && 
           (!classe->mother()->is_upper_class() ) ) 
      {
         out() << indent
         << "  <li class=" << style << ">"
         << reference( classe, "", display ) << comm
         << "<span class=\"parentname\"> (" << reference( classe->mother(), "", classe->mother()->name() ) << ")</span>"
         << "</li>" << endl ;
      }
      else {
         out() << indent << "  <li class=" << style << ">" << reference( classe, "", display ) << comm << "</li>" << endl ;
      }
      class_level = niveau ;
   }
   DOC_Writer::disp_classes_tree( classe,
   niveau,
   list_classes ) ;
}




//--------------------------------------------------------------------
void
DOC_WriterHtml::process_webview( std::string const src,
std::string const dest,
stringVector const& item,
stringVector const& refs,
size_t n )
//--------------------------------------------------------------------
{
   static DOC_WebView* webview = 0 ;
   if( webview==0 ) {
      PEL_Module* mod = PEL_Module::create( PEL_Root::object(), "PEL_Application" ) ;
      mod->add_entry( "concrete_name", PEL_String::create( mod, "webview" ) ) ;
      PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( mod, mod ) ;
      webview = DOC_WebView::create( PEL_Root::object(), exp ) ;
   }
   std::ofstream os( dest.c_str() ) ;
   if( os ) {
      head( os, src, true ) ;
      navbar( os, item, refs, n) ;
      os.close() ;
      stringVector arg( 5 ) ;
      arg(0) = "-g" ;
      arg(1) = "--no-header" ;
      arg(2) = "-o" ;
      arg(3) = dest ;
      arg(4) = src ;
      webview->re_initialize( arg ) ;
      webview->run() ;
      os.open( dest.c_str() , std::ios_base::app ) ;
      os<<"</body>\n</html>"<<endl ;
      os.close() ;
   }
   
}

//--------------------------------------------------------------------
void
DOC_WriterHtml::navbar( std::ostream& os,
stringVector const& items,
stringVector const& refs,
size_t n )
//--------------------------------------------------------------------
{
   PEL_ASSERT( n<items.size() ) ;
   PEL_ASSERT( refs.size()==items.size() ) ;
   
   /** <MODIF> - NAVIGATION BAR */
   os << "\n<!-- ========== START OF NAVBAR ========== -->" << endl ;
   os << "<a name=\"navbar_top\"><!-- --></a>" << endl ;
   os << "<div id=navbar>" << endl ;
   os << "  <a name=\"navbar_top_firstrow\"><!-- --></a>" << endl ;
   os << "  <div id=navbarapp>"<<DOC_Package::application() <<"</div>" << endl ;
   os << "  <div id=navbarmenu>" << endl ;
   for( size_t i=0 ; i<items.size() ; i++ ) {
      if( i==n ) {
         os << "    <span class=selected>" << items(i) << "</span>" << endl ;
      }
      else {
         os << "    <a href=\"" << refs(i) << "\"><span>" << items(i) << "</span></a>" << endl ;
      }
   }
   os << "    <div id=navbarsub>" << endl ;
   os << "      <a href=\"index.html\" target=\"_top\">FRAMES :</a>" << endl ;
   os << "      <a href=\"" << refs(0) << "\" target=\"_top\">: NO FRAMES</a>" << endl ;
   os << "    </div>" << endl ;
   os << "  </div>" << endl ;
   os << "</div>" << endl ;
   os << "<!-- ========== END OF NAVBAR ========== -->" << endl << endl ;
   /** </MODIF> - NAVIGATION BAR */
   
}

//--------------------------------------------------------------------
void
DOC_WriterHtml::process_documentation( DOC_Class const* classe )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::process_documentation" ) ;

   if( classe->package() == 0 ) return ;
   // ----------------------------------

   string pname = classe->package()->top()->name() ;
   
   PEL_List const* elements = classe->owned_components() ;
   std::string inlined ;
   std::string implementation ;
   std::string interface ;
   std::string implementation_dest ;
   std::string interface_dest ;
   std::string inlined_dest ;
   if( classe->publish_source() ) {
      implementation = DOC_Tools::read_file( classe->name()+".cc" ) ;
      interface = DOC_Tools::read_file( classe->name()+".hh" ) ;
      inlined = DOC_Tools::read_file( classe->name()+".icc" ) ;
      if( !implementation.empty() )
         implementation_dest = PEL_Module::basename( implementation )+".html" ;
      if( !interface.empty() )
         interface_dest = PEL_Module::basename( interface )+".html" ;
      if( !inlined.empty() )
         inlined_dest = PEL_Module::basename( inlined )+".html" ;
   }
   
   stringVector item(2), refs(2) ;
   item(0) = "Tree" ; refs(0) = pname + "-tree.html" ;
   item(1) = "Class" ; refs(1) = classe->name() + ".html" ;
   if( !interface_dest.empty() ) {
      item.append( "Header" ) ; refs.append( interface_dest ) ;
   }
   if( !implementation_dest.empty() ) {
      item.append( "Implementation" ) ; refs.append( implementation_dest ) ;
   }
   if( !inlined_dest.empty() ) {
      item.append( "Inlined" ) ; refs.append( inlined_dest ) ;
   }
   navbar( out(), item, refs, 1 ) ;
   if( !interface_dest.empty() ) {
      process_webview( interface, interface_dest, item, refs, 2 ) ;
   }
   if( !implementation_dest.empty() ) {
      process_webview( implementation, implementation_dest, item, refs, 3 ) ;
   }
   if( !inlined_dest.empty() ) {
      process_webview( inlined, inlined_dest, item, refs, 4 ) ;
   }
   
   /** <MODIF> - summary link */
   out() << "<br />\n<div class=link2anchor>SUMMARY :: <a href=\"#method_summary\">METHOD</a></div>\n<br />\n" << endl ;
   out() << "\n<!-- ====== START OF CLASS DATA ======== -->" << endl ;
   
   /** <MODIF> - title */
   string fullpp ;
   DOC_Package const* pp = classe->package() ;
   while( pp != DOC_Package::master() ) {
      fullpp = pp->comment() + " / " + fullpp ;
      pp = pp->father() ;
   }
   out() << "<div id=classtitle>" << endl ;
   out() << "  <span class=classpath>" << fullpp << "</span><br />" << endl ;
   out() << "  <span class=classname>Class  " << classe->name() << "</span>" << endl ;
   out() << "</div>" << endl ;
   out() << "\n<br />" << endl ;
   
   /** <MODIF> - HIERARCHY */
   PEL_List* mms=PEL_List::create(this) ;
   DOC_Class const* mm = classe->mother() ;
   while( mm != 0 ) {
      mms->prepend( const_cast<DOC_Class*>(mm) ) ;
      mm = mm->mother() ;
   }
   size_t i=0 ;
   PEL_ListIterator* class_it = mms->create_iterator(0) ;
   out() << "\n<div id=hierarchy>" << endl ;
   for( class_it->start() ; class_it->is_valid() ; class_it->go_next() ) {
      DOC_Class*cl = static_cast<DOC_Class*>( class_it->item() ) ;
      if (i == 0) {
         out() << "<ul class=rootclass><li>" << reference( cl, "", cl->name() ) << "</li>" << endl ;
      }
      else {
         out() << "  <ul class=parentclass><li>" << reference( cl, "", cl->name() ) << "</li>" << endl ;
      }
      ++i ;
   }
   class_it->destroy() ;
   if( i!=0 ) {
      
      out() << "  <ul class=self><li>" << classe->name() << "</li>" << endl ;
   }
   while (i != 0) {
      out() << "</ul>" ;
      --i ;
   }
   out() << "\n</div>" << endl ;
   
   if( classe->is_virtual() && classe->inherited_classes()->count()>0 ) {
      out() << "\n<!-- ========== DERIVED CLASSES ========== -->" << endl ;
      out() << "<dl><dt><b>Direct Known Derived Classes:</b>\n<dd>" ;
      PEL_List const* inherited_classes = classe->inherited_classes() ;
      PEL_Iterator* child_it = inherited_classes->create_iterator(0) ;
      bool prem = true ;
      for( child_it->start() ; child_it->is_valid() ; child_it->go_next() ) {
         DOC_Class*cl = static_cast<DOC_Class*>( child_it->item() ) ;
         if( !prem ) write( ", " ); else prem=false ;
         
         write( reference( cl,
         "",
         cl->name() ) ) ;
      }
      child_it->destroy() ;
      out() << "</dd>\n</dl>" << endl ;
   }
   hruler() ;
   
   comment(
   DOC_Method::format_comment( *this,
   1,
   classe->comment(),
   '`',
   '\'',
   classe ) ) ;
   hruler() ;
   out() << "<a name=\"method_summary\"></a>" << endl ;
   
   if( elements->count()>0 ) {
      string category = "no init cat" ;
      bool premierDomaine = true ;
      
      string heritage ;
      PEL_Iterator* elm_it = elements->create_iterator(0) ;
      for( elm_it->start() ; elm_it->is_valid() ; elm_it->go_next() ) {
         DOC_ClassItem* el = static_cast<DOC_ClassItem*>( elm_it->item() ) ;
         
         if( classe->has_to_document( el ) ) {
            bool nouvelleDOC_Category = el->category()->name() != category ;
            if( nouvelleDOC_Category ) {
               if( !premierDomaine ) {
                  out() << "</table>" << endl ;
               }
               premierDomaine = false ;
               category = el->category()->name() ;
               string prot = "public" ;
               if( el->protection() == DOC_ClassItem::Protected ) {
                  prot = "protected" ;
               }
               else if( el->protection() == DOC_ClassItem::Private ) {
                  prot = "private" ;
               }
               string cat_title = category ;
               if( cat_title.empty() ) cat_title="Miscellaneous" ;
               start_category( cat_title, prot ) ;
            }
            out() << "  <tr><td colspan=\"2\">" ;
            write( reference( classe, el->bookmark(), el->signature() ) ) ;
            out() << "</td></tr>" << endl;
         }
      }
      out() << "</table>" << endl ;
      
      for( elm_it->start() ; elm_it->is_valid() ; elm_it->go_next() ) {
         DOC_ClassItem* el = static_cast<DOC_ClassItem*>( elm_it->item() ) ;
         if( classe->has_to_document( el ) ) {
            out() << "\n<!-- ========== NEW METHOD ========== -->" << endl ;
            hruler() ;
            write_label( el->bookmark() ) ;
            process_documentation( classe, el ) ;
         }
      }
      elm_it->destroy() ;
   }
}



//--------------------------------------------------------------------
void
DOC_WriterHtml::process_documentation( DOC_Class const* classeAProcess_Documentationr,
DOC_ClassItem * elem )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterHtml::process_documentation II" ) ;
   
   string proto = elem->prototype( *this ) ;
   out() << "\n<div class=method>" << endl ;
   write_prototype( proto ) ;
   comment(
   DOC_Method::format_comment( *this,
   1,
   elem->comment(),
   '`',
   '\'',
   elem->def_class() ) ) ;
   out() << "<br />" << endl ;
   DOC_Method const* method = elem->method() ;
   if( method != 0 ) {
      PEL_List const* preConditions =  method->conditions_list(DOC_Method::pre) ;
      PEL_List const* postConditions =  method->conditions_list(DOC_Method::post) ;
      if ( ( preConditions->count()>0 )||( postConditions->count()>0 ) ) {
         DOC_Method const* referring = method->referring_method() ;
         
         out() << "<div class=conditions>" << endl ;
         if( preConditions->count() > 0 ) {
            if( method->is_precondition() && referring!=0 ) {
               std::string ref = referring->bookmark() ;
               comment( "Is the precondition of method " +
               reference( referring->def_class(),
               ref,
               referring->name()  ) +
               " fulfilled ? True if :" ) ;
            }
            else {
               out() << "  <span class=precondition>Precondition</span>" << endl ;
            }
            out() << "  <ul>" << endl ;
            PEL_Iterator* it = preConditions->create_iterator( 0 ) ;
            for( it->start() ; it->is_valid() ; it->go_next() ) {
               DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
               out() << "    <li>" ;
               write_code( f->referenced_string(
               method->def_class(),
               method,
               this ) ) ;
               out() << "</li>" << endl ;
            }
            it->destroy() ;
            out() << "  </ul>" << endl ;
         }
         if ( ( preConditions->count()>0 )&&( postConditions->count()>0 ) )  out() << "<br />" << endl ;
         if( postConditions->count()>0 ) {
            if( method->is_postcondition()  && referring!=0 ) {
               std::string ref = referring->bookmark() ;
               comment( "Is the postcondition of method " +
               reference( referring->def_class(),
               ref,
               referring->name()  ) +
               " fulfilled ? True if :" ) ;
            }
            else {
               out() << "  <span class=postcondition>Postcondition</span>" << endl ;
            }
            out() << "  <ul>" << endl ;
            PEL_Iterator* it = postConditions->create_iterator( 0 ) ;
            for( it->start() ; it->is_valid() ; it->go_next() ) {
               DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
               out() << "    <li>" ;
               write_code( f->referenced_string(
               method->def_class(),
               method,
               this) ) ;
               out() << "</li>" << endl ;
            }
            out() << "  </ul>" << endl ;
            it->destroy() ;
         }
         out() << "</div>" << endl ;
      }
      
      DOC_Method const* surdefinition = method->overiden_method() ;
      
      out() << "<table class=footer>" << endl ;
      out() << "<tr>" << endl;
      if( classeAProcess_Documentationr!=method->def_class() ) {
         out() << "  <td class=extends>" ;
         write( " Extension from "
         + reference( method->def_class(),
         "",
         method->def_class()->name() ) ) ;
         out() << "</td>" << endl;
      } else if( surdefinition!=0 ) {
         out() << "  <td class=overrides>" ;
         write( " Overridden from "
         + reference( surdefinition->def_class(),
         "",
         surdefinition->def_class()->name() ) ) ;
         out() << "</td>" << endl;
      }
      out() << "  <td class=link2anchor>" ;
      out() << "SUMMARY :: <a href=\"#method_summary\">METHOD</a></td></tr>"
      << endl ;
      out() << "</table>" << endl ;
      out() << "</div>" << endl ;
   }
}
