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

#include <DOC_WriterLaTex.hh>
#include <DOC_Category.hh>
#include <DOC_Function.hh>
#include <DOC_Class.hh>
#include <DOC_Method.hh>
#include <DOC_Package.hh>

#include <PEL_System.hh>

using std::endl ;

//--------------------------------------------------------------------
stringVector DOC_WriterLaTex::menu = stringVector(0) ;
//--------------------------------------------------------------------



//--------------------------------------------------------------------
DOC_WriterLaTex*
DOC_WriterLaTex:: create( PEL_Object* a_owner, std::string const& name ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WriterLaTex:: create" ) ;
   return new DOC_WriterLaTex( a_owner, name ) ;
}



//--------------------------------------------------------------------
DOC_WriterLaTex::DOC_WriterLaTex( PEL_Object* a_owner, std::string const& name ) 
//--------------------------------------------------------------------
      : DOC_Writer( a_owner, name + ".tex" ),
        file( name + ".tex" )
{
   file_principal = name=="index" ;
   
   if( !file_principal )
   {
      menu.append( "\\include{" + name + "}" ) ;
   }
   
   string nouv = file ;
   while( size_t ind = nouv.find( '_' ) <  nouv.length() )
   {
      nouv.erase( ind, ind ) ;
   }
   out() << "\\label{" << nouv << "}"<<std::endl ;
   
}



//--------------------------------------------------------------------
DOC_WriterLaTex::~DOC_WriterLaTex( void ) 
//--------------------------------------------------------------------
{
   if( file_principal )
   {
      for( size_t i=0 ; i<menu.size() ; i++ )
      {
         out() << menu(i) << endl ;
      }
      out() << "\\end{classing}" << endl ;
      out() << "\\end{document}" << endl ;
   }
}



//--------------------------------------------------------------------
string
DOC_WriterLaTex::underlined_text( size_t marqueur,
                          std::string const& text_initial ) const 
//--------------------------------------------------------------------
{
   string mark = "\\markOne" ;
   switch( marqueur )
   {
      case 1 : mark = "\\markOne" ;
         break ;
      case 2 : mark = "\\markTwo" ;
         break ;
      case 10 : mark = "\\markTen" ;
         break ;
   }
   
        
   string ret = mark + "{" + text_initial + "}" ;
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_WriterLaTex::comment( std::string const& com )
//--------------------------------------------------------------------
{
//   out() << "\\begin{verbatim}" << endl << com << "\\end{verbatim}" << endl ;
   out() << com ;
}




//--------------------------------------------------------------------
string 
DOC_WriterLaTex::reference( DOC_Class const* cl,
                        std::string const& label,
                        std::string const& display ) const
//--------------------------------------------------------------------
{
//    string nouv = chap ;
//    while( size_t ind = nouv.find( '_' ) <  nouv.length() )
//    {
//       nouv.erase( ind, ind ) ;
//    }
      
   string ref = underlined_text( 2, display ) ;
   return ref ;
}



//--------------------------------------------------------------------
void
DOC_WriterLaTex::ecrit_code( std::string const& com )
//--------------------------------------------------------------------
{
   out() << "\\begin{alltt}" << endl ;
   string to_print = com ;
   while( !to_print.empty() )
   {
      size_t idx ;
      if( ( idx = to_print.find( "<==>" ) )<to_print.length() )
      {
         to_print.replace( idx, 4, "\\begin{math}\\Longleftrightarrow\\end{math}" ) ;         
      }
      else if( ( idx = to_print.find( "==>" ) )<to_print.length() )
      {
         to_print.replace( idx, 3, "\\begin{math}\\Longrightarrow\\end{math}" ) ;         
      }
      else if( ( idx = to_print.find( "FORALL" ) )<to_print.length() )
      {
         to_print.replace( idx, 6, "\\begin{math}\\forall\\end{math}" ) ;
      }
      else
      {
         out() << true_type( to_print ) ;
         to_print = "" ;
      }
   }
   out() << "\\end{alltt}"  << endl ;
}



//--------------------------------------------------------------------
void
DOC_WriterLaTex::ecrit( std::string const& com )
//--------------------------------------------------------------------
{
   out() << true_type( com ) ;
}

   
//--------------------------------------------------------------------
string
DOC_WriterLaTex::true_type( std::string const& com,
                        bool retour_a_la_line_number )
//--------------------------------------------------------------------
{
   string ret ;
   for( string::const_iterator iter = com.begin() ;
        iter != com.end() ;
        ++iter )
   {
      if( *iter=='_' )
      {
         ret += "\\_" ;
      }
      else if( *iter=='#' )
      {
         ret += "\\#" ;
      }
      else if( *iter=='&' )
      {
         ret += "\\&" ;
      }
      else if( *iter=='<' )
      {
         ret += "$<$" ;
      }
      else if( *iter=='>' )
      {
         ret += "$>$" ;
      }
      else if( retour_a_la_line_number && *iter=='\n' )
      {
         ret += "\\\\\n" ;
      }
      else
      {
         ret += *iter ;
      }
   }
   return ret ;
}

   
//--------------------------------------------------------------------
void
DOC_WriterLaTex::finalize( void )
//--------------------------------------------------------------------
{
   std::string src = "LaTex" ;
   src += PEL_System::path_name_separator() ;
   src += "pelicans_style.tex" ;
   copy( src, "pelicans_style.tex" ) ;
//   copy( "LaTex/pelicans_style.tex", "pelicans_style.tex" ) ;
   
   out() << "\\input pelicans_style" << endl << endl ;

   out() << endl ;
   process_documentation_package( DOC_Package::master() ) ;
   out() << "\\begin{classing}" << endl ;
   
}



//--------------------------------------------------------------------
bool
DOC_WriterLaTex::source_reference( void ) const
//--------------------------------------------------------------------
{
   return false ;
}



//--------------------------------------------------------------------
void
DOC_WriterLaTex::process_documentation( DOC_Class const* classe )
//--------------------------------------------------------------------
{
   out() << "\\begin{class}{" <<
      true_type( classe->name() ) << "}" << endl ;   
   if( classe->mother()!=0 ) 
   {
      out() << "\\ininheritd{" << true_type( classe->mother()->name() ) << "}" << endl ;   
   }
   if( classe->is_virtual() &&
       classe->inherited_classes()->count()>0 )
   {
      out() << "\\begin{derivedClasses}" << endl ;   
      PEL_List const* inherited_classes = classe->inherited_classes() ;
      
      PEL_Iterator* it = inherited_classes->create_iterator(0) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Class*cl = static_cast<DOC_Class*>( it->item() ) ;
         out() << "\\itemDerivedClass{" + true_type( cl->name() ) + "}" << endl ;        
      }
      it->destroy() ;
      out() << "\\end{derivedClasses}" << endl ;   
   }
   out() << "\\begin{classComment}" << endl ;
   comment( DOC_Method::format_comment( *this,
                                          1,
                                          classe->comment(),
                                          '`',
                                          '\'',
                                          classe ) ) ;
   out() << "\\end{classComment}" << endl ;

   
   PEL_List const* elements = classe->owned_components() ;
   
   if( elements->count()>0 )
   {
      out() << "\\begin{classElemList}" << endl ;
      const string defaut_category = "Miscellanous" ;
      string category ;
      bool premierDomaine = true ;
      
      string heritage ;
      PEL_Iterator* it = elements->create_iterator(0) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_ClassItem*elem = static_cast<DOC_ClassItem*>( it->item() ) ;
         if( classe->has_to_document( elem ) )
         {
            string name_cat = elem->category()->name() ;
            if( name_cat.empty() ) name_cat = defaut_category ;
            bool nouvelleDOC_Category = name_cat != category ;
            if( nouvelleDOC_Category )
            {
               if( !premierDomaine )
               {
                  out() << "\\end{category}" << endl ;
               }
               premierDomaine = false ;
               category = name_cat ;
               out() << "\\begin{category}{" << true_type( category ) << "}" << endl ;
            }
            // nouvel_element() ;
            process_documentation( classe, elem ) ;            
         }
      }
      it->destroy() ;
      if( !premierDomaine )
      {
         out() << "\\end{category}" << endl ;
      }
      else
      {
         out() << "\\item {\\huge NO DOCUMENTATION PROVIDED.}" << endl ;
      }
      
      out() << "\\end{classElemList}" << endl ;
   }

   out() << "\\end{class}" << endl ;
}




//--------------------------------------------------------------------
void
DOC_WriterLaTex::process_documentation( DOC_Class const* classeAProcess_Documentationr,
      DOC_ClassItem * elem )  
//--------------------------------------------------------------------
{
   string proto = true_type( elem->prototype( *this ) ) ;
   out() << "\\begin{itemClass}{" << proto << "}" << endl ;
   if( !elem->comment().empty() )
   {
      out() << "\\begin{elemClassComment}" << endl ;
      comment( DOC_Method::format_comment( *this,
                                                  1,
                                                  elem->comment(),
                                                  '`',
                                                  '\'',
                                                 elem->def_class() ) ) ;
      out() << "\\end{elemClassComment}" << endl ;
   }
   else
   {
      out() << true_type( "\n" ) ;
   }
   DOC_Method const* method = elem->method() ;
   if( method != 0 )
   {
       PEL_List const* preConditions =  method->pre_conditions() ;
      PEL_List const* postConditions =  method->post_conditions() ;
      if( preConditions->count()>0 )
      {
         string tit = "Precondition" ;
         if( preConditions->count()>1 ) tit+="s" ;
         
         if( method->is_precondition() )
            tit +=  " implementation for method " +
               true_type( method->short_name() ) ;
         out() << "\\begin{elemClassCondition}{" << tit << "}" << endl  ;
         PEL_Iterator* it = preConditions->create_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
            out() << "\\condition{" ;
            ecrit_code(
               f->referenced_string(
                  method->def_class(),
                  method,
                  this ) ) ;
            out() << "}" << endl ;
         }
         it->destroy() ;
         out() << "\\end{elemClassCondition}" << endl  ;
      }
      if( postConditions->count()>0 )
      {
         string tit = "Postcondition" ;
         if( postConditions->count()>1 ) tit+="s" ;
         
         if( method->is_precondition() )
            tit +=  " implementation for method " +
               true_type( method->short_name() ) ;
         out() << "\\begin{elemClassCondition}{" << tit << "}" << endl  ;
         PEL_Iterator* it = postConditions->create_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Function* f=static_cast<DOC_Function*>( it->item() ) ;
            out() << "\\condition{" ;
            ecrit_code( f->referenced_string(
               method->def_class(),
               method,
               this ) ) ;
            out() << "}" << endl ;
         }
         it->destroy() ;
         out() << "\\end{elemClassCondition}" << endl  ;
      }
      DOC_Method const* surdefinition = method->overiden_method() ;
   
      if( classeAProcess_Documentationr!=method->def_class() )
      {
         out() << "\\extension{ " <<
            true_type( method->def_class()->name() ) << "}" << endl ;
         
      } else if( surdefinition!=0 )
      {
         out() << "\\overriden{ " <<
            true_type( surdefinition->def_class()->name() ) << "}" << endl ;
      }
   }
   out() << "\\end{itemClass}" << endl ;
}

//--------------------------------------------------------------------
void
DOC_WriterLaTex::process_documentation_package( DOC_Package const* pack,
                                   size_t level ) 
//--------------------------------------------------------------------
{   
   if( level==0 )
   {
      out() << "\\begin{packaging}" << endl ;
   }
   else if( level>=1 )
   {
      out() << "\\begin{package}{"<<true_type(pack->comment())<<"}" << endl ;
   }
   DOC_Writer::process_documentation_package( pack, level ) ;
   
   if( level==0 )
   {
      out() << "\\end{packaging}" << endl ;
   }
   else if( level>=1 )
   {
      out() << "\\end{package}" << endl ;
   }
}



//--------------------------------------------------------------------
void
DOC_WriterLaTex::disp_classes_tree( DOC_Class const* classe,
                                        int niveau,
                                        PEL_List* liste_classes )
//--------------------------------------------------------------------
{ 
   string display = classe->name() ;
   if( classe->is_plug_point() ) display = underlined_text(10, classe->name()) ;
   string comm ;
   string comm_court = classe->short_description() ;
   if( !comm_court.empty() )
   {   
      int st = 30 - classe->name().length() ;
      if( st<0 ) st=2 ;
      comm = string( st, '.' ) + comm_court ;
   }
   
   
   out() << "\\packagedClass{" ;
   for(int i=0 ; i<2*niveau ; i++ )
      out() << "\\tabClass " ;
   out() << true_type( display ) << "}" << endl ;

  
   DOC_Writer::disp_classes_tree( classe,
                                     niveau,
                                     liste_classes ) ;
} 
   
