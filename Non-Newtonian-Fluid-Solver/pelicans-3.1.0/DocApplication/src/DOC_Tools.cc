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

#include <iostream>
#include <string>
#include <algorithm>
#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <DOC_Tools.hh>
#include <PEL_Module.hh>
#include <fstream>
#include <dir.h>
#include <stdio.h>

using std::cout ;
using std::cerr ;
using std::endl ;

//------------------------------------------------------------------------  
string DOC_Tools:: buff = "" ;
string DOC_Tools:: file_courant = "" ;
stringVector DOC_Tools:: liste_files(0) ;
stringVector DOC_Tools:: liste_files_lus(0) ;
stringVector DOC_Tools:: absolu_liste_files_lus(0) ;
stringVector DOC_Tools:: liste_repertoires(0) ;
stringVector DOC_Tools:: external_dirs(0) ;
int DOC_Tools:: cc_mode = 0 ;
bool DOC_Tools:: est_cc = false ;
int DOC_Tools:: nbErrors = 0 ;
int DOC_Tools:: lastLineWhereErrorWasFound = -1 ;
int DOC_Tools:: nbLine_Numbers = 1 ;
bool DOC_Tools::verbosity = false ;
bool DOC_Tools::follow = false ;
bool DOC_Tools::doc_private = false ;
bool DOC_Tools::Wno_unresolved = false ;
string DOC_Tools::base_exe = "./" ;
//------------------------------------------------------------------------  


extern void switch_to_buffer( FILE * file ) ;
extern void switch_mode( void ) ;
extern int yyparse( void ) ;
#ifdef YYDEBUG
extern int yydebug ;
#endif
extern int yy_flex_debug ;



//--------------------------------------------------------------------
void
DOC_Tools::parse_explorer( PEL_ModuleExplorer const* exp ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Writer::parse_explorer" ) ;
   
   static bool prem = true ;
   if( prem )
   {
#ifdef YYDEBUG
      yydebug = 0 ;
#endif
      yy_flex_debug = 0 ;
      prem = false ;
   }
   
   if( exp->has_entry( "debug" ) && exp->bool_data( "debug" ) )
   {
#ifdef YYDEBUG
      yydebug = 1 ;
      cerr << "Debugging mode activated" << endl ;
#endif
      yy_flex_debug = 1 ;
   }
   if( exp->has_entry( "include_directories" ) )
   {
      liste_repertoires = exp->stringVector_data( "include_directories" ) ;
   }
   if( exp->has_entry( "external" ) )
   {
      external_dirs = exp->stringVector_data( "external" ) ;
   }
   if( exp->has_entry( "exclude" ) )
   {
      liste_files_lus = exp->stringVector_data( "exclude" ) ;      
   }
   if( exp->has_entry( "verbose" ) )
   {
      verbosity = exp->bool_data( "verbose" ) ;
   }
   if( exp->has_entry( "follow" ) )
   {
      follow = exp->bool_data( "follow" ) ;
   }
   if( exp->has_entry( "private" ) )
   {
      doc_private = exp->bool_data( "private" ) ;
   }
   if( exp->has_entry( "Wno_unresolved" ) )
   {
      Wno_unresolved = exp->bool_data( "Wno_unresolved" ) ;
   }
   liste_repertoires.append( "" ) ;
   liste_repertoires.append( "." ) ;
   if( exp->has_entry( "scanned_files" ) )
   {
      liste_files = exp->stringVector_data( "scanned_files" ) ;
   }
   if( exp->has_entry( "scanned_dir" ) )
   {
      stringVector const& dirs = exp->stringVector_data( "scanned_dir" ) ;
      stringVector ext(0) ;
      ext.append( ".cc" ) ;
      if( exp->has_entry( "ext_filter" ) )
         ext = exp->stringVector_data( "ext_filter" ) ;
      for( size_t i=0 ; i<dirs.size() ; i++ )
      {
         liste_repertoires.append( dirs(i) ) ;
         stringVector res = ls( dirs(i), ext ) ;
         for( size_t j=0 ; j<res.size() ; j++ )
         {
            liste_files.append( res(j) ) ;
         }
      }
   }
   if( liste_files.size()==0 )
   {
      PEL_Error::object()->raise_plain( "Empty scanned file list" ) ;
   }
   else if(message() )
   {
      std::cout << "Scanning following files : \n"
                << liste_files << std::endl ;
   }
   
}



//------------------------------------------------------------------------  
bool
DOC_Tools:: message( void ) 
//------------------------------------------------------------------------  
{
   return verbosity ;
}



//------------------------------------------------------------------------  
bool
DOC_Tools:: private_doc( void ) 
//------------------------------------------------------------------------  
{
   return doc_private ;
}



//------------------------------------------------------------------------  
bool
DOC_Tools:: is_implementation_file( void ) 
//------------------------------------------------------------------------  
{
   return est_cc ;
}

//------------------------------------------------------------------------  
void
DOC_Tools:: toggle_cc_mode( void ) 
//------------------------------------------------------------------------  
{
   cc_mode = 1 - cc_mode ;
}

//------------------------------------------------------------------------  
bool
DOC_Tools:: is_cc_mode( void ) 
//------------------------------------------------------------------------  
{
   return cc_mode==0 ;
}



//------------------------------------------------------------------------  
void
DOC_Tools:: add_new( string const& file ) 
//------------------------------------------------------------------------  
{
   if( message() )
      std::cout << "Adding new file : "<<file << std::endl ;
   if( !liste_files.has(file) )
   {
      liste_files.append( file ) ;
   }
}



//------------------------------------------------------------------------  
void
DOC_Tools:: declare_read( string const& file ) 
//------------------------------------------------------------------------  
{
   if( message() )
      std::cout << "Processed file : "<<file << std::endl ;
   
   liste_files_lus.append( file ) ;
}



//------------------------------------------------------------------------  
std::string
DOC_Tools:: read_file( std::string const& file ) 
//------------------------------------------------------------------------  
{
   std::string result ;
   
   for( size_t i=0 ; i<absolu_liste_files_lus.size() ; i++ )
   {
      if( PEL_Module::basename( absolu_liste_files_lus(i) ) == file ) result = absolu_liste_files_lus(i) ;
   }
   
   return result ;
}



//------------------------------------------------------------------------  
void
DOC_Tools:: evaluate_mode( void ) 
//------------------------------------------------------------------------  
{
   est_cc = true ;
   if( file_courant.find( ".cc" ) >= file_courant.length() &&
       file_courant.find( ".icc" ) >= file_courant.length() )
   {
      est_cc = false ;
   }
   cc_mode = 0 ;
   if( verbosity && est_cc )
   {
      cerr << "Using CC mode" << endl ;
   }
}



//------------------------------------------------------------------------  
std::string
DOC_Tools:: external( std::string const& file ) 
//------------------------------------------------------------------------  
{
   std::string result ;
   for( size_t i=0 ; i<external_dirs.size() ; i++ )
   {
      std::string nam = external_dirs(i)+"/"+file ;
      std::ifstream in( nam.c_str() ) ;
      if( in )
      {
         result = nam ;
         in.close() ;
         break ;
      }
   }
   return result ;
}



//------------------------------------------------------------------------  
bool
DOC_Tools:: read( string const& a_File ) 
//------------------------------------------------------------------------  
{
   if( message() )
      std::cout << "Reading : "<<a_File << std::endl ;
   string file = a_File ;
   
   bool ret = false ;
   if( !liste_files_lus.has( file ) )
   {
      liste_files_lus.append( file ) ;
      int nbLine_NumbersOld = nbLine_Numbers ;
      string fileOld = file_courant ;
      string old_buff ;
      old_buff = buff ;
      string fileToOpen=file ;
      FILE* newFile = 0 ;
      for( size_t i=0 ;
           i<liste_repertoires.size() && newFile==0 ; i++ )
      {
         fileToOpen = liste_repertoires(i) ;
         if( !fileToOpen.empty() ) fileToOpen += "/" ;
         fileToOpen += file ;
         newFile = fopen( fileToOpen.c_str(), "r" ) ;
      }
      if( newFile != 0 )
      {
         ret = true ;
	 if( verbosity && !file_courant.empty() )
	 {
            XWarningE( file_courant, nbLine_Numbers, "Now, reading file " << fileToOpen ) ;
            
	 }
	 file_courant = fileToOpen ;
         absolu_liste_files_lus.append( fileToOpen ) ;
	 nbLine_Numbers = 1 ;
	 evaluate_mode() ;
         switch_to_buffer( newFile ) ;
	 fclose( newFile ) ;
	 if( verbosity )
	 {
	    cerr << "      File " << file << " read " << endl ;
	 } 
      }
      file_courant = fileOld ;
      if( ret && ( !fileOld.empty() ) )
      {
	 evaluate_mode() ;
	 nbLine_Numbers = nbLine_NumbersOld ;
	 buff = old_buff ;
         switch_mode() ;
	 yyparse() ;
      }
   }
   if( !ret )
   {
      if( message() )
         std::cout << a_File << " skipped !" << std::endl ;
   }
   
   return ret ;
}



//------------------------------------------------------------------------  
void
DOC_Tools:: read_files( void ) 
//------------------------------------------------------------------------  
{
   if( message() )
      std::cout << "Begin to process files" << std::endl ;
   for( size_t i=0 ; i<liste_files.size() ; i++ ) 
   {
      read( liste_files(i) ) ;
   }
}


//----------------------------------------------------------------------
int 
DOC_Tools:: record( std::string const& chain ) 
//----------------------------------------------------------------------
{
   size_t cr=chain.find('\n') ;
   size_t cr_last = 0 ;
   for( ;
        cr<chain.length() ;
        cr=chain.find('\n', cr+1 ) )
   {
      nbLine_Numbers++ ;
      cr_last = cr+1 ;
      buff="" ;
   }
   buff += chain.substr( cr_last, chain.length()-cr_last ) ;
   if( chain.find("TODO" ) < chain.length() )
      XWarningE( file_courant, nbLine_Numbers, chain ) ;
   return 0 ;
}
  

//------------------------------------------------------------------------  
void
DOC_Tools:: yyerror( string const& s )
//------------------------------------------------------------------------  
{
   if( nbLine_Numbers != lastLineWhereErrorWasFound )
   {
      if( !s.empty() )
      {
         XWarningE( file_courant, nbLine_Numbers, "peldoc parsing : " << s ) ;
      }
      lastLineWhereErrorWasFound = nbLine_Numbers ;
      std::cerr << "Last line read : " << buff << std::endl ;
      if( nbErrors > maxNbErrors )
      {
         std::cerr << "Too many errors : bailing out" << std::endl ;
         PEL_System::exit( 1 ) ;
      }
      std::cerr << std::endl ;      
   }
}



//------------------------------------------------------------------------  
void
DOC_Tools:: warning( string const& s)
//------------------------------------------------------------------------  
{
   if( verbosity )
   {
      std::cerr << "## Warning in " << file_courant << " L[" << nbLine_Numbers << "] : " << s << endl ;
   }
}



//------------------------------------------------------------------------  
string
DOC_Tools:: text_comment( string const& comment )
//------------------------------------------------------------------------  
{
   string ret ;
   if( comment.length()>2 )
   {
      string::const_iterator iter_deb, iter_end ;
      for( iter_deb = comment.begin() ;
           iter_deb!=comment.end() ;
           ++iter_deb )
      {
         if( !( *iter_deb == '*' ||
                *iter_deb == '/' ||
                *iter_deb == '-' ||
                *iter_deb == ' '
            ) )
         {
            break ;
         }
      }
      size_t deb = iter_deb-comment.begin() ;
      PEL_ASSERT( deb<=comment.length() ) ;

      for( iter_end = comment.end()-1 ;
           iter_end!=comment.begin() ;
           --iter_end )
      {
         if( !( *iter_end == '*' ||
                *iter_end == '/' ||
                *iter_end == '-' ||
                *iter_end == ' '   ) )
         {
            break ;
         }
      }
      int len = iter_end - iter_deb +1 ;
   
      if( len>0 )
      {
         ret = comment.substr( deb, len ) ;
      }
   }
   
   return ret ;
}



//------------------------------------------------------------------------  
string
DOC_Tools:: text_expression_bool( string const& expr )
//------------------------------------------------------------------------  
{
   string ret ;
   size_t deb ;
   size_t fin ;
   if( expr.find( "bool" )<expr.length() )
   {
      deb = expr.find_first_of( '=' ) + 1 ;
      fin = expr.find_last_of( ';' ) - 1 ;
   }
   else
   {
      deb = expr.find_first_of( '(' ) + 1 ;
      fin = expr.find_last_of( ')' ) - 1 ;
   }
   if( deb>fin ) 
   {
      Erreur( "Unable to understand boolean expression " << expr ) ;
   }
   size_t len = fin - deb +1 ;
   
   if( len>0 )
   {
      ret = expr.substr( deb, len ) ;
   }
   return ret ;
}


//------------------------------------------------------------------------  
string
DOC_Tools:: text_check( string const& check )
//------------------------------------------------------------------------  
{
   string ret ;
   size_t deb = check.find_first_of( '(' ) + 1 ;
   size_t fin = check.find_last_of( ')' ) - 1 ;
   while( deb < check.length() && check.at( deb )==' ' )
   {
      deb++ ;
   }
      
   PEL_ASSERT( deb<=fin ) ;
   size_t len = fin - deb +1 ;
   
   if( len>0 )
   {
      ret = check.substr( deb, len ) ;
   }
   return ret ;
}



//------------------------------------------------------------------------  
string
DOC_Tools:: upper( string const& s ) 
//------------------------------------------------------------------------  
{
   string ret ;
   for( string::const_iterator c = s.begin() ;
        c != s.end() ;
        ++c )
   {
      if( *c >='a' && *c <='z' )
      {
         ret += *c + 'A' - 'a' ;
      }
      else
      {
         ret += *c ;
      }
   }
   return ret ;
}



//------------------------------------------------------------------------  
string 
DOC_Tools:: file_include( string const& str )
//------------------------------------------------------------------------  
{
   PEL_ASSERT( str.find( "#include" ) == 0 ) ;
   size_t deb = str.find( '<' ) ;
   if( deb>=str.length() )
   {
      deb = str.find( '"' ) ;
   }
   size_t fin = str.find( '>' ) ;
   if( fin>=str.length() )
   {
      fin = str.find_last_of( '"' ) ;
   }
   
   return str.substr( deb+1, fin-deb-1 ) ;
   
}



//------------------------------------------------------------------------  
string const& 
DOC_Tools:: file( void )
//------------------------------------------------------------------------  
{
   static string file ;
   file = file_courant ;
#ifndef _WIN32
   if( follow )
   {
      char link[1024] ;
      int lien_len ;   
      if( (lien_len=readlink( file.c_str(), link, 1024 )) !=-1 )
      {
         link[lien_len] = '\0' ;
         file = link ;
      }
   }
#endif
   return file ;
}



//------------------------------------------------------------------------  
int
DOC_Tools:: current_line_number( void ) 
//------------------------------------------------------------------------  
{
   return nbLine_Numbers ;
}


//------------------------------------------------------------------------  
stringVector
DOC_Tools:: ls( string const& dirname, stringVector const& exts ) 
//------------------------------------------------------------------------  
{
   PEL_LABEL( "DOC_Tools:: ls" ) ;
   stringVector result(0) ;
   
   DIR*dir = opendir( dirname.c_str() ) ;
   if( dir != 0 )
   {
      dirent* entry ;  
      while( (entry=readdir(dir))!=0 )
      {
         string name( entry->d_name ) ;
         for( size_t j=0 ; j<exts.size() ; j++ )
         {
            string const& ext = exts(j) ;
            size_t idx = name.rfind( ext ) ;
        
            if( idx<name.length() && idx+ext.length()==name.length() )
            {
               result.append( name ) ;
            }
         }
      }
   
      closedir( dir ) ;
   }
   return result ;
}



//------------------------------------------------------------------------  
bool
DOC_Tools:: warn_unreachable_condition( void ) 
//------------------------------------------------------------------------  
{
   return !Wno_unresolved ;
}




