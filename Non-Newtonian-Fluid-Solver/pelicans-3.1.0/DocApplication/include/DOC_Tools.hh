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

#ifndef OUTILS_HH
#define OUTILS_HH

#include <string>
#include <iostream>
#include <stringVector.hh>
#include <PEL_Error.hh>

using std::cerr ;
using std::endl ;
using std::ostream ;
using std::string ;


#define Erreur(x) { std::cerr << "Error : " << x << std::endl ;\
   PEL_Error::object()->raise_plain( "Aborting..." ) ; }

#define XWarningE(file,line,x) { std::cerr << file << ":" << line << ": " << x << std::endl ; }

#define DECLARATION_DEFAUT_SANS_VOID(CLASSE) \
CLASSE operator=(CLASSE const&) ;\
CLASSE( CLASSE const& )

#define DECLARATION_DEFAUT(CLASSE) \
CLASSE(void) ;\
DECLARATION_DEFAUT_SANS_VOID(CLASSE)

#define MOI { return this ; }
#define PASMOI { Erreur( "Pas le bon type" ) ; }

/* Le type des symbols terminaux est declare comme un pointeur
 sur des DOC_Symbol. */
class DOC_Symbol ;
#define YYSTYPE DOC_Symbol *
/* La variable globale yylval est utilisée pour transmettre la valeur
   des symbols */
extern YYSTYPE yylval ;

class PEL_EXPORT DOC_Tools
{
   public://-----------------------------------------------------------------
      static void parse_explorer( PEL_ModuleExplorer const* exp ) ;
      
      static void toggle_cc_mode(void) ;
      static bool is_cc_mode( void ) ;
      static void evaluate_mode( void ) ;
      static bool is_implementation_file( void ) ;
      static void add_new( string const& file ) ;
      static std::string read_file( std::string const& file ) ;
      static void declare_read( string const& file ) ;
      static bool read( string const& file ) ;
      static void read_files( void ) ;
      static bool message( void ) ;
      /* lire : cette procedure stocke dans un buffer la line_number courante
         afin de connaitre en cas d'erreur le dernier mot choisi. */
      static int record( std::string const& chain ) ;
      // Check for `file' in external directories.
      static std::string external( std::string const& file ) ;
      /* ote les balises de comments */
      static string text_comment( string const& comment ) ;
      static string text_check( string const& check ) ;
      static string text_expression_bool( string const& expr ) ;
      /* yyerror : procedure standard de sortie de message d'erreur. */
      static void yyerror( string const& s ) ;
      /* yyerror : Avertissements. */
      static void warning( string const& s ) ;
      static string upper( string const& s ) ;
      static string const& file( void ) ;
      static string file_include( string const& str )  ;
      static int current_line_number( void ) ;
      static bool private_doc( void ) ;      
      static stringVector ls( string const& dirname,
                              stringVector const& exts ) ;
      static bool warn_unreachable_condition( void ) ;
      
   protected://---------------------------------------------------------------
   private://-----------------------------------------------------------------
      DECLARATION_DEFAUT( DOC_Tools ) ;
      enum { maxNbErrors = 5 } ;
         
      static string buff ;
      static string file_courant ;
      static stringVector liste_files ;
      static stringVector liste_files_lus ;
      static stringVector absolu_liste_files_lus ;
      static stringVector liste_repertoires ;
      static stringVector external_dirs ;
      static int cc_mode ;
      static bool est_cc ;
      static bool doc_private ;
      static bool verbosity ;      
      static int nbErrors  ;
      static int lastLineWhereErrorWasFound ;
      static int nbLine_Numbers ;
      static bool follow ;
      static string base_exe ;
      static bool Wno_unresolved ;
      
} ;



#endif // OUTILS_HH
