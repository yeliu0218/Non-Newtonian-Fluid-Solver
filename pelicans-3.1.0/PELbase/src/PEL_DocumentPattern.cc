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

#include <PEL_DocumentPattern.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using std::setw ;

PEL_DocumentPattern const* PEL_DocumentPattern::PROTOTYPE = new PEL_DocumentPattern() ;


//----------------------------------------------------------------------
PEL_DocumentPattern:: PEL_DocumentPattern( void )
//----------------------------------------------------------------------
   : PEL_Application( "document_pattern" )
   , FILENAME("")
   , CSS("")
   , PATTERN(0)
   , NAMES(0)
   , TO_DOCUMENT( 0 )
   , CHECK_SUMS(0)
   , KEYS( 0 )
   , VALS( 0 )
{
}

//----------------------------------------------------------------------
PEL_DocumentPattern:: PEL_DocumentPattern( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , FILENAME("")
   , CSS("")
   , OUT( 0 )
   , PATTERN(0)
   , NAMES(0)
   , TO_DOCUMENT( PEL_Vector::create( this , 0 ) )
   , FORMAT( "html" )
   , OUTPUT_FILENAME( "" )
   , VERBOSE( false )
   , BASE_URL( "file://" )
   , CHECK_SUMS(0)
   , DICO( "" )
   , KEYS( 0 )
   , VALS( 0 )
{
   PEL_LABEL( "PEL_DocumentPattern:: PEL_DocumentPattern" ) ;
   FILENAME = exp->string_data( "pattern_filename" ) ;
   if( exp->has_entry( "output_filename" ) )
      OUTPUT_FILENAME = exp->string_data( "output_filename" ) ;
   if( exp->has_entry( "stylesheet" ) )
      CSS = exp->string_data( "stylesheet" ) ;
   if( exp->has_entry( "base_url" ) )
      BASE_URL = exp->string_data( "base_url" ) ;
   if( exp->has_entry( "format" ) )
   {
      FORMAT = exp->string_data( "format" ) ;
      //exp->test_data_in( "format", "txt,html,latex" ) ;
      exp->test_data_in( "format", "html" ) ;
   }
   if( exp->has_entry( "dictionary" ) )
      DICO = exp->string_data( "dictionary" ) ;
   
   VERBOSE = ( exp->has_entry( "verbose" ) && exp->bool_data( "verbose" ) ) ;
}

//----------------------------------------------------------------------
PEL_DocumentPattern:: PEL_DocumentPattern( PEL_Object* a_owner,
                                           stringVector& args )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , FILENAME("")
   , OUT( 0 )
   , PATTERN(0)
   , NAMES(0)
   , TO_DOCUMENT( PEL_Vector::create( this , 0 ) )
   , FORMAT( "txt" )
   , OUTPUT_FILENAME( "" )
   , CHECK_SUMS(0)
   , DICO( "" )
   , KEYS( 0 )
   , VALS( 0 )
{
   PEL_LABEL( "PEL_DocumentPattern:: PEL_DocumentPattern" ) ;

   for( size_t i=0 ; i<args.size() ; )
   {
      if( args(i)=="-output_filename" && args.size()>i+1 )
      {
         OUTPUT_FILENAME = args(i+1) ;
         args.remove_at(i) ;
         args.remove_at(i) ;
      }
      else if( args(i)=="-stylesheet" && args.size()>i+1 )
      {
         CSS = args(i+1) ;
         args.remove_at(i) ;
         args.remove_at(i) ;
      }

      else if( args(i)=="-dictionary" && args.size()>i+1 )
      {
         DICO = args(i+1) ;
         args.remove_at(i) ;
         args.remove_at(i) ;
      }
      
      else if( args(i)=="-html" )
      {
         FORMAT = "html" ;
         args.remove_at(i) ;
      }
      else
      {
         i++ ;
      }     
   }
   if( args.size() == 0 )
   {
      notify_error_in_arguments() ;
   }
   
   FILENAME = args(0) ;
   args.remove_at(0) ;
}

//----------------------------------------------------------------------
PEL_DocumentPattern:: ~PEL_DocumentPattern( void )
//----------------------------------------------------------------------
{
   if( this == PROTOTYPE )
   {
      PROTOTYPE = 0 ;  
   }
}

//----------------------------------------------------------------------
PEL_DocumentPattern* 
PEL_DocumentPattern:: create_replica( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;

   PEL_DocumentPattern* result = new PEL_DocumentPattern( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_DocumentPattern* 
PEL_DocumentPattern:: create_replica_from_args( PEL_Object* a_owner,
                                      stringVector& args ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: create_replica_from_args" ) ;

   PEL_DocumentPattern* result = new PEL_DocumentPattern( a_owner, args ) ;

   PEL_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: run" ) ;

   if( VERBOSE )
   {  
      PEL::out()
         << "PEL_DocumentPattern : processing documentation of controler file : "
         << FILENAME << std::endl ;
      PEL::out() << " in format " << FORMAT ;
      if( OUTPUT_FILENAME.empty() ) PEL::out() << " to standard output " ;
      else PEL::out() << " to file : " << OUTPUT_FILENAME ;
      PEL::out() << std::endl ;
      if( !DICO.empty() ) PEL::out() << "Using dictionary: " << DICO << std::endl ;
      PEL::out() << "Loading " << FILENAME << " ... " ;
   }
   if( !DICO.empty() ) load_dictionary() ;
   
   PEL_Module const* ctrl = PEL_Module::create( this, "ROOT", FILENAME ) ;
   if( VERBOSE )
   {
      PEL::out() << "done ! " << std::endl ;
   }
   PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( this, ctrl ) ;
   exp = exp->create_subexplorer( exp, "Base" ) ;
   
   PEL_ASSERT( exp->string_data( "revision" ) == "2.0" ) ;
   
   PEL_ModuleExplorer* root = 0 ;
   
   for( exp->start_module_iterator() ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer * sexp = exp->create_subexplorer(this) ;
      if( sexp->name() == "Pattern" )
      {
         PATTERN = sexp ;
      }
      else
      {
         PEL_CHECK( root==0 ) ;
         root = sexp ;
      }
   }

   PEL_ASSERT( root != 0 ) ;
   add( root ) ;

   if( FORMAT == "html" )
   {
      OUT = new OutputInterfaceHtml(OUTPUT_FILENAME,FILENAME,CSS,*this) ;
   }
   else
   {
      PEL_Error::object()->raise_plain(
         "*** PEL_DocumentPattern error:\n"
         "    format: \""+FORMAT+"\" not yet implemented !" ) ;
   }

   PEL_ASSERT( OUT!=0 ) ;
   
   for( size_t i=0 ; i<NAMES.size() ; i++ )
   {
      PEL_ModuleExplorer* module =
         static_cast<PEL_ModuleExplorer*>( TO_DOCUMENT->at( i ) ) ;
      size_t occ =0 ;
      for( size_t j=0 ; j<i ; j++ )
         if( NAMES(j)==NAMES(i) ) occ++ ;
      
      document_module( module,  i, occ ) ;
   }
   
   delete OUT ;
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: print_usage( void ) const
//----------------------------------------------------------------------
{
   PEL::out() << std::endl<<
      "Application to produce documentation about PELICANS datadeck from controler datafile"
              << std::endl << std::endl;
   PEL::out() << "  document_pattern [ -output_filename output_file ] [ -stylesheet stylesheet.css ] "
              << " [ -dictionary dico_file.txt ] [ -html ] <controler filename>" << std::endl ;
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: print_options( void ) const
//----------------------------------------------------------------------
{
  
   PEL::out() << "     -output_filename output_file : name of output file " << std::endl ;
   PEL::out() << "     -stylesheet stylesheet.css : name of alternative stylesheet " << std::endl ;
   PEL::out() << "     -html          : HTML output format (default)" << std::endl ;
}

//----------------------------------------------------------------------
int
PEL_DocumentPattern:: add( PEL_ModuleExplorer* module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: add" ) ;
   size_t result = PEL::bad_index() ;
   
   std::string const& name =  module->name() ;
   // Search for existing module name already registered
   size_t new_check = 0 ;
   if(  NAMES.has( name ) )
   {
      new_check = checksum(module) ;
      // Controls checksum of already registered modules with same name and new one
      for( size_t idx=0 ; idx<NAMES.size() && result==PEL::bad_index() ; idx++ )
      {
         if( NAMES(idx)==name )
         {
            if( CHECK_SUMS(idx)==0 )
            {
               CHECK_SUMS(idx) = checksum( static_cast<PEL_ModuleExplorer const*>( TO_DOCUMENT->at(idx) ) ) ;
            }
            
            size_t check = CHECK_SUMS(idx) ;
            if( new_check==check )
            {
               result = idx ;
            }
         }
      }
   }
   
   // If no registered module matches, store new one
   if( result==PEL::bad_index() )
   {
      NAMES.append( name ) ;      
      TO_DOCUMENT->append( module ) ;
      CHECK_SUMS.append(new_check) ;
      result = NAMES.size()-1 ;
   }
   
   PEL_CHECK_POST( result!=PEL::bad_index() ) ;
   
   return result ;
}

//----------------------------------------------------------------------
size_t
PEL_DocumentPattern:: checksum( PEL_ModuleExplorer const* module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: checksum" ) ;

   std::ostringstream os ;
   module->print(os,0) ;
   std::string str = os.str() ;
   
   const char * data = str.c_str() ;
   size_t result = 0 ;
   while( (*data)!='\0' )
   {
      result += (*data) ;
      data++ ;
   }
   
   return result ;
}


//----------------------------------------------------------------------
void
PEL_DocumentPattern:: document_module( PEL_ModuleExplorer* module,
                                       int order,
                                       size_t occ )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: document_module" ) ;
   
   if( VERBOSE )
   {
      PEL::out() << "Documenting : " << module->absolute_path_name() << std::endl ;
   }
   
   if( module->string_data( "_type" ) != "S" ||
       !module->has_entry( "_access" ) )
   {
      PEL_Error::object()->display_info(
         "*** PEL_DocumentPattern warning:\n"
         "    syntax error in module: \""+module->absolute_path_name()+"\"\n" ) ;
   }
   else
   {
         
      PEL_CHECK( module->string_data( "_type" ) == "S" ) ;
      PEL_CHECK( module->has_entry( "_access" ) ) ;
      std::string help = ( module->has_entry( "_help" ) ? module->string_data("_help") : "" ) ;
      std::string url = ( module->has_entry( "_url" ) ? BASE_URL + module->string_data("_url") : "" ) ;
      std::ostringstream is ;
      stringVector plotter(0) ;
      if( module->has_entry( "_plotter" ) )
      {
         plotter = module->stringVector_data( "_plotter" ) ;
                     
      }
      is << module->string_data("_name") ;
      if( occ>0 ) is << " (" << occ << ")" ;
         
      OUT->start_module(
         is.str(),
         module->string_data("_access"),
         help,
         url,
         order,
         plotter ) ;
      stringVector conditions(1) ;

      conditions(0)="" ;
      
      explore_conditions( module, conditions ) ;
      if( module->has_entry( "_to" ) )
      {
         std::string to = module->string_data(  "_to" ) ;
         explore_conditions( PATTERN->create_subexplorer( this, to ), conditions ) ;
      }
      
      for( size_t i=0 ; i<conditions.size() ; i++ )
      {
         if( i==1 ) OUT->start_list_of_conditions() ;
         std::string cond = conditions(i) ;
         if( !cond.empty() ) OUT->start_condition( cond ) ;
         
         size_t pass_nb = OUT->wished_pass_number() ;
         
         for( PASS=0; PASS<pass_nb ; PASS++ )
         {
            OUT->set_pass(PASS) ;
            explore( module, cond ) ;
            if( module->has_entry( "_to" ) )
            {
               std::string to = module->string_data(  "_to" ) ;
               explore( PATTERN->create_subexplorer( this, to ), cond ) ;
            }
            OUT->finish_pass() ;
         }
         if( !cond.empty() ) OUT->finish_condition( cond ) ;
         if( i==conditions.size()-1 && i>0 ) OUT->finish_list_of_conditions() ;
      }
      
   
      OUT->finish_module(
         module->string_data("_name"),
         module->string_data("_access") ) ;
   }
}
   
//----------------------------------------------------------------------
void
PEL_DocumentPattern:: explore( PEL_ModuleExplorer* module,
                               std::string cond )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: explore" ) ;

   
   if( !cond.empty() )
   {
      for( module->start_module_iterator() ;
           module->is_valid_module() ;
           module->go_next_module() )
      {
         PEL_ModuleExplorer* s = module->create_subexplorer( module ) ;
         std::string const& type = s->string_data( "_type" ) ;
         if( type=="C" )
         {
            std::string const& a_cond = s->string_data( "_if" ) ;
            if( a_cond==cond )
            {
               std::string help = ( s->has_entry( "_help" ) ? s->string_data("_help") : "" ) ;
               std::string url = ( s->has_entry( "_url" ) ? BASE_URL + s->string_data("_url") : "" ) ;
               OUT->add_help_on_condition( help, url ) ;
               explore( s, "" ) ;
            }
         }
      }
   }
   else
   {
      for( size_t j=0 ; j< 4 ; j++ )
      {
         for( module->start_module_iterator() ;
              module->is_valid_module() ;
              module->go_next_module() )
         {
            PEL_ModuleExplorer* s = module->create_subexplorer( module ) ;
            std::string const& type = s->string_data( "_type" ) ;
            if( type=="C" )
            {
            }
            else if( type=="S" )
            {
               if( j==0 )
               {
                  bool show = true ;
                  if( s->has_entry( "_show" ) )
                  {
                     PEL_Data const* data = s->abstract_data( 0, "_show" ) ;
                     if( data->data_type()==PEL_Data::Bool ) show = s->bool_data( "_show" ) ;
                     data->destroy() ;
                  }
                  bool edit = ( s->has_entry( "_edit" ) ? s->bool_data( "_edit" ) : true ) ;
                  
                  if( show && edit && PASS==0 )
                     OUT->add_module( s->string_data( "_name" ), s->string_data( "_access" ), add(s) ) ;
               }
            }
            else if( type=="V" )
            {
               if( j==2 )
               {
                  document_variable( s ) ;
               }  
            }
            else
            {
               if( j==1 )
               {
                  document_entry( s ) ;
               }
            }
         }
      }
   } 
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: explore_conditions( PEL_ModuleExplorer* module,
                                          stringVector& conditions )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: explore_conditions" ) ;
   for( module->start_module_iterator() ;
        module->is_valid_module() ;
        module->go_next_module() )
   {
      PEL_ModuleExplorer* s = module->create_subexplorer( module ) ;
      std::string const& type = s->string_data( "_type" ) ;
      if( type=="C" )
      {
         std::string const& cond = s->string_data( "_if" ) ;
         if( !conditions.has(cond) ) conditions.append(cond) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: document_entry( PEL_ModuleExplorer* module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: document_entry" ) ;
   if( ! ( module->string_data( "_type" ) != "S" &&
           module->string_data( "_type" ) != "C" &&
           module->has_entry( "_access" ) ) )
   {
      PEL_Error::object()->display_info(
         "*** PEL_DocumentPattern warning:\n"
         "    syntax error for description of entry: \""
                           + module->absolute_path_name()+"\"\n" ) ;
   }
   else
   {
   
      std::string def = "" ;
      if( module->has_entry( "_default" ) )
      {
         PEL_Data const* data = module->abstract_data(0,"_default") ;
         def = data->value_as_string(0) ;
         data->destroy() ;
      }
   
      stringVector in(0) ;
      if( module->has_entry( "_in" ) )
      {
         PEL_Data const* data = module->abstract_data(0,"_in") ;
         std::string vec = data->value_as_string(0) ;
         PEL_ASSERT( vec.length()>=4 ) ;
         
         vec = vec.substr(2,vec.length()-4) ; // Remove < and >
         
         in = stringVector(vec,' ') ;
         data->destroy() ;
      }
      
      stringVector help_in(0) ;
      if( module->has_entry( "_help_in" ) )
      {
         help_in = module->stringVector_data("_help_in") ;
      }
   
      std::string sel_in =
         ( module->has_entry( "_select_in" ) ? module->string_data("_select_in") : "" ) ;
   
      std::string where =
         ( module->has_entry( "_where" ) ? module->string_data("_where") : "" ) ;
   
      std::string help =
         ( module->has_entry( "_help" ) ? module->string_data("_help") : "" ) ;

      std::string unit =
         ( module->has_entry( "_unit" ) ? module->string_data("_unit") : "" ) ;

      std::string test =
         ( module->has_entry( "_test" ) ? module->string_data("_test") : "" ) ;

      std::string formula =
         ( module->has_entry( "_formula" ) ? module->string_data("_formula") : "" ) ;

      bool show = true ;
      if( module->has_entry( "_show" ) )
      {
         PEL_Data const* data = module->abstract_data( 0, "_show" ) ;
         if( data->data_type()==PEL_Data::Bool ) show = module->bool_data( "_show" ) ;
         data->destroy() ;
      }                  
      bool edit = ( module->has_entry( "_edit" ) ? module->bool_data( "_edit" ) : true ) ;
      bool unique = ( module->has_entry( "_unique" ) ? module->bool_data( "_unique" ) : false ) ;
      if( show && edit )
      {
         OUT->add_entry( module->string_data("_name"),
                         module->string_data("_type"),
                         module->string_data("_access"),
                         def,
                         in,
                         help_in,
                         sel_in,
                         where,
                         test,
                         unit,
                         help,
                         formula,
                         unique ) ;
      }
   }
}


//----------------------------------------------------------------------
void
PEL_DocumentPattern:: document_variable( PEL_ModuleExplorer* module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: document_variable" ) ;
   PEL_CHECK_PRE( module->string_data( "_type" ) != "S" &&
                  module->string_data( "_type" ) != "C" ) ;
   static std::string const dollar( "$" ) ;
   
   OUT->add_variable( dollar + module->string_data("_name"),
                      module->string_data("_type") ) ;
   
   
}

//----------------------------------------------------------------------
void
PEL_DocumentPattern:: load_dictionary( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: load_dictionary" ) ;
   PEL_CHECK( !DICO.empty() ) ;
   std::ifstream in( DICO.c_str() ) ;
   if( !in )
   {
      PEL_Error::object()->raise_plain(
         "*** PEL_DocumentPattern error:\n"
         "    unable to open dictionary file: "+DICO ) ;
   }
   std::string line ;
   
   while( !in.eof() )
   {
      getline( in, line ) ;
      size_t idx = line.find_first_of( "=" ) ;
      if( idx<line.length() )
      {
         std::string key = line.substr( 0, idx ) ;
         std::string val = line.substr( idx+1, line.length()-idx -1 ) ;
         size_t idx2 = key.find_last_not_of( ' ' ) ;
         if( idx2 < key.length() ) key = key.substr( 0, idx2+1 ) ;
         if( !key.empty() )
         {
            KEYS.append( key ) ;
            VALS.append( val ) ;
            if( VERBOSE ) 
            {
               PEL::out() << " Read in dictionary: " << key << " = " << val << std::endl ;
            }
         }
      }
    }  
}



//----------------------------------------------------------------------
std::string const&
PEL_DocumentPattern:: translation( std::string const& key ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DocumentPattern:: translation" ) ;
   size_t i = KEYS.index_of( key ) ;
   if( i < KEYS.size() ) return VALS(i) ;
   else return key ;
}

//nodoc----------------------------------------------------------------------
PEL_DocumentPattern::OutputInterfaceHtml::OutputInterfaceHtml ( std::string const& filename,
                                                                std::string const& title,
                                                                std::string const& stylesheet,
                                                                PEL_DocumentPattern& document_pattern )
//nodoc----------------------------------------------------------------------
      : OS( &PEL::out() )
      , FIRST_PASS( false )
      , INDEX(0)
      , DOC( document_pattern )
{
   if( !filename.empty() )
   {
      OS = new std::ofstream( filename.c_str() ) ;
      PEL_ASSERT( (*OS) ) ;
   }
   char const* preamble = \
      "body {\n"\
      "  color: black;\n"\
      "  background-color: white;\n"\
      "  font: normal 10pt Arial sans-serif;\n"\
      "}\n"\
      "h1 {\n"\
      "  text-align: left;\n"\
      "  font-size: 16pt;\n"\
      "}\n"\
      "h2 {\n"\
      "  text-align: left;\n"\
      "  font-size: 14pt;\n"\
      "}\n"\
      "h3 {\n"\
      "  text-align: left;\n"\
      "  font-size: 12pt;\n"\
      "  text-decoration: underline;\n"\
      "}\n"\
      "h4 {\n"\
      "  text-align: left;\n"\
      "  font-size: 12pt;\n"\
      "  border:solid 1px;\n"\
      "  background-color:#444;\n"\
      "  color:#FFF;\n"\
      "}\n"\
      "dd { \n"\
      "}\n"\
      "pre,code {\n"\
      "  font-family: monospace \n"\
      "}\n"\
      "a:link    {\n"\
      "  color: #22c;\n"\
      "}\n"\
      "a:visited {color: #22c;\n"\
      "}\n"\
      "a:active  {\n"\
      "  color: #22c;\n"\
      "}\n"\
      "hr {\n"\
      "  height: 1px;\n"\
      "  background-color: #447;\n"\
      "  border: 0px #447 solid;\n"\
      "}\n"\
      ".generic {\n"\
      "  font-family: monospace;\n"\
      "  font-style: italic;\n"\
      "  font-weight: bold;\n"\
      "  color: #000;\n"\
      "}\n"\
      ".mandatory {\n"\
      "  font-family: monospace;\n"\
      "  font-weight: bold;\n"\
      "  color: #000;\n"\
      "}\n"\
      ".optional {\n"\
      "  font-family: monospace;\n"\
      "  color: #111;\n"\
      "}\n"\
      ".defval:before {\n"\
      "  content : '= ';\n"\
      "}\n"\
      ".defval {\n"\
      "  color: #044;\n"\
      "}\n"\
      ".unit:before {\n"\
      "  content : '[';\n"\
      "}\n"\
      ".unit:after {\n"\
      "  content : ']';\n"\
      "}\n"\
      ".unit {\n"\
      "  font-weight: bold;\n"\
      "  color: #333 ;\n"\
      "  padding: 0.5em;\n"\
      "}\n"\
      ".type:before { \n"\
      "  content : '<';\n"\
      "}\n"\
      ".type:after {\n"\
      "  content : '>';\n"\
      "}\n"\
      ".type {\n"\
      "  font-style: italic;\n"\
      "}\n"\
      ".comment {\n"\
      "  font-style: italic;\n"\
      "}\n"\
      ".allowed {font-family:monospace;\n"\
      "  color: green\n"\
      "}\n"\
      ".condition {\n"\
      "  font-weight: bold;\n"\
      "  color: red;\n"\
      "}\n"\
      ".variable {\n"\
      "  font-family: monospace;\n"\
      "  font-style: italic;\n"\
      "  font-weight: bold;\n"\
      "  color: blue;\n"\
      "}\n";
   *(OS) << "<html>" << std::endl ;
   *(OS) << "<head>" << std::endl ;
   *(OS) << "<title>" << PEL_System::basename(title) <<"</title>" << std::endl ;
   if( stylesheet.empty() )
   {
      *(OS) << "<style type=\"text/css\">" << std::endl ;
      *(OS) << preamble << std::endl ;
      *(OS) << "</style>" << std::endl ;
   }
   else
   {
      *(OS) << "<link rel=\"stylesheet\" type=\"text/css\" href=\"" << stylesheet << "\" />" << std::endl ;
   }
   *(OS) << "</head>" << std::endl ;
   *(OS) << "<a name=top>" << std::endl ;
   *(OS) << "<H1>PELICANS data deck structure</H1>" <<
      "This document was generated by PELICANS application <b>document_pattern</b> applied on controler file \""
         << PEL_System::basename(title) << "\"." << std::endl ;
   *(OS) << "<BR>Index of documented modules can be found <a href=#index>at the end of the document</a>." << std::endl ;
   
}

//nodoc----------------------------------------------------------------------
PEL_DocumentPattern::OutputInterfaceHtml::~OutputInterfaceHtml ( void )
//nodoc----------------------------------------------------------------------
{

   hr() ;
   *(OS) << "<HR><a name=index><H1>Index of documented modules</H1></a>" << std::endl ;
   *(OS) << "<ul>" << std::endl ;
   INDEX.sort() ;
   for( size_t i=0 ; i<INDEX.size() ; i++ )
   {
      *(OS) << "<li>" << INDEX(i) << "</li>" << std::endl ;
   }
   *(OS) << "</ul>" << std::endl ;
   
   *(OS) << "</body>" << std::endl ;
   *(OS) << "</html>" << std::endl ;
}


//nodoc----------------------------------------------------------------------
size_t
PEL_DocumentPattern::OutputInterfaceHtml:: wished_pass_number( void ) const
//nodoc----------------------------------------------------------------------
{
   return 2 ;
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::hr(void)
//nodoc----------------------------------------------------------------------
{
   (*OS)<< "<table><tr><td width=\"100%\"><HR><td><a href=\"#top\">[top]</a><td><a href=\"#index\">[index]</a><td><a href=\"javascript:history.back();\">[back]</a><td><a href=\"javascript:history.forward();\">[next]</a></table>\n" ;
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::start_module(
   std::string const& name,
   std::string const& access,
   std::string const& help,
   std::string const& url,
   int index,
   stringVector const& plotter )
//nodoc----------------------------------------------------------------------
{
   hr() ;
   (*OS) << std::endl << "<a name=\"" << index << "\"><H2>Module " << name <<"</H2></a>" << std::endl ;
   if( !help.empty() ) (*OS) << "Definition: <span class=comment>" << convert_to_html(DOC.translation(help)) << "</span>\n" ;
   if( !url.empty() ) (*OS) << "Documentation is available <a href=\""<< url << "\">here</a>" << std::endl ;
   std::ostringstream os ;
   os<<"<a title=" << name << " href=#" << index << ">"<<name<<"</a>";
   if( plotter.size() > 0 )
   {
      (*OS) << "This structure can be displayed by pelguis plotter of class \"" << plotter(0) << "\" with following equivalences between arguments: " <<  std::endl ;
         
      (*OS) << "<table border=\"1\">" ;
            
      (*OS) << "<tr><th>" << plotter(0) << " argument </th>" ;
      (*OS)<<"<th>" << name << " module argument</th>";
      (*OS) << " </tr>\n <span class=allowed>" ;
      bool first = true ;
      for(size_t i=1 ; i<plotter.size() ; i++ )
      {
         size_t idx = plotter(i).find("=") ;
         if( idx >= plotter(i).length() )
         {
            PEL_Error::object()->display_info(
               "*** PEL_DocumentPattern warning:\n"
               "    syntax error in argument equivalence \""+plotter(i)+"\"\n"
               "    of plotter associated to module \""+name+"\"" ) ;
         }
         else
         {
            if( !first ) (*OS) << "<tr>" ;
            first = false ;
            std::string symb = plotter(i).substr(0,idx) ;
            std::string val = plotter(i).substr(idx+1) ;
               
            (*OS) << "<td class=allowed>" << symb << "</td>" ;
            (*OS) << "<td class=allowed>" << val << "</td>"  ;
            (*OS) << "</tr>" ;
         }
      }
      (*OS) << "</table></span>" ;
         
   }
    
   INDEX.append( os.str() ) ;
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::add_module(
   std::string const& name,
   std::string const& access,
   int index )
//nodoc----------------------------------------------------------------------
{
   if( PASS==0 )
   {
      if( FIRST_PASS )
      {
         (*OS) << "<H3>Syntax:</H3><ul>" << std::endl ;
         FIRST_PASS = false ;
      }
      (*OS) << "<li><a href=\"#" << index << "\">" ;
      std::string disp_access = ( access=="optional" || access=="generic" ? access : "mandatory" ) ;
      (*OS) << "<span class=type>Module</span> <span class=" << disp_access << ">" << name << "</span></a></li>" << std::endl;
     
   }
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::finish_module(
   std::string const& name,
   std::string const& access )
//nodoc----------------------------------------------------------------------
{
   (*OS)<<"</dl>"<<std::endl ;   
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::start_list_of_conditions( void )
//nodoc----------------------------------------------------------------------
{

   (*OS) << "<H3>List of conditional structures:</H3><table style='margin-left:3em'>" << std::endl ;   
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::finish_list_of_conditions( void )
//nodoc----------------------------------------------------------------------
{
   (*OS) << "</table>" << std::endl ;   
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::start_condition(
   std::string const& condition )
//nodoc----------------------------------------------------------------------
{

   (*OS) << "<tr><td><H4>On condition: " << condition << " </H4>" << std::endl ;   
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::add_help_on_condition(
   std::string const& help, std::string const& url )
//nodoc----------------------------------------------------------------------
{

   if( PASS==0 )
   {
      if( !help.empty() ) (*OS) << "Alternative description: " << convert_to_html(DOC.translation(help)) << std::endl ;   
      if( !url.empty() ) (*OS) << "Documentation is available <a href=\"" << url << "\">here</a>" << std::endl ;
   }
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::finish_condition(
   std::string const& condition )
//nodoc----------------------------------------------------------------------
{
   (*OS) << "</dd></dt>" << std::endl ;   
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::add_entry(
   std::string const& name,
   std::string const& type,
   std::string const& access,
   std::string const& def,
   stringVector const& in,
   stringVector const& help_in,
   std::string const& sel_in,
   std::string const& where,
   std::string const& test,
   std::string const& unit,
   std::string const& help,
   std::string const& formula,
   bool unique )
//nodoc----------------------------------------------------------------------
{
   bool has_doc = !help.empty()
      || in.size()>0
      || !sel_in.empty()
      || !test.empty() ;
   
   if( FIRST_PASS )
   {
      if(PASS==0)
      {
         (*OS) << "<H3>Syntax:</H3><ul>" << std::endl ;
         FIRST_PASS = false ;
      }
      if( PASS==1 && has_doc )
      {
         (*OS)<<"<H3>Description:</H3><dl>"<<std::endl ;
         FIRST_PASS = false ;
      }
   }
   
   if( PASS==0 )
   {
      (*OS) << "<li><span class=type>" << type << "</span> " ;
      std::string disp_access = ( access=="optional" || access=="generic" ? access : "mandatory" ) ;
      (*OS) << "<span class=" << disp_access << ">" << name << "</span>" ;
      
      if( !def.empty() )  (*OS) << " <span class=defval>" << def <<"</span>" ;
      if( !unit.empty() )  (*OS) << " <span class=unit>" << unit <<"</span>" ;
      (*OS)<<"</li>"<<std::endl;
   }
   
   if( PASS==1 && has_doc )
   {
      (*OS) << "<dt>" << name  ;
      if( !help.empty() ) (*OS) << "<dd>Definition: <span class=comment>" << convert_to_html(DOC.translation(help)) << "</span></dd>" ;
      if( in.size()>0 || ( !sel_in.empty() ) || ( !test.empty() ) )
      {
         (*OS) << "<dd>Validity:" << std::endl ;
         if( in.size()>0 )
         {
            (*OS) << "<table border=\"1\">" ;
            bool has_comment = help_in.size() >0 ;
            
            (*OS) << "<tr><th>Valid choices</th>" ;
            if(has_comment) (*OS)<<"<th>Comment</th>";
            (*OS) << " </tr>\n <span class=allowed>" ;
            for(size_t i=0 ; i<in.size() ; i++ )
            {
               if( i>0 ) (*OS) << "<tr>" ;
               (*OS) << "<td class=allowed>" << in(i) << "</td>" ;
               if( help_in.size() > i ) (*OS) << "<td class=allowed>" << convert_to_html(DOC.translation(help_in(i))) << "</td>"  ;
               (*OS) << "</tr>" ;
            }
            (*OS) << "</table></span>" ;
         }
         if( !sel_in.empty() )
         {
            (*OS) << " <span class=condition>SELECT (" << sel_in ;
      
            if( !where.empty() ) (*OS) << ") and (" << where ;
            (*OS) <<")</span>";
         }
         if( !test.empty() )
         {
            (*OS) << " <span class=condition>" << test << "</span>" ;
         }
         if( unique )
         {
            (*OS) << " <span class=condition>uniqueness : all values must be distincts</span>" ;
         }
         (*OS)<<"</dd>"<<std::endl;
         
      }
      (*OS) << "</dt>" << std::endl  ;
   }
}

//nodoc----------------------------------------------------------------------
std::string const&
PEL_DocumentPattern::OutputInterfaceHtml::compact_type( std::string const& a_type ) 
//nodoc----------------------------------------------------------------------
{
   static std::map<std::string,std::string> compact ;
   
   if( compact.size()==0 ) 
   {
      compact["Double"] = "D" ;
      compact["Int"] = "I" ;
      compact["Bool"] = "B" ;
      compact["String"] = "S" ;
      compact["DoubleVector"] = "D[]" ;
      compact["IntVector"] = "I[]" ;
      compact["BoolVector"] = "B[]" ;
      compact["StringVector"] = "S[]" ;
      compact["DoubleArray2D"] = "D[[]]" ;
      compact["IntArray2D"] = "I[[]]" ;      
   }
   return compact[a_type] ;
}

      
//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml::add_variable(
   std::string const& name,
   std::string const& type )
//nodoc----------------------------------------------------------------------
{
   if( PASS==0 )
   {
      (*OS) <<  "<li>Running time variable: <span class=variable>"
            << name << "</span></li>" << std::endl ;
   }
}

//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml:: set_pass( size_t pass ) 
//nodoc----------------------------------------------------------------------
{
   PEL_CHECK_PRE( pass<=wished_pass_number() ) ;
   PASS=pass;
   FIRST_PASS = true ;
   
}


//nodoc----------------------------------------------------------------------
void
PEL_DocumentPattern::OutputInterfaceHtml:: finish_pass( void ) 
//nodoc----------------------------------------------------------------------
{
   if( !FIRST_PASS )
   {
      if( PASS==0 )
      {
         (*OS) << "</ul>" << std::endl ;
      }
      if( PASS==1 )
      {
         (*OS) << "</dl>" << std::endl ;
      }
   }
}

//nodoc----------------------------------------------------------------------
std::string
PEL_DocumentPattern::OutputInterfaceHtml:: convert_to_html( std::string const& str ) 
//nodoc----------------------------------------------------------------------
{
   std::string result ;
   // Convert UTF-16 \u escaped to html
   for( size_t i=0 ; i<str.length() ; i++ )
   {
      if( i<str.length()-6 && str[i]=='\\' && str[i+1]=='u' )
      {
         result += "&#x" ;
         result += str.substr( i+2, 4 ) ;
         result += ";" ;
         i += 5 ;
      }
      
      else
      {
         result += str[i] ;
      }
   }
   return result ;
}
