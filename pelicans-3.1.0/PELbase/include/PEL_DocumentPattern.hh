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

#ifndef PEL_DocumentPattern_HH
#define PEL_DocumentPattern_HH

#include <PEL_Application.hh>
#include <string>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <fstream>
class PEL_ModuleExplorer ;
class PEL_Vector ;

/*
  Document a pattern model.
  
   COMMAND LINE :
   Syntaxes recognized in command line mode are :
   <exe> -A document_pattern [ -output_filename output_file ]
                             [ -stylesheet stylesheet.css ] 
                             [ -dictionary dico_file ] 
                             <controler filename>

   DATA DECK :
   Data files recognized are :
   
MODULE PEL_Application
   concrete_name = "document_pattern"
   pattern_filename = <controler file to document>
   [ output_filename = <name of output file> ]
   [ base_url = <basename for url resolving> ]
   [ styleshhet = <stylesheet.css> ]
   [ dictionary = <dico_file> ]
END MODULE PEL_Application
   
*/

class PEL_EXPORT PEL_DocumentPattern : public PEL_Application
{

   public: //-----------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   //-- Plug in
      
      virtual ~PEL_DocumentPattern( void ) ; 
      
      PEL_DocumentPattern( void ) ;
      
      PEL_DocumentPattern( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) ;
      
      PEL_DocumentPattern( PEL_Object* a_owner,
                           stringVector& args ) ;
      
      virtual PEL_DocumentPattern* create_replica( 
                           PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) const ;
      
      virtual PEL_DocumentPattern* create_replica_from_args( 
                           PEL_Object* a_owner,
                           stringVector& args ) const ;
      
   //-- Command line
      
      virtual void print_usage( void ) const ;
      
      virtual void print_options( void ) const ;
      
   private: //----------------------------------------------------------

      PEL_DocumentPattern( PEL_DocumentPattern const& other ) ;
      PEL_DocumentPattern& operator=( PEL_DocumentPattern const& other ) ;

      void document_module( PEL_ModuleExplorer* module,
                            int order,
                            size_t occ ) ;
      void document_entry( PEL_ModuleExplorer* module ) ;
      void document_variable( PEL_ModuleExplorer* module ) ;
      void explore( PEL_ModuleExplorer* module, std::string cond ) ;
      void explore_conditions( PEL_ModuleExplorer* module, stringVector& conditions ) ;

      
      // HTML writer
      class OutputInterfaceHtml 
      {
         public: //--------------------------------------------------------
            
            OutputInterfaceHtml( std::string const& filename,
                                 std::string const& title,
                                 std::string const& stylesheet,
                                 PEL_DocumentPattern& document_pattern ) ;
            virtual ~OutputInterfaceHtml( void ) ;
            virtual void start_module( std::string const& name,
                                       std::string const& access,
                                       std::string const& help,
                                       std::string const& url,
                                       int index,
                                       stringVector const& plotter ) ;
            virtual void finish_module( std::string const& name,
                                        std::string const& access ) ;
            virtual void start_list_of_conditions( void ) ;
            virtual void finish_list_of_conditions( void ) ;
            virtual void start_condition( std::string const& condition ) ;
            virtual void add_help_on_condition(
               std::string const& help, std::string const& url ) ;
            virtual void finish_condition( std::string const& condition ) ;
            virtual void add_module(  
               std::string const& name, std::string const& access,  int index );
             virtual void add_entry(  
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
               bool unique ) ;
             virtual void add_variable( 
               std::string const& name,
               std::string const& type ) ;
            virtual size_t wished_pass_number( void ) const ;
            virtual void set_pass( size_t pass ) ;
            virtual void finish_pass( void ) ;
            void display_conventions( void );
            
        protected: //--------------------------------------------------------
            
        private: //----------------------------------------------------------
            
            OutputInterfaceHtml( void ) ;
            OutputInterfaceHtml( OutputInterfaceHtml const& other ) ;
            OutputInterfaceHtml& operator=( OutputInterfaceHtml const& other ) ;
            static std::string const& compact_type( std::string const& a_type ) ;
            void hr( void ) ;
            static std::string convert_to_html( std::string const& str ) ;
            
            std::ostream * OS ;
            size_t PASS ;
            bool FIRST_PASS ;
            stringVector INDEX ;
            PEL_DocumentPattern const& DOC ;
      } ;

      friend class OutputInterfaceHtml ;
      int add( PEL_ModuleExplorer* module ) ;
      static size_t checksum( PEL_ModuleExplorer const* module ) ;
      void load_dictionary( void ) ;
      std::string const& translation( std::string const& key ) const ;
      
   //-- Class attribute
      
      static PEL_DocumentPattern const* PROTOTYPE ;

   //-- Attribute
      
      std::string FILENAME ;
      std::string CSS ;
      OutputInterfaceHtml * OUT ;
      PEL_ModuleExplorer const* PATTERN ;  
      stringVector NAMES ;
      PEL_Vector* TO_DOCUMENT ;
      std::string FORMAT ;
      std::string OUTPUT_FILENAME ;
      bool VERBOSE ;
      std::string BASE_URL ;
      size_t_vector CHECK_SUMS ;
      size_t PASS ;
      std::string DICO ;
      stringVector KEYS ;
      stringVector VALS ;
} ;

#endif



