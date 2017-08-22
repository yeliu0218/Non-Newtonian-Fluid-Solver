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

#ifndef DOC_WriterHtml_HH
#define DOC_WriterHtml_HH


#include <DOC_Package.hh>
#include <DOC_Writer.hh>
#include <PEL_List.hh>
#include <string>
#include <stringVector.hh>

// DOC_Writer to format HTML.
// To be published, comment must exhib this word : PUBLISHED.

class PEL_EXPORT DOC_WriterHtml : public DOC_Writer
{
   public://----------------------------------------------------------

   //-- Creation
      static DOC_WriterHtml* create( PEL_Object* a_owner,
                                 std::string const& str ) ;
      
   //-- Documentation processing
      virtual void finalize( void ) ;
      virtual string underlined_text( size_t marqueur,
                                   std::string const& text_initial ) const ;
      virtual std::string reference( DOC_Class const* cl,
                                std::string const& label,
                                std::string const& display ) const ;
      virtual bool source_reference( void ) const ;
      virtual void process_documentation( DOC_Class const* classe ) ;
      
   //-- Auto-test
      void dummy_method( void ) const ;
      
   protected://-------------------------------------------------------
      
   private://---------------------------------------------------------
   //-- Package parsing
      virtual void process_documentation_package( DOC_Package const* pack,
                                      size_t level=0 ) ;
      virtual void disp_classes_tree( DOC_Class const* classe,
                                          int niveau,
                                          PEL_List* list_classes ) ;
      
      virtual void start_package_doc( DOC_Package const* pack,
                                      size_t level ) ;
      DOC_WriterHtml( PEL_Object* a_owner, std::string const& str ) ;
      DOC_WriterHtml( void ) ;
      ~DOC_WriterHtml( void ) ;
      DOC_WriterHtml& operator=( DOC_WriterHtml const& r ) ;
      DOC_WriterHtml( DOC_WriterHtml const& r ) ;
      
      void process_documentation( DOC_Class const* classeAProcess_Documentationr,
                           DOC_ClassItem * elem ) ;

		//-- HTML Elements
		void list( std::string const& id ) ;
      void new_element( void ) ;
      void write_label( std::string const& label ) ;
      void hruler( void ) ;
      void start_category( std::string const& name, 
                            std::string const& protection ) ;
      void title( std::string const& tit, int niveau ) ;
      void comment( std::string const& com ) ;
      void write_code( std::string const& com ) ;
      void write_prototype( std::string const& com ) ;
      void write( std::string const& com ) ;
      static std::string expand( std::string const& com, std::string const& html_cmd ) ;
      static void process_webview( std::string const src,
                                   std::string const dest,
                                   stringVector const& item,
                                   stringVector const& refs,
                                   size_t n ) ;
      static void head( std::ostream& ovf, std::string const& title, bool set_body ) ;
      static void navbar( std::ostream& os,
                          stringVector const& items,
                          stringVector const& refs,
                          size_t n ) ;
      
      enum some_fonts { font_text = 4 , font_title_list = 3 } ;
      
      static stringVector menu ;
      std::string file ;
      bool file_principal ;
      size_t package_level ;
      size_t class_level ;
} ;

#ifndef OUTLINE
#include <DOC_WriterHtml.icc>
#endif

#endif // DOC_WriterHtml_HH
