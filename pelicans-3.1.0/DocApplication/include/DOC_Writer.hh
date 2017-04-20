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

#ifndef DOC_Writer_HH
#define DOC_Writer_HH


#include <fstream>
#include <string>

#include <PEL_Object.hh>
#include <PEL_List.hh>

#include <DOC_Tools.hh>

class DOC_Class ;
class DOC_ClassItem ;
class DOC_Method ;
class DOC_Package ;

// Abstract class realizing documentation.
// For each class, a new instance is created and method `process_documentation'
// is then called.
// When all classes are processed, a last instance is created and method
// `finalize' is called.
// FRAMEWORK INSTANTIATION :
// A writer inherited class must implement DOC_Writer pure virtual method.
// It can use `out()' stream to process output.
// It can also use DOC_Package specific purpose method to parse packaging tree.

class PEL_EXPORT DOC_Writer : public PEL_Object
{
   public://----------------------------------------------------------

   //-- Creation
      enum Format { DocumentationHtml,
                    DocumentationLaTex,
                    DocumentationJS }  ;

      static DOC_Writer* create( PEL_Object* a_owner, string const& name ) ;

      // Parse input datafile.
      static void parse_explorer( PEL_ModuleExplorer const* exp ) ;
      
   //-- Documentation processing
         
      // convert `text_initial' string using `marqueur' indicator
      virtual string underlined_text( size_t marqueur,
                                  string const& text_initial ) const = 0 ;      
      // cross-referenced string for a class `cl' and an optional
      // item of it `label'
      // `display' is the string which must be seen by reader
      virtual string reference( DOC_Class const* cl,
                                string const& label,
                                string const& display ) const = 0  ;

       // Is croos-referencing enable ?
      virtual bool source_reference( void ) const = 0 ;
      
      // Process `classe' documentation.
      virtual void process_documentation( DOC_Class const* classe ) = 0 ;
      
      // Finalize documentation process by realizing index, if any..
      virtual void finalize( void ) = 0 ;
      
   protected://-------------------------------------------------------

   //-- Package parsing(8)
      
      void doc_class_from( DOC_Package const* pack,
                                         int level ) ;
      // Parse package structure.
      // As instance, for further package structure :
      // MASTER
      //   PACK1
      //    class1
      //      inherit_class1
      //    PACK2
      //      class2
      //  following methods will be called :
      //  `::doc_class_from'( MASTER, 0 )
      //  `::doc_class_from'( PACK1, 1 )
      //  `::disp_classes_tree'( class1, 0, <class1,inherit_class1,class2> )
      //  `::disp_classes_tree'( inherit_class1, 1, <inherit_class1,class2> )
      //  `::process_documentation_package'( PACK2, 2 )
      //  `::disp_classes_tree'( class2, 0, <class1,inherit_class1,class2> ).
      virtual void process_documentation_package( DOC_Package const* pack,
                                      size_t level=0 ) ;
      virtual void disp_classes_tree( DOC_Class const* classe,
                                          int niveau,
                                          PEL_List* liste_classes ) ;
      
      
  //-- Creation
      DOC_Writer( PEL_Object* a_owner, string const& chaine ) ;
      virtual ~DOC_Writer( void ) ;

  //-- Tools(9.0)
      // output stream
      ostream& out( void ) const ;
      // convert integer to string
      string itoa( int i ) const ;
      // Copy $PELICANSHOME/`source' to `target'.
      void copy( string const& source,
                   string const& target ) ;
      
      
   private://---------------------------------------------------------

      DOC_Writer& operator=( DOC_Writer const& red ) ;
      DOC_Writer(DOC_Writer  const& red ) ;
      DOC_Writer( void ) ;
      std::ofstream leStream ;
      
      static Format format_de_sortie ;
      static bool interactif ;
      
} ;

#endif // DOC_Writer_HH
