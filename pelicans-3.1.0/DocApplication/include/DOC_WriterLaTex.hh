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

#ifndef DOC_WriterLaTex_HH
#define DOC_WriterLaTex_HH


#include <DOC_Writer.hh>
#include <PEL_List.hh>
#include <string>
#include <stringVector.hh>

// DOC_Writer to format laTex.

class PEL_EXPORT DOC_WriterLaTex : public DOC_Writer
{
   public://----------------------------------------------------------

   //-- Creation
      static DOC_WriterLaTex* create( PEL_Object* a_owner,
                                   std::string const& name ) ;
   //-- Documentation processing
      virtual std::string underlined_text( size_t marqueur,
                                   std::string const& text_initial ) const ;
      virtual std::string reference( DOC_Class const* cl,
                                     std::string const& label,
                                     std::string const& display ) const ;
      
      virtual void finalize( void ) ;
      virtual bool source_reference( void ) const ;
      virtual void process_documentation( DOC_Class const* classe ) ;
      
     
   protected://-------------------------------------------------------
   private://---------------------------------------------------------      
   //-- Package parsing
      void process_documentation( DOC_Class const* classeAProcess_Documentationr,
                      DOC_ClassItem * elem )  ;
      virtual void process_documentation_package( DOC_Package const* pack,
                                      size_t level=0 ) ;
      virtual void disp_classes_tree( DOC_Class const* classe,
                                          int niveau,
                                          PEL_List* liste_classes ) ;
      DOC_WriterLaTex( PEL_Object* a_owner, std::string const& name ) ;
      DOC_WriterLaTex( void ) ;
      ~DOC_WriterLaTex( void ) ;
      
      DOC_WriterLaTex& operator=( DOC_WriterLaTex const& r ) ;
      DOC_WriterLaTex( DOC_WriterLaTex const& r ) ;
      void comment( std::string const& com ) ;
      void ecrit_code( std::string const& com ) ;
      void ecrit( std::string const& com ) ;
      
      static std::string true_type( std::string const& com,
                               bool retour_a_la_line_number = true ) ;
      static stringVector menu ;
      std::string file ;
      std::string para ;
      bool file_principal ;
} ;

#endif // DOC_WriterLaTex_HH
