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

#ifndef TOOL_WEB_VIEW_HH
#define TOOL_WEB_VIEW_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;

// Convert source file to HTML format.

class PEL_EXPORT DOC_WebView : public PEL_Application
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization(0.1)

      static DOC_WebView* create( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) ;
   
      void re_initialize( stringVector& args ) ;

   //-- Program core execution(0.2)

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~DOC_WebView( void ) ;
      DOC_WebView( DOC_WebView const& other ) ;
      DOC_WebView& operator=( DOC_WebView const& other ) ;

      DOC_WebView( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

      DOC_WebView( PEL_Object* a_owner, stringVector& args ) ;

   //-- Plug in(11.0)

      DOC_WebView( void ) ;

      virtual DOC_WebView* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

      virtual DOC_WebView* create_replica_from_args( 
                                       PEL_Object* a_owner,
				       stringVector& args ) const ;
   //-- Class attributes

      static DOC_WebView const* PROTOTYPE ;

   //-- Attributes

      stringVector MY_ARGS ;      
} ;

#endif
