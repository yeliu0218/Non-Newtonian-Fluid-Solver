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

#include <DOC_WebView.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>

#include <webview>

DOC_WebView const* DOC_WebView::PROTOTYPE = new DOC_WebView() ;

//----------------------------------------------------------------------------
DOC_WebView:: DOC_WebView( void )
//----------------------------------------------------------------------------
   : PEL_Application( "webview" )
   , MY_ARGS(0)
{
}

//----------------------------------------------------------------------------
DOC_WebView*
DOC_WebView:: create( PEL_Object* a_owner, 
                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WebView:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   DOC_WebView* result = new DOC_WebView( a_owner, exp ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
DOC_WebView*
DOC_WebView:: create_replica( PEL_Object* a_owner, 
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WebView:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   DOC_WebView* result = new DOC_WebView( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
DOC_WebView:: DOC_WebView( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , MY_ARGS( 0 )
{
   if( exp->has_entry( "command_line_arguments" ) )
   {
      MY_ARGS = exp->stringVector_data( "command_line_arguments" ) ;
   }
   
}

//----------------------------------------------------------------------------
DOC_WebView*
DOC_WebView:: create_replica_from_args( PEL_Object* a_owner,
                                          stringVector& args ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WebView:: create_replica_from_args" ) ;

   DOC_WebView* result = new DOC_WebView( a_owner, args ) ;

   PEL_CHECK( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
DOC_WebView:: DOC_WebView( PEL_Object* a_owner,
                             stringVector& args ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , MY_ARGS(0)
{
   re_initialize( args ) ;
}

//----------------------------------------------------------------------------
DOC_WebView:: ~DOC_WebView( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
DOC_WebView:: re_initialize( stringVector& args ) 
//----------------------------------------------------------------------------
{
   if( args.size()>0 )
   {
      MY_ARGS = args ;
      args.re_initialize(0) ;
   }
}

//----------------------------------------------------------------------------
void
DOC_WebView:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_WebView:: run" ) ;
   stringVector argv(0) ;
   argv.append( "webview" ) ;
   for( size_t i=0 ; i<MY_ARGS.size() ; i++ )
      argv.append( MY_ARGS(i) ) ;
   
   int argc = argv.size() ;
   WebView::webview_main( argc, argv ) ;
}
