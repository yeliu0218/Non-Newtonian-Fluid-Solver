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

#include <PEL_ExternalAPI.hh>

#include <PEL_Error.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_assertions.hh>

#include <stringVector.hh>

#include <iostream>

//----------------------------------------------------------------------
PEL_ExternalAPI:: PEL_ExternalAPI( std::string const& a_name,
                                   size_t a_priority_level )
//----------------------------------------------------------------------
   : PEL_Object( 0 )
   , MY_NAME( a_name )
   , MY_PRIORITY( a_priority_level )
{
   PEL_LABEL( "PEL_ExternalAPI:: PEL_ExternalAPI" ) ;

   if( a_priority_level>9 )
   {
      PEL_Error::object()->raise_plain(
         "PEL_ExternalAPI error :\n"
         "   external API of name : \n"+a_name+"\"\n"
         "   is defined with a priority level greater than 9" ) ;
   }
   
   plugins_map()->register_item( a_name, this ) ;
   plugins_names().append( a_name ) ;
}

//----------------------------------------------------------------------
PEL_ExternalAPI:: ~PEL_ExternalAPI( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_ExternalAPI:: initialize_all_APIs( int& argc, char**& argv )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExternalAPI:: initialize_all_APIs" ) ;

   size_t const nb_APIs = plugins_names().size() ;
   for( size_t level=9 ; level<10 ; --level )
   {
      for( size_t i=0 ; i<nb_APIs ; ++i )
      {
         std::string const& name = plugins_names()(i) ;
         PEL_ExternalAPI* api =
            static_cast<PEL_ExternalAPI*>( plugins_map()->item( name ) ) ;
         if( api->MY_PRIORITY==level )
         {
            api->initialize( argc, argv ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
PEL_ExternalAPI:: terminate_all_APIs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExternalAPI:: terminate_all_APIs" ) ;
   
   size_t const nb_APIs = plugins_names().size() ;
   for( size_t level=0 ; level<10 ; ++level )
   {
      for( size_t i=nb_APIs-1 ; i<nb_APIs ; --i )
      {
         std::string const& name = plugins_names()(i) ;
         if( plugins_map()->has( name ) )
         {
            PEL_ExternalAPI* api =
               static_cast<PEL_ExternalAPI*>( plugins_map()->item( name ) ) ;
            if( api->MY_PRIORITY==level )
            {
               plugins_map()->unregister_item( name ) ;
               api->destroy() ;
            }
         }
      }
   }
   plugins_map()->destroy() ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_ExternalAPI:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
         PEL_ObjectRegister::create( 0, "PEL_ExternalAPI descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector&
PEL_ExternalAPI:: plugins_names( void )
//----------------------------------------------------------------------
{
   static stringVector result(0) ;
   return( result ) ;
}
