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

#include <AP_Greeting.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>
#include <fstream>

using std::cout ;
using std::endl ;

AP_Greeting const* AP_Greeting::PROTOTYPE = new AP_Greeting() ;

//----------------------------------------------------------------------------
AP_Greeting:: AP_Greeting( void )
//----------------------------------------------------------------------------
   : PEL_Application( "AP_Greeting" )
{
}

//----------------------------------------------------------------------------
AP_Greeting:: ~AP_Greeting( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
AP_Greeting*
AP_Greeting:: create_replica( PEL_Object* a_owner, 
                              PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_Greeting:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_Greeting* result = new AP_Greeting( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
AP_Greeting:: AP_Greeting( PEL_Object* a_owner, 
                           PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , NAME( exp->string_data( "name" ) )
   , HAS_CARD( false )
{
   if( exp->has_module( "greeting_card" ) )
   {
      PEL_ModuleExplorer const* se = 
                         exp->create_subexplorer( 0, "greeting_card" ) ;
      HAS_CARD = true ;
      INST = se->string_data( "institution" ) ;
      TOWN = se->string_data( "town" ) ;
      ZIP  = se->int_data( "zip" ) ;
      se->destroy() ;
   }
}

//----------------------------------------------------------------------------
void
AP_Greeting:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_Greeting:: run" ) ;
   
   std::cout << endl << "Hello " << NAME << std::endl ;
   if( HAS_CARD )
   {
      std::string ff = NAME + ".card" ;
      std::ofstream os( ff.c_str() ) ;
      if( !os )
         PEL_Error::object()->raise_plain( "file opening failed" ) ; 
      std::cout << "Editing greeting card: " << ff << std::endl ;
      os << "TO: " << NAME << std::endl ;
      os << "    " << INST << std::endl ;
      os << "    " << ZIP  << " " << TOWN << std::endl ;
      os << "THANK YOU!" << std::endl ;
      os.close() ;
   }
}
