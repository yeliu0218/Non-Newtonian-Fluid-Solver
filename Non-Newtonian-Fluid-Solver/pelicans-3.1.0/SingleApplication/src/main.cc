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

#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <stringVector.hh>

#include <iostream>
#include <fstream>
#include <string>

using std::cout ; 
using std::endl ;
using std::ofstream ;
using std::string ;

//---------------------------------------------------------------------------
void print_greeting_card( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "greeting" ) ;

   string name = ee->string_data( "name" ) ;
   cout << endl << "Hello " << name << endl ;

   string ff = name + ".card" ;
   ofstream os( ff.c_str() ) ;
   if( !os )
   {
      PEL_Error::object()->raise_plain( "Unable to open \"greeting.txt\"" ) ; 
   }

   PEL_ModuleExplorer const* se = ee->create_subexplorer( ee, "id_card" ) ;
   cout << "Editing your id card into file : " << ff << endl ;
   os << "Institution : " << se->string_data( "institution" ) << endl ;
   os << "Name  : " << name << endl ;
   os << "Town  : " << se->string_data( "town" ) << endl ;
   os << "Zip   : " << se->int_data( "zip" ) << endl ;
   os << "email : " << se->string_data( "email" ) << endl ;
   os.close() ;
   
   ee->destroy() ;
}

//---------------------------------------------------------------------------
int main( int argc, char* argv[] )
//---------------------------------------------------------------------------
{
   stringVector args( 0 ) ;

   PEL_Exec::initialize( argc, argv, args ) ;

   try
   {
      if( args.size() == 1 )
      {
         string fname = args( 0 ) ;

         PEL_Module const* mm = 
                     PEL_Module::create( PEL_Root::object(), "MAIN", fname ) ;

         PEL_ModuleExplorer* exp = PEL_ModuleExplorer::create( 0, mm ) ;

         print_greeting_card( exp ) ;

         exp->destroy() ;
      }
      else
      {
         cout << "invalid command line" << endl ;
         PEL_Exec::set_exit_code( 1 ) ;
      }
   }
   catch( PEL_Exceptions::Error )
   {
      PEL_Exec::set_exit_code( 1 ) ;
   }

   PEL_Exec::terminate() ;
   return( PEL_Exec::exit_code() ) ;
}
