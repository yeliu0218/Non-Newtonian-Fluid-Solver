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

#include <PEL_System_TEST.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_Context.hh>
#include <PEL_Map.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intVector.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <stringVector.hh>


#include <fstream>
#include <iostream>

using std::string ;
using std::cerr ;
using std::endl ;

//-------------------------------------------------------------------------
PEL_System_TEST*
PEL_System_TEST::unique_instance = new PEL_System_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
PEL_System_TEST:: PEL_System_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_System", "PEL_System_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_System_TEST:: ~PEL_System_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_System_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_System_TEST:: process_one_test" ) ;
   
   if( exp->has_module( "path_comparison" ) )
   {
      string init_dir = PEL_System::working_directory() ;
      
      PEL_ModuleExplorer* se = 
         exp->create_subexplorer( 0, "path_comparison" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module()  )
      {
         PEL_ModuleExplorer const* eee = se->create_subexplorer( 0 ) ;
         string nn1 = eee->string_data( "path_1" ) ;
         string nn2 = eee->string_data( "path_2" ) ;
         
         bool ok = PEL_System::changedir( nn1 ) ;
         if( !ok ) display_chdir_failed( nn1 ) ;
         string wd1 = PEL_System::working_directory() ;
         
         ok = ok && PEL_System::changedir( init_dir ) ;
         if( !ok ) display_chdir_failed( init_dir ) ;

         ok = ok && PEL_System::changedir( nn2 ) ;
         if( !ok ) display_chdir_failed( nn2 ) ;
         string wd2 = PEL_System::working_directory() ;
                  
         ok = ok && ( wd1 == wd2 ) ;
         if( !ok )
         {
            PEL::out() << "cd " << nn1 << " & pwd : " << wd1 << endl ;
            PEL::out() << "cd " << nn2 << " & pwd : " << wd2 << endl ;
         }
         
         ok = ok && PEL_System::changedir( init_dir ) ;
         if( !ok ) display_chdir_failed( init_dir ) ;
         
         notify_one_test_result( eee->name() +
                                 ": changedir & working_directory", ok ) ;
         eee->destroy() ;
      }
      se->destroy() ;
   }
   if( exp->has_module( "run" ) )
   {
      PEL_ModuleExplorer* se = 
         exp->create_subexplorer( 0, "run" ) ;
      
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module()  )
      {
         PEL_ModuleExplorer const* eee = se->create_subexplorer( se ) ;
         string cmd = eee->string_data( "command" ) ;
         string dir = eee->string_data( "directory" ) ;
         string output = eee->string_data( "output" ) ;

         bool ok = PEL_System::mkdir( dir ) ;
         notify_one_test_result( eee->name() +
                                 ": mkdir "+dir, ok ) ;
         if( ok ) 
         {
            
            ok = PEL_System::run( cmd, dir, output ) == 0 ;
            notify_one_test_result( eee->name() +
                                    ": run "+cmd, ok ) ;
            if( ok ) 
            {
               string file = dir + PEL_System::path_name_separator()
                  + output ;
               ok = PEL_System::can_read( file ) ;
               notify_one_test_result( eee->name() +
                                       ": verify "+output, ok ) ;
            }
         }    
         
      }
      
      se->destroy() ;
   }

   static bool first = true ;
   if( first )
   {
      double t0 = PEL_System::epoch_time() ;
      PEL_System::sleep(5 * 1000) ;
      double t1 = PEL_System::epoch_time() ;
      bool ok = t1 - t0 >= 5 ;
      if( !ok ) out() << " t0 = " << t0 << "t1 = " << t1 << std::endl ;
      
      notify_one_test_result("epoch_time & sleep ", ok ) ;
      first = false ;

      std::string sep ;
      sep += PEL_System::path_name_separator() ;
      
      std::string path = sep + sep + "toto" ;
      ok = PEL_System::absolute_path(path)==sep+"toto" ;
      if( !ok ) out() << "sep = " << sep << " path = " << path <<
         " PEL_System::absolute_path(path) = " << PEL_System::absolute_path(path) << std::endl ;
      
      notify_one_test_result(" absolute_path ", ok ) ;

      std::ofstream toto("t o(to") ;
      toto << "titi" ;
      toto.close() ;

      ok = PEL_System::copy( "t o(to", "toto" ) ;
      if(ok) 
      {   
         std::ifstream toto_in("toto") ;
         std::string line ;
         toto_in >> line ;
         ok = ok && line=="titi" ;
         if( !ok ) out() << "line = " << line << std::endl ;
         toto_in.close() ;
      }
      
      notify_one_test_result(" copy ", ok ) ;

      std::string bad = sep + "toto" + sep + "trucmuch" ;
      PEL::out() << "bad file : " << bad << std::endl ;
      
      ok = PEL_System::can_write( "toto" ) &&
         !PEL_System::can_write( bad ) ;
      
      notify_one_test_result(" can_write", ok ) ;
      
   }
   
}

//---------------------------------------------------------------------------
void
PEL_System_TEST:: display_chdir_failed( std::string const& dir )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_System_TEST:: display_chdir_failed" ) ;
   
   PEL::out() << "unable to change directory to :" << dir << endl ;
}


