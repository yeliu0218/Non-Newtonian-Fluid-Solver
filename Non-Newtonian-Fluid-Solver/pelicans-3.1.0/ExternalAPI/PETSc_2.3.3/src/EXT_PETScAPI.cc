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

#include <EXT_PETScAPI.hh>

#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Root.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Variable.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <sstream>
#include <iostream>

#if( ! ( PETSC_VERSION_MAJOR    == 2 && \
         PETSC_VERSION_MINOR    == 3 && \
         PETSC_VERSION_SUBMINOR == 3 ) )
 "Bad version of PETSC ( Version 2.3.3 should be used )" ;
#endif

EXT_PETScAPI* EXT_PETScAPI:: SINGLETON = new EXT_PETScAPI() ;
PEL_Timer* EXT_PETScAPI:: timer = 0 ;

//----------------------------------------------------------------------
EXT_PETScAPI:: EXT_PETScAPI( void )
//----------------------------------------------------------------------
   : PEL_ExternalAPI( "EXT_PETScAPI", 0 )
{
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_PETSc" ), val ) ;

   std::ostringstream ver ;
   ver << PETSC_VERSION_MAJOR << "." 
       << PETSC_VERSION_MINOR << "." 
       << PETSC_VERSION_SUBMINOR ;   
   PEL_String* rev = PEL_String::create( 0, ver.str() ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "SS_PETSc_REV" ), rev ) ;
}

//----------------------------------------------------------------------
EXT_PETScAPI:: ~EXT_PETScAPI( void )
//----------------------------------------------------------------------
{
   PetscFinalize() ;
}

//----------------------------------------------------------------------
void
EXT_PETScAPI:: initialize( int& argc, char **& argv )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScAPI:: initialize" ) ;
   const char* help = "EXT_PETScAPI:: initialize" ;

   char** my_argv = new char* [ argc ] ;
   int my_argc=0 ;
   char** new_argv = new char* [ argc ] ;
   int new_argc = 0 ;
   
   if( argc>0 )
   {
      my_argv[my_argc++] = argv[0] ;
   }
   
   for( int i=0 ; i<argc ; ++i )
   {
      std::string str = argv[i] ;
      if(  str == "-Xpetsc" )
      {
         i++ ;
         PEL_ASSERT( i < argc ) ;
         std::string str1 = argv[i] ;
         if( str1 == "-trace" )
         {
            timer = PEL_Timer::create( PEL_Root::object() ) ;
         }
         else
         {
            my_argv[my_argc++] = argv[i] ;
         }  
      }
      else
      {
	 new_argv[new_argc++] = argv[i] ;
      }
   }
   if( my_argc > 1 || timer!=0 )
   {
      argc = new_argc ;
      argv = new_argv ;
   }
   else
   {
      delete [] new_argv ;
   }
   if( my_argc > 1 )
   {
      std::cout << "PETSc init : " << std::endl ;
      for( int i=1 ; i<my_argc ; ++i ) 
      {
         std::cout << "    " << my_argv[i] << std::endl ;
      }
   }
   PetscInitialize( &my_argc, &my_argv, PETSC_NULL, help ) ;
   delete [] my_argv ;
}


//----------------------------------------------------------------------
bool
EXT_PETScAPI:: parse_options( PEL_ModuleExplorer const* exp,
                              bool verbose )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScAPI:: parse_options" ) ;
   bool result =  exp->has_module( "options" ) ;
   
   if( result )
   {
      PEL_ModuleExplorer * sexp = exp->create_subexplorer( 0, "options" )  ;
      for( sexp->start_entry_iterator() ;
           sexp->is_valid_entry()   ;
           sexp->go_next_entry() )
      {
         std::string const& name = "-"+sexp->keyword() ;
         PEL_Data * data = sexp->data( 0 ) ;
         std::string const& val = data->to_string() ;
         if( verbose )
            std::cout << "EXT_PETScAPI - setting option: " 
                      << name << " " << val << std::endl ;
         
         if( val.empty() )
         {
            PetscOptionsSetValue( name.c_str(), PETSC_NULL ) ;
         }
         else
         {
            PetscOptionsSetValue( name.c_str(), val.c_str() ) ;
         }
         
         data->destroy() ;
      }
      sexp->destroy() ;
   }
   return result ;
}


//----------------------------------------------------------------------
void
EXT_PETScAPI:: going_to_do(  char const* action )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScAPI:: going_to_do" ) ;
   if( timer!=0 )
   {
      PEL::out()<<"["<<PEL_System::used_memory()/1024/1024<<"Mo]"
                <<"Start -> "<<action<<std::endl ;
      if( timer->is_running() )
      {
         timer->stop() ;
         PEL::out() << "*** Warning : imbricated times" << std::endl
                    << "    PETSc timer stopped at cumulative time : (s) " << timer->time() << std::endl
                    << "    Restart PETSc timer..." << std::endl ;
      }
      timer->start() ;
   }
}

//----------------------------------------------------------------------
void
EXT_PETScAPI:: verify( char const * action, int result )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScAPI:: verify" ) ;
   if( timer!=0 )
   {
      PEL::out() << "[" << PEL_System::used_memory()/1024/1024<<"Mo]"
                 << action << " <- End." ;
      if( timer->is_running() )
      {
         timer->stop() ;
         PEL::out() << " PETSc cumulative time : (s) " << timer->time() ;
      }
      PEL::out() << std::endl << std::endl ;
   }
   if( result!=0 )
   {
      std::string mess = "Internal Petsc error encountered in " ;
      mess += action ;
      PEL_Error::object()->raise_internal( mess ) ;
   }
}

