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

#include <PEL_CoupledApplications.hh>

#include <PEL.hh>
#include <PEL_Context.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <stringVector.hh>
#include <size_t_vector.hh>

#include <fstream>
#include <sstream>
#include <iostream>
#include <fcntl.h>

   
PEL_CoupledApplications const* 
PEL_CoupledApplications::PROTOTYPE = new PEL_CoupledApplications() ;

stringVector PEL_CoupledApplications::CODENAMES( 0 ) ;

std::string const coupling( "coupling.ppid" ) ;
std::string const lock( "lock" ) ;
int WAIT_MILLISECOND = 10 ;
bool TRACE = false ;
int fd = -1 ;

class PEL_OneCoupledAppli : public PEL_Object
{
   public :
      
      PEL_OneCoupledAppli( PEL_Object* a_owner ) ;
     ~PEL_OneCoupledAppli( void ) ;

      std::string NAME ;
      std::string EXE ;
      std::string DATAFILE ;
      stringVector EXTRAS ;
      stringVector MPI_OPTIONS ;
      stringVector MPI_MACHINES ;
      std::string MPI_MACHINEFILE ;
      std::string MPI_TMP_FILE ;
} ;
      
//----------------------------------------------------------------------
PEL_CoupledApplications:: PEL_CoupledApplications( void )
//----------------------------------------------------------------------
   : PEL_Application( "PEL_CoupledApplications" )
   , SILENT( false )
   , APPLIS( 0 )
   , PROCESSES( 0 )
   , NB_PROCESS( 0 )
   , MPI_RUN( "" )
{
}

//----------------------------------------------------------------------
PEL_CoupledApplications* 
PEL_CoupledApplications:: create_replica( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_CoupledApplications* result = new PEL_CoupledApplications( a_owner, 
                                                                  exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_CoupledApplications:: PEL_CoupledApplications( 
                                               PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , SILENT( false )
   , APPLIS( PEL_Vector::create( this, 0 ) )
   , PROCESSES( 0 )
   , NB_PROCESS( 0 )
   , MPI_RUN( "" )
{
   PEL_LABEL( "PEL_CoupledApplications:: PEL_CoupledApplications" ) ;
   
   if( exp->has_entry( "mpirun" ) )
   {
      MPI_RUN = exp->string_data( "mpirun" ) ;
      exp->test_file( "mpirun", "read" ) ;
   }
   else
   {
      PEL_Context const* exe_ctx = PEL_Exec::execution_context() ;
      PEL_Variable const* with_mpi =  PEL_Variable::object( "BS_with_MPI" ) ;
      if( exe_ctx->has_variable( with_mpi ) &&
          exe_ctx->value( with_mpi )->to_bool() )
      {
         PEL_Variable const* mpi_run = PEL_Variable::object( "SS_MPI_RUN" ) ;
         if( !exe_ctx->has_variable( mpi_run ) )
         {
            PEL_Error::object()->raise_internal(
               "*** PEL_RunTest error:\n"
               "    variable SS_MPI_RUN not set" ) ;
         }
         MPI_RUN = exe_ctx->value( mpi_run )->to_string() ;
      }  
   }

   PEL_ModuleExplorer* it =
                       exp->create_subexplorer( 0, "list_of_coupled_codes" ) ;
   it->start_module_iterator() ;
   for( ; it->is_valid_module() ; it->go_next_module() ) 
   {
      PEL_ModuleExplorer* a_code = it->create_subexplorer(0) ;
      
      PEL_OneCoupledAppli* appli = new PEL_OneCoupledAppli( APPLIS ) ;
      APPLIS->append( appli ) ;

      appli->NAME = a_code->string_data( "name" )  ;
      
      if( a_code->has_entry( "executable" ) )
         appli->EXE = a_code->string_data( "executable" ) ;
      
      appli->DATAFILE = a_code->string_data( "datafile" ) ;
      
      if( a_code->has_entry( "extra_command" ) )
         appli->EXTRAS = a_code->stringVector_data( "extra_command" ) ;
   
      if( a_code->has_entry( "mpi_options" ) )
      {
         appli->MPI_OPTIONS = a_code->stringVector_data( "mpi_options" ) ;
         if( MPI_RUN.empty() )
         {
            PEL_Error::object()->raise_plain(
               "*** PEL_CoupledApplications error:\n"
               "    MPI options encountered whereas MPI is not detected" ) ;
         }
         if( a_code->has_entry( "mpi_machines" ) )
            appli->MPI_MACHINES = a_code->stringVector_data( "mpi_machines" ) ;
         if( a_code->has_entry( "mpi_machinefile" ) )
            appli->MPI_MACHINEFILE = a_code->string_data( "mpi_machinefile" ) ;
         if( a_code->has_entry( "mpi_machines" ) &&
             a_code->has_entry( "mpi_machinefile" ) )
         {
            PEL_Error::object()->raise_plain(
               "*** PEL_CoupledApplications error:\n"
               "    entries \"mpi_machines\" and \"mpi_machinefile\" are incompatible" ) ;
         }
      }

      a_code->destroy() ;
   }
   it->destroy() ;

   PEL_ASSERT( APPLIS->count() > 0 ) ;
   
   SILENT = exp->has_entry( "verbose" ) && !exp->bool_data( "verbose" ) ;
   TRACE = !SILENT ;
}

//----------------------------------------------------------------------
PEL_CoupledApplications:: ~PEL_CoupledApplications( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: ~PEL_CoupledApplications" ) ;
   
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   for( size_t i=0 ; i<NB_PROCESS ; i++ )
   {
      delete PROCESSES[i] ;
   }
   delete [] PROCESSES ;
   if( !TRACE ) delete_temporary_files() ;
   if( PEL_System::can_read( coupling ) && !PEL_System::erase( coupling ) )
   {
      PEL::out() << "Unable to clear " << coupling << std::endl ;
   }
}

//----------------------------------------------------------------------
bool
PEL_CoupledApplications:: has_coupled_applications( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications::has_coupled_applications" ) ;

   static bool first = true ;
   static bool COUPLED_APPLICATIONS = false ;
   if( first )
   {
      first = false ;
      wait_for_lock() ;
      std::ifstream in( coupling.c_str() ) ;
      COUPLED_APPLICATIONS = in ;
      release_lock() ;
   }
   return( COUPLED_APPLICATIONS ) ;
}

//----------------------------------------------------------------------
bool
PEL_CoupledApplications:: has_application( std::string const& appli_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: has_application" ) ;
   PEL_CHECK_PRE( ! appli_name.empty() ) ;

   bool result = has_coupled_applications() ;
   if( result )
   {
      init_on_demand() ;
      result = CODENAMES.has( appli_name ) ;
   }

   PEL_CHECK_POST( IMPLIES( !has_coupled_applications(), result == false ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer *
PEL_CoupledApplications:: receive( PEL_Object* a_owner,
                                   std::string const& me,
                                   std::string const& src )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: receive" ) ;
   PEL_CHECK_PRE( has_coupled_applications() ) ;
   PEL_CHECK_PRE( !me.empty() && has_application( me ) ) ;
   PEL_CHECK_PRE( !src.empty() && has_application( src ) ) ;
   PEL_CHECK_PRE( PEL_Exec::communicator()->rank() == 0 ) ;

   init_on_demand() ;

   log( me+" waiting data from "+src ) ;
   
   std::string exchange = name_of_exchange( src, me ) ;
   
   PEL_ModuleExplorer* result = 0 ;
   while( result == 0 )
   {
      wait_for_lock() ;
      
      if( exist( exchange ) )
      {
         PEL_Module* read = PEL_Module::create( 0, "Root", exchange ) ;
         PEL_ModuleIterator* it = read->create_module_iterator( 0 ) ;
         it->start() ;
   
         PEL_ASSERT( it->is_valid() ) ;
         PEL_Module* child = it->item() ;
         log( "Receive : ", child ) ;
         result = PEL_ModuleExplorer::create( a_owner, child ) ;
         read->set_owner( result ) ;
         release( exchange ) ;
         
         it->destroy() ; it=0 ;
      }

      release_lock() ;

      if( result == 0 ) 
      {
         PEL_System::sleep( WAIT_MILLISECOND ) ;
         // log( "waiting" ) ;
      }
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: send( PEL_Module const* module,
                                std::string const& me,
                                std::string const& dest )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications::send" ) ;
   PEL_CHECK_PRE( has_coupled_applications() ) ;
   PEL_CHECK_PRE( !me.empty() && has_application( me ) ) ;
   PEL_CHECK_PRE( !dest.empty() && has_application( dest ) ) ;
   PEL_CHECK_PRE( PEL_Exec::communicator()->rank() == 0 ) ;

   init_on_demand() ;
   
   log( me + " sending data to " + dest, module ) ;
   
   std::string exchange = name_of_exchange( me, dest ) ;
   while( exist( exchange ) )
   {
      PEL_System::sleep(WAIT_MILLISECOND) ;
   }

   wait_for_lock() ;

   std::ofstream os( exchange.c_str() ) ;

   os.precision( 14 ) ;
   module->print( os, 0 ) ;
   os.flush() ;

   os.close() ;
   PEL_ASSERT( exist( exchange ) ) ;

   release_lock() ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: run" ) ;

   delete_temporary_files() ;
   
   start_coupled_execution() ;

   PEL_System::Process::wait_for_child_processes( NB_PROCESS, PROCESSES ) ;

}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: start_coupled_execution( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: start_coupled_execution" ) ;

   // The file "coupling", to be shared with the processes associated
   // to all the codes, will be created => a lock is requested
   PEL_ASSERT( ask_for_lock() ) ;

   size_t ppid = PEL_System::process_id() ;
   std::ofstream out( coupling.c_str() ) ;
   out << ppid << " " << (!SILENT) << std::endl ;
   out.flush() ;
   
   PROCESSES = new PEL_System::Process*[APPLIS->count()] ;

   for( size_t i=0 ; i<APPLIS->count() ; i++ ) 
   {
      PEL_OneCoupledAppli* appli =
                     static_cast<PEL_OneCoupledAppli*>( APPLIS->at(i) ) ;
      
      stringVector args(0) ;
      if( appli->MPI_OPTIONS.size()>0 )
      {
         args.append( MPI_RUN ) ;
         if( appli->MPI_MACHINES.size()>0 )
         {
            appli->MPI_TMP_FILE =
               PEL_System::absolute_path( "machines_"+appli->NAME+".tmp" ) ;
            std::ofstream ff( appli->MPI_TMP_FILE.c_str() ) ;
            if( !ff ) PEL_Error::object()->raise_file_handling(
                                             appli->MPI_TMP_FILE, "open" ) ;
            for( size_t j=0 ; j<appli->MPI_MACHINES.size() ; ++j )
               ff << appli->MPI_MACHINES(j) << std::endl ;
            args.append( "-machinefile" ) ;
            args.append( appli->MPI_TMP_FILE ) ;
         }
         else if( !appli->MPI_MACHINEFILE.empty() )
         {
            args.append( "-machinefile" ) ;
            args.append( appli->MPI_MACHINEFILE ) ;
         }
         for( size_t j=0 ; j<appli->MPI_OPTIONS.size() ; ++j )
            args.append( appli->MPI_OPTIONS(j) ) ;
      }
      args.append( appli->EXE ) ;
      args.append( appli->DATAFILE ) ;
      args.append( "-o" ) ;
      args.append( "resu_"+appli->NAME ) ;
      for( size_t j=0 ; j<appli->EXTRAS.size() ; ++j )
         args.append( appli->EXTRAS(j) ) ;
      
      PEL::out() << "Starting : " << args <<std::endl ;

      PROCESSES[i] = new PEL_System::Process( args ) ;
      std::ostringstream os2 ;
      os2 << "process " << i << " named " << appli->NAME ; 
      PROCESSES[i]->set_name( os2.str() ) ;
      out << appli->NAME << " " << PROCESSES[i]->id() << std::endl ;
      PEL::out() << "Process : " << PROCESSES[i]->id() << std::endl ;
   }
   out.close() ;

   release_lock() ;

   NB_PROCESS = APPLIS->count() ;

   PEL_CHECK_POST( has_coupled_applications() ) ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: delete_temporary_files( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: delete_temporary_files" ) ;

   if( APPLIS != 0 )
   {
      for( size_t i=0 ; i<APPLIS->count() ; i++ ) 
      {
         PEL_OneCoupledAppli* appli_i =
            static_cast<PEL_OneCoupledAppli*>( APPLIS->at(i) ) ;
         for( size_t j=0 ; j<APPLIS->count() ; j++ ) 
         {
            PEL_OneCoupledAppli const* appli_j =
               static_cast<PEL_OneCoupledAppli const*>( APPLIS->at(j) ) ;
            if( j!=i )
            {
               std::string exchange =
                           name_of_exchange( appli_i->NAME, appli_j->NAME ) ;
               if( PEL_System::can_read( exchange ) && 
                   !PEL_System::erase( exchange ) )
               {
                  PEL::out() << "Unable to clear " << exchange
                             << std::endl ;
               }
            }
         }
         if( ! appli_i->MPI_TMP_FILE.empty() )
         {
            if( PEL_System::can_read( appli_i->MPI_TMP_FILE ) &&
                !PEL_System::erase( appli_i->MPI_TMP_FILE ) )
            {
               PEL::out() << "Unable to clear " << appli_i->MPI_TMP_FILE
                          << std::endl ;
            }
            appli_i->MPI_TMP_FILE = "" ;
         }
      }
   }
   if( PEL_System::can_read( lock ) && !PEL_System::erase( lock ) )
   {
      PEL::out() << "Unable to clear " << lock << std::endl ;
   }
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: init_on_demand( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications::init_on_demand" ) ;

   if( CODENAMES.size() == 0 ) // executé la première fois pour chaque process
   {
      wait_for_lock() ;
      std::ifstream in( coupling.c_str() ) ;
      if( !in ) 
      {
         PEL_Error::object()->raise_plain(
            "PEL_CoupledApplications : Unable to open main coupling repertory "
            + coupling ) ;
      }
   
      int ppid ;
   
      in >> ppid >> TRACE ;
      while( in.good() ) 
      {
         std::string name ;
         int pid ;

         in >> name >> pid ;
		
         CODENAMES.append( name ) ;
      }
      release_lock() ;
   }
}

//----------------------------------------------------------------------
bool
PEL_CoupledApplications::exist( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications::exist" ) ;
   return PEL_System::can_read( file ) ;
}

//----------------------------------------------------------------------
bool
PEL_CoupledApplications:: ask_for_lock( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: ask_for_lock" ) ;
   PEL_ASSERT( fd<0 ) ;
   bool result = false ;
   
   // O_RDONLY : Open for reading only.
   // If O_CREAT and O_EXCL are set, open() will fail if the file exists.   
   if( (fd = PEL_System::open_file_descriptor( lock, O_RDONLY | O_CREAT | O_EXCL )) >= 0 )
   {
      result = true ;
      // log("-L- get lock");
   }
   else 
   {
      // log("-L- can't get lock");
   }
   
   return result ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: wait_for_lock( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: wait_for_lock" ) ;

   while( !ask_for_lock() ) PEL_System::sleep( WAIT_MILLISECOND ) ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: release_lock( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications::release_lock" ) ;

   PEL_ASSERT( fd>=0 ) ; 
   PEL_ASSERT( PEL_System::close_file_descriptor( fd ) == 0 )  ;
   fd = -1 ;
   PEL_System::erase( lock ) ; 

}

//----------------------------------------------------------------------
std::string
PEL_CoupledApplications:: name_of_exchange( std::string const& from,
                                            std::string const& to )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: name_of_exchange" ) ;
   return( ".#CA" + from + "_" + to ) ;
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: log( std::string const& msg,
                               PEL_Module const* extra )
//----------------------------------------------------------------------
{
   if( TRACE ) 
   {
      PEL::out() << std::endl << "*******************************" 
                 << std::endl ;
      PEL::out() << msg << std::endl ;
      if( extra!=0 ) 
      {
         extra->print( PEL::out() , 1 ) ;
      }
      PEL::out() << "*******************************"  
                 << std::endl << std::endl ;
      PEL::out().flush() ;
   }
}

//----------------------------------------------------------------------
void
PEL_CoupledApplications:: release( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CoupledApplications:: release" ) ;

   if( exist(file) ) 
   {
      PEL_System::erase( file ) ;
   }
   PEL_ASSERT( !exist( file ) ) ;
}

//internal-------------------------------------------------------------------
PEL_OneCoupledAppli:: PEL_OneCoupledAppli( PEL_Object* a_owner )
//internal-------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NAME()
   , EXE( PEL_Exec::name_of_exe() )
   , DATAFILE()
   , EXTRAS( stringVector( "-v" ) )
   , MPI_OPTIONS(0)
   , MPI_MACHINES(0)
   , MPI_MACHINEFILE()
   , MPI_TMP_FILE()
{
}

//internal-------------------------------------------------------------------
PEL_OneCoupledAppli:: ~PEL_OneCoupledAppli( void )
//internal-------------------------------------------------------------------
{
}
