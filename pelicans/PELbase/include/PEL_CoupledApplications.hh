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

#ifndef PEL_COUPLED_APPLICATIONS_HH
#define PEL_COUPLED_APPLICATIONS_HH

#include <PEL_Application.hh>
#include <PEL_System.hh>
#include <intVector.hh>
#include <string>
#include <stringVector.hh>

class PEL_Vector ;
class PEL_ModuleExplorer ;

/*
Coupled concurrent applications with smart exchange of data.

Each instance is an application that is split into elementary components
called codes, such that:
   - codes are executed as separate processes;
   - codes communicate by passing messages (databases of the PELICANS
     Hierarchical Data System).

Unix and Windows operating systems are handled.

INSTANCE DELIVERY AND INITIALIZATION
------------------------------------
  
Instances are delivered by `create' whose second argument is associated 
to a Module Hierarchy whose structure, in the particuliar example
of a coupled application running two codes, is sketched below :

   MODULE PEL_Application
      concrete_name = "PEL_CoupledApplications"
      verbose = true // Optional entry
      MODULE list_of_coupled_codes
         MODULE Code1
            //  executable = join( "..", "bin", "exe" ) Optional entry
            datafile = join( this_file_dir(), "data1.pel" )
            extra_command = < "-Call" >
            mpi_machinefile = vector( host_name() )        // mpi machines
            mpi_options = < "-np", "2" >                   // mpi options
            name = "code1"
         END MODULE Code1
         MODULE Code2
            //  executable = join( "..", "bin", "exe" ) Optional entry
            datafile = join( this_file_dir(), "data2.pel" )
            name = "code2"
            extra_command = < "-Call" >
         END MODULE Code2
      END MODULE list_of_coupled_codes
    END MODULE PEL_Application   

The entry of keyword "verbose" may be used to trace message passing between 
the processes.
If the entry of keyword "executable" is missing for a given code, the 
associated process will be launched using the same executable as that
used for the current coupled application.

MESSAGE PASSING
---------------

Message passing between the processes is performed by sending and receiving 
databases of the PELICANS Hierarchical Data System via respectively
the `send' and `receive' methods.
The reponsibility of the synchronization between the various processes is 
delegated the processes themselves.

TERMINATION
-----------

When execution terminates, the exit status of each process is verified.
If a process returns a failure code, all the remaining processes are killed
and PEL_CoupledApplications returns a failure exit status to system.
*/

class PEL_EXPORT PEL_CoupledApplications : public PEL_Application
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_CoupledApplications* create( PEL_Object* a_owner,
                                              PEL_ModuleExplorer const* exp ) ;

      static bool has_coupled_applications( void ) ;

   //-- Utilities : message passing between applications (0.1)

      static bool has_application( std::string const& appli_name ) ;
      
      // navigator to interrogate the database expected from the code
      // called `src'
      static PEL_ModuleExplorer* receive( PEL_Object* a_owner,
                                          std::string const& me,
                                          std::string const& src ) ;
      
      // Send the message `module' to the code called `dest'.
      static void send( PEL_Module const* module,
                        std::string const& me,
                        std::string const& dest ) ;
      
   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_CoupledApplications( void ) ; 
      PEL_CoupledApplications( PEL_CoupledApplications const& other ) ;
      PEL_CoupledApplications& operator=( 
                               PEL_CoupledApplications const& other ) ;

      PEL_CoupledApplications( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) ;
      PEL_CoupledApplications( PEL_Object* a_owner,
                               stringVector& args ) ;
      
   //-- Plug in
      
      PEL_CoupledApplications( void ) ;
      
      virtual PEL_CoupledApplications* create_replica( 
                                     PEL_Object* a_owner,
				     PEL_ModuleExplorer const* exp ) const ;

   //-- Internals of the PEL_CoupledApplication instance

      /*
      the access to shared ressources between processes (here files) is
      restricted with a lock mechanism (a lock file is used somehow as a 
      semaphore) :
         - When a process wants to access to a shared ressource, he must ask
           for a lock.
         - As soon as a lock is delivered to that process, no other process
           can access to the ressource.
         - After finishing with that ressource, the lock is removed.
      */

      // Create a process for each code, and build a file that contains
      // the correspondance between the code names and the processus ids.
      void start_coupled_execution( void ) ;

      void delete_temporary_files( void ) ;

   //-- Internals for the current code process

      // Initialize the static members of the current code process
      // by reading the file created by `start_coupled_execution'.
      static void init_on_demand( void ) ;

      static bool exist( std::string const& file ) ;

      // Is is possible to create a lock file ? If yes, do create it.
      static bool ask_for_lock( void ) ;
      
      // Wait until it is possible to create a lock file (and create it).
      static void wait_for_lock( void ) ;

      // Remove the lock file
      static void release_lock( void ) ;

      // name of the file that will be used for the message passing
      // from process `from' to process 'who'
      static std::string name_of_exchange( std::string const& from,
                                           std::string const& to ) ;

      static void log( std::string const& msg,
                       PEL_Module const* extra = 0 ) ;

      static void release( std::string const& file ) ;
    
   //-- Class attributes
      
      static PEL_CoupledApplications const* PROTOTYPE ;
      static stringVector CODENAMES ;

   //-- Attribute
      
      bool SILENT ;
      PEL_Vector* APPLIS ;
      PEL_System::Process** PROCESSES ;
      size_t NB_PROCESS ;
      std::string MPI_RUN ;
} ;

#endif



