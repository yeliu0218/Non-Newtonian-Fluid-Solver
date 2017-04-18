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

#ifndef PEL_SYSTEM_HH
#define PEL_SYSTEM_HH

#include <PEL_Object.hh>

#include <iosfwd>
class stringVector ;

/*
interfaces to operating system features that are available on
all operating systems
*/

class PEL_EXPORT PEL_System : public PEL_Object
{
   public: //-----------------------------------------------------------------
      
      // host name
      static std::string const& host_name( void ) ;
      
      // process id
      static int process_id( void ) ;

      // character used to form path name description on native system
      static char path_name_separator( void ) ;

      // Does `str' matches `pattern' description ?
      static bool matches(  std::string const& str,
                            std::string const & pattern ) ;
      
      // search for `filename' entry in `directory'
      static void find( std::string const& directory,
                        std::string const & filename,
                        stringVector & result,
                        bool recurse=true,
                        std::string const & prefix="" ) ;

      // execute command `cmd' in `directory' and save output in `output_file'
      static int run( stringVector const& cmd, 
                      std::string const& directory,
                      std::string const& output_file ) ;
      
      // build directory
      static bool mkdir( std::string const& directory ) ;

      static bool changedir( std::string const& directory ) ;
      
      static bool copy( std::string const& src, std::string const& dest ) ;

      // remove the file called `filename'
      static bool erase( std::string const& filename ) ;

      // working directory
      static std::string const& working_directory( void ) ;

      // kernel name
      static std::string const& sysname( void ) ;

      // currently amount of used TIME
      static double user_time( void ) ;

      // absolute TIME elapsed from 1/1/1970 in s
      static double epoch_time( void ) ;

      // delay execution of the calling process for (at least) 
      // `msec' milliseconds
      static void sleep( size_t msec ) ;

      // currently amount of used memory
      static size_t used_memory( void ) ;
      
      // Name of compiler used to build PELICANS library.
      static std::string const& compiler_name( void ) ;

      // Stripped non-directory suffix from file name `path_and_name'.
      // Directory separator is given by `separator'.
      static std::string dirname( std::string const& path_and_name,
                                  char separator = path_name_separator() ) ;
      
      // Return `path_and_name' with any leading directory components removed.
      // Directory separator is given by `separator'.
      static std::string basename( std::string const& path_and_name,
                                   char separator = path_name_separator() ) ;

      // Return absolute path corresponding to relative or absolute `filename'.
      static std::string absolute_path( std::string const& filename ) ;
      
      // Active exception treatment.
      static void exception_trapping( void ) ;

      // Is big endian encoding format supported on the current machine ?
      static bool big_endian_encoding( void ) ;

      // new recommended block size for a block of memory that needs
      // to be made larger (`current_size' is the actual size of the block
      // and `expected_size' is the new needed size)
      static size_t new_block_size( size_t current_size, 
                                    size_t expected_size ) ;

      // file access on reading
      static bool can_read( std::string const& filename ) ;
      
      // file access on writing
      static bool can_write( std::string const& filename ) ;

      // close file descriptor
      static int close_file_descriptor(int fildes) ;

      // open and possibly create a file or device
      static int open_file_descriptor(std::string filename, int flags ) ;

      // process management
      class Process 
      {
         public :
            
            Process( stringVector const& args ) ;
            
            size_t id( void ) const ;
            
            static void wait_for_child_processes( size_t nb, Process ** child ) ;
            
            ~Process(void) ;
            
            void set_name( std::string a_name ) ;
            
            std::string const& name( void ) const ;
            
         private : 
            
            bool running ;
            std::string my_name ;
            size_t pid ;
            void * sdef ;
      } ;
      
      // exit current process
      static void exit( int exit_status ) ;

      // environment variable content (empty string if not defined)
      static std::string getenv( std::string const& str ) ;

   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------

      PEL_System( void ) ;
     ~PEL_System( void ) ;
      PEL_System( PEL_System const& other ) ;
      PEL_System& operator=( PEL_System const& other ) ;

      static std::string const& string( std::string const& cmd ) ;
      
      static PEL_System const* SINGLETON ;
      
} ;

#endif
