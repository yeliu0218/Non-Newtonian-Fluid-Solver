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

#ifndef PEL_OBJECT_WRITER_HH
#define PEL_OBJECT_WRITER_HH

#include <PEL_Object.hh>

#include <stack>
#include <string>

class PEL_Data ;
class PEL_ModuleExplorer ;

/*
Servers used to store objects so that they can be retrieved with
associated `PEL_ObjectReader::' instances.

Objects are stored in files according to some options
specified in the Hierarchical Data Structure provided 
at creation of `self'. That data structure is briefly described below.

The entry of keyword "output_format" defines the format of the
saving files. There are two possibilities:
   - "text": human readable but not exact (truncated values)
   - "hybrid": parts of the data remain readable but 
               double or integer values are stored in binary format 
               (file with "bin" extension)
                
Several saving strategies are available :

   - all the saved cycles are stored in the same file

     example :
     
        MODULE PEL_ObjectWriter
           type = "all_cycles_in_one_file"
           file_name = join( getcwd(), "saving.pel" )
           output_format = "hybrid"
        END MODULE PEL_ObjectWriter

        A text file named "saving.pel" is created to store all the cycles
        (a companion binary file named "saving.pel.bin" is also created 
        to store the double and integer values).

   - each saved cycle is stored in a separate file (one cycle per file)

     example :
     
        MODULE PEL_ObjectWriter
           type = "cycles_in_separate_files"
           files_basename = join( getcwd(), "saving" )
           output_format = "hybrid"
        END MODULE PEL_ObjectWriter

        A sequence of text files named "saving.00001.pel", 
        "saving.00002.pel",... is created to store respectively 
        the first cycle, the second cycle,...
        (a sequence of companion binary files named "saving.00001.pel.bin", 
        "saving.00002.pel.bin",... is also created to stored
         the double and integer values).

   - only the last two cycles are stored

     example :
     
        MODULE PEL_ObjectWriter
           type = "last_two_cycles"
           file_name_0 = join( getcwd(), "saving_0.pel" )
           file_name_1 = join( getcwd(), "saving_1.pel" )
           output_format = "hybrid"
        END MODULE PEL_ObjectWriter

        The text files "saving_0.pel" and "saving_1.pel" are created to store
        the last two cycles ; the time of last modification of these files 
        identifies that of the more recent saving.
        (the companion binary files named "saving_0.pel.bin" and 
        "saving_1.pel.bin" are also created to store with the double and 
        integer values).


See `PEL_ObjectReader::' and `PEL_ApplicationRestorer::' for restoring
objects stored with `PEL_ObjectWriter::' objects.

PUBLISHED
*/

class PEL_EXPORT PEL_ObjectWriter : public PEL_Object
{
   public: //----------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PEL_ObjectWriter* create( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp,
                                       PEL_ModuleExplorer const* header_exp ) ;

   //-- Cycles(1.)

      // Start a new cycle.
      void start_cycle( void ) ;

      // Terminate the current cycle. `::finalize_object' must have been
      // called as many times as `::start_new_object'. If not, a fatal error 
      // is raised.
      void terminate_cycle( void ) ;

      // Is there a cycle that is started and not terminated ?
      bool has_an_opened_cycle( void ) const ;

      // cycle number
      size_t cycle_number( void ) const ;

   //-- Object storing(2.)

      // Notify that the storage of a new object is starting so that all 
      // subsequent calls to `::add_entry' are relative that object,
      // until `::finalize_object' or `::start_new_object' are called.
      void start_new_object( std::string const& class_name ) ;

      // nonzero number associated to the object being currently stored if any,
      // 0 otherwize
      size_t current_object_number( void ) const ;

      void add_entry( std::string const& keyword, PEL_Data* data ) ;

      // Notify that the storage of the current object is completed.
      void finalize_object( void ) ;

   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      // all_cycles :      save all the cycles in one file
      // per_one_cycle :   save all the cycles, but one per file
      // last_two_cycles : save only the two last cycles
      enum PEL_ObjectWriterType 
      {
         all_cycles,
         per_one_cycle,
         last_two_cycles
      } ;
         

      PEL_ObjectWriter( void ) ;
     ~PEL_ObjectWriter( void ) ;
      PEL_ObjectWriter( PEL_ObjectWriter const& other ) ;
      PEL_ObjectWriter& operator=( PEL_ObjectWriter const& other ) ;

      PEL_ObjectWriter( PEL_Object* a_owner,
                        PEL_ObjectWriterType const writer_type,
                        PEL_ModuleExplorer const* exp,
                        PEL_ModuleExplorer const* header_exp ) ;

   //-- Internals

      void initialize_saving_file( void ) ;
      void set_file_name( void ) ;
      void write_header( void ) const ;
      void write_communicator( void ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_ObjectWriterType const TYPE ;

      std::string const OFILE_FORMAT ;
      PEL_ModuleExplorer const* const HEADER_EXP ;

      // File names :
      std::string OFILE_NAME ;
      std::string OFILE_NAME1 ;
      std::string OFILE_NAME2 ;
      
      size_t iCYCLE ;
      size_t NB_OBJECTS ;
      std::stack< PEL_Module* > MODS ;
} ;

#endif
