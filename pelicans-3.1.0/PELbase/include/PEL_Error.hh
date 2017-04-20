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

#ifndef PEL_ERROR_HH
#define PEL_ERROR_HH

#include <PEL_Object.hh>
#include <PEL_Data.hh>

#include <iosfwd>
#include <string>

class PEL_Map ;
class PEL_ModuleExplorer ;
class stringVector ;

//---------------------------------------------------------------------------
//   Facilities for adapting the error handling mechanism
//---------------------------------------------------------------------------

class PEL_EXPORT PEL_Error
{
   public: //----------------------------------------------------------------

      static PEL_Error* object( void ) ;

      void raise_plain( std::string message ) ;

      void raise_internal( std::string message ) ;

      void raise_read_syntax_error( std::string const& file,
                                    int line,
                                    std::string const& last,
                                    std::string const& nature ) ;

      void display_info( std::string message ) ;
      
   //-- Incomplete implementation

      void raise_not_tested( std::string file, int line, std::string test ) ;

      void raise_not_implemented( PEL_Object const* oo,
                                  std::string method ) ;

   //-- Simulation of virtual methods with respect to many objects

      void raise_bad_types( PEL_Object const* oo,
                            std::string method,
                            PEL_Object const* argument ) ;

      void raise_bad_types( PEL_Object const* oo,
                            std::string method,
                            PEL_Object const* argument_1, 
                            PEL_Object const* argument_2 ) ;

   //-- Registration

      void raise_missing_keyword( PEL_ModuleExplorer const* exp,
                                  std::string keyword ) ;

      void raise_missing_module( PEL_ModuleExplorer const* exp,
                                 std::string path_and_name ) ;
      
      void raise_module_error( PEL_ModuleExplorer const* exp,
                               std::string error ) ;

      void raise_bad_data_type( PEL_ModuleExplorer const* exp,
                                std::string keyword,
				PEL_Data::Type query_kind ) ;

      void raise_bad_file( PEL_ModuleExplorer const* exp,
                           std::string const& filename,
                           std::string const& access ) ;

      void raise_not_evaluable( PEL_ModuleExplorer const* exp,
                                std::string keyword,
                                stringVector const& undefined_variables ) ;

      void raise_bad_data_value( PEL_ModuleExplorer const* exp,
                                 std::string keyword,
				 std::string allowed_values ) ;
      
      void raise_data_error( PEL_ModuleExplorer const* exp,
                             std::string keyword,
                             std::string error ) ;

   //-- Contract

      void raise_precondition_violation( std::string test ) ;

      void raise_postcondition_violation( std::string file,
                                          int line,
                                          std::string test ) ;
      void raise_invariant_violation( std::string file,
                                      int line,
                                      std::string test ) ;
      void raise_assertion_violation( std::string file,
                                      int line,
                                      std::string test ) ;

      void raise_bad_object( std::string message,
                             PEL_Object const* a_object ) ;
      
   //-- File manipulation

      void raise_file_handling( std::string file,
                                std::string operation ) ;

   //-- Warnings

      void trace( std::string const& event ) ;
      
      void display_new_syntax( PEL_ModuleExplorer const* older_module,
                               PEL_ModuleExplorer const* result_module ) const ;
      
      static void notify( PEL_Module* list,
                          std::string const& message,
                          std::string const& where,
                          stringVector const& choices )  ;
      
      static PEL_Module* notify( PEL_Module* list,
                                  std::string const& message,
                                  std::string const& where )  ;
      
      void display_data_checking( PEL_ModuleExplorer* report ) const ;
      
   //-- Abnormal termination
      
      // Exits on error and return exit_status to shell.
      static void exit( size_t status = 1 ) ;
      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

     ~PEL_Error( void ) ;
      PEL_Error( PEL_Error const& other ) ;
      PEL_Error const& operator=( PEL_Error const& other ) ;

      PEL_Error( void ) ;

      void print_invocated_methods( void ) ;
      void print_client_and_server( void ) ;
      
      void begin( void ) ;
      void end( void ) ;
      std::string const& hline( void ) const ;
      
      static void exit_with_internal_error( void ) ;
      
      //---------------------------------------------------------------------
      //   ATTRIBUTES
      //---------------------------------------------------------------------
      std::ostream& os ;

} ;

#endif
