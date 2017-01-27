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

#ifndef PEL_FILE_TO_MODULE_HH
#define PEL_FILE_TO_MODULE_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT PEL_FileToModule : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_FileToModule const* object( std::string const& format ) ;
      
      static bool has( std::string const& format ) ;
      
      static void find_file_format( std::string const& a_file_name,
                                    std::string& a_format ) ;
      
      static std::string const& list_of_formats( void ) ;

   //-- Characteristics

      std::string const& format( void ) const ;
      
      std::string const& default_motif( void ) const ;

   //-- Module building
      
      virtual PEL_Module* create_from_file( 
                                 PEL_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const = 0 ;
      
   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PEL_FileToModule( void ) ;

      PEL_FileToModule( std::string const& a_format,
                        std::string const& a_default_motif ) ;

   //-- Preconditions, Postconditions, Invariant

      bool create_from_file_PRE( PEL_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const ;
      
      bool create_from_file_POST( PEL_Module const* result,
                                  PEL_Object* a_owner,
                                  std::string const& module_name,
                                  std::string const& file_name ) const ;
      
   private: //----------------------------------------------------------

      PEL_FileToModule( void ) ;
      PEL_FileToModule( PEL_FileToModule const& other ) ;
      PEL_FileToModule& operator=( PEL_FileToModule const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
      
      static std::string& formats( void ) ;
      
   //-- Attributes

      std::string MY_FORMAT ;
      std::string MY_MOTIF ;
} ;

#endif
