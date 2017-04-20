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

#ifndef PEL_NATIVE_FILE_TO_MODULE_HH
#define PDE_NATIVE_FILE_TO_MODULE_HH

#include <PEL_FileToModule.hh>

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT PEL_NativeFileToModule : public PEL_FileToModule
{
   public: //-----------------------------------------------------------

   //-- Module building

      virtual PEL_Module* create_from_file( 
                                 PEL_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PEL_NativeFileToModule( void ) ;
     ~PEL_NativeFileToModule( void ) ;
      PEL_NativeFileToModule( PEL_NativeFileToModule const& other ) ;
      PEL_NativeFileToModule& operator=( 
                              PEL_NativeFileToModule const& other ) ;
      
      PEL_NativeFileToModule( std::string const& a_format,
                              std::string const& a_default_motif ) ;

   //-- Internals
      
      PEL_Module* create_from_multicolumns( 
                              PEL_Object* a_owner,
                              std::string const& name,
                              std::string const& file_name,
                              std::string const& separator ) const ;
      
   //-- Class attributes

      static PEL_NativeFileToModule* PEL_OBJ ;
      static PEL_NativeFileToModule* GENE_OBJ ;
      static PEL_NativeFileToModule* CSV_OBJ ;
} ;

#endif


