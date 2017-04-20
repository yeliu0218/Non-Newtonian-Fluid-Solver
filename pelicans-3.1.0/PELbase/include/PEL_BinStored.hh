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

#ifndef PEL_BIN_STORED_HH
#define PEL_BIN_STORED_HH

#include <PEL_TransferExp.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <boolVector.hh>

class PEL_Vector ;

// Tools to perform binary storage.

class PEL_EXPORT PEL_BinStored : public PEL_TransferExp
{
   public: //----------------------------------------------------------

   //-- Access
      
      virtual PEL_Data::Type data_type( void ) const ;

  //-- Binary file related access

      // Is type supported for binary storage through PEL_BinStored ?
      static bool is_type_supported( PEL_Data::Type a_type ) ;

      // Initialize file for writting
      static void init_binary_file( std::string const& file_name ) ;

      // Add `data' to binary file `file_name' and return
      // refering expression.
      static PEL_BinStored const* create_reference(
                    PEL_Object* a_owner,
                    PEL_Data const* a_data,
                    std::string const& file_name,
                    bool local_reference_file_name = false ) ;
      
      // Is `file_name' corresponding to a valid binary file ?
      static bool is_valid_binary_file( std::string const& file_name ) ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      struct BinaryRecord 
      {
            PEL_Data::Type type ;
            size_t length ;
            size_t number ;
            size_t foo[8] ;
      } ;
      
     ~PEL_BinStored( void ) ;
      PEL_BinStored( PEL_BinStored const& other ) ;
      PEL_BinStored& operator=( PEL_BinStored const& other ) ;
      
      PEL_BinStored( PEL_Object* a_owner,
                     PEL_Sequence const* argument_list ) ;
      
      // Retrieve data from binary file `file_name' whose identifier is
      // `record_number'.
      static PEL_Data const* restore_from_binary(
         PEL_Object * a_owner,
         std::string const& file_name,
         size_t record_number ) ;

      // Bad record number
      static const size_t bad_record ;

      // Last record number of file named `file_name' or bad_record if
      // file is an empty valid binary file.
      static size_t last_record_number( std::string const& file_name ) ;
      
   //-- Plug in

      PEL_BinStored( void ) ;
      
      virtual PEL_BinStored* create_replica(
         PEL_Object * a_owner,
         PEL_Sequence const* argument_list ) const ;
      
   //-- Identification

      virtual std::string const& usage( void ) const ;
      
   //-- Arguments

      virtual bool valid_arguments(
                           PEL_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      PEL_Data const* data( PEL_Context const* ct ) const ;
      
   //-- Class attributes
      
      static PEL_BinStored const* static_OpComponent ;
      static size_t preamble_length ;
      static size_t last_record ;
      static std::string last_file ;
      
   //-- Attributes
      
      mutable PEL_Data const* my_data ;
};

#endif
