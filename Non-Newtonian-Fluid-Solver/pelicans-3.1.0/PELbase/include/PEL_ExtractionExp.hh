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

#ifndef PEL_EXTRACTION_EXP_HH
#define PEL_EXTRACTION_EXP_HH

#include <PEL_TransferExp.hh>
#include <string>

class PEL_ModuleExplorer ;

/*
Expressions extracting data from an attached data structure that has
to be primarily specified by calling `::initialize'( `mod' ).
mod must not be modified until having called `::reset' method.

Before calling `::initialize', none of the expressions implemented here
can be used.

---
name     : has_data
argument : String
type     : Bool

The returned value is true if there exists, in the attached data structure,
an entry whose keyword match the argument.

Example:
  if( `has_data'( "/TOTO/my_entry" ) )
  MODULE titi
     ...
  END MODULE titi

---
name      : extracted_data
arguments : String, optional second argument
type      : everything

The returned value is the result of the evaluation of the data
whose keyword matches the first argument (in the attached data structure),
if such an entry exists (in the attached data structure). If not, the
second argument can be used to define a default value.

Example 1:

   toto = `extracted_data'( "/TOTO/my_entry" )
   
Example 2:

   toto = `extracted_data'( "/TOTO/my_entry", 3. )
   (toto is "/TOTO/my_entry" in `mod' if any, 3. elsewhere)
      
   titi = `extracted_data'( "/TOTO/my_entry", < "a" "b" > )
   (titi is "/TOTO/my_entry" in `mod' if any, < "a" "b" > elsewhere)

Example 3:

   For implementation reasons, the second argument becomes mandatory in
   "if" constructions, even when it is not used.

   if( has_data( "/TOTO/my_entry" ) )
   MODULE titi
      toto = `extracted_data'( "/TOTO/my_entry", 0. )
   END MODULE titi
   
---
name     : has_module
argument : String
type     : Bool

The returned value is true if there exists, in the attached data structure,
a module whose keyword match the argument.

Example:
  if( `has_module'( "/TOTO" ) )
  MODULE titi
     toto = `extracted_data'( "/TOTO/my_entry", 3. )
     ...
  END MODULE titi
  
---
name     : extracted_module
argument : String, String, optional third argument
type     : String

The returned value is the name of a temporary file containing the
MODULE called according to the first argument in the attached data structure.
The extracted module is renamed as the second argument of the function.

Example1:

    MODULE titi
       #include( `extracted_module'( "/TOTO", "TITI" ) )
    END MODULE titi

    The module of name "TOTO" is extracted from the database,
    and included as module "TITI".

Example2:

   For implementation reasons, a special syntax (with a dummy third argument)
   is required in "if" constructions.

   if( has_module( "/TOTO/module" ) )
   MODULE titi
      #include( `extracted_module'( "/TOTO/module", "module", "" ) )
   END MODULE titi
   
PUBLISHED
*/

class PEL_EXPORT PEL_ExtractionExp : public PEL_TransferExp
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static void initialize( PEL_Module const* mod ) ;
      static void reset( void ) ;

      static bool is_initialized( void ) ;
      static PEL_Module const* data_base( void ) ;
      
   //-- Context
      
      virtual void declare( PEL_List* lst ) const ;
      virtual bool context_has_required_variables( 
                                               PEL_Context const* ct ) const ;
      
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------
      
      PEL_ExtractionExp( void ) ;
     ~PEL_ExtractionExp( void ) ;
      PEL_ExtractionExp( PEL_ExtractionExp const& other ) ;
      PEL_ExtractionExp& operator=( PEL_ExtractionExp const& other ) ;
      
      enum ExtractionExp{ has_data, ext_data, has_mod, ext_mod } ;
      
      PEL_ExtractionExp( PEL_Object* a_owner,
                         ExtractionExp op,
                         std::string const& a_name,
                         PEL_Sequence const* argument_list ) ;

   //-- Plug in

      PEL_ExtractionExp( std::string const& a_name, ExtractionExp op ) ;

      virtual PEL_ExtractionExp* create_replica( 
                      PEL_Object* a_owner,
                      PEL_Sequence const* argument_list ) const ;

      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      virtual PEL_Data const* data( PEL_Context const* ct ) const ;
      
   //-- Private methods

      static std::string const& temporary_file( void ) ;

      static std::string const& data_name( std::string const& exp_name,
                                           PEL_Data const* d ) ;
      
      void extract_module( std::string const& file_name,
                           std::string const& d_name,
                           std::string const& m_name ) const ;

   //-- Class attributes
      
      static PEL_Module const* DB_MOD ;

      static PEL_ExtractionExp* PROTO_HAS_DATA ;
      static PEL_ExtractionExp* PROTO_DATA ;
      static PEL_ExtractionExp* PROTO_HAS_MOD ;
      static PEL_ExtractionExp* PROTO_EXTRACTED_MODULE ;
 
   //-- Attributes

      ExtractionExp const OP ;
      std::string TEMP_FILE_NAME ;
      PEL_Data const* SRC ;
} ;

#endif
