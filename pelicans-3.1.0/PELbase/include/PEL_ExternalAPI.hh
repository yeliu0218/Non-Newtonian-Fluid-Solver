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

#ifndef PEL_EXTERNAL_API_HH
#define PEL_EXTERNAL_API_HH

#include <PEL_Object.hh>

class PEL_ObjectRegister ;
class stringVector ;

/*
External applications, performing their specific initialization
and termination.

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass, say MyExtAPI.
   2. Choose a name for MyAppli, say "my_api", 
      and a priority level, say pl.
   3. Declare all constructors private.
   4. Define an instance to be registered :
      4.1 Implement a default constructor that initializes the
          `PEL_ExternalAPI::' subobject by calling
               `PEL_ExternalAPI( std::string const&, size_t )'
          with "my_api" and pl as argument.
          Example of pseudo-code :
          | MyExtAPI:: MyExtAPI( void ) : PEL_ExternalAPI( "my_api", pl )" ) {}
      4.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyExtAPI.hh) :
             | static MyExtAPI const* SINGLETON ;
             definition (in the implementation file, eg MyExtAPI.cc) :
             | MyExtAPI const* MyExtAPI::SINGLETON = new MyExtAPI() ;'
   5. Implement `::initialize' that does the initialization required
      by the external API at hand.
   6. Implement a private destructor that does the deinitialization required
      by the external API at hand. 

PUBLISHED*/

class PEL_EXPORT PEL_ExternalAPI : public PEL_Object
{

   public: //-----------------------------------------------------------
      
   //-- Management of all registered instances

      // Call `::initialize' for all registered instances in the
      // reverse order of their priority level (instances with
      // higher priority are initialized first).
      static void initialize_all_APIs( int& argc, char**& argv ) ;

      // Call `PEL_Object::destroy' for all registered instances in the
      // order of their priority level (instances with
      // lower priority are terminanted first).
      static void terminate_all_APIs( void ) ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PEL_ExternalAPI( void ) ;

      // registration of `self', calling it `a_name' and
      // setting its priority level to `a_priority_level'.
      PEL_ExternalAPI( std::string const& a_name, size_t a_priority_level ) ;

   //-- Current instance management

      // Initialize `self'.
      virtual void initialize( int& argc, char**& argv ) = 0 ;

   private: //----------------------------------------------------------

      PEL_ExternalAPI( void ) ;
      PEL_ExternalAPI( PEL_ExternalAPI const& other ) ;
      PEL_ExternalAPI& operator=( PEL_ExternalAPI const& other ) ;
      
      static PEL_ObjectRegister* plugins_map( void ) ;
      static stringVector& plugins_names( void ) ;

   //-- Attributes

      std::string MY_NAME ;
      size_t MY_PRIORITY ;
} ;

#endif



