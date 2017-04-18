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

#ifndef PEL_MODULE_EXPANDER_HH
#define PEL_MODULE_EXPANDER_HH

#include <PEL_Application.hh>
#include <string>

class PEL_Module ;
class PEL_ModuleExplorer ;

/*
Builders of data structures expanded from two parts : the first one is
a skeleton data structure with missing pieces, the second one is a
complementary data structure providing these missing pieces.

The missing pieces of the skeleton are entries whose data expressions
implemented by the `PEL_ExtractionExp::' class.

PUBLISHED
*/

class PEL_EXPORT PEL_ModuleExpander : public PEL_Application
{
   public: //-----------------------------------------------------------

      static PEL_ModuleExpander* create( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) ;
      
   //-- Program core execution

      virtual void run( void ) ;

   //-- Expander

      static PEL_Module const* create_expanded_module(
                          PEL_Object* a_owner,
                          PEL_Module const* input_mod,
                          std::string const& skeleton_file_name ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_ModuleExpander( void ) ; 
      PEL_ModuleExpander( PEL_ModuleExpander const& other ) ;
      PEL_ModuleExpander& operator=( PEL_ModuleExpander const& other ) ;

      PEL_ModuleExpander( PEL_Object* a_owner,
                          std::string const& skeleton_file,
                          std::string const& input_file,
                          std::string const& expanded_file,
                          std::string const& submod_name ) ;
      
   //-- Plug in
      
      PEL_ModuleExpander( void ) ;
      
      virtual PEL_ModuleExpander* create_replica( 
                            PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) const ;

      virtual PEL_ModuleExpander* create_replica_from_args(
                            PEL_Object* a_owner,
                            stringVector& args ) const ;
    //-- Command line

      virtual void print_usage( void ) const ;
      virtual void print_operands( void ) const ;

   //-- Expander

      static PEL_Module* create_skeleton_module(
                   PEL_Object* a_owner, std::string const& file_name ) ;

   //-- Class attribute
      
      static PEL_ModuleExpander const* PROTOTYPE ;
      
   //-- Attribute
      
      std::string const BASE ;
      std::string const OUTPUT ;
      std::string const INPUT ;
      std::string SUBMOD ;

} ;

#endif



