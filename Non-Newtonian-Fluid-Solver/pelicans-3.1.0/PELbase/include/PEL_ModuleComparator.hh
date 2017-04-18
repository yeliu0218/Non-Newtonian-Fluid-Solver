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

#ifndef PEL_MODULE_COMPARATOR_HH
#define PEL_MODULE_COMPARATOR_HH

#include <PEL_Object.hh>

#include <string>

class PEL_Module;
class PEL_ModuleExplorer ;
class PEL_Context;
class PEL_Data;

class PEL_EXPORT PEL_ModuleComparator : public PEL_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Creates and return an instance containing neither modules nor entries.
      static PEL_ModuleComparator* create( PEL_Object* a_owner, 
                                           PEL_ModuleExplorer const* exp ) ;

      int compare( PEL_Module const* m1, 
                   PEL_Module const* m2, 
                   PEL_Module* result ) ;

   protected: //-------------------------------------------------------

   private: //-------------------------------------------------------

      PEL_ModuleComparator( void ) ;
     ~PEL_ModuleComparator( void ) ;
      PEL_ModuleComparator( PEL_ModuleComparator const& other ) ;
      PEL_ModuleComparator const& operator=( 
                            PEL_ModuleComparator const& other ) ;

      PEL_ModuleComparator( PEL_Object* a_owner, 
                            PEL_ModuleExplorer const* exp ) ;

      bool is_valid_module( std::string const& name ) ;
      bool is_valid_data( std::string const& name, std::string const& abs_path_name ) ;
      bool is_verbose( void ) ;
      int internalCompare( PEL_Module const * m1, 
                           PEL_Module const * m2, 
                           PEL_Module* result, 
                           std::string const& path) ;
      int compare( std::string const& dataname, 
                   PEL_Module const* m1, 
                   PEL_Module const* m2, 
                   PEL_Module* result) ;
      int compare( double v1, double v2, double& status ) ;
      int compare( bool v1, bool v2 ) ;
      int compare( int v1, int v2 ) ;
      int compare( std::string const& v1, std::string const& v2 ) ;

   //-- Attributes

      std::string left;
      std::string right;
      bool verbose;
      class stringVector *valid_module;
      class stringVector *valid_data;
      class stringVector *ignore_data;
      double MY_DBL_EPS ;
      double MY_DBL_MIN ;
      
};

#endif
