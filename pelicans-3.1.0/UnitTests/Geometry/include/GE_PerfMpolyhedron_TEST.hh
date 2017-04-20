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

#ifndef GE_PERF_M_POLYHEDRON_TEST_HH
#define GE_PERF_M_POLYHEDRON_TEST_HH

#include <PEL_ObjectTest.hh>

class PEL_ModuleExplorer ;
class PDE_DomainAndFields ;
class PDE_LocalFE ;

class PEL_EXPORT GE_PerfMpolyhedron_TEST : public PEL_ObjectTest
{
   public: //-----------------------------------------------------------------

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      GE_PerfMpolyhedron_TEST( void ) ;
     ~GE_PerfMpolyhedron_TEST( void ) ;
      GE_PerfMpolyhedron_TEST( GE_PerfMpolyhedron_TEST const& other ) ;
      GE_PerfMpolyhedron_TEST& operator=( 
                            GE_PerfMpolyhedron_TEST const& other ) ;

   //-- Elementary tests management

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

   //-- Internals

      void do_mesh_go_next( PDE_LocalFE* fe,
                            PEL_ModuleExplorer const* exp,
                            std::string const& mesg ) ;

      void do_mesh_go_i_th( PDE_LocalFE* fe,
                            PEL_ModuleExplorer const* exp,
                            std::string const& mesg ) ;

      void do_mesh_iteration_on_IPs( PDE_LocalFE* fe,
                                     PEL_ModuleExplorer const* exp,
                                     std::string const& mesg ) ;

   //-- Class attributes

      static GE_PerfMpolyhedron_TEST const* REGISTRATOR ;

   //-- Attributes

      std::string TEST_NAME ;
} ;

#endif
