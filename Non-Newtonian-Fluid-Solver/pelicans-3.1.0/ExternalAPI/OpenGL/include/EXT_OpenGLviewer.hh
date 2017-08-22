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

#ifndef EXT_OpenGLviewer_HH
#define EXT_OpenGLviewer_HH

#include <PEL_DataOnMeshingWriter.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

#include <PEL.hh>
#include <PEL_Application.hh>

#include <sys/types.h>
#include <sys/stat.h>

#include <EXT_OpenGLwindow.hh>

/*
Visualization of PELICANS data file with OpenGL windows.
*/

class EXT_OpenGLviewer : public PEL_Application
{
   public: //-----------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;
      virtual void go_to( EXT_OpenGLwindow::CycleMove relative_cycle_number ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~EXT_OpenGLviewer(void ) ;
      EXT_OpenGLviewer( EXT_OpenGLviewer const& other ) ;
      EXT_OpenGLviewer& operator=( EXT_OpenGLviewer const& other ) ;

      EXT_OpenGLviewer( PEL_Object* a_owner,
                        stringVector& args ) ;

   //-- Plug in

      EXT_OpenGLviewer( void ) ;

      virtual EXT_OpenGLviewer* create_replica(
                                   PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const ;

      virtual EXT_OpenGLviewer* create_replica_from_args(
                                        PEL_Object* a_owner,
                                        stringVector& args ) const ;

      PEL_ModuleExplorer const* read_cycle( PEL_Object * a_owner,
                                            size_t cycle_number ) ;


   //-- Write

      void write_cycle( PEL_ModuleExplorer const* exp ) ;

      void write_grid( PEL_ModuleExplorer const* exp ) ;

      void write_field( int cycle_number,
                        PEL_ModuleExplorer const* exp ) ;

      void write_integration_domain( PEL_ModuleExplorer const* exp )  ;

   //-- Command line

      virtual void print_usage( void ) const ;
      virtual void print_operands( void ) const ;
      virtual void print_exit_status( void ) const ;

   //-- Class attributes

      static EXT_OpenGLviewer const* PROTOTYPE ;

   //-- Attributes

      std::string filename ;
      EXT_OpenGLwindow * ogl ;
      doubleArray2D vertices ;
      intArray2D connectivity ;
      bool first ;
      off_t last_modif ;
      size_t last_pos ;
      size_t last_displayed_cycle ;
      size_t last_requested_cycle ;


};

#endif
