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

#ifndef EXT_PETSc_API_HH
#define EXT_PETSc_API_HH

#include <PEL_ExternalAPI.hh>
#include <PEL_Timer.hh>
#include <PEL.hh>
#include <mpi.h>
class PEL_ModuleExplorer ;
// Some defs to allow PETSc tracelogs.
#define PETSC_USE_DEBUG 1
#define PETSC_USE_LOG 1
#define PETSC_USE_STACK 1
#include <mpi.h>

extern "C"
{
#include <petscao.h>
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscversion.h"
}

/*
PETSc applications, performing their specific initialization
and termination.

PUBLISHED
*/

#define PETSC_IMPLEMENTATION 3
class EXT_PETScAPI : public PEL_ExternalAPI
{

   public: //-----------------------------------------------------------

      static bool parse_options( PEL_ModuleExplorer const* exp,
                                 bool verbose ) ;
      static void going_to_do( char const* action ) ;
      static void verify( char const* action, int result ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScAPI( void ) ;
     ~EXT_PETScAPI( void ) ;
      EXT_PETScAPI( EXT_PETScAPI const& other ) ;
      EXT_PETScAPI& operator=( EXT_PETScAPI const& other ) ;

      static PEL_Timer* timer ;
   //-- Current instance management

      virtual void initialize( int& argc, char **& argv ) ;

   //-- Class attributes

      static EXT_PETScAPI* SINGLETON ;
} ;

#define PETSc_do(X) { EXT_PETScAPI::going_to_do( #X ) ; { PEL_Marker pspy( #X ) ; EXT_PETScAPI::verify( #X, X ) ; } }

#endif



