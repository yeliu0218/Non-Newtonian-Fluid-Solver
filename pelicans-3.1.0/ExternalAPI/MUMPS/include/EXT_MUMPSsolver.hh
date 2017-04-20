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

#ifndef EXT_MUMPSsolver_HH
#define EXT_MUMPSsolver_HH

#include <doubleVector.hh>
#include <intVector.hh>

#include <LA_Solver.hh>

#include <EXT_MUMPS_API.hh>

/*
 * Wrappers around MUMPS direct solver.
 * Data deck instanciation :
 
 MODULE LA_Solver
    concrete_name = "EXT_MUMPSsolver"
    [ verbose = <boolean> ]
    [ MUMPS_verbosity = -1, 0, 1, 2, 3, 4 ]         // Activate MUMPS (very) verbose mode
    [ out_of_core = <boolean> ]                     // Allow MUMPS to use disk to store part of assembled matrix

    [ icntl1 =  <int> ]
    ...
    [ icntl40 = <int> ]

    [ cntl1 = <double> ]
    ...
    [ cntl6 = <double> ]
    
 END MODULE LA_Solver

 For details about icntl and cntl control parameters, let read MUMPS User's guide.
 
*/

class EXT_MUMPSsolver : public LA_Solver
{
   public: //-----------------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual EXT_MUMPSsolver* create_clone( PEL_Object* a_owner ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~EXT_MUMPSsolver( void ) ;
      EXT_MUMPSsolver( EXT_MUMPSsolver const& other ) ;
      EXT_MUMPSsolver& operator=( EXT_MUMPSsolver const& other ) ;

      EXT_MUMPSsolver( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;
      EXT_MUMPSsolver( PEL_Object* a_owner, EXT_MUMPSsolver const* other ) ;

      bool do_the_job( int job ) ;
      
      int& icntl( size_t i ) ;

      int& infog( size_t i ) ;
      
      int& info( size_t i ) ;
      
      double& cntl( size_t i ) ;
      
   //-- Plug in

      EXT_MUMPSsolver( void ) ;

      virtual EXT_MUMPSsolver* create_replica( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution
      
      virtual void set_matrix_self( LA_Matrix const* mat, bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
      virtual void unset_matrix_self( void ) ;

      void initialize( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Class attributes

      static EXT_MUMPSsolver const* PROTOTYPE ;
      
    //-- Attributes
      
      PEL_ModuleExplorer const* const EXP ;
      DMUMPS_STRUC_C SOLVER ;
      intVector EXTRA_ICNTL ;
      doubleVector EXTRA_CNTL ;
} ;

#endif
