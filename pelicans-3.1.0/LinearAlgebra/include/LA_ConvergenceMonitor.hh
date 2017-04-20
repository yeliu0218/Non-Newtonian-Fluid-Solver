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

#ifndef LA_CONVERGENCE_MONITOR_HH
#define LA_CONVERGENCE_MONITOR_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class LA_ConvergenceTest ;
class LA_IterativeSolver ;
class LA_Matrix ;
class LA_Preconditioner ;
class LA_Vector ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT LA_ConvergenceMonitor : public PEL_Object
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      static LA_ConvergenceMonitor* make( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) ;
            
      virtual LA_ConvergenceMonitor* create_clone( 
                                            PEL_Object* a_owner ) const = 0 ;
      
   //-- Displays
      
      virtual void display_at_entry( LA_Matrix const* A, 
                                     LA_Vector const* b,
                                     LA_Preconditioner* prec,
                                     bool zero_initial_guess, 
                                     LA_Vector const* x,
                                     LA_ConvergenceTest const* cvgt ) = 0 ;
      
      virtual void monitor( size_t iter, double r_norm ) = 0 ;
      
      virtual void display_at_exit( LA_Matrix const* A, 
                                    LA_Vector const* b,
                                    LA_Preconditioner* prec,
                                    LA_Vector const* x,
                                    LA_ConvergenceTest const* cvgt ) = 0 ;
      
   protected: //--------------------------------------------------------
      
   //-- Plug in
      
      virtual ~LA_ConvergenceMonitor( void ) ;
      
      LA_ConvergenceMonitor( std::string const& a_name ) ;

      LA_ConvergenceMonitor( PEL_Object* a_owner ) ;
      
      LA_ConvergenceMonitor( PEL_Object* a_owner,
                             LA_ConvergenceMonitor const* other ) ;
      
      virtual LA_ConvergenceMonitor* create_replica(
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) const = 0 ;
      
      bool is_a_prototype( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      bool display_at_entry_PRE( LA_Matrix const* A, 
                                 LA_Vector const* b,
                                 LA_Preconditioner* prec,
                                 bool zero_initial_guess, 
                                 LA_Vector const* x,
                                 LA_ConvergenceTest const* cvgt ) const;
            
      bool display_at_exit_PRE( LA_Matrix const* A, 
                                LA_Vector const* b,
                                LA_Preconditioner* prec,
                                LA_Vector const* x,
                                LA_ConvergenceTest const* cvgt ) const ;
      
      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      
      bool create_replica_POST( LA_ConvergenceMonitor const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;
      
   private: //----------------------------------------------------------
      
      LA_ConvergenceMonitor( void ) ;
      LA_ConvergenceMonitor( LA_ConvergenceMonitor const& other ) ;
      LA_ConvergenceMonitor& operator=( LA_ConvergenceMonitor const& other ) ;
            
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const IS_PROTO ;
} ;

#endif
