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

#ifndef PEL_TIMER_HH
#define PEL_TIMER_HH

#include <PEL_Object.hh>

#include <ctime>

/*
Timers for measuring and reporting the elapsed time passed
between start and stop events.

PUBLISHED
*/

class PEL_EXPORT PEL_Timer : public PEL_Object
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_Timer* create( PEL_Object * a_owner ) ;
      
   //-- Commands(1.)
      
      // Start or restart self.
      void start( void ) ;
      
      // Stop self.
      void stop( void ) ;

      // Stop self and nullify cumulative time.
      void reset( void ) ;

   //-- Access
      
      // cumulative user time spent between 
      // the `::start' and the `::stop' method calls
      double time( void ) const ;

      // cumulative wall clock time spent between 
      // the `::start' and the `::stop' method calls
      double elapsed_time( void ) const ;

      // Is self running ?
      bool is_running( void ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      static void print_time( double a_time, 
                              std::ostream& os, size_t indent_width ) ;

   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      PEL_Timer( void ) ;
     ~PEL_Timer( void ) ;
      PEL_Timer( PEL_Timer const& other ) ;
      PEL_Timer& operator=( PEL_Timer const& other ) ;

      PEL_Timer( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      double CUMUL_TIME ;
      double CURRENT_TIME ;
      bool RUNNING ;
      double CUMUL_ELAPSED_TIME ;
      double CURRENT_ELAPSED_TIME ;
      
} ;


#endif

