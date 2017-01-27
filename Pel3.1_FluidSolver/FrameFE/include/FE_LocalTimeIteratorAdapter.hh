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

#ifndef FE_LocalTimeIteratorAdapter_HH
#define FE_LocalTimeIteratorAdapter_HH

#include <PEL_Object.hh>

class PDE_DomainAndFields ;
class PDE_SetOfDomains ;

class PDE_SetOfDiscreteFields ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class FE_SetOfParameters ;
class FE_TimeIterator ;

/**@brief (Almost) identical to FE_TimeIteratorAdapter 
 *
 * Difference: The routine 
 * restart_iteration_with_new_time_step( double dt )
 * and the destructor 
 * are public instead of protected.
 * 
PUBLISHED*/

class PEL_EXPORT FE_LocalTimeIteratorAdapter : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return the concrete instance of `FE_LocalTimeIteratorAdapter::'
      // linked to `t_it' and of name `exp'->string_data( "concrete_name" )
      // for single domain applications.
      static FE_LocalTimeIteratorAdapter* make(
                                FE_TimeIterator* t_it,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer const* exp ) ;
      
      // Create and return the concrete instance of `FE_LocalTimeIteratorAdapter::'
      // linked to `t_it' and of name `exp'->string_data( "concrete_name" )
      // for multi domains applications.
      static FE_LocalTimeIteratorAdapter* make(
                                FE_TimeIterator* t_it,
                                PDE_SetOfDomains const* sdoms,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer const* exp ) ;

      // Create and return the default `FE_LocalTimeIteratorAdapter::'
      // linked to `t_it' for single domain applications.
      static FE_LocalTimeIteratorAdapter* make_default(
                                FE_TimeIterator* t_it,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms ) ;
      
      // Create and return the default `FE_LocalTimeIteratorAdapter::'
      // linked to `t_it' for multi domains applications.
      static FE_LocalTimeIteratorAdapter* make_default(
                                FE_TimeIterator* t_it,
                                PDE_SetOfDomains const* sdoms,
                                FE_SetOfParameters const* prms ) ;

   //-- Associated FE_TimeIterator object(1.0)

      FE_TimeIterator const* time_iterator( void ) const ;
      
   //-- Adaptation of the associated FE_TimeIterator object(2.0)

      // Initialize for a new time step.
      void initialize_time_step( void ) ;

      // Adapt new time iterator.
      void adapt_time_iterator( void ) ;
      
      // IMPLEMENTATION : raise error and terminate program.
      virtual void set_time_iteration_failed( void ) ;

   //-- Proposed modifications to the associated FE_TimeIterator object(2.2)

      // Propose `dt' for the next time step of `::time_iterator()' 
      // (`::time_iterator()' is not modified and the proposed modification 
      // will be effective only after calling `::adapt_time_iterator()').
      void propose_next_time_step( double dt ) ;
      
      // value currently proposed for the next time step of `::time_iterator()'
      double next_time_step( void ) const ;

      // Propose to terminate current time iterations
      // (`::time_iterator()' is not modified and the proposed modification 
      // will be effective only after calling `::adapt_time_iterator()').
      void propose_to_finish_iterations( void ) ;
      
      // Are the current time iterations terminated ?
      bool iterations_are_finished( void ) const ;
      
      // Restart time iteration with `dt' as new time step :
      // should be called by `::set_time_iteration_failed'().
      void restart_iteration_with_new_time_step( double dt ) ;

   //-- Persistence

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;
      
      virtual ~FE_LocalTimeIteratorAdapter( void ) ;

   protected: //--------------------------------------------------------

   //-- Plug in

      FE_LocalTimeIteratorAdapter( std::string const& a_concrete_name ) ;

      FE_LocalTimeIteratorAdapter( FE_TimeIterator* t_it ) ;
      
      virtual FE_LocalTimeIteratorAdapter* create_replica( 
                                     FE_TimeIterator* t_it,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp ) const ;
      virtual FE_LocalTimeIteratorAdapter* create_replica( 
                                     FE_TimeIterator* t_it,
                                     PDE_SetOfDomains const* sdoms,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp ) const ;
      
      bool is_a_prototype( void ) const ;

   //-- Default adaptator

      FE_LocalTimeIteratorAdapter( void ) ;
      virtual FE_LocalTimeIteratorAdapter* create_replica(
                                     FE_TimeIterator* t_it,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms ) const ;
      virtual FE_LocalTimeIteratorAdapter* create_replica(
                                     FE_TimeIterator* t_it,
                                     PDE_SetOfDomains const* sdoms,
                                     FE_SetOfParameters const* prms ) const ;
      
   //-- Adaptation of the associated FE_TimeIterator object

      // Inner initialization called by `::initialize_time_step'().
      // IMPLEMENTATION : do nothing.
      virtual void initialize_inner( void ) ;
      
      /*
      Indicate the desired parameters for the next time iteration.
      On exit, the meaning of the arguments is the following :
         - `restart'=true: restart current iteration with 
                          `next_dt' as time step
         - `finished'=true: time iterations are finished
         - `finished'=false: perform next iteration with time step `next_dt'.
      IMPLEMENTATION : arguments are unchanged.
      */
      virtual void define_parameters_for_next_iteration(
                           bool& finished, bool& restart, double& next_dt ) ;
      
      // new time step from the one proposed by the client
      // IMPLEMENTATION : `dt' is returned without modification
      virtual double next_time_step_from_proposed_one( double dt ) const ;


   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool set_time_iteration_failed_PRE( void ) const ;
      virtual bool set_time_iteration_failed_POST( void ) const ;

      virtual bool initialize_inner_PRE( void ) const ;
      virtual bool initialize_inner_POST( void ) const ;
      
      virtual bool define_parameters_for_next_iteration_PRE( 
                       bool finished, bool restart, double next_dt ) const ;
      virtual bool define_parameters_for_next_iteration_POST(
                       bool finished, bool restart, double next_dt ) const ;
      
      virtual bool next_time_step_from_proposed_one_PRE( double dt ) const ;
      virtual bool next_time_step_from_proposed_one_POST(
                                          double result, double dt ) const ;
      
      virtual bool create_replica_PRE(
                               FE_TimeIterator const* t_it,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST(
                               FE_LocalTimeIteratorAdapter const* result,
                               FE_TimeIterator const* t_it,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      
      virtual bool create_replica_PRE(
                               FE_TimeIterator const* t_it,
                               PDE_SetOfDomains const* sdoms,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST(
                               FE_LocalTimeIteratorAdapter const* result,
                               FE_TimeIterator const* t_it,
                               PDE_SetOfDomains const* sdoms,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      
      virtual bool create_replica_PRE(
                               FE_TimeIterator const* t_it,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms ) const ;
      virtual bool create_replica_POST(
                               FE_LocalTimeIteratorAdapter const* result,
                               FE_TimeIterator const* t_it,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms ) const ;
      
      virtual bool create_replica_PRE(
                               FE_TimeIterator const* t_it,
                               PDE_SetOfDomains const* sdoms,
                               FE_SetOfParameters const* prms ) const ;
      virtual bool create_replica_POST(
                               FE_LocalTimeIteratorAdapter const* result,
                               FE_TimeIterator const* t_it,
                               PDE_SetOfDomains const* sdoms,
                               FE_SetOfParameters const* prms ) const ;

      
   private: //----------------------------------------------------------
      
      FE_LocalTimeIteratorAdapter( FE_LocalTimeIteratorAdapter const& other ) ;
      FE_LocalTimeIteratorAdapter& operator=(
                              FE_LocalTimeIteratorAdapter const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Static attributes

      static FE_LocalTimeIteratorAdapter const* DEFAULT_PROTOTYPE ;

   //-- Attributes

      bool const IS_PROTO ;
      FE_TimeIterator* const T_IT ;
      double DT ;
      bool FINISHED ;
} ;

#endif
