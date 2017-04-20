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

#ifndef FE_ONE_STEP_ITERATION_HH
#define FE_ONE_STEP_ITERATION_HH

#include <PEL_Object.hh>

#include <PDE_SetOfBCs.hh>

#include <map>
#include <set>

class PEL_Communicator ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;
class PEL_ObjectWriter ;
class PEL_Timer ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalFE ;
class PDE_ResultSaver ;
class PDE_SetOfDomains ;
class PDE_SystemNumbering ;

class FE_Parameter ;
class FE_SetOfParameters ;
class FE_TimeIteratorAdapter ;
class FE_TimeIterator ;

/*
Objects decomposing the steps performed in `FE_StepByStepProgression::run'
into substeps possibly involving inner iterations. The substeps are intented
to perform some computations related to the solution of a system of partial
differential equations. The major components required to perform this task
can be delivered by an associated `PDE_DomainAndFields::' instance.

FRAMEWORK INSTANTIATION

   CASE 1 : derivation of a concrete subclass

   1. Derive a concrete subclass, say MySt.
   2. Choose a name for MyAppli, say "my_st".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `FE_OneStepIteration::' subobject by calling
               `FE_OneStepIteration( std::string const& )'
          with "my_st" as argument.
          Example of pseudo-code :
          | MySt:: MySt( void ) : FE_OneStepIteration( "my_st" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MySt.hh) :
             | static MySt const* PROTOTYPE ;
             definition (in the implementation file, eg MySt.cc) :
             | MySt const* MySt::PROTOTYPE = new MySt() ;'
   6. Implement a private constructor that initializes the
      `FE_OneStepIteration::' subobject by calling
             `FE_OneStepIteration( PEL_Object*, PDE_DomainAndFields const*, PEL_ModuleExplorer const* )'.
      Example of pseudo-code :
      | MySt:: MySt( PEL_Object* a_owner,
      |              PDE_DomainAndFields const* dom,
      |              FE_SetOfParameters const* prms,
      |              PEL_ModuleExplorer const* exp )
      |    : FE_OneStepIteration( a_owner, dom, exp ), ...
      | { ... }
   7. Implement the `::create_replica' method that allocates an object
      of type `MySt' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | MySt* MySt::create_replica( PEL_Object* a_owner,
      |                             PDE_DomainAndFields const* dom,
      |                             FE_SetOfParameters const* prms,
      |                             PEL_ModuleExplorer* exp ) const
      | {
      |    PEL_LABEL( "MySt::create_replica" ) ;
      |    PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;
      |    MySt* result = new MySt( a_owner, dom, prms, exp ) ;
      |    PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ;
      |    return result ;
      | }
   8. Implement the pure virtual methods, and possibly overrid the non
      pure virtual methods.


   CASE 2 : derivation of an abstract subclass

   1. Derive an abstract subclass, say MySt.
   2. Implement a protected virtual destructor.
   3. Implement a protected constructor that initializes the
      `FE_OneStepIteration::' subobject by calling
               `FE_OneStepIteration( std::string const& )'
      Example of pseudo-code :
      | MySt:: MySt( std::string const& name )
      |    : FE_OneStepIteration( name ) {}
      This constructor is devoted to be used by the concrete subclasses
      of MySt for the registration of their prototype.
   4. Implement a protected constructor that initializes the
      `PEL_Application::' subobject by calling
            `FE_OneStepIteration( PEL_Object*, PDE_DomainAndFields const*, PEL_ModuleExplorer const* )'.
      Example of pseudo-code :
      | MySt:: MySt( PEL_Object* a_owner,
      |              PDE_DomainAndFields const* dom,
      |              PEL_ModuleExplorer const* exp )
      |    : FE_OneStepIteration( a_owner, dom, exp ), ...
      | { ... }
      This constructor is devoted to be used to initialize the MySt
      base class subobject when creating objects of concrete subclasses
      of MySt (such creations are performed in the `create_replica::'
      method whose implementation is deferred into those concrete subclasses).

PUBLISHED
*/

class PEL_EXPORT FE_OneStepIteration : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance according to the data attainable by
      // `exp', devoted to perform computation associated to a system
      // of partial differential equations whose numerical solution requires
      // objects delivered by `dom'.
      static FE_OneStepIteration* make( PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) ;

      static FE_OneStepIteration* make( PEL_Object* a_owner,
                                        PDE_SetOfDomains const* sdoms,
                                        FE_SetOfParameters const* prms,
                                        PEL_ModuleExplorer* exp ) ;

   //-- Substeps of the step by step progression

      // Before starting time marching in `FE_StepByStepProgression::run',
      // perform initial computations.
      // IMPLEMENTATION : do nothing.
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      // Within a step of the time marching procedure in
      // `FE_StepByStepProgression::run', perform an initial stage
      // just before start of the inner iterations.
      // IMPLEMENTATION : do nothing.
      virtual void do_before_inner_iterations_stage(
                                           FE_TimeIterator const* t_it ) ;

      virtual void do_inner_iterations_stage( FE_TimeIterator const* t_it ) ;

      bool inner_iterations_stage_failed( void ) const ;

      // Within a step of the time marching procedure in
      // `FE_StepByStepProgression::run', perform an inner iteration stage.
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) = 0 ;

      // Are inner iterations performed in `::do_one_inner_iteration'
      // completed ?
      // IMPLEMENTATION : true
      virtual bool inner_iterations_are_completed(
                                     FE_TimeIterator const* t_it ) const ;

      // Within a step of the time marching procedure in
      // `FE_StepByStepProgression::run', perform a final stage just
      // after completion of the inner iterations.
      // IMPLEMENTATION : do nothing.
      virtual void do_after_inner_iterations_stage(
                                     FE_TimeIterator const* t_it ) ;

      virtual void do_after_time_adaptation(
                                     FE_TimeIterator const* t_it ) ;

      // After completion of the time marching procedure in
      // `FE_StepByStepProgression::run', perform final computations.
      // IMPLEMENTATION : do nothing.
      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      static void reset_standard_times( void ) ;

      static void print_standard_times( std::ostream& os,
                                        size_t indent_width ) ;

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Time iterator modification

      virtual void adapt_time_iterator( FE_TimeIteratorAdapter* t_adapter ) ;

      virtual void notify_inner_iterations_stage_failure( void ) ;

      virtual void reset_after_inner_iterations_stage_failure( void ) ;

   //-- Savings for post-processing

      // Use `rs' to save data other than time and fields at the current time
      // step (attainable through `t_it') in the time marching procedure
      // of `FE_StepByStepProgression::run'.
      // IMPLEMENTATION : do nothing.
      virtual void do_additional_savings( FE_TimeIterator const* t_it,
                                          PDE_ResultSaver* rs ) ;

   //-- Persistence

      // Extend `list' so that it contains all objects required by the
      // storage and retrieval mechanisms concerning the `FE_OneStepIteration::'
      // base class subobject, then call `::add_storable_objects'.
      void register_storable_objects( PEL_ListIdentity* list ) ;

      // Extend `list' so that it contains all objects required by the
      // storage and retrieval mechanisms that are not part of the
      // `FE_OneStepIteration::' base class subobject.
      // IMPLEMENTATION : do nothing, i.e. leave `list' unchanged.
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      static std::string const& indent( void ) ;

      static void increase_indent( void ) ;

      static void decrease_indent( void ) ;

   protected: //--------------------------------------------------------------

      virtual ~FE_OneStepIteration( void ) ;

      // Construction of an instance whose owner is `a_owner'.
      // On exit, `self' still refers to the object identified by `dom'
      FE_OneStepIteration( PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           PEL_ModuleExplorer const* exp ) ;

      FE_OneStepIteration( PEL_Object* a_owner,
                           PDE_SetOfDomains const* sdoms,
                           PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      // for prototype registration only
      FE_OneStepIteration( std::string const& name ) ;

      virtual FE_OneStepIteration* create_replica(
	                                PEL_Object* a_owner,
	                                PDE_DomainAndFields const* dom,
	                                FE_SetOfParameters const* prms,
	                                PEL_ModuleExplorer* exp ) const = 0 ;

      virtual FE_OneStepIteration* create_replica(
                                   PEL_Object* a_owner,
                                   PDE_SetOfDomains const* sdoms,
                                   FE_SetOfParameters const* prms,
	                                PEL_ModuleExplorer* exp ) const ;

      bool is_a_prototype( void ) const ;

      void configure_multilevel_preconditioner(
                                   PEL_ModuleExplorer const* exp,
                                   std::string const& keyword,
                                   PDE_DomainAndFields const* dom,
                                   PDE_SystemNumbering const* nmb ) const ;

   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields(
                                                  FE_TimeIterator const* t_it,
                                                  PDE_ResultSaver* rs ) ;

   //-- Distributed processing

      static PEL_Communicator const* communicator( void ) ;

   //-- Timer and screen output

      size_t verbose_level( void ) const ;

      void start_assembling_timer( bool silent = false ) const ;
      void stop_assembling_timer( bool silent = false ) const ;

      void start_solving_timer( bool silent = false ) const ;
      void stop_solving_timer( bool silent = false ) const ;

      void start_total_timer( std::string const& mesg,
                              bool silent = false ) const ;
      void stop_total_timer( bool silent = false ) const ;

   //-- Errors

      static void check_param_nb_components( FE_Parameter const* prm,
                                             std::string const& keyword,
                                             size_t required_nb_cmps ) ;

      static void check_field_nb_components( PDE_DiscreteField const* ff,
                                             size_t required_nb_cmps ) ;

      void check_field_storage_depth( PDE_DiscreteField const* ff,
                                      size_t requested_max_level ) const ;

      void check_nb_local_nodes( PDE_DiscreteField const* ff,
                                 PDE_LocalFE const* fe,
                                 size_t required_nb_nodes ) const ;

      void check_boundary_condition(
                          PDE_SetOfBCs const* bcs,
                          GE_Color const* color,
                          PDE_DiscreteField const* field,
                          size_t ic = PDE_SetOfBCs::all_components ) const ;

      void raise_bad_BC_type(
                          std::string const& type,
                          std::string const& allowed_types,
                          PDE_DiscreteField const* field,
                          size_t ic = PDE_SetOfBCs::all_components ) const ;

   //-- Preconditions, Postconditions, Invariant

      bool do_before_time_stepping_PRE( FE_TimeIterator const* t_it ) const ;

      bool do_before_inner_iterations_stage_PRE(
                                       FE_TimeIterator const* t_it ) const ;
      bool do_before_inner_iterations_stage_POST(
                                       FE_TimeIterator const* t_it ) const ;

      bool do_inner_iterations_stage_PRE(
                                       FE_TimeIterator const* t_it ) const ;

      bool do_one_inner_iteration_PRE( FE_TimeIterator const* t_it ) const ;

      bool inner_iterations_are_completed_PRE(
                                       FE_TimeIterator const* t_it ) const ;

      bool do_after_inner_iterations_stage_PRE(
                                       FE_TimeIterator const* t_it ) const ;

      bool do_after_time_adaptation_PRE(
                                       FE_TimeIterator const* t_it ) const ;

      bool notify_inner_iterations_stage_failure_PRE( void ) const ;
      bool notify_inner_iterations_stage_failure_POST( void ) const ;

      bool reset_after_inner_iterations_stage_failure_PRE( void ) const ;
      bool reset_after_inner_iterations_stage_failure_POST( void ) const ;

      bool adapt_time_iterator_PRE(
                           FE_TimeIteratorAdapter const* t_adapter ) const ;

      bool do_additional_savings_PRE( FE_TimeIterator const* t_it,
                                      PDE_ResultSaver* rs ) const ;

      bool save_other_than_time_and_fields_PRE(
                                       FE_TimeIterator const* t_it,
                                       PDE_ResultSaver* rs ) const ;

      bool add_storable_objects_PRE( PEL_ListIdentity* list ) const ;

      bool create_replica_PRE( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( FE_OneStepIteration const* result,
                                PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer const* exp ) const ;

      bool create_replica_PRE( PEL_Object* a_owner,
                               PDE_SetOfDomains const* sdoms,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( FE_OneStepIteration const* result,
                                PEL_Object* a_owner,
                                PDE_SetOfDomains const* sdoms,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer const* exp ) const ;

   private: //----------------------------------------------------------------

      FE_OneStepIteration( void ) ;
      FE_OneStepIteration( FE_OneStepIteration const& other ) ;
      FE_OneStepIteration& operator=( FE_OneStepIteration const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

      void start_internal_timer( std::map< std::string, PEL_Timer* >& timers,
                                 std::string const& display,
                                 bool do_display ) const ;
      void stop_internal_timer(
                         std::map< std::string, PEL_Timer* > const& timers,
                         bool do_display ) const ;

   //-- Class attributes

      static std::string INDENT ;
      static std::map< std::string, PEL_Timer* > ASS_TIMER ;
      static std::map< std::string, PEL_Timer* > SOL_TIMER ;
      static std::map< std::string, PEL_Timer* > TOT_TIMER ;

   //-- Attributes

      bool IS_PROTO ;
      std::set< PDE_DomainAndFields const* > DOMS ;
      bool FAILURE ;
      size_t VERBOSE_LEVEL ;

} ;

#endif
