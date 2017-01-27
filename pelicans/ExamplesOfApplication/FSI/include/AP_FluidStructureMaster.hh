#ifndef AP_FLUID_STRUCTURE_MASTER_HH
#define AP_FLUID_STRUCTURE_MASTER_HH

#include <PEL_Application.hh>

#include <doubleVector.hh>

class PEL_ModuleExplorer ;
class PEL_Timer ;
class PEL_Vector ;

class LA_Vector ;

class GE_Color ;
class GE_SetOfPoints ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_GridMover ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_ResultSaver ;
class PDE_SetOfDomains ;

class FE_OneStepIteration ;
class FE_SetOfParameters ;
class FE_TimeIterator ;


/*
Applications involving a time marching procedure in which unknowns are
computed at each step from their values at the previous steps.

Each instance is associated to
   1. an instance of FE_TimeIterator that monitors the various steps
      of the time marching procedure ;
   2. an instance of `FE_OneStepIteration::' (constructed from the above 
      object) to which the progress of each step is delegated.

PUBLISHED 
*/

class AP_FluidStructureMaster : public PEL_Application
{
    public: //----------------------------------------------------------

    //-- Instance delivery and initialization

      static AP_FluidStructureMaster* create( 
                                         PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) ;

    //-- Instance characteristics

      FE_TimeIterator const* time_iterator( void ) const ;

      FE_SetOfParameters const* set_of_parameters( void ) const ;

    //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~AP_FluidStructureMaster( void ) ;
      AP_FluidStructureMaster( AP_FluidStructureMaster const& other ) ;
      AP_FluidStructureMaster& operator=( 
                               AP_FluidStructureMaster const& other ) ;

      AP_FluidStructureMaster( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_FluidStructureMaster( void ) ;

      virtual AP_FluidStructureMaster* create_replica( 
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      enum CouplingType { Weak, FixedPoint, FixedPointAitken } ;

      void update_FGD_at_interf( void ) ;

      void check_displacement_continuity( void ) ;

      void update_FV_at_interf( double dt ) ;

      void change_fluid_grid_position( void ) const ;

      bool check_convergence( void ) const ;

      bool check_closeness( double d_solid, double d_grid ) const ;

      double relaxation_coefficient( void ) const ;

      void save_for_restart( void ) ;

      void save_for_post_processing( void ) ;

      void do_save_for_post( PDE_ResultSaver* rs ) ;

      double next_time( doubleVector const& dates ) const ;

      void check_times( std::string const& list_name,
                        doubleVector const& dates ) const ;
      
      static bool greater_or_equal( double t1, double t2 ) ;

      static void print_memory_info( void ) ;

      void check_field_nb_components( PDE_DiscreteField const* ff,
                                      size_t required_nb_cmps ) const ;

      void check_field_storage_depth( PDE_DiscreteField const* ff,
                                      size_t requested_level ) const ;

   //-- Persistence
      
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static AP_FluidStructureMaster const* PROTOTYPE ;

      static size_t graphics_level ;

   //-- Attributes

      FE_TimeIterator* TIME_IT ;
      PDE_SetOfDomains* SDOMS ;
      PDE_DomainAndFields* F_DOM ;
      PDE_DomainAndFields* S_DOM ;
      FE_SetOfParameters* PRMS ;

      FE_OneStepIteration* FGRID ;
      FE_OneStepIteration* FLUID ;
      FE_OneStepIteration* LOADC ;
      FE_OneStepIteration* SOLID ;

      PDE_DiscreteField* FGD ;
      PDE_DiscreteField* FV ;
      PDE_DiscreteField* FP ;
      PDE_DiscreteField* SD ;
      PDE_DiscreteField* SV ;

      size_t L_NEW ;
      size_t L_EXPLICIT ;
      size_t L_EXPLICIT_EXPLICIT ;
      size_t L_OLD ;
      size_t BC_FV_INTERF ;

      PEL_Vector* VERTS_INI ;
      PEL_Vector* VERTS_NEW ;

      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;
      GE_Color const* INTERF_COL ;

      CouplingType TT ;
      size_t N_d_INTERF ;
      LA_Vector* DELTA_d_INTERF ;
      LA_Vector* DELTA_d_INTERF_OLD ;
      double OMEGA ;
      double TOL_EPS ;
      double TOL_MIN ;
      bool SAVE_ALL ;

      doubleVector graphics_times ;
      double graphics_next_time ;

      doubleVector saver_times ;
      double saver_next_time ;

      PEL_Timer* overall ;
      bool SAVEFG ;

      double DBL_EPS_GRID_POS ;
      double DBL_MIN_GRID_POS ;

      bool VERBOSE ;
} ;


#endif

