#ifndef AS_NavierStokes_HH
#define AS_NavierStokes_HH

#include <FE_OneStepIteration.hh>

#include <vector>
class PEL_ContextSimple ;
class PEL_DoubleVector ;
class GE_Point ;
class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;

class FE_LocalBCsBuilder ;
class ML_NavierStokesSystem ;
class FE_Parameter ;

/** \class AS_NavierStokes
    \brief
 *  Adjusted version of
 *  pelicans-3.1.0/ExamplesOfApplication/FSI/src/AP_FSINavierStokes
 *
 *  \note For the implementation structure see
 *  <a href="../../NavierStokes/doc/ImplStructure_AS_NavierStokes.html"> ImplStructure_AS_NavierStokes.html </a>
 *
PUBLISHED
*/

class AS_NavierStokes : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
   ///@name Substeps of the step by step progression
   //@{
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;
      virtual void do_before_inner_iterations_stage( FE_TimeIterator const* t_it ) ;
      virtual void do_inner_iterations_stage( FE_TimeIterator const* t_it );
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      virtual bool inner_iterations_are_completed( FE_TimeIterator const* t_it ) const;

      double compute_L2_norm( PDE_DiscreteField const * FIELD, size_t level1, size_t level2) const;
      double compute_flowrate( PDE_DiscreteField const * FIELD, size_t level1 ) const;


      void update_pressuredrop();
   //@}

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   ///@name Constructors, Destructors
   //@{
     ~AS_NavierStokes( void ) ;
      AS_NavierStokes( AS_NavierStokes const& other ) ;
      AS_NavierStokes& operator=( AS_NavierStokes const& other ) ;

      AS_NavierStokes( PEL_Object* a_owner,
                          PDE_DomainAndFields const* dom,
                          FE_SetOfParameters const* prms,
                          PEL_ModuleExplorer const* exp ) ;
   //@}

   //-- Plug in
   ///@name Plug In
   //@{
      AS_NavierStokes( void ) ;

      virtual AS_NavierStokes* create_replica(
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const ;
   //@}

   //-- Discrete system building
   ///@name Discrete System
   //@{
      void loop_on_cells_Outer( FE_TimeIterator const* t_it ) ;
      void loop_on_bounds_Outer( FE_TimeIterator const* t_it ) ;

      void loop_on_cells_Inner( FE_TimeIterator const* t_it ) ;
      void loop_on_bounds_Inner( FE_TimeIterator const* t_it ) ;

      void setup_Outer( FE_TimeIterator const* t_it ) ;
      void setup_Inner( FE_TimeIterator const* t_it ) ;
      void setup_Flowrate( FE_TimeIterator const* t_it ) ;
   //@}

   //-- Class attributes

      static AS_NavierStokes const* PROTOTYPE ;

   //-- Attributes
   //@name Discrete Fields
   //@{
      PDE_DiscreteField* UU ;   ///< @brief Velocity discrete field
      size_t L_UU ;             ///< @brief Level of explicite velocity
      size_t L_UPDATE_UU ;      ///< @brief Level of velocity term to update

      PDE_DiscreteField* PP ;   ///< @brief Pressure discrete field
      size_t L_PP ;             ///< @brief Level of explicite pressure
      size_t L_UPDATE_PP ;      ///< @brief Level of pressure to update

      PDE_DiscreteField* STRESS ;   ///< @brief  \f$\gamma \f$ discrete field
      PDE_DiscreteField* GAMMADOT ;   ///< @brief \f$\dot\gamma \f$ discrete field
      PDE_DiscreteField* GAMMA ;   ///< @brief  \f$\gamma \f$ discrete field

      std::vector< PDE_DiscreteField const* > AAs ; ///< @brief forms for the advective term to calculate \f$ u^\star \f$
      std::vector< size_t > L_AAs ;				    ///< @brief Levels of the advective terms
      std::vector< FE_Parameter* > COEF_AAs ;		///< @brief Coefficient of the advective terms -> \f$ \gamma_j\f $
   //@}
      bool ADV ;				///< @brief Advective term true?false

      size_t ORDER ;			///< @brief Order of the time discretisation

	///@name Parameters
	//@{
      FE_Parameter* ALPHA ;		///< @brief Advection parameter 0?1
      FE_Parameter* RE ;		///< @brief Reynolds number
      FE_Parameter* RHSU ;		///< @brief RHS of the momentum equation
      FE_Parameter* KAPPA;		///< @brief nondimensional viscosity term
      FE_Parameter* BN;         ///< @brief Bingham number
      FE_Parameter* HB_n;         ///< @brief Herschel Bulkley coefficient (n), if n=1 we have Bingham
      FE_Parameter* E_FIELD;    ///< @brief electrical field E=-grad(phi)
	//@}
      bool LAPL_UU ;

   ///@name Nonstandard Boundary Conditions:
   //@{
      PDE_SetOfBCs const* BCs ;
      bool BC_STRESS ;
      FE_Parameter* BC_STRESS_VALUE;
      FE_LocalBCsBuilder* LOCAL_BC ;
   //@}

   ///@name Discrete system builder
   //@{
      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      ML_NavierStokesSystem* GLOBAL_EQ ;
   //@}

    ///@name Piccard Iteration Parameters
	//@{
    bool Viscoplastic;              ///< @brief Viscoplastic Problem ?
    std::string ViscoplasticType;   ///< @brief Type of the nonlinearity
    int ViscoplasticMaxiter;        ///< @brief Maximum number of inner iterations
    double ViscoplasticResidual;    ///< @brief Residual of the inner iterations
    double ViscoplasticAug;         ///< @brief Augmentation parameter (r)
    int ViscoplasticIter;           ///< @brief Nonlinear iteration counter
    bool FACTORIZE;                       ///< @brief Usage of the Umfpack solver ?
    size_t L_ViscoplasticU;         ///< @brief Storage level of the velocity field where the old solution is stored for convergence analysis.
	//@}
	///@name Flowrate Parameters:
	//@{
	bool Flowrate;				///< @brief Flowrate ?
	size_t FlowrateCP;			///< @brief Component to be used.
	size_t FlowrateMaxiter;		///< @brief Maximal number of iterations
	double FlowrateTarget;		///< @brief Target flowrate
	double FlowrateResidual;	///< @brief Residual
	std::vector<double> FlowratedPdL;	///< @brief Vector containing the pressuredrops
	std::vector<double>  FlowrateFr;	///< @brief Vector containing the flowrates.
	double FlowrateScale;		///< @brief Scale for the flowrate calculation
	std::string FlowrateMethod;	///< @brief Method for flowrate calculation
	//@}
	PEL_Communicator const* const COM ; ///< @brief PEL communicator
	PEL_ContextSimple* CONTEXT ;///< @brief Context used to determine the value of the stress BC
	PEL_DoubleVector* COORDS; 	///< @brief Coordinate coefficients for calculating the norms.

    bool ElectricField; /// @brief electric field ?
    double FUNC_E; /// @brief f(|E|) used with electrical field

    double param_euclidian_norm( FE_TimeIterator const* t_it , FE_Parameter const * param) const;
    double param_yield_stress_function( double const param_norm) const ;

} ;

#endif
