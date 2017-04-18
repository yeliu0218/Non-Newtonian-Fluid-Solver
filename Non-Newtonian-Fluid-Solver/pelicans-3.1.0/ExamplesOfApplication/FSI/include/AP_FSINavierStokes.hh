#ifndef AP_FSI_NAVIER_STOKES_HH
#define AP_FSI_NAVIER_STOKES_HH

#include <FE_OneStepIteration.hh>

#include <vector>

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;

class FE_LocalBCsBuilder ;
class AP_FSINavierStokesSystem ;
class FE_Parameter ;

/*
PUBLISHED
*/

class AP_FSINavierStokes : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_inner_iterations_stage( 
                                              FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_FSINavierStokes( void ) ;
      AP_FSINavierStokes( AP_FSINavierStokes const& other ) ;
      AP_FSINavierStokes& operator=( AP_FSINavierStokes const& other ) ;

      AP_FSINavierStokes( PEL_Object* a_owner, 
                          PDE_DomainAndFields const* dom,
                          FE_SetOfParameters const* prms,
                          PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_FSINavierStokes( void ) ;

      virtual AP_FSINavierStokes* create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;

   //-- Class attributes

      static AP_FSINavierStokes const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* UU ;
      size_t L_UU ;
      size_t L_UPDATE_UU ;

      PDE_DiscreteField* PP ;
      size_t L_PP ;
      size_t L_UPDATE_PP ;

      std::vector< PDE_DiscreteField const* > AAs ;
      std::vector< size_t > L_AAs ;
      std::vector< FE_Parameter* > COEF_AAs ;
      bool ADV ;

      size_t ORDER ;

      FE_Parameter* ALPHA ;
      FE_Parameter* MU ;
      FE_Parameter* RHSU ;
      bool LAPL_UU ;

      PDE_SetOfBCs const* BCs ;

      FE_LocalBCsBuilder* LOCAL_BC ;

      PDE_LocalEquation* ELEMENT_EQ ;

      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      AP_FSINavierStokesSystem* GLOBAL_EQ ;
} ;

#endif
