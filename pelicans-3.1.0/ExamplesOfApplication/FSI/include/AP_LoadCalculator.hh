#ifndef AP_LOAD_CALCULATOR_HH
#define AP_LOAD_CALCULATOR_HH

#include <FE_OneStepIteration.hh>

#include <map>

class GE_Color ;
class GE_QRprovider ;

class LA_PelMatrix ;
class LA_SeqVector ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;

class FE_Parameter ;

/*
PUBLISHED
*/

class AP_LoadCalculator : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_LoadCalculator* object( std::string const& a_name ) ;

   //-- Instance characteristics

      std::string const& name( void ) const ;

      PDE_DiscreteField const* fluid_velocity( void ) const ;

      PDE_DiscreteField const* structure_displacement( void ) const ;

   //-- Substeps of the step by step progression
      
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Results derived from the computations

      double fluid_load( size_t n, size_t ic ) const ;

      double solid_load( size_t n, size_t ic ) const ;

      void update_load( PDE_LinkDOF2Unknown const* link,
                        LA_SeqVector* load ) const ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_LoadCalculator( void ) ;
      AP_LoadCalculator( AP_LoadCalculator const& other ) ;
      AP_LoadCalculator& operator=(
                         AP_LoadCalculator const& other ) ;

      AP_LoadCalculator( PEL_Object* a_owner, 
		         PDE_SetOfDomains const* sdoms,
		         FE_SetOfParameters const* prms,
		         PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_LoadCalculator( void ) ;

      virtual AP_LoadCalculator* create_replica( 
                                       PEL_Object* a_owner,
				       PDE_DomainAndFields const* dom,
				       FE_SetOfParameters const* prms,
				       PEL_ModuleExplorer* exp ) const ;

      virtual AP_LoadCalculator* create_replica( 
                                       PEL_Object* a_owner,
				       PDE_SetOfDomains const* sdoms,
				       FE_SetOfParameters const* prms,
				       PEL_ModuleExplorer* exp ) const ;
   //-- Internals

      size_t idx( PDE_DiscreteField const* ff,
                  size_t n, size_t ic ) const ;

   //-- Class attributes

      static std::map< std::string, AP_LoadCalculator* > OBJS ;

      static AP_LoadCalculator const* PROTOTYPE ;

   //-- Attributes

      std::string NAME ;

      PDE_DiscreteField* SD ;
      PDE_DiscreteField const* FV ;
      size_t L_FV ;
      PDE_DiscreteField const* FP ;
      size_t L_FP ;
      FE_Parameter* MU ;

      size_t NB_COMPS ;
      size_t NB_DIMS ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEbound* bFE ;

      GE_Color const* SURFCOL ;

      LA_PelMatrix* S2F ;
      double MY_DBL_EPS ;
      double MY_DBL_MIN ;
      LA_SeqVector* F_LOAD ;
      LA_SeqVector* S_LOAD ;
} ;

#endif
