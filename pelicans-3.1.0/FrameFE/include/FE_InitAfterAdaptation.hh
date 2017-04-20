#ifndef FE_INIT_AFTER_ADAPTATION_HH
#define FE_INIT_AFTER_ADAPTATION_HH

#include <FE_OneStepIteration.hh>
#include <vector>

class GE_QRprovider ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEcell ;
class PDE_SystemNumbering ;

/*
PUBLISHED
*/

class PEL_EXPORT FE_InitAfterAdaptation : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~FE_InitAfterAdaptation( void ) ;
      FE_InitAfterAdaptation( FE_InitAfterAdaptation const& other ) ;
      FE_InitAfterAdaptation& operator=( 
                              FE_InitAfterAdaptation const& other ) ;

      FE_InitAfterAdaptation( PEL_Object* a_owner, 
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_InitAfterAdaptation( void ) ;

      virtual FE_InitAfterAdaptation* create_replica( 
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_cells( FE_TimeIterator const* t_it, size_t i ) ;

   //-- Class attributes

      static FE_InitAfterAdaptation const* PROTOTYPE ;

   //-- Attributes

      std::vector< PDE_DiscreteField* > FIELDs ;
      std::vector< size_t > L_FIELDs ;
      std::vector< PDE_DiscreteField const* > FIELD_EXPs ;
      std::vector< size_t > L_FIELD_EXPs ;
      std::vector< std::string > TYPEs ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      
      std::vector< PDE_SystemNumbering* > NMBs ;
      
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* U ;
      LA_SeqVector* U_LOC ;

      LA_Solver* SOLVER ;
} ;

#endif
