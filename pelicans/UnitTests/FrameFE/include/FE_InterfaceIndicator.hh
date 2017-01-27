#ifndef FE_INTERFACE_INDICATOR_HH
#define FE_INTERFACE_INDICATOR_HH

#include <PDE_AdaptationIndicator.hh>

#include <doubleVector.hh>

class PEL_Context ;
class PEL_DataWithContext ;
class PEL_DoubleVector ;
class PEL_Int ;
class PEL_ModuleExplorer ;

class GE_Mpolyhedron ;
class GE_Point ;
class GE_QRprovider ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_ReferenceElement ;
class PDE_SetOfBCs ;

class FE_InterfaceIndicator : public PDE_AdaptationIndicator
{
   public: //-----------------------------------------------------------

   //-- Indicator calculation

      virtual void reset( void ) ;

      virtual void build( void ) ;

   //-- Indicator results

      virtual double cell_indicator( size_t cell_id ) const ;

      virtual bool to_be_refined( double bf_indicator,
                                  GE_Mpolyhedron const* poly,
                                  PDE_ReferenceElement const* elm,
                                  size_t local_node ) const ;

      virtual bool to_be_unrefined( double bf_indicator,
                                    GE_Mpolyhedron const* poly,
                                    PDE_ReferenceElement const* elm,
                                    size_t local_node ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_InterfaceIndicator( void ) ;
      FE_InterfaceIndicator( FE_InterfaceIndicator const& other ) ;
      FE_InterfaceIndicator& operator=( 
                             FE_InterfaceIndicator const& other ) ;

      FE_InterfaceIndicator( PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             PEL_ModuleExplorer const* exp,
                             size_t a_verbose_level ) ;

   //-- Plug in

      FE_InterfaceIndicator( void ) ;

      virtual FE_InterfaceIndicator* create_replica( 
                                            PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            PEL_ModuleExplorer const* exp,
                                            size_t a_verbose_level ) const ;

   //-- Class attributes

      static FE_InterfaceIndicator const* PROTOTYPE ;

   //-- Attributes

      PEL_Context* CTX ;
      mutable PEL_DoubleVector* COORDS ;
      mutable PEL_Int* ITER ;
      PEL_DataWithContext* PHF ;
      PDE_LocalFEcell* cFE ;
      GE_QRprovider const* QRP ;
      doubleVector CELL_INTERF ;
      double H_INTERF ;
      double BF_MIN_REFI ;
      double BF_MAX_REFI ;
      double BF_MIN_UNREFI ;
      double BF_MAX_UNREFI ;
      size_t ICALL ;
} ;

#endif
