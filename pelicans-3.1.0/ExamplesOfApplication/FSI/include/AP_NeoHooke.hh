#ifndef AP_NEO_HOOKE_HH
#define AP_NEO_HOOKE_HH

#include <AP_ConstitutiveLaw.hh>

class AP_NeoHooke : public AP_ConstitutiveLaw
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_NeoHooke* create( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) ;

   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const ;
      
   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_NeoHooke( void ) ;
     ~AP_NeoHooke( void ) ;
      AP_NeoHooke( AP_NeoHooke const& other ) ;
      AP_NeoHooke& operator=( AP_NeoHooke const& other ) ;

      AP_NeoHooke( PEL_Object* a_owner,
                   PEL_ModuleExplorer const* exp ) ;

   //-- Attributes

      double C1 ;
} ;

#endif
