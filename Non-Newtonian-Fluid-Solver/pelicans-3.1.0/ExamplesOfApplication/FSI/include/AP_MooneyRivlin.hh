#ifndef AP_MOONEY_RIVLIN_HH
#define AP_MOONEY_RIVLIN_HH

#include <AP_ConstitutiveLaw.hh>

class AP_MooneyRivlin : public AP_ConstitutiveLaw
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_MooneyRivlin* create( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) ;

   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const ;
            
   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_MooneyRivlin( void ) ;
     ~AP_MooneyRivlin( void ) ;
      AP_MooneyRivlin( AP_MooneyRivlin const& other ) ;
      AP_MooneyRivlin& operator=( AP_MooneyRivlin const& other ) ;

      AP_MooneyRivlin( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Attributes

      double C1 ;
      double C2 ;
} ;

#endif
