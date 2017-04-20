#ifndef AP_LAW_L2_TEST_HH
#define AP_LAW_L2_TEST_HH

#include <AP_ConstitutiveLaw.hh>

class AP_LawL2_TEST : public AP_ConstitutiveLaw
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_LawL2_TEST* create( PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) ;

   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const ;
      
   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_LawL2_TEST( void ) ;
     ~AP_LawL2_TEST( void ) ;
      AP_LawL2_TEST( AP_LawL2_TEST const& other ) ;
      AP_LawL2_TEST& operator=( AP_LawL2_TEST const& other ) ;

      AP_LawL2_TEST( PEL_Object* a_owner,
                     PEL_ModuleExplorer const* exp ) ;

   //-- Attributes

} ;

#endif
