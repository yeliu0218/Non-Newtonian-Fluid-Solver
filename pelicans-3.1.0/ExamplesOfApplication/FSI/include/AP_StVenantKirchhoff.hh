#ifndef AP_ST_VENANT_KIRCHHOFF_HH
#define AP_ST_VENANT_KIRCHHOFF_HH

#include <AP_ConstitutiveLaw.hh>

class AP_StVenantKirchhoff : public AP_ConstitutiveLaw
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_StVenantKirchhoff* create( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp ) ;

   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const ;

   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_StVenantKirchhoff( void ) ;
     ~AP_StVenantKirchhoff( void ) ;
      AP_StVenantKirchhoff( AP_StVenantKirchhoff const& other ) ;
      AP_StVenantKirchhoff& operator=( AP_StVenantKirchhoff const& other ) ;

      AP_StVenantKirchhoff( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) ;

   //-- Attributes

      double YOUNG ;
      double POISSON ;
} ;

#endif
