#ifndef AP_CONSTITUTIVE_LAW_HH
#define AP_CONSTITUTIVE_LAW_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class doubleArray2D ;
class doubleArray4D ;

class AP_KinematicState ;

class AP_ConstitutiveLaw : public PEL_Object
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_ConstitutiveLaw* create( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) ;

   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const = 0 ;

   protected: //------------------------------------------------------

      virtual ~AP_ConstitutiveLaw( void ) ;

      AP_ConstitutiveLaw( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      bool update_S_dSdE_PRE( AP_KinematicState const* st,
                              doubleArray2D& Pio2,
                              doubleArray4D& dPio2dGL ) const ;

   private: //--------------------------------------------------------

      AP_ConstitutiveLaw( void ) ;
      AP_ConstitutiveLaw( AP_ConstitutiveLaw const& other ) ;
      AP_ConstitutiveLaw& operator=( AP_ConstitutiveLaw const& other ) ;

   //-- Attributes

} ;

#endif
