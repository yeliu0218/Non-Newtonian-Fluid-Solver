#ifndef AP_CONSTITUTIVE_LAW_TEST_HH
#define AP_CONSTITUTIVE_LAW_TEST_HH

#include <PEL_ObjectTest.hh>

/*
PUBLISHED
*/

class AP_ConstitutiveLaw_TEST : public PEL_ObjectTest
{
   public: //-----------------------------------------------------------------

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      AP_ConstitutiveLaw_TEST( void ) ;
     ~AP_ConstitutiveLaw_TEST( void ) ;
      AP_ConstitutiveLaw_TEST( AP_ConstitutiveLaw_TEST const& other ) ;
      AP_ConstitutiveLaw_TEST& operator=( 
                              AP_ConstitutiveLaw_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void display_error( double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static AP_ConstitutiveLaw_TEST const* PROTOTYPE ;

   //-- Attributes

} ;

#endif
