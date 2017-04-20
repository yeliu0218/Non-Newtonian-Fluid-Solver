#ifndef AP_LOAD_CALCULATOR_TEST_HH
#define AP_LOAD_CALCULATOR_TEST_HH

#include <PEL_ObjectTest.hh>

/*
PUBLISHED
*/

class AP_LoadCalculator_TEST : public PEL_ObjectTest
{
   public: //-----------------------------------------------------------------

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      AP_LoadCalculator_TEST( void ) ;
     ~AP_LoadCalculator_TEST( void ) ;
      AP_LoadCalculator_TEST( AP_LoadCalculator_TEST const& other ) ;
      AP_LoadCalculator_TEST& operator=( 
                              AP_LoadCalculator_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void display_error( double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static AP_LoadCalculator_TEST const* PROTOTYPE ;

   //-- Attributes

} ;

#endif
