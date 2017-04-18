#ifndef AP_KINEMATIC_STATE_TEST_HH
#define AP_KINEMATIC_STATE_TEST_HH

#include <PEL_ObjectTest.hh>

class doubleArray2D ;

/*
PUBLISHED
*/

class AP_KinematicState_TEST : public PEL_ObjectTest
{
   public: //-----------------------------------------------------------------

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      AP_KinematicState_TEST( void ) ;
     ~AP_KinematicState_TEST( void ) ;
      AP_KinematicState_TEST( AP_KinematicState_TEST const& other ) ;
      AP_KinematicState_TEST& operator=( 
                              AP_KinematicState_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void display_error( double xx_1, double xx_2 ) const ;

      bool ok_rCG_inv_rCG( doubleArray2D const& rcg,
                           doubleArray2D const& inv_rcg,
                           size_t dim ) const ;

      bool ok_rCG_GL( doubleArray2D const& rcg,
                      doubleArray2D const& gl,
                      size_t dim ) const ;

      bool ok_F_inv_F( doubleArray2D const& F,
                      doubleArray2D const& inv_F,
                      size_t dim ) const ;

   //-- Class attributes

      static AP_KinematicState_TEST const* PROTOTYPE ;

   //-- Attributes

      double D_EPS ;
      double D_MIN ;

} ;

#endif
