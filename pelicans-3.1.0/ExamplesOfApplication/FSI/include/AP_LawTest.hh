#ifndef AP_LAW_TEST_HH
#define AP_LAW_TEST_HH

#include <PEL_Application.hh>

class doubleVector ;
class doubleArray2D ;
class doubleArray4D ;
class AP_ConstitutiveLaw ;

class AP_LawTest : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_LawTest( void ) ;
      AP_LawTest( AP_LawTest const& other ) ;
      AP_LawTest& operator=( AP_LawTest const& other ) ;

      AP_LawTest( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_LawTest( void ) ;

      virtual AP_LawTest* create_replica( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static AP_LawTest const* PROTOTYPE ;

   //-- Attributes

      AP_ConstitutiveLaw* LAW ;

      std::string OUTPUT_FILE ;
} ;

#endif
