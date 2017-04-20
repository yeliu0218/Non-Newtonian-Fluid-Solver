#ifndef PEL_RUN_TEST_DEMO_1_HH
#define PEL_RUN_TEST_DEMO_1_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;

class PEL_RunTest_DEMO1 : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~PEL_RunTest_DEMO1( void ) ;
      PEL_RunTest_DEMO1( PEL_RunTest_DEMO1 const& other ) ;
      PEL_RunTest_DEMO1& operator=( PEL_RunTest_DEMO1 const& other ) ;

      PEL_RunTest_DEMO1( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_RunTest_DEMO1( void ) ;

      virtual PEL_RunTest_DEMO1* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static PEL_RunTest_DEMO1 const* PROTOTYPE ;

   //-- Attributes

      std::string OUT_NAME ;
} ;

#endif
