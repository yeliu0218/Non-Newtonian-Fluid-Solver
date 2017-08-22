#ifndef GE_COLOR_DEMO_1_HH
#define GE_COLOR_DEMO_1_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;

class GE_Color_DEMO1 : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~GE_Color_DEMO1( void ) ;
      GE_Color_DEMO1( GE_Color_DEMO1 const& other ) ;
      GE_Color_DEMO1& operator=( GE_Color_DEMO1 const& other ) ;

      GE_Color_DEMO1( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      GE_Color_DEMO1( void ) ;

      virtual GE_Color_DEMO1* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static GE_Color_DEMO1 const* PROTOTYPE ;

   //-- Attributes

      std::string OUT_NAME ;
} ;

#endif
