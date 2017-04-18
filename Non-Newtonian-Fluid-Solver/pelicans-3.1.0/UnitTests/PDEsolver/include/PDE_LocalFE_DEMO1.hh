#ifndef PDE_LOCAL_FE_DEMO_1_HH
#define PDE_LOCAL_FE_DEMO_1_HH

#include <PEL_Application.hh>

class PEL_ModuleExplorer ;

class PDE_DomainAndFields ;

class PDE_LocalFE_DEMO1 : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~PDE_LocalFE_DEMO1( void ) ;
      PDE_LocalFE_DEMO1( PDE_LocalFE_DEMO1 const& other ) ;
      PDE_LocalFE_DEMO1& operator=( PDE_LocalFE_DEMO1 const& other ) ;

      PDE_LocalFE_DEMO1( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PDE_LocalFE_DEMO1( void ) ;

      virtual PDE_LocalFE_DEMO1* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      void do_1( std::ostream& os ) ;
      void do_2( std::ostream& os ) ;
      void do_3( std::ostream& os ) ;
      void do_4( std::ostream& os ) ;
      void do_5( std::ostream& os ) ;
      void do_6( std::ostream& os ) ;
      void do_7( std::ostream& os ) ;
      void do_8( std::ostream& os ) ;

   //-- Class attributes

      static PDE_LocalFE_DEMO1 const* PROTOTYPE ;

   //-- Attributes

      PDE_DomainAndFields const* DOM ;
      std::string OUT_NAME ;
} ;

#endif
