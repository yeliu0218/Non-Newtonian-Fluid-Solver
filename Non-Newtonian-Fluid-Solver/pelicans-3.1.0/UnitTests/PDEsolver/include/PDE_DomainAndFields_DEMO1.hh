#ifndef PDE_DOMAIN_AND_FIELDS_DEMO_1_HH
#define PDE_DOMAIN_AND_FIELDS_DEMO_1_HH

#include <PEL_Application.hh>

class PDE_DomainAndFields ;
class PEL_ModuleExplorer ;

class PDE_DomainAndFields_DEMO1 : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~PDE_DomainAndFields_DEMO1( void ) ;
      PDE_DomainAndFields_DEMO1( PDE_DomainAndFields_DEMO1 const& other ) ;
      PDE_DomainAndFields_DEMO1& operator=( 
                               PDE_DomainAndFields_DEMO1 const& other ) ;

      PDE_DomainAndFields_DEMO1( PEL_Object* a_owner, 
                               PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PDE_DomainAndFields_DEMO1( void ) ;

      virtual PDE_DomainAndFields_DEMO1* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static PDE_DomainAndFields_DEMO1 const* PROTOTYPE ;

   //-- Attributes

      PDE_DomainAndFields const* DOM ;
      std::string OUT_NAME ;
} ;

#endif
