#ifndef AP_1D_CHECK_HH
#define AP_1D_CHECK_HH

#include <PEL_Application.hh>

class FE_TimeIterator ;

class AP_1Dcheck : public PEL_Application
{
   public: //-----------------------------------------------------------------

      virtual void run( void ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_1Dcheck( void ) ;
      AP_1Dcheck( AP_1Dcheck const& other ) ;
      AP_1Dcheck& operator=( AP_1Dcheck const& other ) ;

      AP_1Dcheck( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_1Dcheck( void ) ;

      virtual AP_1Dcheck* create_replica( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      void run_cartesian( std::ostream& os ) ;

      void run_axisymetrical( std::ostream& os ) ;

      double relaxation_coefficient( double A_old, double A ) const ;

   //-- Class attributes

      static AP_1Dcheck const* PROTOTYPE ;

   //-- Attributes

      std::string GEOMETRY ;
      double a ;
      double c0 ;
      double b ;
      double RHO ;
      double VISC ; 
      double YOUNG ;
      double POIS ;
      double tauN ;
      double OMEGA ;
      double TOLERANCE ;
      double MY_EPS ;
      double MY_MIN ;
      size_t TIME_ORDER ;
      FE_TimeIterator* TIME_IT ;
      std::string OUTPUT_FILE ;
} ;

#endif
