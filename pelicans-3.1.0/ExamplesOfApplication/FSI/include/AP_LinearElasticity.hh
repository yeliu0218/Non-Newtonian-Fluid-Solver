#ifndef AP_LINEAR_ELASTICITY_HH
#define AP_LINEAR_ELASTICITY_HH

#include <AP_ConstitutiveLaw.hh>

class AP_LinearElasticity : public AP_ConstitutiveLaw
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_LinearElasticity* create( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) ;


   //-- Stress and elasticiy tensors

      virtual void update_S_dSdE( AP_KinematicState const* st,
                                  doubleArray2D& Pio2,
                                  doubleArray4D& dPio2dGL ) const ;

   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_LinearElasticity( void ) ;
     ~AP_LinearElasticity( void ) ;
      AP_LinearElasticity( AP_LinearElasticity const& other ) ;
      AP_LinearElasticity& operator=( AP_LinearElasticity const& other ) ;

      AP_LinearElasticity( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp ) ;

   //-- Utilities
    
      static double kronecker( size_t i, size_t j ) ;

      static double I( size_t i, size_t j, size_t k, size_t l ) ;

      static double J( size_t i, size_t j, size_t k, size_t l ) ;

      static double K( size_t i, size_t j, size_t k, size_t l ) ; 

      static double CauchyTensor( doubleArray2D const& grad_disp, 
                                  size_t i, size_t j ) ;

      double compute_rigidity_matrix( size_t i, size_t j, size_t k, 
                                      size_t l, double young,
                                      double poisson ) const ;

   //-- Attributes

      double YOUNG ;
      double POISSON ;
} ;

#endif
