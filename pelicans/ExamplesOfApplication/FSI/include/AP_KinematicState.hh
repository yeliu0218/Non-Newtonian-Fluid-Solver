#ifndef AP_KINEMATIC_STATE_HH
#define AP_KINEMATIC_STATE_HH

#include <PEL_Object.hh>

#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

class PDE_DiscreteField ;
class PDE_LocalFE ;

class AP_KinematicState : public PEL_Object
{
   public: //---------------------------------------------------------

   //-- Instance delivery and initialization

      static AP_KinematicState* create( PEL_Object* a_owner,
                                        size_t nb_disp_dims ) ;

   //-- Object configuration

      void set_state( doubleArray2D const& g_disp ) ;

      void set_state( PDE_DiscreteField const* uu,
                      size_t level,
                      PDE_LocalFE const* fe ) ;

   //-- Kinematic data

      size_t nb_displacement_dimensions( void ) const ;

      doubleArray2D const& grad_disp( void ) const ;

      // the deformation gradient
      doubleArray2D const& DG( void ) const ;

      // the determinant of the deformation gradient
      double detF( void ) const ;

      // the inverse of the deformation gradient
      doubleArray2D const& inv_F( void ) const ;

      // the right Cauchy-Green tensor
      doubleArray2D const& rCG( void ) const ;

      // the inverse of the right Cauchy-Green tensor
      doubleArray2D const& inv_rCG( void ) const ;

      // the Green-Lagrange tensor
      doubleArray2D const& GL( void ) const ;

   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      AP_KinematicState( void ) ;
     ~AP_KinematicState( void ) ;
      AP_KinematicState( AP_KinematicState const& other ) ;
      AP_KinematicState& operator=( AP_KinematicState const& other ) ;

      AP_KinematicState( PEL_Object* a_owner, size_t nb_disp_dims ) ;

   //-- Internals

      void compute_F_J( void ) const ;

      void compute_Fm1( void ) const ;

      void compute_C( void ) const ;

      void compute_Cm1( void ) const ;

      void compute_E( void ) const ;

      static double determinant( doubleArray2D const& mat,
                                 size_t mat_size ) ;

   //-- Attributes

      size_t NB_DISP_DIMS ;

      doubleArray2D grad_DISP ;

      mutable bool CALC_F_J ;
      mutable doubleArray2D F ;
      mutable double J ;

      mutable bool CALC_Fm1 ;
      mutable doubleArray2D Fm1 ;

      mutable bool CALC_C ;
      mutable doubleArray2D C ;

      mutable bool CALC_Cm1 ;
      mutable doubleArray2D Cm1 ;

      mutable bool CALC_E ;
      mutable doubleArray2D E ;
} ;

#endif
