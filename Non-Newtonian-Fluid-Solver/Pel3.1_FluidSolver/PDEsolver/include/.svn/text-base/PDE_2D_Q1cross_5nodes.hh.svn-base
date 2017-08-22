#ifndef PDE_2D_Q1cross_5nodes_HH
#define PDE_2D_Q1cross_5nodes_HH

#include <PDE_ReferenceElement.hh>

/*
Lagrange reference finite element with the following characteristics :
   Dimension         : 2
   Geometrical shape : quadrilateral
   Polynomial degree : 1
   number of nodes   : 5

   criss-cross grid
*/

class PDE_2D_Q1cross_5nodes : public PDE_ReferenceElement
{

   public: //-----------------------------------------------------------

   //-- Basis functions

      virtual double N_local( size_t node,
                              GE_Point const* pt_ref ) const ;

      virtual double dN_local( size_t node,
                               size_t a,
                               GE_Point const* pt_ref ) const ;

      virtual double d2N_local( size_t node,
                                size_t a,
                                size_t b,
                                GE_Point const* pt_ref ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_2D_Q1cross_5nodes( void ) ;
     ~PDE_2D_Q1cross_5nodes( void ) ;
      PDE_2D_Q1cross_5nodes( PDE_2D_Q1cross_5nodes const& other ) ;
      PDE_2D_Q1cross_5nodes& operator=( PDE_2D_Q1cross_5nodes const& other ) ;

   //-- Class attributes

      static PDE_2D_Q1cross_5nodes* uniqueInstance ;

} ;

#endif


