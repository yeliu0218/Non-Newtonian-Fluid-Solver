#ifndef PDE_2D_P1DC_3_NODES_HH
#define PDE_2D_P1DC_3_NODES_HH

#include <PDE_ReferenceElement.hh>

/** @brief The Class PDE_2D_P0isoP2_4nodes

Triangular finite element with linear discontinuous basis functions.

This construction leads to the following characteristics :
   Dimension         : 2
   Geometrical shape : triangle
   Polynomial degree : 1
   number of nodes   : 3

@author A. Wachs - Particulate flow project 2007-2009 */

class PDE_2D_P1DC_3nodes : public PDE_ReferenceElement
{

   public: //-----------------------------------------------------------

   //-- Basis functions

      virtual double N_local( size_t node, GE_Point const* pt_ref ) const ;

      virtual double dN_local( size_t node,
                               size_t a,
                               GE_Point const* pt_ref ) const ;

      virtual double d2N_local( size_t node,
                                size_t a, 
                                size_t b,
                                GE_Point const* pt_ref ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_2D_P1DC_3nodes( void ) ;
     ~PDE_2D_P1DC_3nodes( void ) ;
      PDE_2D_P1DC_3nodes( PDE_2D_P1DC_3nodes const& other ) ;
      PDE_2D_P1DC_3nodes const& operator=( PDE_2D_P1DC_3nodes const& other ) ;

   //-- Class attributes

      static PDE_2D_P1DC_3nodes* uniqueInstance ;

} ;

#endif
