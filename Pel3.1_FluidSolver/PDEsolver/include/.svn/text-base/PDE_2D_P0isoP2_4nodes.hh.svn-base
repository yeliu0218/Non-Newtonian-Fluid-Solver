#ifndef PDE_2D_P0_ISO_P2_4_NODES_HH
#define PDE_2D_P0_ISO_P2_4_NODES_HH

#include <PDE_ReferenceElement.hh>

/** @brief The Class PDE_2D_P0isoP2_4nodes

Triangular finite element divided into 4 sub-triangles on which the basis 
function is constant (polynomial of order 0).

This construction leads to the following characteristics :
   Dimension         : 2
   Geometrical shape : triangle
   Polynomial degree : 0
   number of nodes   : 4

@author A. Wachs - Particulate flow project 2007-2009 */

class PDE_2D_P0isoP2_4nodes : public PDE_ReferenceElement
{

   public: //-----------------------------------------------------------

   //-- Basis functions

      /** @name Basis functions */
      //@{
      /** @brief Return the basis function at a given point
      @param node basis function number
      @param pt_ref pointer to the point */ 
      virtual double N_local( size_t node, GE_Point const* pt_ref ) const ;

      /** @brief Return the derivative of a basis function at a given point
      @param node basis function number
      @param a derivative direction
      @param pt_ref pointer to the point */      
      virtual double dN_local( size_t node,
                               size_t a,
                               GE_Point const* pt_ref ) const ;

      /** @brief Return the second derivative of a basis function at a given 
      point
      @param node basis function number
      @param a first derivative direction
      @param b second derivative direction      
      @param pt_ref pointer to the point */      
      virtual double d2N_local( size_t node,
                                size_t a, 
                                size_t b,
                                GE_Point const* pt_ref ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */
      PDE_2D_P0isoP2_4nodes( void ) ;
      
      /** @brief Destructor */
      ~PDE_2D_P0isoP2_4nodes( void ) ;
      
      /** @brief Copy constructor 
      @param other the element to be copied */ 
      PDE_2D_P0isoP2_4nodes( PDE_2D_P0isoP2_4nodes const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side element */ 
      PDE_2D_P0isoP2_4nodes& operator=( 
                             PDE_2D_P0isoP2_4nodes const& other ) ;
      //@}      

   //-- Class attributes

      static PDE_2D_P0isoP2_4nodes const* REGISTRATOR ;
} ;

#endif
