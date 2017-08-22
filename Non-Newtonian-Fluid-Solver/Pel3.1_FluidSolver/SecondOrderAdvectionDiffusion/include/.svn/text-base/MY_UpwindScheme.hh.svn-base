#ifndef IS_UPWIND_SCHEME_HH
#define IS_UPWIND_SCHEME_HH

#include <MY_AdvectiveScheme.hh>

/*

Upwind convective scheme.

The convective flux of a scalar `Phi' from a mesh K to a mesh L is :

    `FluxK->L' is the "upwind" flux defined as :

        if `vK->L'>0 :  `FuK->L' = `mK_L'*`rhoK_L'*`vK->L'*`PhiK'
        if `vK->L'<0 :  `FuK->L' = `mK_L'*`rhoK_L'*`vK->L'*`PhiL'

With :

    `mK_L' : area of the side K_L (m)
    
    `rhoK_L' : density of the side K_L (kg/m3)
    
    `vK->L' : normal velocity on the side from K to L (m/s)
    
Rems :

    1/ Without source terms :
    
          - the resulting matrix is a M_Matrix ;

          - the maximum principle is valid

    2/ `rhoK_L' is choosen according to the discretization used to solve the
       mass balance equation.

PUBLISHED
*/

class MY_UpwindScheme : public MY_AdvectiveScheme 
{
   public: //-----------------------------------------------------------------
      
   //-- Convective assembling
      
      virtual void assemble_convection_on_side( FE_TimeIterator const* t_it,
						PDE_CursorFEside* fe,
						PDE_LocalEquation* leq ) ;
       
      virtual void assemble_convection_on_bound(  FE_TimeIterator const* t_it,
						  PDE_LocalFEbound* fe,
						  double const bound_value,
						  PDE_LocalEquation* leq ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------
      
   //-- Constructors
      
      MY_UpwindScheme( void ) ;
     ~MY_UpwindScheme( void ) ;
      MY_UpwindScheme( MY_UpwindScheme const& other ) ;
      MY_UpwindScheme& operator=( MY_UpwindScheme const& other ) ;

      MY_UpwindScheme( PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms, 
                       PDE_DiscreteField const* field,
                       PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      MY_UpwindScheme( std::string const& a_name ) ;

      virtual MY_AdvectiveScheme* create_replica(
                       PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms, 
                       PDE_DiscreteField const* field,
                       PEL_ModuleExplorer const* exp ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Class attributes

      static MY_UpwindScheme const* PROTOTYPE ;
      
} ;

#endif
