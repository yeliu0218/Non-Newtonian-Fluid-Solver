#ifndef MY_MUSCL_SCHEME_HH
#define MY_MUSCL_SCHEME_HH

#include <MY_AdvectiveScheme.hh>

class MY_MUSCL_DataStructure ;

/*

MUSCL convective scheme for a scalar `Phi'.

Two steps are necessary to obtain the fluxes at the interface
`K|L' of cells `K' and `L' :

- First step : reconstruction of `Phi' at the cell interfaces

   `Phi_{K|L}^-' = `Phi_K' + 1/2 `Psi({Delta Phi_K^+}/{{Delta Phi_K^-}})' * `Delta Phi_K^-'
   `Phi_{K|L}^+' = `Phi_L' - 1/2 `Psi({Delta Phi_L^+}/{{Delta Phi_L^-}})' * `Delta Phi_K^+'

   where `Delta Phi_K^+' = `Delta Phi_L^-' = `Phi_L' - `Phi_K'

- Second step : 

    `FluxK->L' is the "upwind" flux defined as :

        if `vK->L'>0 :  `FluxK->L' = `mK_L'*`rhoK_L'*`vK->L'*`Phi_{K|L}^-'
        if `vK->L'<0 :  `FluxK->L' = `mK_L'*`rhoK_L'*`vK->L'*`Phi_{K|L}^+'

With :

    `mK_L' : area of the side K_L (m)
    
    `rhoK_L' : density of the side K_L (kg/m3)
    
    `vK->L' : normal velocity on the side from K to L (m/s)
    
PUBLISHED
*/

class MY_MUSCL_Scheme : public MY_AdvectiveScheme 
{
   public: //-----------------------------------------------------------------
      
   //-- Convective assembling
      
      virtual void assemble_convection_on_side(
                           FE_TimeIterator const* t_it,
                           PDE_CursorFEside* fe,
                           PDE_LocalEquation* leq ) ;
       
      virtual void assemble_convection_on_bound(
                           FE_TimeIterator const* t_it,
                           PDE_LocalFEbound* fe,
                           double const bound_value,
                           PDE_LocalEquation* leq ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------
      
   //-- Constructors
      
      MY_MUSCL_Scheme( void ) ;
     ~MY_MUSCL_Scheme( void ) ;
      MY_MUSCL_Scheme( MY_MUSCL_Scheme const& other ) ;
      MY_MUSCL_Scheme& operator=( MY_MUSCL_Scheme const& other ) ;

      MY_MUSCL_Scheme( PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms, 
                       PDE_DiscreteField const* field,
                       PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      MY_MUSCL_Scheme( std::string const& a_name ) ;

      virtual MY_AdvectiveScheme* create_replica(
                       PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms, 
                       PDE_DiscreteField const* field,
                       PEL_ModuleExplorer const* exp ) const ;

   //-- Implementation
      virtual double limit_slope_ratio(
                           double numerator_slope,
                           double denominator_slope,
                           double smallest_nonzero_den=1.E-10 ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Class attributes

      static MY_MUSCL_Scheme const* PROTOTYPE ;
      
   //-- Attributes
 
      MY_MUSCL_DataStructure const* MDS ; 
} ;

#endif
