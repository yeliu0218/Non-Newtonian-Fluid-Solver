#ifndef FE_TensorFormAssembling_HH
#define FE_TensorFormAssembling_HH

#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_CrossProcessUnknownNumbering.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


class PDE_DiscreteField ;
class PDE_LocalFEcell ;
class PDE_LocalEquation ;


/** @brief The Class FE_TensorFormAssembling.

Some matrix formula that are not available in Pelicans.

@author A. Wachs - Particulate flow project 2007-2009 */

class FE_TensorFormAssembling
{
   public: //-----------------------------------------------------------------

   //-- Matrix assembling routines

      /** @name Matrix assembling routines */
      //@{
      /** @brief Assemble the elementary matrix at the level of the element
      for the term div(lambda) where lambda is a tensor stored as a vector
      @param leq elementary matrix
      @param fe local finite element
      @param coef parameter */
      static void add_grad_row_col( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe,
		double coef ) ;
      //by Leo
      static void add_grad_row_col_A( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe,
		double coef ) ;
      static void add_grad_row_col_Leo( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe,
		double coef ) ;
      static void add_row_col( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe,
		double coef ) ;

      /** @brief Assemble the elementary mass matrix at the level of the element
      using a lumping formula
      @param leq elementary matrix
      @param fe local finite element
      @param coef parameter */
      static void add_lumped_row_col( PDE_LocalEquation* leq,
                                      PDE_LocalFEcell const* fe,
                                      double coef ) ;
      //@}

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */
      ~FE_TensorFormAssembling( void ) {}

      /** @brief Copy constructor */
      FE_TensorFormAssembling( FE_TensorFormAssembling const& other ) {}

      /** @brief Constructor without argument */
      FE_TensorFormAssembling( void ) {}
      //@}

} ;

#endif
