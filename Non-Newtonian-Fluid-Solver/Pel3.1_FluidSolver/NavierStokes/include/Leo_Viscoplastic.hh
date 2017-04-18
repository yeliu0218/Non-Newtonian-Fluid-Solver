#ifndef Leo_Viscoplastic_HH
#define Leo_Viscoplastic_HH


#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_CrossProcessUnknownNumbering.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class FE_Parameter;
class FE_TimeIterator;
class PDE_DiscreteField ;
class PDE_LocalFEcell ;
class PDE_LocalEquation ;
class doubleVector ;


/** @brief The Class Leo_Viscoplastic.
 * Server for the resolution of a viscoplastic problem by an Augmented Lagrangian
 * algorithm. This is the non-dimensionalised version.
 * The Bingham number can be a general FE_Parameter class
 *
 * Conventions:
 *  - Anthony calculates with \f$ D(u):=\frac{1}{2} u_{i,j}+u_{j,i}\f$ form solid mechanic text books while I use
 *    \f$ \dot \gamma_{ij}(u) = u_{i,j}+u_{j,i} \f$, the usual convention form fluid mechanics text books.
 * 	  This change manifests in factors of \f$ \frac{1}{2} \f$ appearing and disappearing during the algorithm.
 *  - Norm of a tensor: Second Invariant
 * @author A. Wachs - Particulate flow project 2007-2009
 * @author A. Putz - CFL project
 *
 *
PUBLISHED
*/

class Leo_Viscoplastic
{
public: //-----------------------------------------------------------------

    static void compute_strain_rate_tensor( PDE_LocalFEcell* cFE,
                                            PDE_DiscreteField const* UU,
                                            PDE_DiscreteField* gammadot);
    static void update_Lagrange_multiplier( PDE_LocalFEcell* cFE,
                                            PDE_DiscreteField* LAMBDA,
                                            PDE_DiscreteField const* gammadot,
                                            PDE_DiscreteField const* gamma,
                                            double const& aug_param);
    static void update_strain_rate_tensor_d( PDE_LocalFEcell* cFE,
            PDE_DiscreteField const* LAMBDA,
            PDE_DiscreteField const* DUU,
            PDE_DiscreteField* d,
            double const& aug_param,
            FE_Parameter const* kappa,
            FE_Parameter const* Bn,
            double func_E,
			FE_TimeIterator const* t_it);

    static void update_strain_rate_tensor_d_HB( PDE_LocalFEcell* cFE,
            PDE_DiscreteField const* LAMBDA,
            PDE_DiscreteField const* DUU,
            PDE_DiscreteField* d,
            double const& aug_param,
            FE_Parameter const* kappa,
            FE_Parameter const* Bn,
            FE_Parameter const* HB_n,
            double func_E,
			FE_TimeIterator const* t_it);

    static double compute_norm_Dminusd( PDE_LocalFEcell* cFE,
                                        PDE_DiscreteField const* DUU,
                                        PDE_DiscreteField const* d);


    //-- Utilities

    /** @name Utilities */
    //@{
    static pair<size_t,double> strain_rate_component_number(
        size_t const &dimension,
        size_t const& compIdx, size_t const &direcIdx );
    static double Euclidian_norm( doubleVector const& tensor );
    //@}

protected: //--------------------------------------------------------------

private: //----------------------------------------------------------------

    //-- Constructors & Destructor

    /** @name Constructors & Destructor */
    //@{
    /** @brief Destructor */
    ~Leo_Viscoplastic( void ) {}

    /** @brief Copy constructor */
    Leo_Viscoplastic( Leo_Viscoplastic const& other ) {}

    /** @brief Constructor without argument */
    Leo_Viscoplastic( void ) {}
    //@}

} ;

#endif
