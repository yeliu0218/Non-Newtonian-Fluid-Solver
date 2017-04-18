#ifndef UT_Viscoelastic_HH
#define UT_Viscoelastic_HH


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
class LA_SeqVector ;
class LA_SymmetricMatrix ;
class PDE_DiscreteField ;
class PDE_LocalFEcell ;
class PDE_LocalEquation ;
class doubleVector ;


/** @brief The Class UT_Viscoelastic.
 *
PUBLISHED
*/

class UT_Viscoelastic
{
public: //-----------------------------------------------------------------


    static double compute_norm_DminusGAMMADOT( PDE_LocalFEcell* cFE,
                                        PDE_DiscreteField const* DD,
                                        PDE_DiscreteField const* GAMMADOT);

    static void elastic_model_PPT(LA_SeqVector* c_minus_I,
								  LA_SeqVector* const eigenvals,
								  double const relax_time,
								  double const eps ) ;

protected: //--------------------------------------------------------------

private: //----------------------------------------------------------------

    //-- Constructors & Destructor

    /** @name Constructors & Destructor */
    //@{
    /** @brief Destructor */
    ~UT_Viscoelastic( void ) {}

    /** @brief Copy constructor */
    UT_Viscoelastic( UT_Viscoelastic const& other ) {}

    /** @brief Constructor without argument */
    UT_Viscoelastic( void ) {}
    //@}

} ;

#endif
