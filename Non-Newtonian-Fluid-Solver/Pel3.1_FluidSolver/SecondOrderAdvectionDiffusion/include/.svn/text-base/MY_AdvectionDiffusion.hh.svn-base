#ifndef MY_ADVECTION_DIFFUSION_HH
#define MY_ADVECTION_DIFFUSION_HH

#include <FE_OneStepIteration.hh>

#include <vector>
#include <doubleVector.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PEL_ContextSimple ;
class PEL_DoubleVector ;
class PEL_Double ;

class GE_QRprovider ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

class FE_Parameter ;
class FE_SetOfParameters ;
class FE_LocalTimeIteratorAdapter;
class FE_TimeIterator;

class MY_AdvectiveScheme ;


/** \class MY_AdvectionDiffusion
    \brief
    Solves the advection diffusion equation:
    \f{eqnarray*}
        \frac{\partial \rho C}{\partial t} + \nabla\cdot \left( {\bf u}\ C\right) + D\Delta C  = f       
    \f}

    <table border="1">
    <td>
    <tr> AD_unknown_field = \f$ C \f$  </tr>
    <tr> AD_coeff_unsteady = \f$ \rho \f$        </tr>
    <tr> AD_coeff_diffusion = \f$ D \f$          </tr>
    <tr> AD_param_advective_velocity = \f$ {\bf u} \f$ </tr>
    <tr> AD_param_source = \f$ f \f$             </tr></td>
    </table>
    
  Corresponding HDS:
\code
MODULE FE_OneStepIteration#concentration
            concrete_name = "MY_AdvectionDiffusion"
            
            AD_unknown_field="CC"     
            AD_coeff_diffusion="diffusivity"
            AD_coeff_unsteady ="density"    
            AD_param_source="source"

            MODULE convective_scheme
               concrete_name = "MY_MUSCL_Scheme"
               AD_param_advective_velocity="velocity" 
               AD_coeff_unsteady ="density"
            END MODULE convective_scheme

            MODULE LA_Matrix
            END MODULE LA_Matrix
            
            MODULE LA_Solver
            END MODULE LA_Solver
            
         END MODULE FE_OneStepIteration#concentration
\endcode
    
*/

class MY_AdvectionDiffusion : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

    //-- Substeps of the step by step progression
        ///@name Substeps of the step by step progression
        //@{
        virtual void do_before_inner_iterations_stage( FE_TimeIterator const* t_it );
        virtual void do_inner_iterations_stage( FE_TimeIterator const* t_it );
        virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
        
        virtual void reset_discrete_problem( FE_TimeIterator const* t_it ) ;
	    //}@

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

       ~MY_AdvectionDiffusion( void ) ;
        MY_AdvectionDiffusion( MY_AdvectionDiffusion const& other ) ;
        MY_AdvectionDiffusion& operator=( MY_AdvectionDiffusion const& other ) ;

        MY_AdvectionDiffusion( PEL_Object* a_owner,
		PDE_DomainAndFields const* dom,
		FE_SetOfParameters const* prms,
		PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

        MY_AdvectionDiffusion( void ) ;

        virtual MY_AdvectionDiffusion* create_replica( PEL_Object* a_owner,
				      PDE_DomainAndFields const* dom,
				      FE_SetOfParameters const* prms,
				      PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

        void loop_on_cells( FE_TimeIterator const* t_it ) ;
        void loop_on_sides( FE_TimeIterator const* t_it ) ;
        void loop_on_bounds( FE_TimeIterator const* t_it ) ;
      
        void estimate_unknowns( void ) ;
        void update_fields( void ) ;	      

  //-- Class attributes

        static MY_AdvectionDiffusion const* PROTOTYPE ;

   //-- Attributes

      // Field : concentration
        PDE_DiscreteField* CC ;
        size_t L_UPDATE ;

      // Parameters
        FE_Parameter* DENSITY ; ///< @brief AD_coeff_unsteady     
        FE_Parameter* DIFFU ;   ///< @brief AD_coeff_diffusion
        FE_Parameter* SOURCE ;  ///< @brief AD_param_source 

      // Local equation assembling
        PDE_SetOfBCs const* BCs ;
        PEL_ContextSimple* CONTEXT ;
        PEL_DoubleVector* COORDS ;
	    PEL_Double* TT ;
        PDE_LocalEquation* ELEMENT_EQ ;
        PDE_LocalFEbound* bFE ;
        PDE_LocalFEcell* cFE ;
        PDE_CursorFEside* sFE ;

      // Global assembling
        PDE_SystemNumbering* NMB ;

      // Problem parameters :
        double const EPSILON ;

      // Convective scheme
        MY_AdvectiveScheme* AS ;


	    ///@name Courant Number
        //@{
        bool useCOURANT;              ///< @brief use Courant number: true or false
        PDE_SetOfDiscreteFields const* const FIELDS ;
        size_t LEVEL0 ;
        stringVector FIELDS_TABLE ;
        double COURANT;               ///< @brief Courant number
        //@}
        FE_LocalTimeIteratorAdapter* t_local;
        FE_TimeIterator* t_global;

        PEL_Object* A_OWNER;
        PDE_DomainAndFields const* DOM;
        FE_SetOfParameters const* PRMS;
        PEL_ModuleExplorer const* EXP;

        bool AD_inner_iterations_are_completed( FE_TimeIterator * t_global, FE_TimeIterator const* t_it, double epsilon ) const;
        double time_step_by_Courant( double dt_new, FE_TimeIterator const* t_it ) ;
        double smallest_vertices( void ) ;
        doubleVector gDT;
        doubleVector valOutlet; 
        doubleVector Xoutlet, Youtlet; 
        
        LA_Matrix* A ;
        LA_Vector* F ;
        LA_Vector* X ;
        LA_SeqVector* X_LOC ;
        LA_Solver* SOLVER ;
} ;

#endif
