#ifndef MY_ADVECTION_DIFFUSION_HH
#define MY_ADVECTION_DIFFUSION_HH

#include <FE_OneStepIteration.hh>

#include <vector>
#include <doubleVector.hh>
#include <intVector.hh>

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


/** \class MY_PLIC3_Advection
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
            concrete_name = "MY_PLIC3_Advection"
            
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

class MY_PLIC3_Advection : public FE_OneStepIteration
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

       ~MY_PLIC3_Advection( void ) ;
        MY_PLIC3_Advection( MY_PLIC3_Advection const& other ) ;
        MY_PLIC3_Advection& operator=( MY_PLIC3_Advection const& other ) ;

        MY_PLIC3_Advection( PEL_Object* a_owner,
		PDE_DomainAndFields const* dom,
		FE_SetOfParameters const* prms,
		PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

        MY_PLIC3_Advection( void ) ;

        virtual MY_PLIC3_Advection* create_replica( PEL_Object* a_owner,
				      PDE_DomainAndFields const* dom,
				      FE_SetOfParameters const* prms,
				      PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building
        void setup_connectivity();
        void compute_gradient_C(size_t ic);
        void loop_on_cellsXX( FE_TimeIterator const* t_it, size_t ic ) ;
        void loop_on_cellsZZ( FE_TimeIterator const* t_it, size_t ic ) ;
        void regularize_C(FE_TimeIterator const* t_it);
        //void estimate_unknowns( void ) ;
        void update_fields( void ) ;
        double vol2( double mx, double mz, double alfa, double b) ;
        void print_mesh( FE_TimeIterator const* t_it )	;      

  //-- Class attributes

        static MY_PLIC3_Advection const* PROTOTYPE ;

   //-- Attributes

      // Field : concentration
        PDE_DiscreteField* CC ;
        PDE_DiscreteField* UU ;
        size_t L_UPDATE ;
        
      // Local equation assembling
        PDE_SetOfBCs const* BCs ;
        PEL_ContextSimple* CONTEXT ;
        PEL_DoubleVector* COORDS ;
	    PEL_Double* TT ;
        PDE_LocalFEcell* cFE ;
        PDE_CursorFEside* sFE ;
        PDE_LocalFEbound* bFE;
      // Problem parameters :
        double const EPSILON ;
        //size_t NX;
        //size_t NY;
        //doubleVector Xvector;
        //doubleVector Yvector;
        //doubleVector Xcoord;
        intVector UP;
        intVector DOWN;
        intVector LEFT;
        intVector RIGHT;
        doubleVector GRADx;
        doubleVector GRADz;
        doubleVector C1;
        doubleVector C2;
        doubleVector C3;
        double Cmax;
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
        void update_interface(double diff, size_t ic) ;


} ;

#endif
