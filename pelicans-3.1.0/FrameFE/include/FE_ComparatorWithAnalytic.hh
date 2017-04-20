/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#ifndef FE_COMPARE_WITH_ANALYTIC_HH
#define FE_COMPARE_WITH_ANALYTIC_HH

#include <FE_OneStepIteration.hh>

#include <doubleVector.hh>
#include <intVector.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

class doubleVector ;

class GE_Mpolyhedron ;
class GE_Point ;
class GE_QRprovider ;

class PEL_Communicator ;
class PEL_DataWithContext ;
class PEL_Double ;
class PEL_DoubleVector ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;

/*
Applications involving the L_inf, Lp, H1 or H1_0 norm comparison between a 
computed solution and a given analytical one.

The computed norms are defined with the `norms' table in the data deck.
In the next section, S referes to the analytical solution, Sh to the 
projection of the analytical solution in the finite element space, 
and Xh to the solution computed by the program.

1/ L_inf norm :

    - `Infinity_error_norm' :
          max(S(i)-Xh(i)) for all the nodes i of the finite element 
                          discretization
    - `Infinity_solution_norm' :
          max(S(i)) for all the nodes i of the finite element discretization

2/ Lp norm (with a decent value for p):

    - `Lp_solution_norm' :
          [ int{ S(x)^p } ]^(1/p), integration for a given
                                   quadrature rule provider
    - `Lp_interpolation_error_norm' :
          [ int{ (S(x)-Sh(x))^p } ]^(1/p), integration for a given 
                                          quadrature rule provider
    - `Lp_error_norm' :
          [ int{ (S(x)-Xh(x))^p } ]^(1/p), integration for a given 
                                           quadrature rule provider
    - `Lp_error_D_norm' :
          [ int{ (Sh(x)-Xh(x))^p } ]^(1/p), finite volume discrete Lp norm, 
                   obtained through a one point integration at the finite 
                   volume center 

3/ H1_0 semi norm :

    - `H1_0_solution_norm' or `sW12_solution_norm' :
          sqrt[ int{ (dS(x)/dx)^2 } ], integration for a given 
                                       quadrature rule provider
    - `H1_0_interpolation_error_norm' or `sW12_interpolation_error_norm' :
          sqrt[ int{ (dS(x)/dx-dSh(x)/dx)^2 } ], integration for a given 
                                                 quadrature rule provider
    - `H1_0_error_norm' or `sW12_error_norm' :
          sqrt[ int{ (dS(x)/dx-dXh(x)/dx)^2 } ], integration for a given 
                                                 quadrature rule provider
    - `H1_D_error_D_norm' :
          discrete H1 semi-norm of the error Sh-Xh (corresponding to Neumann
          boundary conditions)
          reference : Eymard et al., Finite Volume Methods, in Handbook of
                      Numerical Analysis VII, Editors Ciarlet and Lions,
                      definition 10.2
          meaningful for discrete fields that are constant on each cell
    - `H1_Dirichlet_D_error_D_norm' :
          discrete H1_gamma semi-norm of the error Sh-Xh (corresponding to
          Dirichlet boundary conditions on a subset gamma of the frontier)
          reference : Eymard et al., Finite Volume Methods, in Handbook of
                      Numerical Analysis VII, Editors Ciarlet and Lions,
                      definition 9.3
          meaningful for discrete fields that are constant on each cell

4/ H1 norm :

    - `H1_solution_norm' or `W12_solution_norm' :
          sqrt[ int{ S(x)^2 + (dS(x)/dx)^2) } ],
                integration for a given quadrature rule provider
    - `H1_interpolation_error_norm' or `W12_interpolation_error_norm' :
          sqrt[ int{ (S(x)-Sh(x))^2) + (dS(x)/dx-dSh(x)/dx)^2 } ],
                integration for a given quadrature rule provider
    - `H1_error_norm' or `W12_error_norm' :
          sqrt[ int{ (S(x)-Xh(x))^2 + (dS(x)/dx-dXh(x)/dx)^2 } ],
                integration for a given quadrature rule provider
          
5/ W1p semi-norms (with a decent value for p) :

   - `sW1p_solution_norm'
          [ int{ (dS(x)/dx)^p } ]^(1/p), integration for a given 
                                         quadrature rule provider
   - `sW1p_interpolation_error_norm'
          [ int{ (dS(x)/dx-dSh(x)/dx)^p } ]^(1/p), integration for a given 
                                                   quadrature rule provider
   - `sW1p_error_norm'
          [ int{ (dS(x)/dx-dXh(x)/dx)^p } ]^(1/p), integration for a given 
                                                   integration rule provider
   
6/ W1p norms (with a decent value for p) :

   - `W1p_solution_norm'
          [ int{ S(x)^p + (dS(x)/dx)^p) } ]^(1/p),
            integration for a given quadrature rule provider
   - `W1p_interpolation_error_norm'
          [ int{ (S(x)-Sh(x))^p) + (dS(x)/dx-dSh(x)/dx)^p } ]^(1/p),
            integration for a given quadrature rule provider
   - `W1p_error_norm'
          [ int{ (S(x)-Xh(x))^p + (dS(x)/dx-dXh(x)/dx)^p } ]^(1/p),
            integration for a given quadrature rule provider
   
Rem1 : The norm values, the time step, and the meshing size are stored in the
       saving file (key `norm_saving_names' for the norms) ;
Rem2 : for field which are defined to within a constant (example of the 
       pressure in some cases), this constant is defined 
       (key `nullify_integral' in the data deck)
       such as :
                         int{ S(x) }=0   and   int{ Xh(x) }=0

Example :

MODULE PEL_Application

   $DS_NU   = $DS_VISC/$DS_RHO
   $DS_expT  = exp( -2.*sqr( pi() )*$DS_NU*$DS_T )
   $DS_exp2T  = exp( -4.*sqr( pi() )*$DS_NU*$DS_T )
   $DS_x = component( $DV_X, 0 )
   $DS_y = component( $DV_X, 1 )
   $DS_Vx = -cos( pi()*$DS_x )*sin( pi()*$DS_y )*$DS_expT
   $DS_Vy =  sin( pi()*$DS_x )*cos( pi()*$DS_y )*$DS_expT
   $DS_dVxdx =  pi()*sin( pi()*$DS_x )*sin( pi()*$DS_y )*$DS_expT
   $DS_dVxdy = -pi()*cos( pi()*$DS_x )*cos( pi()*$DS_y )*$DS_expT
   $DS_dVydx =  pi()*cos( pi()*$DS_x )*cos( pi()*$DS_y )*$DS_expT
   $DS_dVydy = -pi()*sin( pi()*$DS_x )*sin( pi()*$DS_y )*$DS_expT
   $DS_P =  -0.25*( cos( 2.*pi()*$DS_x )+cos( 2.*pi()*$DS_y ) )*$DS_exp2T

   ...

   MODULE FE_OneStepIteration

      concrete_name = "FE_SplitSystem"

      MODULE list_of_FE_OneStepIteration

         ...

         MODULE FE_OneStepIteration#velocity_solution_storage
            concrete_name = "FE_ComparatorWithAnalytic"
            field = "velocity"
            level = 0
            solution = vector( $DS_Vx, $DS_Vy )
            d_solution = array( vector( $DS_dVxdx, $DS_dVxdy ),
                                vector( $DS_dVydx, $DS_dVydy ) )
            quadrature_rule_provider = "GE_QRprovider_5"
            norms = < "L2_error_norm" "H1_error_norm" >
            norm_saving_names = < "EVL2" "EVH1" >
         END MODULE FE_OneStepIteration#velocity_solution_storage
         
         MODULE FE_OneStepIteration#pressure_solution_storage
            concrete_name = "FE_ComparatorWithAnalytic"
            field = "pressure"
            level = 0
            solution = vector( $DS_P )
            quadrature_rule_provider = "GE_QRprovider_1"
            norms = < "L2_error_norm" >
            norm_saving_names = < "EPL2" >
            nullify_integral = true
         END MODULE FE_OneStepIteration#pressure_solution_storage

      END MODULE list_of_FE_OneStepIteration           

      ...
      
   END MODULE FE_OneStepIteration

   ...
END MODULE PEL_Application

Output :
--------

By default, the computed norms are saved via the PDE_ResultSaver object
of the associated `PDE_DomainAndFields' instance.
Those norms may also be saved in a text column file, whose name is the 
data of keyword "output_file".

By default, the saving times correspond to the 
`::save_other_than_time_and_fields' times. It is possible to specify 
different saving times with the data of keyword "saving_times".

Example :

MODULE PEL_Application
   ...
   MODULE FE_OneStepIteration
      concrete_name = "FE_SplitSystem"
      ...
      MODULE list_of_FE_OneStepIteration
         ...
         MODULE FE_OneStepIteration#ca
            concrete_name = "FE_ComparatorWithAnalytic"
            ...
            output_file = "error_norms.txt"
            saving_times = < 0.01 0.02 >
         END MODULE FE_OneStepIteration#ca
      END MODULE list_of_FE_OneStepIteration           
   END MODULE FE_OneStepIteration
END MODULE PEL_Application

PUBLISHED
*/

class PEL_EXPORT FE_ComparatorWithAnalytic : public FE_OneStepIteration
{ 

   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual void do_after_time_adaptation( FE_TimeIterator const* t_it ) ;
      
   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields(
                                            FE_TimeIterator const* t_it, 
                                            PDE_ResultSaver* rs ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      enum NormType
      {
         Invalid,
         
         sol_Linf,
         err_Linf,
         
         sol_Lp,
         eip_Lp,
         err_Lp,
         
         sol_sW1p,
         eip_sW1p,
         err_sW1p,
         
         sol_W1p,
         err_W1p,
         eip_W1p,
         
         err_LpD,
         err_H1D,
         err_H1DirD
      } ;
      
     ~FE_ComparatorWithAnalytic( void ) ;
      FE_ComparatorWithAnalytic( FE_ComparatorWithAnalytic const& other ) ;
      FE_ComparatorWithAnalytic& operator=(
                              FE_ComparatorWithAnalytic const& other ) ;
      
      FE_ComparatorWithAnalytic( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer* exp ) ;
      
   //-- Plug in

      FE_ComparatorWithAnalytic( void ) ;
      
      virtual FE_ComparatorWithAnalytic* create_replica(
                         PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         FE_SetOfParameters const* prms,
                         PEL_ModuleExplorer* exp ) const ;
   //-- Internals

      virtual void save_norms( FE_TimeIterator const* t_it, 
                               PDE_ResultSaver* rs ) ;
      
      virtual void compute_norms( FE_TimeIterator const* t_it,
                                  doubleVector& norm_values ) ;
      
      // compute the infinity norm of
      // * the solution : `sol_norm'
      // * the approximation error : `err_norm'.
      void calc_Linf( double sol_shift, double app_shift,
                      double& sol_norm, double& err_norm ) ;
      
      // compute the p-th power of the Lp norm of
      // * the solution : `sol_norm'
      // * the interpolation error for the solution : `eip_norm'
      // * the approximation error : `err_norm'.
      void calc_Lp(  double sol_shift, double app_shift, double exponent,
                     double& sol_norm, 
                     double& eip_norm, 
                     double& err_norm ) ;
      
      // compute the p-th power of the W1p seminorm of
      // * the solution : `sol_snorm'
      // * the interpolation error for the solution : `eip_snorm'
      // * the approximation error : `err_snorm'.
      void calc_sW1p( double exponent,
                     double& sol_snorm, 
                     double& eip_snorm, 
                     double& err_snorm ) ;

      // compute the p-th power of the Lp discrete norm of
      // * the approximation error : `err_norm'.
      void calc_LpD( double sol_shift, double app_shift, double exponent,
                     double& err_norm ) ;
      
      void calc_H1D( double& err_norm ) ;
      
      void calc_H1DirD( double& err_norm ) ;
      
      void compute_shift( double& sol_shift, double& app_shift ) ;

      double max_inter_vertices_distance_of_cells( void ) ;
      double max_equivalent_ball_diameter_of_cells( void ) ;

      /* center used for Lp discrete norm
         (finite volume center if exists and is located in the polyhedron,
         cell barycenter elsewhere) */
      GE_Point const* D_norm_center( GE_Mpolyhedron const* poly ) const ;

      static void read_norm( std::string const& norm_name, 
                             NormType& nt, double& exponent,
                             PEL_ModuleExplorer const* exp ) ;
      
      static void parse_norm_name( std::string const& norm_name,
                                   std::string const& start,
                                   double& pp, std::string& end ) ;
      
  //-- Exponents 
      
      void add_exponent( double pp ) ;
      
      size_t nb_exponents( void ) const ;
      
      size_t idx_of_exponent( double pp ) const ;
      
  //-- Static attributes

      static FE_ComparatorWithAnalytic const* PROTOTYPE ;

  //-- Class attributes

      // Field :
      PDE_DiscreteField const* const FIELD ;
      size_t const LEVEL ;

      // Iterators :
      PDE_LocalFEcell* const cFE ;
      PDE_CursorFEside* const sFE ;
      PDE_LocalFEbound* const bFE ;

      GE_QRprovider const* QRP ;
      PDE_SetOfBCs const* BCs ;

      // Norms :
      // * EXPOS( i ) is the exponent associated to NORMS( i )
      // * each element of IDX_EXPOS is different
      intVector NORMS ;
      doubleVector EXPOS ;
      doubleVector IDX_EXPOS ;
      stringVector const NORM_NAMES ;
      stringVector const NORM_SAVE_NAMES ;
      size_t MAX_LENN ;

      // Context :
      PEL_DoubleVector* COORDS ;
      PEL_Double* TT ;
      
      // Solution :
      PEL_DataWithContext const* SOL ;
      PEL_DataWithContext const* DSOL ;

      // Preprocessing :
      bool const NULLIFY_INTEGRAL ;

      // D norm :
      double const DIST_REF ;

      // Communicator :
      PEL_Communicator const* const COM ;

      // Output
      std::string FILE_OUT ;
      size_t SAVING_NUMBER ;
      double NEXT_SAVING_TIME ;
      doubleVector SAVING_TIMES ;
} ;

#endif
